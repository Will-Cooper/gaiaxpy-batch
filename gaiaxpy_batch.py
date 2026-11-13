# coding: utf-8
"""
Batch downloading of Gaia DR3 XP spectra, calibrated (or converted) with GaiaXPy.

GaiaXPy sends every source ID it is given to the Gaia archive in a single DataLink request, logging in each time,
and the archive limits each request to 50 (anonymous) or 5000 (logged in) source IDs.
This script downloads the XP continuous spectra with astroquery in chunks, in parallel threads that each log in once,
retries server-side failures with increasing waits, and passes the downloaded coefficients to GaiaXPy.
The spectra are saved to csv and (optionally) per-object fits/txt files.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed, wait as futures_wait
from contextlib import redirect_stdout
from getpass import getpass
from io import StringIO
import logging
import os
import re
import shutil
import socket
import sys
import tempfile
import threading
from typing import Callable, Dict, List, Optional, Tuple, Union

from astropy.io import fits
from astropy.table import Table
from astroquery import log as astroquery_log
from astroquery.gaia import GaiaClass
import gaiaxpy
from gaiaxpy import calibrate, convert, pwl_to_wl, wl_to_pwl
from gaiaxpy.core.server import data_release, gaia_server
import numpy as np
import pandas as pd
from tqdm import tqdm

__version__ = '2.0.2'  # the single source of the version, read by pyproject.toml

log = logging.getLogger('gaiaxpy_batch')

SAMPLING_FUNCTIONS = ('linspace', 'arange', 'geomspace')
DEFAULT_PWL = np.linspace(0, 60, 600)  # default pseudo-wavelength sampling of gaiaxpy.convert
WL_LIMITS = (330., 1050.)  # absolute wavelength [nm] limits accepted by gaiaxpy.calibrate
PWL_LIMITS = (-10., 70.)  # pseudo-wavelength limits accepted by gaiaxpy.convert

# Gaia archive user quotas: https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas
ANONYMOUS_CHUNKSIZE = 50  # DataLink sources per request, anonymous
LOGGED_IN_CHUNKSIZE = 1000  # the logged in quota is 5000, but smaller requests lose less when one fails
LOGGED_IN_MAX_WORKERS = 5  # active jobs allowed per registered user
RETRIES = 3  # retries of a request after a server-side failure
RETRY_DELAY = 30.  # seconds to wait before the first retry, doubled for each one after
TIMEOUT = 900.  # seconds without a response before a connection is given up
ABORT_AFTER = 3  # chunks in a row failing every retry before the archive is taken to be unavailable

ID_PATTERN = re.compile(r'^\s*(?:Gaia\s*)?(?:DR3\s*)?(\d{1,19})\s*$', re.IGNORECASE)  # 123 or Gaia DR3 123
UNSAFE_PATTERN = re.compile(r'[^\w.+-]+')  # characters replaced by _ in file names
LIMIT_PATTERN = re.compile(r'Maximum number of Source Ids reached \(max: (\d+)')  # archive's request size error
STATUS_PATTERN = re.compile(r'\b(?:HTTP Error|Error|status code)\s*:?\s*(\d{3})\b', re.IGNORECASE)
FATAL_ERRORS = (AttributeError, ImportError, NameError, TypeError)  # bugs, not archive failures: no point retrying

# reasons given in the failed csv
NO_XP = 'no XP spectrum'
DOWNLOAD_FAILED = 'download failed'
NOT_ATTEMPTED = 'not attempted'  # the downloads were stopped (archive unavailable, or interrupted)
CALIBRATION_FAILED = 'calibration failed'
INVALID_ID = 'invalid source ID'
INCOMPLETE = (DOWNLOAD_FAILED, NOT_ATTEMPTED)  # reasons worth trying again later


def parse_sampling(func: str, start: float, end: float, step: float) -> np.ndarray:
    """Wavelength sampling [nm] from linspace/geomspace (step = number of points) or arange (step in nm)"""
    if func not in SAMPLING_FUNCTIONS:
        raise ValueError(f'Sampling function must be one of {SAMPLING_FUNCTIONS}, not {func!r}')
    start, end, step = float(start), float(end), float(step)
    if func == 'arange':
        return np.arange(start, end, step)
    return getattr(np, func)(start, end, int(step))


def sysargs(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Command line arguments, named as the keywords of batch"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--filename', required=True, metavar='Filename',
                        help='Input CSV file containing the Gaia DR3 Source IDs')
    parser.add_argument('-s', '--sampling', nargs=4, metavar=('numpy function', 'start', 'end', 'step'),
                        help='Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120" '
                             f'(functions: {", ".join(SAMPLING_FUNCTIONS)})')
    parser.add_argument('-t', '--truncate', action='store_true',
                        help='Use only the basis functions GaiaXPy recommends for each source (smoother, less noisy '
                             'spectra)')
    parser.add_argument('-o', '--outputstyle', choices=('fits', 'txt'),
                        help='Output of produced spectra (default: None)')
    parser.add_argument('-i', '--idcolname', default='source_id', metavar='Source ID column name',
                        help='Name of the column relating to DR3 Source ID')
    parser.add_argument('-n', '--namecol', metavar='Object name column name',
                        help='Name of the column relating to the object name (for saving spectra, '
                             'default: GaiaDR3_<source_id>)')
    parser.add_argument('-u', '--uncalibrated', action='store_true',
                        help='Get the internally calibrated BP or RP spectrum (gaiaxpy.convert, electron/s) '
                             'instead of the calibrated one (gaiaxpy.calibrate, W m-2 nm-1)')
    parser.add_argument('-x', '--xp', type=str.upper, choices=('BP', 'RP'),
                        help='With -u, which spectrum: BP or RP (default: RP)')
    parser.add_argument('-c', '--chunksize', type=int,
                        help=f'Number of source IDs per archive request (default: {ANONYMOUS_CHUNKSIZE} anonymous, '
                             f'{LOGGED_IN_CHUNKSIZE} logged in; the archive allows up to 50 and 5000)')
    parser.add_argument('-w', '--workers', default=4, type=int,
                        help=f'Number of parallel archive requests (default: 4, at most {LOGGED_IN_MAX_WORKERS} '
                             f'when logged in)')
    parser.add_argument('--retries', default=RETRIES, type=int,
                        help=f'Retries of a request after a server-side failure (default: {RETRIES})')
    parser.add_argument('--retry-delay', default=RETRY_DELAY, type=float,
                        help=f'Seconds to wait before the first retry, doubled for each one after '
                             f'(default: {RETRY_DELAY:.0f})')
    parser.add_argument('--timeout', default=TIMEOUT, type=float,
                        help=f'Seconds without a response before a connection is given up (default: {TIMEOUT:.0f})')
    parser.add_argument('-a', '--anonymous', action='store_true',
                        help='Query the archive without logging in (otherwise uses $GAIA_USERNAME and '
                             '$GAIA_PASSWORD, or prompts)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print details of every request, login and failure')
    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__} (GaiaXPy {gaiaxpy.__version__})')
    args = parser.parse_args(argv)
    if (args.chunksize is not None and args.chunksize < 1) or args.workers < 1 or args.retries < 0 \
            or args.retry_delay < 0 or args.timeout <= 0:
        parser.error('-c and -w must be at least 1, --retries and --retry-delay at least 0, --timeout above 0')
    if args.sampling is not None:
        args.sampling = parse_sampling(*args.sampling)
    return args


def str2float(s: str) -> np.ndarray:
    """An array of the output csv's flux or flux_error column, from its text"""
    return np.array(s.strip('[]').split(), dtype=float)


def get_credentials(username: Optional[str] = None, password: Optional[str] = None,
                    anonymous: bool = False) -> Tuple[Optional[str], Optional[str]]:
    """Gaia archive username and password from the arguments, then $GAIA_USERNAME/$GAIA_PASSWORD, then a prompt"""
    if anonymous:
        return None, None
    try:
        username = username or os.environ.get('GAIA_USERNAME') or input('Gaia archive username: ')
        password = password or os.environ.get('GAIA_PASSWORD') or getpass('Gaia archive password: ')
    except (EOFError, OSError):  # no terminal to ask on, e.g. run with nohup or from a script
        raise ValueError('No Gaia archive credentials: set GAIA_USERNAME and GAIA_PASSWORD, '
                         'or use -a to download anonymously') from None
    return username, password


def get_sampling(sampling: Optional[np.ndarray] = None, uncalibrated: bool = False,
                 xp: Optional[str] = 'RP') -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Validates the wavelength sampling [nm] and returns (wavelengths in GaiaXPy's output order, sampling for GaiaXPy):
    nm for calibrate (both None for GaiaXPy's default), ascending pseudo-wavelength for convert
    """
    if sampling is not None and not len(sampling):
        raise ValueError('The sampling has no points, e.g. arange with a step that is not positive')
    if not uncalibrated:
        if sampling is None:
            return None, None
        wave = np.asarray(sampling, dtype=float)
        if np.any(np.diff(wave) <= 0) or wave[0] < WL_LIMITS[0] or wave[-1] > WL_LIMITS[1]:
            raise ValueError(f'Sampling must be ascending and within {WL_LIMITS[0]:.0f}-{WL_LIMITS[1]:.0f} nm')
        return wave, wave
    xp = (xp or 'RP').upper()
    if xp not in ('BP', 'RP'):
        raise ValueError(f'The spectrum must be BP or RP, not {xp!r}')
    if sampling is None:
        pwl = DEFAULT_PWL.copy()
        wave = pwl_to_wl(xp, pwl)
    else:
        wave = np.asarray(sampling, dtype=float)
        pwl = wl_to_pwl(xp, wave)
        order = np.argsort(pwl)  # BP pseudo-wavelength decreases with wavelength, GaiaXPy needs it ascending
        pwl, wave = pwl[order], wave[order]
    if pwl[0] < PWL_LIMITS[0] or pwl[-1] > PWL_LIMITS[1]:
        raise ValueError(f'Sampling {wave.min():.0f}-{wave.max():.0f} nm is outside the {xp} pseudo-wavelength '
                         f'range {PWL_LIMITS}')
    return wave, pwl


def calibrate_wrap(coeffs: pd.DataFrame, **kwargs) -> Tuple[pd.DataFrame, np.ndarray]:
    """GaiaXPy calibrate (or convert, for kwargs xp) of downloaded spectra: (source_id/flux/flux_error, sampling)"""
    sampling = kwargs.get('sampling', None)
    truncate = kwargs.get('truncate', False)
    if not kwargs.get('uncalibrated', False):
        dfout, sampling = calibrate(coeffs, sampling=sampling, truncation=truncate, save_file=False)
    else:
        dfout, sampling = convert(coeffs, sampling=sampling, truncation=truncate, save_file=False)
        dfout = dfout.loc[dfout.xp == kwargs['xp']].drop(columns='xp')
    return dfout[['source_id', 'flux', 'flux_error']], sampling


class ArchiveError(Exception):
    """An archive request failed on every retry"""


class ArchiveUnavailable(ArchiveError):
    """The archive appears to be down, or the login failed: stop downloading"""


def http_status(e: Exception) -> Optional[int]:
    """HTTP status code in an archive exception's message (astroquery puts it only there), if any"""
    match = STATUS_PATTERN.search(str(e))
    return int(match.group(1)) if match else None


def is_transient(e: Exception) -> bool:
    """
    Whether an archive failure is probably temporary and worth retrying: server errors (HTTP 5xx), timeouts,
    lost sessions and broken responses (e.g. 'File is not a zip file'), rather than a bad request (other HTTP 4xx)
    """
    status = http_status(e)
    return status is None or status >= 500 or status in (401, 403, 408, 429)


def _short(e: Exception) -> str:
    """The first lines of an exception, for log messages, with HTML error pages reduced to their title and message"""
    text = str(e)
    if '<html' in text.lower():  # e.g. Tomcat's "HTTP Status 401 – Unauthorized ... Message Credentials expiration"
        found = [re.search(pattern, text, re.S | re.I) for pattern in (r'<title>(.*?)</title>',
                                                                       r'<b>Message</b>\s*(.*?)</p>')]
        text = ': '.join(match.group(1).strip() for match in found if match) or text.split('<')[0].strip()
    lines = [line for line in text.splitlines() if line.strip()]
    return f'{type(e).__name__}: {" ".join(lines[:2])}'[:200]


class Archive:
    """
    Downloads XP continuous spectra from the Gaia archive with astroquery, using one session per thread
    (each logged in once; username None for anonymous), and retrying server-side failures with increasing waits
    """

    def __init__(self, username: Optional[str] = None, password: Optional[str] = None, retries: int = RETRIES,
                 delay: float = RETRY_DELAY):
        if retries < 0 or delay < 0:
            raise ValueError(f'retries ({retries}) and delay ({delay}) cannot be negative')
        self.username = username
        self.password = password
        self.retries = retries
        self.delay = delay
        self.stop = threading.Event()  # set when the archive is taken to be unavailable
        self._local = threading.local()

    def session(self) -> GaiaClass:
        """The astroquery session of this thread, created (and logged in) on first use"""
        gaia = getattr(self._local, 'gaia', None)
        if gaia is None:
            gaia = GaiaClass(gaia_tap_server=gaia_server, gaia_data_server=gaia_server, show_server_messages=False)
            if self.username:
                try:
                    self.retry(lambda: self._login(gaia), 'Login')
                except ArchiveError as e:
                    raise ArchiveUnavailable(f'Could not log in to the Gaia archive as {self.username}: '
                                             f'check the username and password, or the archive status') from e
            self._local.gaia = gaia
        return gaia

    def _login(self, gaia: GaiaClass):
        """Logs a session in, raising if it did not work (astroquery logs login errors rather than raising them)"""
        gaia.login(user=self.username, password=self.password)
        tap = getattr(gaia, '_TapPlus__isLoggedIn', None)  # astroquery internals, None if they change
        data = getattr(getattr(gaia, '_GaiaClass__gaiadata', None), '_TapPlus__isLoggedIn', None)
        if tap is False or data is False:
            raise ConnectionError('Gaia archive login failed')

    def retry(self, func: Callable, what: str):
        """Calls func (described by what, for log messages), retrying while it fails in a way that seems temporary"""
        for attempt in range(self.retries + 1):
            if self.stop.is_set():
                raise ArchiveUnavailable('Stopped, the Gaia archive appears to be unavailable')
            try:
                return func()
            except FATAL_ERRORS + (ArchiveError,):
                raise
            except Exception as e:
                if LIMIT_PATTERN.search(str(e)) or not is_transient(e):
                    raise  # handled by the caller
                if attempt == self.retries:
                    raise ArchiveError(f'{what} failed {attempt + 1} times, lastly with {_short(e)}') from e
                wait = self.delay * 2 ** attempt
                if http_status(e) in (401, 403):  # the session has expired (after about 30 min), log in again now
                    self._local.gaia = None
                    wait = 0.
                log.warning(f'{what} failed ({_short(e)}), retrying in {wait:.0f} s ({attempt + 1}/{self.retries})')
                self.stop.wait(wait)

    def load(self, idlist: List[int]) -> Optional[pd.DataFrame]:
        """XP continuous spectra of the source IDs, downloaded as GaiaXPy does; None if none of them have any"""
        def attempt():
            result = self.session().load_data(ids=idlist, format='csv', data_release=data_release,
                                              data_structure='DATAMODEL_GAIA', retrieval_type='XP_CONTINUOUS',
                                              avoid_datatype_check=True)
            keys = [key for key in result if 'continuous' in key.lower()]
            return result[keys[0]][0].to_pandas() if keys else None
        return self.retry(attempt, f'Download of {len(idlist)} sources starting {idlist[0]}')


Results = Tuple[List[pd.DataFrame], List[np.ndarray], Dict[int, str]]  # spectra, their samplings, failed IDs


def _merge(results: List[Results]) -> Results:
    """Combines the results of several pieces of a chunk"""
    frames, samplings, failed = [], [], {}
    for pieceframes, piecesamplings, piecefailed in results:
        frames.extend(pieceframes)
        samplings.extend(piecesamplings)
        failed.update(piecefailed)
    return frames, samplings, failed


def _calibrate_rows(coeffs: pd.DataFrame, **kwargs) -> Results:
    """Calibrates downloaded spectra, splitting them in half and retrying each half on failure (no new download)"""
    try:
        dfout, sampling = calibrate_wrap(coeffs, **kwargs)
    except FATAL_ERRORS:
        raise
    except Exception as e:
        if len(coeffs) == 1:
            xid = int(coeffs.source_id.iloc[0])
            log.debug(f'Calibration of {xid} failed: {_short(e)}')
            return [], [], {xid: CALIBRATION_FAILED}
        half = (len(coeffs) + 1) // 2
        return _merge([_calibrate_rows(coeffs.iloc[:half], **kwargs), _calibrate_rows(coeffs.iloc[half:], **kwargs)])
    # sources GaiaXPy returned no spectrum for: e.g. the BP spectrum of a source with only an RP one is None
    usable = dfout.flux.map(_usable).astype(bool)
    dfout = dfout[usable]
    returned = set(dfout.source_id.astype(np.int64))
    reason = f'no {kwargs["xp"]} spectrum' if kwargs.get('uncalibrated') else CALIBRATION_FAILED
    failed = {int(xid): reason for xid in coeffs.source_id if int(xid) not in returned}
    return ([dfout, ], [sampling, ], failed) if len(dfout) else ([], [], failed)


def _usable(flux) -> bool:
    """Whether a spectrum returned by GaiaXPy has any values (it is None for a missing band)"""
    if flux is None:
        return False
    arr = np.asarray(flux, dtype=float)
    return arr.ndim == 1 and bool(np.isfinite(arr).any())


def _process(idlist: List[int], archive: Archive, state: dict, pbar: tqdm, **kwargs) -> Results:
    """
    Downloads and calibrates a chunk of source IDs.
    If the archive says the chunk is too large, the chunk is split to the archive's limit and the limit is kept for
    the other chunks. A chunk failing with a bad request is split in half and each half tried again.
    Chunks that fail every retry are given up, and after several in a row the archive is taken to be unavailable.
    state is shared between threads: the archive's limit ('limit'), chunks in a row failing every retry
    ('exhausted'), why the downloads were stopped ('abort') and a lock ('lock')
    """
    def pieces(size: int) -> Results:
        return _merge([_process(idlist[i:i + size], archive, state, pbar, **kwargs)
                       for i in range(0, len(idlist), size)])

    maxids = state['limit']
    if maxids is not None and len(idlist) > maxids:
        return pieces(maxids)
    try:
        coeffs = archive.load(idlist)
    except FATAL_ERRORS:
        raise
    except ArchiveUnavailable as e:
        with state['lock']:
            if not archive.stop.is_set():
                archive.stop.set()
                state['abort'] = str(e)
        pbar.update(len(idlist))
        return [], [], dict.fromkeys(idlist, NOT_ATTEMPTED)
    except ArchiveError as e:
        log.error(str(e))
        with state['lock']:
            state['exhausted'] += 1
            if state['exhausted'] >= ABORT_AFTER and not archive.stop.is_set():
                archive.stop.set()
                state['abort'] = (f'{ABORT_AFTER} requests in a row failed every retry, '
                                  f'the Gaia archive appears to be unavailable')
        pbar.update(len(idlist))
        return [], [], dict.fromkeys(idlist, DOWNLOAD_FAILED)
    except Exception as e:
        toomany = LIMIT_PATTERN.search(str(e))
        if toomany and int(toomany.group(1)) < len(idlist):
            maxids = int(toomany.group(1))
            with state['lock']:
                if state['limit'] != maxids:
                    state['limit'] = maxids
                    log.warning(f'The Gaia archive allows {maxids} source IDs per request, using chunks of {maxids}')
            return pieces(maxids)
        # a bad request, perhaps because of one of the source IDs
        log.debug(f'Download of {len(idlist)} sources starting {idlist[0]} failed: {_short(e)}')
        if len(idlist) == 1:
            pbar.update(1)
            return [], [], {idlist[0]: DOWNLOAD_FAILED}
        return pieces((len(idlist) + 1) // 2)
    with state['lock']:
        state['exhausted'] = 0
    if coeffs is None or not len(coeffs):
        pbar.update(len(idlist))
        return [], [], dict.fromkeys(idlist, NO_XP)
    found = set(coeffs.source_id.astype(np.int64))
    frames, samplings, failed = _calibrate_rows(coeffs, **kwargs)
    failed.update({xid: NO_XP for xid in idlist if xid not in found})
    pbar.update(len(idlist))
    return frames, samplings, failed


def download_xp(xpids: List[int], chunksize: Optional[int] = None, workers: int = 4,
                username: Optional[str] = None, password: Optional[str] = None, retries: int = RETRIES,
                retry_delay: Optional[float] = None, timeout: float = TIMEOUT,
                **kwargs) -> Tuple[pd.DataFrame, np.ndarray, Dict[int, str]]:
    """
    Downloads XP spectra for a list of source IDs in chunks, in parallel threads, and calibrates them with GaiaXPy
    (kwargs: sampling, truncate, uncalibrated, xp; retry_delay None for RETRY_DELAY).
    Returns (source_id/flux/flux_error, ascending wavelengths [nm], {source ID: reason} for those without a spectrum).
    While downloading, the working directory is a temporary one (astroquery writes files to it), the socket timeout
    is timeout, and stdout is hidden (GaiaXPy and astroquery print progress bars and raw error pages to it)
    """
    xpids = list(dict.fromkeys(int(xid) for xid in xpids))  # unique, keeping order
    if not xpids:
        raise ValueError('No source IDs given')
    kwargs['xp'] = (kwargs.get('xp') or 'RP').upper()
    wave, gaiaxpy_sampling = get_sampling(kwargs.get('sampling', None), kwargs.get('uncalibrated', False),
                                          kwargs['xp'])
    kwargs['sampling'] = gaiaxpy_sampling
    if username and workers > LOGGED_IN_MAX_WORKERS:
        log.warning(f'The Gaia archive allows {LOGGED_IN_MAX_WORKERS} active jobs per user, '
                    f'using {LOGGED_IN_MAX_WORKERS} workers')
        workers = LOGGED_IN_MAX_WORKERS
    if chunksize is None:  # the default, made smaller if needed to give every worker a chunk
        chunksize = min(LOGGED_IN_CHUNKSIZE if username else ANONYMOUS_CHUNKSIZE, -(-len(xpids) // workers))
    xchunk = [xpids[i:i + chunksize] for i in range(0, len(xpids), chunksize)]
    workers = min(workers, len(xchunk))
    log.info(f'{len(xpids)} sources in {len(xchunk)} chunks of up to {chunksize}, using {workers} workers, '
             f'{"logged in as " + username if username else "anonymously"}')
    archive = Archive(username, password, retries, RETRY_DELAY if retry_delay is None else retry_delay)
    state = dict(limit=None, exhausted=0, abort=None, lock=threading.Lock())
    res, samplings, failed = [], [], {}
    oldtimeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(timeout)  # astroquery connections have no timeout
    # astroquery writes (and deletes, unless killed) a temp_<time> directory in the working directory for every
    # request, so work in a temporary directory instead of the user's
    oldcwd = os.getcwd()
    workdir = tempfile.mkdtemp(prefix='gaiaxpy_batch_')
    os.chdir(workdir)
    collected = set()

    def collect(future):
        chunkres, chunksamp, chunkfailed = future.result()
        res.extend(chunkres)
        samplings.extend(chunksamp)
        failed.update(chunkfailed)
        collected.add(future)

    try:
        with redirect_stdout(StringIO()), \
                tqdm(total=len(xpids), desc='Downloading XP spectra', unit='source', file=sys.stderr) as pbar:
            pool = ThreadPoolExecutor(max_workers=workers)
            futures = [pool.submit(_process, chunk, archive, state, pbar, **kwargs) for chunk in xchunk]
            try:
                for future in as_completed(futures):
                    collect(future)
            except KeyboardInterrupt:  # stop starting requests, keep what has finished
                archive.stop.set()
                state['abort'] = 'Interrupted, the remaining sources were not attempted'
                pool.shutdown(wait=False, cancel_futures=True)
                running = [future for future in futures if future.running()]
                log.warning(f'Interrupted, waiting for {len(running)} requests in progress to finish '
                            f'(Ctrl-C again to quit without saving)')
                try:  # cancelled futures never count as done for wait()
                    futures_wait([future for future in futures if not future.cancelled()])
                except KeyboardInterrupt:
                    os._exit(130)
                for future in futures:
                    if future not in collected and not future.cancelled():
                        collect(future)
            finally:
                pool.shutdown(wait=False)
    finally:
        socket.setdefaulttimeout(oldtimeout)
        os.chdir(oldcwd)
        shutil.rmtree(workdir, ignore_errors=True)
    if state['abort']:
        log.error(state['abort'])
    if not res:
        if state['abort']:
            raise ArchiveUnavailable(state['abort'])
        reasons = list(failed.values())
        summary = ', '.join(f'{reasons.count(reason)} {reason}' for reason in sorted(set(reasons)))
        if any(reason in INCOMPLETE for reason in reasons):
            raise ArchiveError(f'No XP spectra downloaded, {summary}')
        raise ValueError(f'No XP spectra found, {summary}')
    stackdf = pd.concat(res, ignore_index=True).drop_duplicates(subset='source_id')  # in case the archive repeats
    if wave is None:  # calibrate with default sampling, only known now
        wave = np.asarray(samplings[0], dtype=float)
    order = np.argsort(wave)
    if np.any(order != np.arange(len(wave))):  # make the wavelengths ascending (BP in uncalibrated mode)
        wave = wave[order]
        for col in ('flux', 'flux_error'):
            stackdf[col] = [np.asarray(arr)[order] for arr in stackdf[col]]
    found = set(stackdf.source_id.astype(np.int64))
    unknown = NOT_ATTEMPTED if state['abort'] else DOWNLOAD_FAILED  # e.g. chunks cancelled by an interruption
    failed = {xid: failed.get(xid, unknown) for xid in xpids if xid not in found}  # input order
    return stackdf, wave, failed


def output_names(fname: str) -> Tuple[str, str, str]:
    """The full, cut (only sources with spectra) and failed csv names for an input csv"""
    stem = os.path.splitext(fname)[0] + '_XP_spectra'
    return stem + '.csv', stem + '_cut.csv', stem + '_failed.csv'


def read_ids(values: pd.Series) -> pd.Series:
    """
    Gaia DR3 source IDs (Int64, missing where invalid) from a column read as strings (floats lose the last digits
    of 19 digit IDs), as plain integers or names such as 'Gaia DR3 123'
    """
    def parse(value):
        match = None if pd.isna(value) else ID_PATTERN.match(str(value))
        return int(match.group(1)) if match and int(match.group(1)) < 2 ** 63 else pd.NA
    return pd.Series([parse(value) for value in values], index=values.index, dtype='Int64')


def validate_inputs(fname: str, idcolname: str = 'source_id', namecol: Optional[str] = None,
                    outputstyle: Optional[str] = None, sampling: Optional[np.ndarray] = None,
                    uncalibrated: bool = False, xp: Optional[str] = None):
    """Checks the input csv and options before anything is downloaded, raising a ValueError or OSError if wrong"""
    if not os.path.isfile(fname):
        raise FileNotFoundError(f'Input file {fname} not found')
    if os.path.splitext(fname)[1].lower() != '.csv':
        raise ValueError(f'Input file {fname} must be a csv (comma separated values) file ending .csv')
    columns = list(pd.read_csv(fname, nrows=0).columns)
    for column, flag in ((idcolname, '-i'), (namecol, '-n')):
        if column is not None and column not in columns:
            raise ValueError(f'Column {column!r} not found in {fname}, which has the columns {columns}; '
                             f'choose the column with {flag}')
    if outputstyle not in (None, 'fits', 'txt'):
        raise ValueError(f'The output style must be fits or txt, not {outputstyle!r}')
    if xp is not None and not uncalibrated:
        log.warning(f'-x {xp} only applies to uncalibrated (-u) spectra, calibrated spectra combine BP and RP')
    get_sampling(sampling, uncalibrated, xp)
    outdir = os.path.dirname(os.path.abspath(fname))
    if not os.access(outdir, os.W_OK):
        raise PermissionError(f'Cannot write the outputs to {outdir}, copy the input csv somewhere writable')


def object_names(ids: pd.Series, names: Optional[pd.Series] = None) -> List[str]:
    """
    File names (without extension) for the per-object spectra: the object name made safe for a file name,
    or GaiaDR3_<source_id> without one (or names None). Names shared by different sources get _<source_id> added
    """
    files = []
    for i, xid in enumerate(ids):
        name = None if names is None else names.iloc[i]
        name = '' if name is None or pd.isna(name) else UNSAFE_PATTERN.sub('_', str(name).strip()).strip('._')
        files.append(name or f'GaiaDR3_{int(xid)}')
    counts = Counter(files)
    duplicates = sorted(name for name, count in counts.items() if count > 1)
    if duplicates:
        log.warning(f'{len(duplicates)} names are shared by several sources (e.g. {duplicates[0]}), '
                    f'their files are named <name>_<source_id>')
    return [f'{name}_{int(xid)}' if counts[name] > 1 else name for name, xid in zip(files, ids)]


def save(df: pd.DataFrame, fname: str, wave: np.ndarray, failed: Optional[Dict[Union[int, str], str]] = None,
         **kwargs):
    """
    Saves the input csv merged with the spectra (df) to the full and cut csvs, failed ({source ID: reason}) to the
    failed csv, and with kwargs outputstyle fits or txt, one file per object (kwargs idcolname, namecol,
    uncalibrated and xp as for batch)
    """
    outputstyle = kwargs.get('outputstyle', None)
    idcolname = kwargs.get('idcolname', 'source_id')
    namecol = kwargs.get('namecol', None)
    uncalibrated = kwargs.get('uncalibrated', False)
    fnameout, fnamecut, fnamefail = output_names(fname)
    fpath = os.path.join(os.path.dirname(os.path.abspath(fname)), 'outputspectra')

    dfcut = df.dropna(subset=['flux', ])
    with np.printoptions(threshold=sys.maxsize, linewidth=sys.maxsize):  # full arrays, not numpy's summary
        df.to_csv(fnameout, index=False)
        dfcut.to_csv(fnamecut, index=False)
    if failed:
        pd.DataFrame({'source_id': list(failed.keys()), 'reason': list(failed.values())}).to_csv(fnamefail,
                                                                                               index=False)
    elif os.path.exists(fnamefail):
        os.remove(fnamefail)

    if outputstyle is None:
        return
    fluxunit = 'electron/s' if uncalibrated else 'W m-2 nm-1'
    objects = dfcut.drop_duplicates(subset=idcolname)  # a source listed twice is saved once
    names = object_names(objects[idcolname], objects[namecol] if namecol and namecol != idcolname else None)
    os.makedirs(fpath, exist_ok=True)
    for (_, row), objname in tqdm(zip(objects.iterrows(), names), total=len(objects), desc='Saving spectra to files',
                                  file=sys.stderr):
        objpath = os.path.join(fpath, objname)
        arr = np.array((wave, row.flux, row.flux_error)).T
        if outputstyle == 'txt':
            np.savetxt(objpath + '.txt', arr, header=f'wave[nm] flux[{fluxunit}] fluxerror[{fluxunit}]')
        else:
            t = Table(arr, names=('wave', 'flux', 'fluxerror'),
                      units=('nm', None, None) if uncalibrated else ('nm', fluxunit, fluxunit))
            t.meta['SOURCEID'] = int(row[idcolname])
            t.meta['OBJECT'] = str(row[namecol]) if namecol and pd.notna(row[namecol]) else objname
            t.meta['XPMODE'] = f'convert {kwargs["xp"]}' if uncalibrated else 'calibrate'
            t.meta['GAIAXPY'] = gaiaxpy.__version__
            hdu = fits.table_to_hdu(t)
            if uncalibrated:  # electron is not a FITS standard unit, which astropy would silently drop
                hdu.header['TUNIT2'] = hdu.header['TUNIT3'] = fluxunit
            hdu.writeto(objpath + '.fits', overwrite=True)


def batch(fname: str, **kwargs) -> pd.DataFrame:
    """
    Downloads, calibrates and saves the XP spectra of the source IDs in a csv; kwargs are the long command line
    options (sampling as an array). Returns the csv (source ID column as integers) with flux and flux_error added;
    its attrs['incomplete'] is the number of sources that could not be downloaded (worth trying again later)
    """
    fname = os.path.abspath(fname)
    idcolname = kwargs.get('idcolname', 'source_id')
    validate_inputs(fname, idcolname, kwargs.get('namecol', None), kwargs.get('outputstyle', None),
                    kwargs.get('sampling', None), kwargs.get('uncalibrated', False), kwargs.get('xp', None))
    kwargs['xp'] = (kwargs.get('xp') or 'RP').upper()
    kwargs['username'], kwargs['password'] = get_credentials(kwargs.get('username', None),
                                                             kwargs.get('password', None),
                                                             kwargs.pop('anonymous', False))
    df = pd.read_csv(fname, dtype={idcolname: str})  # as strings, floats would lose digits of 19 digit IDs
    ids = read_ids(df[idcolname])
    invalid = {'' if pd.isna(value) else str(value): INVALID_ID for value in df[idcolname][ids.isna()]}
    if invalid:
        examples = ', '.join(repr(value) for value in list(invalid)[:3])
        hint = (' (written in scientific notation, e.g. by a spreadsheet, so digits have been lost)'
                if any('e+' in value.lower() for value in invalid) else '')
        nrows = int(ids.isna().sum())
        log.warning(f'{nrows} row{"s have" if nrows > 1 else " has"} no valid Gaia DR3 source ID{hint}, '
                    f'e.g. {examples}')
    df[idcolname] = ids
    if ids.isna().all():
        raise ValueError(f'No valid Gaia DR3 source IDs in column {idcolname!r} of {fname}')
    stackdf, wave, failed = download_xp(ids.dropna().tolist(), **kwargs)
    stackdf = stackdf.astype({'source_id': 'Int64'}).rename(columns={'source_id': idcolname})
    combdf = pd.merge(df, stackdf, on=idcolname, how='left')
    failed = {**failed, **invalid}
    reasons = list(failed.values())
    summary = ''.join(f', {reasons.count(reason)} {reason}' for reason in sorted(set(reasons)))
    log.info(f'{len(stackdf)} of {len(stackdf) + len(failed)} sources have XP spectra{summary}')
    partial = int(sum(np.isnan(np.asarray(flux, dtype=float)).any() for flux in stackdf.flux))
    if partial:  # e.g. calibrated spectra of sources with only an RP spectrum are NaN below about 640 nm
        log.warning(f'{partial} of the spectra are partial, NaN where a band is missing (e.g. below about 640 nm '
                    f'for sources with only an RP spectrum, as for many very red sources): check for NaNs before '
                    f'integrating or fitting them')
    save(combdf, fname, wave, failed, **kwargs)
    combdf.attrs['incomplete'] = sum(reason in INCOMPLETE for reason in reasons)
    if combdf.attrs['incomplete']:
        log.error(f'{combdf.attrs["incomplete"]} sources could not be downloaded, try them again later with '
                  f'-f {output_names(fname)[2]} (the reason column says which)')
    return combdf


class ConsoleFormatter(logging.Formatter):
    """Log messages as 'LEVEL: message', but plain 'message' for information"""

    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        return message if record.levelno == logging.INFO else f'{record.levelname}: {message}'


def main(argv: Optional[List[str]] = None):
    """
    Command line entry point. Exits with status 0 if everything was downloaded (or has no XP spectrum), 1 if some
    sources could not be downloaded, and 2 if nothing could be downloaded or the inputs are wrong
    """
    kwargs = vars(sysargs(argv))
    fname, verbose = kwargs.pop('filename'), kwargs.pop('verbose')
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(ConsoleFormatter())
    log.handlers = [handler]
    log.setLevel('DEBUG' if verbose else 'INFO')
    log.propagate = False
    astroquery_log.setLevel('INFO' if verbose else 'ERROR')
    astroquery_log.propagate = False  # astroquery has its own handler, don't print its messages twice
    try:
        combdf = batch(fname, **kwargs)
    except (ArchiveError, OSError, ValueError) as e:
        log.error(str(e) if isinstance(e, ArchiveError) and 'downloaded' in str(e) else f'Nothing downloaded: {e}')
        sys.exit(2)
    if combdf.attrs.get('incomplete'):
        sys.exit(1)


if __name__ == '__main__':
    main()
