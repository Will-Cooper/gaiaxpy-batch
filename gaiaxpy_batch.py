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
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import nullcontext, redirect_stdout
from getpass import getpass
from io import StringIO
import logging
import os
import re
import socket
import sys
import threading
from typing import Callable, Dict, List, Optional, Tuple

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

__version__ = '2.0.0'

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

LIMIT_PATTERN = re.compile(r'Maximum number of Source Ids reached \(max: (\d+)')  # archive's request size error
STATUS_PATTERN = re.compile(r'\b(?:HTTP Error|Error|status code)\s*:?\s*(\d{3})\b', re.IGNORECASE)
FATAL_ERRORS = (AttributeError, ImportError, NameError, TypeError)  # bugs, not archive failures: no point retrying

# reasons given in the failed csv
NO_XP = 'no XP spectrum'
DOWNLOAD_FAILED = 'download failed'
NOT_ATTEMPTED = 'not attempted, archive unavailable'
CALIBRATION_FAILED = 'calibration failed'
INCOMPLETE = (DOWNLOAD_FAILED, NOT_ATTEMPTED)  # reasons worth trying again later


def parse_sampling(func: str, start: float, end: float, step: float) -> np.ndarray:
    """
    Creates a wavelength sampling from a whitelisted numpy function, e.g. ('linspace', 600, 1050, 120)

    Parameters
    ----------
    func: str
        One of linspace, arange or geomspace
    start: float
        The first wavelength [nm]
    end: float
        The last wavelength [nm]
    step: float
        The number of points (linspace, geomspace) or the step size [nm] (arange)

    Returns
    -------
    _: np.ndarray
        The wavelength sampling [nm]
    """
    if func not in SAMPLING_FUNCTIONS:
        raise ValueError(f'Sampling function must be one of {SAMPLING_FUNCTIONS}, not {func!r}')
    start, end, step = float(start), float(end), float(step)
    if func == 'arange':
        return np.arange(start, end, step)
    return getattr(np, func)(start, end, int(step))


def sysargs(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """
    These are the system arguments given after calling this python script

    Parameters
    ----------
    argv: Optional[List[str]]
        Arguments to parse, defaults to sys.argv

    Returns
    -------
    _args
        The different argument parameters, can be grabbed via their long names (e.g. _args.filename)
    """
    _args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    _args.add_argument('-f', '--filename', dest='filename', required=True, type=str, metavar='Filename',
                       help='Input CSV file containing the Gaia DR3 Source IDs')
    _args.add_argument('-s', '--sampling', dest='sampling', default=None, type=str, nargs=4,
                       metavar=('numpy function', 'start', 'end', 'step'),
                       help='Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120" '
                            f'(functions: {", ".join(SAMPLING_FUNCTIONS)})')
    _args.add_argument('-t', '--truncate', dest='truncate', help='Truncate set of bases?', default=False,
                       action='store_true')
    _args.add_argument('-o', '--outputstyle', dest='outputstyle', choices=('fits', 'txt'), default=None,
                       help='Output of produced spectra (default: None)')
    _args.add_argument('-i', '--idcolname', dest='idcolname', help='Name of the column relating to DR3 Source ID',
                       default='source_id', metavar='Source ID column name')
    _args.add_argument('-n', '--namecol', dest='namecol', default=None, metavar='Object name column name',
                       help='Name of the column relating to the object name (for saving spectra, '
                            'default: GaiaDR3_<source_id>)')
    _args.add_argument('-u', '--uncalibrated', dest='uncalibrated', action='store_true', default=False,
                       help='Switch to uncalibrated spectra mode?')
    _args.add_argument('-x', '--xp', dest='xp', default='RP', type=str.upper, choices=('BP', 'RP'),
                       help='If in uncalibrated mode, RP or BP? (default: RP)')
    _args.add_argument('-c', '--chunksize', dest='chunksize', default=None, type=int,
                       help=f'Number of source IDs per archive request (default: {ANONYMOUS_CHUNKSIZE} anonymous, '
                            f'{LOGGED_IN_CHUNKSIZE} logged in; the archive allows up to 50 and 5000)')
    _args.add_argument('-w', '--workers', dest='workers', default=4, type=int,
                       help=f'Number of parallel archive requests (default: 4, at most {LOGGED_IN_MAX_WORKERS} '
                            f'when logged in)')
    _args.add_argument('--retries', dest='retries', default=RETRIES, type=int,
                       help=f'Retries of a request after a server-side failure (default: {RETRIES})')
    _args.add_argument('--retry-delay', dest='retry_delay', default=RETRY_DELAY, type=float,
                       help=f'Seconds to wait before the first retry, doubled for each one after '
                            f'(default: {RETRY_DELAY:.0f})')
    _args.add_argument('--timeout', dest='timeout', default=TIMEOUT, type=float,
                       help=f'Seconds without a response before a connection is given up (default: {TIMEOUT:.0f})')
    _args.add_argument('-a', '--anonymous', dest='anonymous', default=False, action='store_true',
                       help='Query the archive without logging in (otherwise uses $GAIA_USERNAME and '
                            '$GAIA_PASSWORD, or prompts)')
    _args.add_argument('-v', '--verbose', dest='verbose', default=False, action='store_true',
                       help='Print failure errors?')
    _args.add_argument('--version', action='version', version=f'%(prog)s {__version__} (GaiaXPy {gaiaxpy.__version__})')
    _args = _args.parse_args(argv)
    if _args.sampling is not None:
        _args.sampling = parse_sampling(*_args.sampling)
    return _args


def str2float(s: str) -> np.ndarray:
    """
    Convenience function for converting the flux/ error cells of the output csv back into arrays

    Parameters
    ----------
    s: str
        What should be an array of floats, being interpreted as a string

    Returns
    -------
    _: np.ndarray
        Array of floats as parsed from string
    """
    return np.fromiter(filter(None, s[1:-1].replace('\n', ' ').split(' ')), float)


def get_credentials(username: Optional[str] = None, password: Optional[str] = None,
                    anonymous: bool = False) -> Tuple[Optional[str], Optional[str]]:
    """
    Gets the Gaia archive credentials from the arguments, then the environment, then a prompt

    Parameters
    ----------
    username: Optional[str]
        Cosmos username
    password: Optional[str]
        Cosmos password
    anonymous: bool
        Switch to skip logging in

    Returns
    -------
    _: Tuple[Optional[str], Optional[str]]
        The username and password
    """
    if anonymous:
        return None, None
    username = username or os.environ.get('GAIA_USERNAME') or input('Gaia archive username: ')
    password = password or os.environ.get('GAIA_PASSWORD') or getpass('Gaia archive password: ')
    return username, password


def get_sampling(sampling: Optional[np.ndarray] = None, uncalibrated: bool = False,
                 xp: str = 'RP') -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Validates the wavelength sampling and, in uncalibrated mode, converts it to pseudo-wavelength

    Parameters
    ----------
    sampling: Optional[np.ndarray]
        Wavelength [absolute nm] sampling, None for the GaiaXPy default
    uncalibrated: bool
        Switch to uncalibrated (gaiaxpy.convert) mode
    xp: str
        BP or RP, for uncalibrated mode

    Returns
    -------
    wave: Optional[np.ndarray]
        The wavelength [nm] of each sample, in the order GaiaXPy returns them (None if only known after calibrate)
    gaiaxpy_sampling: Optional[np.ndarray]
        The sampling to pass to GaiaXPy (absolute nm for calibrate, ascending pseudo-wavelength for convert)
    """
    if not uncalibrated:
        if sampling is None:
            return None, None
        wave = np.asarray(sampling, dtype=float)
        if np.any(np.diff(wave) <= 0) or wave[0] < WL_LIMITS[0] or wave[-1] > WL_LIMITS[1]:
            raise ValueError(f'Sampling must be ascending and within {WL_LIMITS[0]:.0f}-{WL_LIMITS[1]:.0f} nm')
        return wave, wave
    xp = xp.upper()
    if xp not in ('BP', 'RP'):
        raise ValueError('Needs to be BP or RP for dispersion')
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
    """
    Wraps the GaiaXPy calibrate or convert function with downloaded XP continuous spectra

    Parameters
    ----------
    coeffs: pd.DataFrame
        The XP continuous spectra (basis function coefficients) as downloaded from the archive
    kwargs
        Containing sampling information, truncation switch and mode switch

    Returns
    -------
    dfout: pd.DataFrame
        The source_id, flux and flux_error of each spectrum
    sampling: np.ndarray
        The sampling GaiaXPy used
    """
    sampling = kwargs.get('sampling', None)
    truncate = kwargs.get('truncate', False)
    uncalibrated = kwargs.get('uncalibrated', False)
    xp = kwargs.get('xp', 'RP').upper()
    if not uncalibrated:
        dfout, sampling = calibrate(coeffs, sampling=sampling, truncation=truncate, save_file=False)
    else:
        dfout, sampling = convert(coeffs, sampling=sampling, truncation=truncate, save_file=False)
        dfout = dfout.loc[dfout.xp == xp].drop(columns='xp')
    return dfout[['source_id', 'flux', 'flux_error']], sampling


class ArchiveError(Exception):
    """An archive request failed on every retry"""


class ArchiveUnavailable(ArchiveError):
    """The archive appears to be down, or the login failed: stop downloading"""


def http_status(e: Exception) -> Optional[int]:
    """
    Gets the HTTP status code of an exception, if it has one

    Parameters
    ----------
    e: Exception
        The exception raised by an archive request

    Returns
    -------
    _: Optional[int]
        The HTTP status code
    """
    status = getattr(getattr(e, 'response', None), 'status_code', None)
    if status is None:
        match = STATUS_PATTERN.search(str(e))
        status = int(match.group(1)) if match else None
    return status


def is_transient(e: Exception) -> bool:
    """
    Whether an archive failure is probably temporary and worth retrying: server errors (HTTP 5xx), timeouts,
    lost sessions and broken responses (e.g. 'File is not a zip file'), rather than a bad request (other HTTP 4xx)

    Parameters
    ----------
    e: Exception
        The exception raised by an archive request

    Returns
    -------
    _: bool
        True if the request should be retried
    """
    status = http_status(e)
    return status is None or status >= 500 or status in (401, 403, 408, 429)


def _short(e: Exception) -> str:
    """The first line of an exception, for log messages"""
    lines = [line for line in str(e).splitlines() if line.strip()]
    return f'{type(e).__name__}: {" ".join(lines[:2])}'[:200]


class Archive:
    """
    Downloads XP continuous spectra from the Gaia archive with astroquery, using one session per thread
    (each logged in once), and retrying server-side failures with increasing waits
    """

    def __init__(self, username: Optional[str] = None, password: Optional[str] = None, retries: int = RETRIES,
                 delay: float = RETRY_DELAY):
        """
        Parameters
        ----------
        username: Optional[str]
            Cosmos username, None for anonymous
        password: Optional[str]
            Cosmos password
        retries: int
            Retries of a request after a server-side failure
        delay: float
            Seconds to wait before the first retry, doubled for each one after
        """
        self.username = username
        self.password = password
        self.retries = max(int(retries), 0)
        self.delay = max(float(delay), 0.)
        self.stop = threading.Event()  # set when the archive is taken to be unavailable
        self._local = threading.local()

    def session(self) -> GaiaClass:
        """
        The astroquery session of this thread, created (and logged in) on first use

        Returns
        -------
        gaia: GaiaClass
            The session
        """
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
        """
        Logs a session in, raising if it did not work (astroquery logs login errors rather than raising them)

        Parameters
        ----------
        gaia: GaiaClass
            The session
        """
        gaia.login(user=self.username, password=self.password)
        tap = getattr(gaia, '_TapPlus__isLoggedIn', None)  # astroquery internals, None if they change
        data = getattr(getattr(gaia, '_GaiaClass__gaiadata', None), '_TapPlus__isLoggedIn', None)
        if tap is False or data is False:
            raise ConnectionError('Gaia archive login failed')

    def retry(self, func: Callable, what: str):
        """
        Calls a function, retrying it with increasing waits while it fails in a way that is probably temporary

        Parameters
        ----------
        func: Callable
            The request
        what: str
            What the request is, for log messages

        Returns
        -------
        _
            Whatever the function returns
        """
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
                if http_status(e) in (401, 403):
                    self._local.gaia = None  # the session has expired, log in again
                wait = self.delay * 2 ** attempt
                log.warning(f'{what} failed ({_short(e)}), retrying in {wait:.0f} s ({attempt + 1}/{self.retries})')
                self.stop.wait(wait)

    def load(self, idlist: List[int]) -> Optional[pd.DataFrame]:
        """
        Downloads the XP continuous spectra of a list of source IDs, as GaiaXPy does itself

        Parameters
        ----------
        idlist: List[int]
            List of Gaia DR3 Source IDs

        Returns
        -------
        _: Optional[pd.DataFrame]
            The XP continuous spectra, None if none of the sources have any
        """
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
    """
    Calibrates downloaded spectra, splitting them in half and retrying each half on failure (locally, no download)

    Parameters
    ----------
    coeffs: pd.DataFrame
        The XP continuous spectra as downloaded from the archive
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    _: Results
        The spectra, their samplings and the source IDs that failed, with the reason
    """
    try:
        dfout, sampling = calibrate_wrap(coeffs, **kwargs)
    except FATAL_ERRORS:
        raise
    except Exception as e:
        if len(coeffs) == 1:
            xid = int(coeffs.source_id.iloc[0])
            log.info(f'Calibration of {xid} failed: {_short(e)}')
            return [], [], {xid: CALIBRATION_FAILED}
        half = (len(coeffs) + 1) // 2
        return _merge([_calibrate_rows(coeffs.iloc[:half], **kwargs), _calibrate_rows(coeffs.iloc[half:], **kwargs)])
    return [dfout, ], [sampling, ], {}


def _process(idlist: List[int], archive: Archive, state: dict, pbar: tqdm, **kwargs) -> Results:
    """
    Downloads and calibrates a chunk of source IDs.
    If the archive says the chunk is too large, the chunk is split to the archive's limit and the limit is kept for
    the other chunks. A chunk failing with a bad request is split in half and each half tried again.
    Chunks that fail every retry are given up, and after several in a row the archive is taken to be unavailable

    Parameters
    ----------
    idlist: List[int]
        List of Gaia DR3 Source IDs
    archive: Archive
        The archive connections
    state: dict
        Shared between threads: the archive's limit ('limit'), chunks in a row failing every retry ('exhausted'),
        why the downloads were stopped ('abort') and a lock ('lock')
    pbar: tqdm
        The progress bar, updated with the number of source IDs finished
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    _: Results
        The spectra, their samplings and the source IDs that failed, with the reason
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
        log.info(f'Download of {len(idlist)} sources starting {idlist[0]} failed: {_short(e)}')
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
                username: Optional[str] = None, password: Optional[str] = None, retries: Optional[int] = None,
                retry_delay: Optional[float] = None, timeout: Optional[float] = None,
                **kwargs) -> Tuple[pd.DataFrame, np.ndarray, Dict[int, str]]:
    """
    Downloads XP spectra for a list of source IDs in chunks, in parallel threads, and calibrates them with GaiaXPy

    Parameters
    ----------
    xpids: List[int]
        List of Gaia DR3 Source IDs
    chunksize: Optional[int]
        Number of source IDs per archive request (default: 50 anonymous, 1000 logged in)
    workers: int
        Number of parallel archive requests (at most 5 when logged in)
    username: Optional[str]
        Cosmos username, None for anonymous
    password: Optional[str]
        Cosmos password
    retries: Optional[int]
        Retries of a request after a server-side failure (default: 3)
    retry_delay: Optional[float]
        Seconds to wait before the first retry, doubled for each one after (default: 30)
    timeout: Optional[float]
        Seconds without a response before a connection is given up (default: 900)
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    stackdf: pd.DataFrame
        The source_id, flux and flux_error of each spectrum found
    wave: np.ndarray
        The wavelength [nm] of each sample, ascending
    failed: Dict[int, str]
        The source IDs with no spectrum, with the reason
    """
    xpids = list(dict.fromkeys(int(xid) for xid in xpids))  # unique, keeping order
    if not xpids:
        raise ValueError('No source IDs given')
    wave, gaiaxpy_sampling = get_sampling(kwargs.get('sampling', None), kwargs.get('uncalibrated', False),
                                          kwargs.get('xp', 'RP'))
    kwargs['sampling'] = gaiaxpy_sampling
    if chunksize is None:
        chunksize = LOGGED_IN_CHUNKSIZE if username else ANONYMOUS_CHUNKSIZE
    chunksize = max(int(chunksize), 1)
    if username and workers > LOGGED_IN_MAX_WORKERS:
        log.warning(f'The Gaia archive allows {LOGGED_IN_MAX_WORKERS} active jobs per user, '
                    f'using {LOGGED_IN_MAX_WORKERS} workers')
        workers = LOGGED_IN_MAX_WORKERS
    xchunk = [xpids[i:i + chunksize] for i in range(0, len(xpids), chunksize)]  # chunk up the list of IDs
    workers = max(min(int(workers), len(xchunk)), 1)
    log.info(f'{len(xpids)} sources in {len(xchunk)} chunks of up to {chunksize}, using {workers} workers, '
             f'{"logged in as " + username if username else "anonymously"}')
    archive = Archive(username, password, RETRIES if retries is None else retries,
                      RETRY_DELAY if retry_delay is None else retry_delay)
    state = dict(limit=None, exhausted=0, abort=None, lock=threading.Lock())
    res, samplings, failed = [], [], {}
    oldtimeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(TIMEOUT if timeout is None else timeout)  # astroquery connections have no timeout
    try:
        with tqdm(total=len(xpids), desc='Downloading XP spectra', unit='source', file=sys.stderr) as pbar:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_process, chunk, archive, state, pbar, **kwargs) for chunk in xchunk]
                for future in as_completed(futures):
                    chunkres, chunksamp, chunkfailed = future.result()
                    res.extend(chunkres)
                    samplings.extend(chunksamp)
                    failed.update(chunkfailed)
    finally:
        socket.setdefaulttimeout(oldtimeout)
    if state['abort']:
        log.error(state['abort'])
    res = [caldf for caldf in res if len(caldf)]  # make a list of the successful calibrations
    if not res:
        if state['abort']:
            raise ArchiveUnavailable(state['abort'])
        reasons = list(failed.values())
        summary = ', '.join(f'{reasons.count(reason)} {reason}' for reason in sorted(set(reasons)))
        if any(reason in INCOMPLETE for reason in reasons):
            raise ArchiveError(f'No XP spectra downloaded, {summary}')
        raise ValueError(f'No XP spectra found, {summary}')
    stackdf = pd.concat(res, ignore_index=True).drop_duplicates(subset='source_id')
    if wave is None:  # calibrate with default sampling, only known now
        wave = np.asarray(samplings[0], dtype=float)
    order = np.argsort(wave)
    if np.any(order != np.arange(len(wave))):  # make the wavelengths ascending (BP in uncalibrated mode)
        wave = wave[order]
        for col in ('flux', 'flux_error'):
            stackdf[col] = [np.asarray(arr)[order] for arr in stackdf[col]]
    found = set(stackdf.source_id.astype(np.int64))
    failed = {xid: failed.get(xid, DOWNLOAD_FAILED) for xid in xpids if xid not in found}  # input order
    return stackdf, wave, failed


def save(df: pd.DataFrame, fname: str, wave: np.ndarray, failed: Optional[Dict[int, str]] = None, **kwargs):
    """
    Saves the found XP spectra, always to two csvs (full and cut - to those with spectra only)

    Parameters
    ----------
    df: pd.DataFrame
        The input csv, merged with the spectra
    fname: str
        The input csv name
    wave: np.ndarray
        The wavelength [nm] of each sample
    failed: Optional[Dict[int, str]]
        The source IDs with no spectrum, with the reason, saved to a third csv
    kwargs
        Contains the name of the source ID column plus name of object name column (e.g. shortname)
        and the outputstyle (fits or txt)
    """
    outputstyle = kwargs.get('outputstyle', None)
    idcolname = kwargs.get('idcolname', 'source_id')
    namecol = kwargs.get('namecol', None) or idcolname
    uncalibrated = kwargs.get('uncalibrated', False)
    xp = kwargs.get('xp', 'RP').upper()
    fnameout = fname.replace('.csv', '_XP_spectra.csv')
    fnamecut = fnameout.replace('.csv', '_cut.csv')
    fnamefail = fnameout.replace('.csv', '_failed.csv')
    fpath = os.path.join(os.path.dirname(os.path.abspath(fname)), 'outputspectra')

    # saving to files, writing the full arrays rather than numpy's summarised repr
    dfcut = df.dropna(subset=['flux', ])  # drop null rows to make cut table
    with np.printoptions(threshold=sys.maxsize, linewidth=sys.maxsize):
        df.to_csv(fnameout, index=False)
        dfcut.to_csv(fnamecut, index=False)
    if failed:
        pd.DataFrame({'source_id': list(failed.keys()), 'reason': list(failed.values())}).to_csv(fnamefail,
                                                                                               index=False)
    elif os.path.exists(fnamefail):
        os.remove(fnamefail)

    if outputstyle is None or outputstyle not in ('fits', 'txt'):
        return
    fluxunit = 'electron/s' if uncalibrated else 'W m-2 nm-1'
    dfcutgrp = dfcut.groupby(namecol)
    dfcutgrplen = dfcut[namecol].nunique()
    os.makedirs(fpath, exist_ok=True)
    for objname, objgrp in tqdm(dfcutgrp, total=dfcutgrplen, desc='Saving spectra to files', file=sys.stderr):
        row = objgrp.iloc[0]
        objname = f'GaiaDR3_{int(row.source_id)}' if namecol == idcolname else str(objname)
        objpath = os.path.join(fpath, objname)
        arr = np.array((wave, row.flux, row.flux_error)).T
        if outputstyle == 'txt':
            np.savetxt(objpath + '.txt', arr, header=f'wave[nm] flux[{fluxunit}] fluxerror[{fluxunit}]')
        elif outputstyle == 'fits':
            t = Table(arr, names=('wave', 'flux', 'fluxerror'),
                      units=('nm', None, None) if uncalibrated else ('nm', fluxunit, fluxunit))
            t.meta['SOURCEID'] = int(row.source_id)
            t.meta['OBJECT'] = objname
            t.meta['XPMODE'] = f'convert {xp}' if uncalibrated else 'calibrate'
            t.meta['GAIAXPY'] = gaiaxpy.__version__
            hdu = fits.table_to_hdu(t)
            if uncalibrated:  # electron is not a FITS standard unit, which astropy would silently drop
                hdu.header['TUNIT2'] = hdu.header['TUNIT3'] = fluxunit
            hdu.writeto(objpath + '.fits', overwrite=True)
    return


def batch(fname: str, **kwargs) -> pd.DataFrame:
    """
    Batch function to either be run through main, or imported

    Parameters
    ----------
    fname: str
        The input csv name
    kwargs
        Keyword arguments for the calibrate function (sampling, truncate, uncalibrated, xp), the download
        (username, password, anonymous, chunksize, workers, retries, retry_delay, timeout),
        column names (idcolname, namecol) and output styling (outputstyle)

    Returns
    -------
    combdf: pd.DataFrame
        The full csv as given, with newly added flux and flux_error columns.
        combdf.attrs['incomplete'] is the number of sources that could not be downloaded (worth trying again later)
    """
    if not fname.endswith('.csv'):
        raise ValueError('Input file must be a csv')
    idcolname = kwargs.get('idcolname', 'source_id')
    kwargs['username'], kwargs['password'] = get_credentials(kwargs.get('username', None),
                                                             kwargs.get('password', None),
                                                             kwargs.get('anonymous', False))
    df = pd.read_csv(fname)
    df[idcolname] = pd.to_numeric(df[idcolname], errors='coerce').astype('Int64')
    xpids = df[idcolname].dropna().tolist()
    stackdf, wave, failed = download_xp(xpids, **kwargs)
    stackdf = stackdf.astype({'source_id': 'Int64'})
    if idcolname == 'source_id':
        combdf = pd.merge(df, stackdf, on='source_id', how='left')
    else:
        combdf = pd.merge(df, stackdf, left_on=idcolname, right_on='source_id', how='left',
                          suffixes=('', '_xp'))
    reasons = list(failed.values())
    summary = ''.join(f', {reasons.count(reason)} {reason}' for reason in sorted(set(reasons)))
    log.warning(f'{len(stackdf)} of {len(stackdf) + len(failed)} sources have XP spectra{summary}')
    save(combdf, fname, wave, failed, **kwargs)
    combdf.attrs['incomplete'] = sum(reason in INCOMPLETE for reason in reasons)
    if combdf.attrs['incomplete']:
        fnamefail = fname.replace('.csv', '_XP_spectra_failed.csv')
        log.error(f'{combdf.attrs["incomplete"]} sources could not be downloaded, try them again later with '
                  f'-f {fnamefail} (the reason column says which)')
    return combdf


def main(argv: Optional[List[str]] = None):
    """
    Main function, runs only when file run from command line using system arguments
    """
    args = sysargs(argv)
    kwargs = dict(sampling=args.sampling, truncate=args.truncate,
                  outputstyle=args.outputstyle, idcolname=args.idcolname, namecol=args.namecol,
                  uncalibrated=args.uncalibrated, xp=args.xp, chunksize=args.chunksize, workers=args.workers,
                  anonymous=args.anonymous, retries=args.retries, retry_delay=args.retry_delay, timeout=args.timeout)
    kwargs['username'], kwargs['password'] = get_credentials(anonymous=args.anonymous)  # prompt before hiding stdout
    logging.basicConfig(format='%(levelname)s: %(message)s', level='INFO' if args.verbose else 'WARNING')
    astroquery_log.setLevel('INFO' if args.verbose else 'ERROR')
    # GaiaXPy prints its own messages and progress bars to stdout, hide them unless verbose
    with nullcontext() if args.verbose else redirect_stdout(StringIO()):
        try:
            combdf = batch(args.filename, **kwargs)
        except ArchiveError as e:
            log.error(f'Nothing downloaded: {e}')
            sys.exit(2)
    if combdf.attrs.get('incomplete'):
        sys.exit(1)
    return


if __name__ == '__main__':
    main()
