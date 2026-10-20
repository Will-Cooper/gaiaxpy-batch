# coding: utf-8
"""
Batch downloading of Gaia DR3 XP spectra through GaiaXPy.

GaiaXPy sends every source ID it is given to the Gaia archive in a single DataLink request.
This script splits a large list of source IDs into chunks, downloads the chunks in parallel threads,
retries failed chunks by splitting them in half, and saves the spectra to csv and (optionally) per-object
fits/txt files.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import nullcontext, redirect_stdout
from getpass import getpass
from io import StringIO
import logging
import os
import re
import sys
from typing import List, Optional, Tuple

from astropy.table import Table
from astroquery import log as astroquery_log
import gaiaxpy
from gaiaxpy import calibrate, convert, pwl_to_wl, wl_to_pwl
import numpy as np
import pandas as pd
from tqdm import tqdm

__version__ = '2.0.0'

log = logging.getLogger('gaiaxpy_batch')

SAMPLING_FUNCTIONS = ('linspace', 'arange', 'geomspace')
DEFAULT_PWL = np.linspace(0, 60, 600)  # default pseudo-wavelength sampling of gaiaxpy.convert
WL_LIMITS = (330., 1050.)  # absolute wavelength [nm] limits accepted by gaiaxpy.calibrate
PWL_LIMITS = (-10., 70.)  # pseudo-wavelength limits accepted by gaiaxpy.convert
NO_DATA_MESSAGE = 'No continuous BP/RP data found'  # raised by GaiaXPy when no source in a request has XP spectra
LIMIT_PATTERN = re.compile(r'Maximum number of Source Ids reached \(max: (\d+)')  # archive's request size error
FATAL_ERRORS =(AttributeError, ImportError, NameError, TypeError)  # bugs, not archive failures: no point retrying


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
    _args.add_argument('-c', '--chunksize', dest='chunksize', default=50, type=int,
                       help='Number of source IDs per archive request (default: 50, the anonymous archive limit)')
    _args.add_argument('-w', '--workers', dest='workers', default=4, type=int,
                       help='Number of parallel archive requests (default: 4)')
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


def calibrate_wrap(idlist: List[int], **kwargs) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Wraps the GaiaXPy calibrate or convert function with a list of source IDs

    Parameters
    ----------
    idlist: List[int]
        List of Gaia DR3 Source IDs
    kwargs
        Containing sampling information, truncation switch, mode switch, username and password information

    Returns
    -------
    dfout: pd.DataFrame
        The source_id, flux and flux_error of each spectrum
    sampling: np.ndarray
        The sampling GaiaXPy used
    """
    sampling = kwargs.get('sampling', None)
    truncate = kwargs.get('truncate', False)
    username = kwargs.get('username', None)
    password = kwargs.get('password', None)
    uncalibrated = kwargs.get('uncalibrated', False)
    xp = kwargs.get('xp', 'RP').upper()
    if not uncalibrated:
        dfout, sampling = calibrate(idlist, sampling=sampling, truncation=truncate, username=username,
                                    password=password, save_file=False)
    else:
        dfout, sampling = convert(idlist, sampling=sampling, truncation=truncate, username=username,
                                  password=password, save_file=False)
        dfout = dfout.loc[dfout.xp == xp].drop(columns='xp')
    return dfout[['source_id', 'flux', 'flux_error']], sampling


def _split(idlist: List[int], pbar: tqdm, limit: dict, size: int,
           **kwargs) -> Tuple[List[pd.DataFrame], List[np.ndarray]]:
    """
    Splits a chunk of source IDs into pieces of a given size and tries each piece

    Parameters
    ----------
    idlist: List[int]
        List of Gaia DR3 Source IDs
    pbar: tqdm
        The progress bar, updated with the number of source IDs finished
    limit: dict
        The archive's maximum number of source IDs per request (key 'max'), shared between threads
    size: int
        The size of each piece
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    res: List[pd.DataFrame]
        The successful spectra
    samplings: List[np.ndarray]
        The sampling of each success
    """
    res, samplings = [], []
    for i in range(0, len(idlist), size):
        pieceres, piecesamp = _bisect(idlist[i:i + size], pbar, limit, **kwargs)
        res.extend(pieceres)
        samplings.extend(piecesamp)
    return res, samplings


def _bisect(idlist: List[int], pbar: tqdm, limit: dict, **kwargs) -> Tuple[List[pd.DataFrame], List[np.ndarray]]:
    """
    Tries a chunk of source IDs, splitting it in half and retrying each half on failure.
    If the archive says the chunk is too large, the chunk is split to the archive's limit instead,
    and the limit is kept for the other chunks

    Parameters
    ----------
    idlist: List[int]
        List of Gaia DR3 Source IDs
    pbar: tqdm
        The progress bar, updated with the number of source IDs finished
    limit: dict
        The archive's maximum number of source IDs per request (key 'max'), shared between threads
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    res: List[pd.DataFrame]
        The successful spectra
    samplings: List[np.ndarray]
        The sampling of each success
    """
    maxids = limit.get('max')
    if maxids is not None and len(idlist) > maxids:
        return _split(idlist, pbar, limit, maxids, **kwargs)
    try:
        dfout, sampling = calibrate_wrap(idlist, **kwargs)
    except FATAL_ERRORS:
        raise
    except Exception as e:
        if NO_DATA_MESSAGE in str(e):  # none of these sources have XP spectra, no point splitting
            log.info(f'No XP spectra for {len(idlist)} sources starting {idlist[0]}')
            pbar.update(len(idlist))
            return [], []
        toomany = LIMIT_PATTERN.search(str(e))
        if toomany and int(toomany.group(1)) < len(idlist):
            maxids = int(toomany.group(1))
            if limit.get('max') != maxids:
                limit['max'] = maxids
                log.warning(f'The Gaia archive allows {maxids} source IDs per request, using chunks of {maxids}')
            return _split(idlist, pbar, limit, maxids, **kwargs)
        log.info(f'Chunk of {len(idlist)} sources starting {idlist[0]} failed: {e!r}')
        if len(idlist) == 1:
            pbar.update(1)
            return [], []
        return _split(idlist, pbar, limit, (len(idlist) + 1) // 2, **kwargs)
    pbar.update(len(idlist))
    return [dfout, ], [sampling, ]


def download_xp(xpids: List[int], chunksize: int = 50, workers: int = 4,
                **kwargs) -> Tuple[pd.DataFrame, np.ndarray, List[int]]:
    """
    Downloads XP spectra for a list of source IDs in chunks, in parallel threads

    Parameters
    ----------
    xpids: List[int]
        List of Gaia DR3 Source IDs
    chunksize: int
        Number of source IDs per archive request
    workers: int
        Number of parallel archive requests
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    stackdf: pd.DataFrame
        The source_id, flux and flux_error of each spectrum found
    wave: np.ndarray
        The wavelength [nm] of each sample, ascending
    failed: List[int]
        The source IDs with no spectrum returned
    """
    xpids = list(dict.fromkeys(int(xid) for xid in xpids))  # unique, keeping order
    if not xpids:
        raise ValueError('No source IDs given')
    wave, gaiaxpy_sampling = get_sampling(kwargs.get('sampling', None), kwargs.get('uncalibrated', False),
                                          kwargs.get('xp', 'RP'))
    kwargs['sampling'] = gaiaxpy_sampling
    chunksize = max(int(chunksize), 1)
    xchunk = [xpids[i:i + chunksize] for i in range(0, len(xpids), chunksize)]  # chunk up the list of IDs
    workers = max(min(int(workers), len(xchunk)), 1)
    log.info(f'{len(xpids)} sources in {len(xchunk)} chunks of up to {chunksize}, using {workers} workers')
    res, samplings = [], []
    limit = {}  # filled in when the archive reports its maximum number of source IDs per request
    with tqdm(total=len(xpids), desc='Downloading XP spectra', unit='source', file=sys.stderr) as pbar:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_bisect, chunk, pbar, limit, **kwargs) for chunk in xchunk]
            for future in as_completed(futures):
                chunkres, chunksamp = future.result()
                res.extend(chunkres)
                samplings.extend(chunksamp)
    res = [caldf for caldf in res if len(caldf)]  # make a list of the successful calibrations
    if not res:
        raise ValueError('All batches failed')
    stackdf = pd.concat(res, ignore_index=True).drop_duplicates(subset='source_id')
    if wave is None:  # calibrate with default sampling, only known now
        wave = np.asarray(samplings[0], dtype=float)
    order = np.argsort(wave)
    if np.any(order != np.arange(len(wave))):  # make the wavelengths ascending (BP in uncalibrated mode)
        wave = wave[order]
        for col in ('flux', 'flux_error'):
            stackdf[col] = [np.asarray(arr)[order] for arr in stackdf[col]]
    found = set(stackdf.source_id.astype(np.int64))
    failed = [xid for xid in xpids if xid not in found]
    return stackdf, wave, failed


def save(df: pd.DataFrame, fname: str, wave: np.ndarray, failed: Optional[List[int]] = None, **kwargs):
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
    failed: Optional[List[int]]
        The source IDs with no spectrum returned, saved to a third csv
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
        pd.DataFrame({'source_id': failed}).to_csv(fnamefail, index=False)
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
            t = Table(arr, names=('wave', 'flux', 'fluxerror'), units=('nm', fluxunit, fluxunit))
            t.meta['SOURCEID'] = int(row.source_id)
            t.meta['OBJECT'] = objname
            t.meta['XPMODE'] = f'convert {xp}' if uncalibrated else 'calibrate'
            t.meta['GAIAXPY'] = gaiaxpy.__version__
            t.write(objpath + '.fits', overwrite=True)
    return


def batch(fname: str, **kwargs) -> pd.DataFrame:
    """
    Batch function to either be run through main, or imported

    Parameters
    ----------
    fname: str
        The input csv name
    kwargs
        Keyword arguments for the calibrate function (sampling, truncate, uncalibrated, xp, username, password,
        anonymous), the download (chunksize, workers), column names (idcolname, namecol) and output styling
        (outputstyle)

    Returns
    -------
    combdf: pd.DataFrame
        The full csv as given, with newly added flux and flux_error columns
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
    stackdf, wave, failed = download_xp(xpids, kwargs.pop('chunksize', 50), kwargs.pop('workers', 4), **kwargs)
    stackdf = stackdf.astype({'source_id': 'Int64'})
    if idcolname == 'source_id':
        combdf = pd.merge(df, stackdf, on='source_id', how='left')
    else:
        combdf = pd.merge(df, stackdf, left_on=idcolname, right_on='source_id', how='left',
                          suffixes=('', '_xp'))
    log.warning(f'{len(stackdf)} of {len(stackdf) + len(failed)} sources have XP spectra')
    save(combdf, fname, wave, failed, **kwargs)
    return combdf


def main(argv: Optional[List[str]] = None):
    """
    Main function, runs only when file run from command line using system arguments
    """
    args = sysargs(argv)
    kwargs = dict(sampling=args.sampling, truncate=args.truncate,
                  outputstyle=args.outputstyle, idcolname=args.idcolname, namecol=args.namecol,
                  uncalibrated=args.uncalibrated, xp=args.xp, chunksize=args.chunksize, workers=args.workers,
                  anonymous=args.anonymous)
    kwargs['username'], kwargs['password'] = get_credentials(anonymous=args.anonymous)  # prompt before hiding stdout
    logging.basicConfig(format='%(levelname)s: %(message)s', level='INFO' if args.verbose else 'WARNING')
    astroquery_log.setLevel('INFO' if args.verbose else 'ERROR')
    # GaiaXPy prints its own messages and progress bars to stdout, hide them unless verbose
    with nullcontext() if args.verbose else redirect_stdout(StringIO()):
        batch(args.filename, **kwargs)
    return


if __name__ == '__main__':
    main()
