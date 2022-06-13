# coding: utf-8
from contextlib import redirect_stdout
from io import StringIO
import logging
from functools import partial
import multiprocessing as mp
import os
import time
from typing import Optional, List

from astropy.table import Table
from gaiaxpy import calibrate
import pandas as pd
from tqdm import tqdm

from utils import *


def calibrate_wrap(idlist: List[int], **kwargs) -> Optional[pd.DataFrame]:
    """
    Wraps the GaiaXPy calibrate function with a list of source IDs

    Parameters
    ----------
    idlist: List[int]
        List of Gaia DR3 Source IDs
    kwargs
        Containing sampling information, truncation switch, verbosity switch, username and password information

    Returns
    -------
    _: pd.DataFrame
    """
    sampling = kwargs.get('sampling', None)
    truncate = kwargs.get('truncate', False)
    username = kwargs.get('username', None)
    password = kwargs.get('password', None)
    verbose = kwargs.get('verbose', False)
    lock.acquire()
    time.sleep(1)
    lock.release()
    try:
        dfout, _ = calibrate(idlist, sampling=sampling, truncation=truncate, username=username, password=password,
                             save_file=False)
    except (BaseException, TypeError, IndexError) as e:
        if verbose:
            print(repr(e))
        return None
    return dfout


def allcalib(xpids: List[int], nearest: int, **kwargs) -> Optional[pd.DataFrame]:
    """
    Depending on number of source IDs, will either thread in multiprocess pool or serialise

    Parameters
    ----------
    xpids: List[int]
        List of Gaia DR3 Source IDs
    nearest: int
        The nearest order of magnitude to try and chunk the list up by
    kwargs
        To be passed to the calibrate function

    Returns
    -------
    stackdf
        The combined dataframes from each thread having run through calibrate
    """
    cpucount = mp.cpu_count() - 1 or 1  # always leave a core free
    chunksize = int(nearest * np.ceil(len(xpids) // cpucount / nearest))  # size of each chunk to nearest order mag
    verbose = kwargs.get('verbose', False)
    if chunksize < cpucount:  # serialise if not many sources
        res = [calibrate_wrap([xid, ], **kwargs) for i, xid in enumerate(xpids)]
    else:
        xchunk = [xpids[i:i+chunksize] for i in range(0, len(xpids), chunksize)]  # chunk up the list of IDs
        if verbose:
            print(f'Chunksize {chunksize} for {len(xchunk)} threads')
        mp.freeze_support()
        with mp.Pool() as pool:  # wrap into a process pool
            res = pool.map(partial(calibrate_wrap, **kwargs), xchunk)
    res = [caldf for caldf in res if caldf is not None]  # make a list of the successful calibrations
    if len(res):
        stackdf = pd.concat(res, ignore_index=True)
        return stackdf
    return None


def download_xp(fname: str, **kwargs):
    """
    Iterative process for downloading XP spectra

    Parameters
    ----------
    fname: str
        The input csv name
    kwargs
        Contains the name of the source ID column plus kwargs to be passed to the calibrate function

    Returns
    -------
    combdf
        The full csv as given, with newly added flux and fluxerror columns
    """
    with tqdm(total=100, desc='Calibrating in steps') as pbar:
        idcolname = kwargs.get('idcolname', 'source_id')
        # downloaded from archive
        df = pd.read_csv(fname)
        xpidsstart = df[idcolname].values.tolist()
        pbar.update(20)

        # checking in steps of nearest 1000
        stackall = allcalib(xpidsstart, 1000, **kwargs)
        combdf = pd.merge(df, stackall, left_on=idcolname,
                          right_on='source_id', how='left') if stackall is not None else pd.DataFrame()
        if not len(combdf):
            raise ValueError('All batches failed')
        xpidsfail1000 = combdf[idcolname][combdf.flux.isna()].values.tolist()
        pbar.update(20)

        if len(xpidsfail1000) > 1000:
            # checking in steps of nearest 100
            stackdf100 = allcalib(xpidsfail1000, 100, **kwargs)
            stackall = pd.concat([stackall, stackdf100], ignore_index=True)
            combdf = pd.merge(df, stackall, left_on=idcolname,
                              right_on='source_id', how='left') if stackall is not None else combdf
            xpidsfail100 = combdf[idcolname][combdf.flux.isna()].values.tolist()
        else:
            xpidsfail100 = xpidsfail1000
        pbar.update(20)

        if len(xpidsfail100) > 100:
            # checking in steps of nearest 10
            stackdf10 = allcalib(xpidsfail100, 10, **kwargs)
            stackall = pd.concat([stackall, stackdf10], ignore_index=True)
            combdf = pd.merge(df, stackall, left_on=idcolname,
                              right_on='source_id', how='left') if stackall is not None else combdf
            xpidsfail10 = combdf[idcolname][combdf.flux.isna()].values.tolist()
        else:
            xpidsfail10 = xpidsfail100
        pbar.update(20)

        if len(xpidsfail10):
            # checking each object indivudally
            stackdf1 = allcalib(xpidsfail10, 1, **kwargs)
            stackall = pd.concat([stackall, stackdf1], ignore_index=True)
        combdf = pd.merge(df, stackall, left_on=idcolname,
                          right_on='source_id', how='left') if stackall is not None else combdf
        pbar.update(20)
    return combdf


def save(df: pd.DataFrame, fname: str, **kwargs):
    """
    Saves the found XP spectra, always to two csvs (full and cut - to those with spectra only)

    Parameters
    ----------
    df: pd.DataFrame
        The input csv
    fname: str
        The input csv name
    kwargs
        Contains the name of the source ID column plus name of object name column (e.g. shortname)
        and the outputstyle (fits or txt)
    """
    outputstyle = kwargs.get('outputstyle', None)
    idcolname = kwargs.get('idcolname', 'source_id')
    fnameout = fname.replace('.csv', '_XP_spectra.csv')
    fnamecut = fnameout.replace('.csv', '_cut.csv')
    fpath = os.path.dirname(os.path.abspath(fname)) + '/outputspectra/'
    namecol = kwargs.get('namecol', idcolname)

    # saving to files
    dfcut = df.dropna(subset=['flux', ])  # drop null rows to make cut table
    df.to_csv(fnameout, index=False)
    dfcut.to_csv(fnamecut, index=False)

    if outputstyle is None or outputstyle not in ('fits', 'txt'):
        return
    dfcutgrp = dfcut.groupby(namecol)
    dfcutgrplen = dfcut[namecol].nunique()
    wave: np.ndarray = kwargs.get('sampling', np.arange(336., 1201., 2.))  # default arange from help pages
    # https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Calibrator%20tutorial.html
    if not os.path.exists(fpath):
        os.mkdir(fpath)
    for i, (objname, objgrp) in tqdm(enumerate(dfcutgrp), total=dfcutgrplen, desc='Saving spectra to files'):
        if namecol == idcolname:
            objname = 'GaiaDR3_' + objname
        objpath = fpath + str(objname)
        row = objgrp.iloc[0]
        flux = row.flux
        fluxerr = row.flux_error
        arr = np.array((wave, flux, fluxerr)).T
        if outputstyle == 'txt':
            np.savetxt(objpath + '.txt', arr)
        elif outputstyle == 'fits':
            t = Table(arr, names=('wave', 'flux', 'fluxerror'))
            t.write(objpath + '.fits', overwrite=True)
    return


def batch(fname: str, **kwargs):
    """
    Batch function to either be run through main, or imported

    Parameters
    ----------
    fname: str
        The input csv name
    kwargs
        Keyword arguments for the calibrate function plus column names and output styling
    """
    if not fname.endswith('.csv'):
        raise ValueError('Input file must be a csv')
    combdf = download_xp(fname, **kwargs)
    save(combdf, fname, **kwargs)
    return


def main():
    """
    Main function, runs only when file run from command line using system arguments
    """
    args = sysargs()
    fname = args.filename
    kwargs = dict(sampling=args.sampling, username=args.username, password=args.password, truncate=args.truncate,
                  outputstyle=args.outputstyle, idcolname=args.idcolname, namecol=args.namecol, verbose=args.verbose)
    loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
    if args.verbose:
        [logger.setLevel('INFO') for logger in loggers]
        batch(fname, **kwargs)
    else:
        fio = StringIO()
        [logger.setLevel('ERROR') for logger in loggers]
        with redirect_stdout(fio):
            batch(fname, **kwargs)
    return


lock = mp.Lock()
if __name__ == '__main__':
    main()
