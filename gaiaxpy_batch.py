# coding: utf-8
from functools import partial
import multiprocessing as mp
import os
from re import search
import time
from typing import Optional

from astropy.table import Table
from gaiaxpy import calibrate
import pandas as pd
from tqdm import tqdm

from utils import *


def calibrate_wrap(idlist: list, **kwargs) -> Optional[pd.DataFrame]:
    i = kwargs.get('i', None)
    sampling = kwargs.get('sampling', None)
    truncate = kwargs.get('truncate', False)
    username = kwargs.get('username', None)
    password = kwargs.get('password', None)
    try:
        procnum = search('\d+', mp.current_process().name).group()
    except AttributeError:
        procnum = i
    time.sleep(0.01 * int(procnum))
    try:
        return calibrate(idlist, sampling=sampling, truncation=truncate, save_file=False,
                         username=username, password=password)[0]
    except (BaseException, TypeError, IndexError) as e:
        print(repr(e))
        return None


def allcalib(xpids: list, nearest: int, **kwargs) -> Optional[pd.DataFrame]:
    cpucount = mp.cpu_count() - 1 or 1
    chunksize = int(nearest * np.ceil(len(xpids) // cpucount / nearest))
    if chunksize < cpucount:
        res = [calibrate_wrap([xid, ], i=i, **kwargs) for i, xid in enumerate(xpids)]
    else:
        xchunk = [xpids[i:i+chunksize] for i in range(0, len(xpids), chunksize)]
        print(f'Chunksize {chunksize} for {len(xchunk)} threads')
        with mp.Pool() as pool:
            res = pool.map(partial(calibrate_wrap, **kwargs), xchunk)
    res = [caldf for caldf in res if caldf is not None]
    if len(res):
        stackdf = pd.concat(res, ignore_index=True)
        return stackdf
    return None


def download_xp(fname: str, idcolname: str, **kwargs):
    # downloaded from archive
    df = pd.read_csv(fname)
    xpidsstart = df[idcolname].values.tolist()
    # checking in steps of nearest 1000
    stackall = allcalib(xpidsstart, 1000, **kwargs)
    combdf = pd.merge(df, stackall, left_on=idcolname,
                      right_on='source_id', how='left') if stackall is not None else pd.DataFrame()
    if not len(combdf):
        raise ValueError('All batches failed')
    xpidsfail1000 = combdf.source_id[combdf.flux.isna()].values.tolist()
    if len(xpidsfail1000) > 1000:
        # checking in steps of nearest 100
        stackdf100 = allcalib(xpidsfail1000, 100, **kwargs)
        stackall = pd.concat([stackall, stackdf100], ignore_index=True)
        combdf = pd.merge(df, stackall, left_on=idcolname,
                          right_on='source_id', how='left') if stackall is not None else combdf
        xpidsfail100 = combdf.source_id[combdf.flux.isna()].values.tolist()
    else:
        xpidsfail100 = xpidsfail1000
    if len(xpidsfail100) > 100:
        # checking in steps of nearest 10
        stackdf10 = allcalib(xpidsfail100, 10, **kwargs)
        stackall = pd.concat([stackall, stackdf10], ignore_index=True)
        combdf = pd.merge(df, stackall, left_on=idcolname,
                          right_on='source_id', how='left') if stackall is not None else combdf
        xpidsfail10 = combdf.source_id[combdf.flux.isna()].values.tolist()
    else:
        xpidsfail10 = xpidsfail100
    if len(xpidsfail10):
        stackdf1 = allcalib(xpidsfail10, 1, **kwargs)
        stackall = pd.concat([stackall, stackdf1], ignore_index=True)
    combdf = pd.merge(df, stackall, left_on=idcolname,
                      right_on='source_id', how='left') if stackall is not None else combdf
    return combdf


def save(df: pd.DataFrame, fname: str, outputstyle: str, idcolname: str, **kwargs):
    fnameout = fname.replace('.csv', '_XP_spectra.csv')
    fnamecut = fnameout.replace('.csv', '_cut.csv')
    fpath = os.path.dirname(os.path.abspath(fname)) + '/outputspectra/'
    namecol = kwargs.get('namecol', idcolname)
    dfcut = df.dropna(subset=['flux', ])
    df.to_csv(fnameout, index=False)
    dfcut.to_csv(fnamecut, index=False)
    if outputstyle is None or outputstyle not in ('fits', 'txt'):
        return
    dfcutgrp = dfcut.groupby(namecol)
    dfcutgrplen = dfcut[namecol].nunique()
    wave: np.ndarray = kwargs.get('sampling', np.arange(336., 1201., 2.))
    if not os.path.exists(fpath):
        os.mkdir(fpath)
    for i, (objname, objgrp) in tqdm(enumerate(dfcutgrp), total=dfcutgrplen, desc='Saving spectra as txt'):
        if namecol == idcolname:
            objname = 'GaiaDR3_' + objname
        objpath = fpath + str(objname)
        row = objgrp.iloc[0]
        flux = str2arr(row.flux)
        fluxerr = str2arr(row.error)
        arr = np.array((wave, flux, fluxerr)).T
        if outputstyle == 'txt':
            np.savetxt(objpath + '.txt', arr)
        elif outputstyle == 'fits':
            t = Table(arr, names=('wave', 'flux', 'fluxerror'))
            t.write(objpath + 'fits', overwrite=True)
    return


def batch(fname: str, **kwargs):
    outputstyle = kwargs.get('outputstyle', None)
    idcolname = kwargs.get('idcolname', 'source_id')
    combdf = download_xp(fname, idcolname, **kwargs)
    save(combdf, fname, outputstyle, idcolname, **kwargs)
    return


def main():
    args = sysargs()
    fname = args.filename
    kwargs = dict(sampling=args.sampling, username=args.username, password=args.password, truncate=args.truncate,
                  outputstyle=args.outputstyle, idcolname=args.idcolname, namecol=args.namecol)
    batch(fname, **kwargs)
    return


if __name__ == '__main__':
    main()
