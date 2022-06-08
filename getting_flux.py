# coding: utf-8
from functools import partial
import multiprocessing as mp
import os
from re import search
import time
from typing import Optional

import numpy as np
import pandas as pd
from gaiaxpy import calibrate


def calibrate_wrap(idlist: list, **kwargs) -> Optional[pd.DataFrame]:
    i = kwargs.get('i', None)
    sampling = kwargs.get('sampling', np.linspace(600, 1050, 120))
    try:
        procnum = search('\d+', mp.current_process().name).group()
    except AttributeError:
        procnum = i
    time.sleep(0.01 * int(procnum))
    try:
        return calibrate(idlist, sampling=sampling, save_file=False)[0]
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


def download_rp(fname: str, fnameout: str, **kwargs):
    # downloaded from archive
    df = pd.read_csv(fname)
    xpidsstart = df.source_id.values.tolist()
    # checking in steps of nearest 1000
    stackall = allcalib(xpidsstart, 1000, **kwargs)
    combdf = pd.merge(df, stackall, on='source_id', how='left') if stackall is not None else pd.DataFrame()
    if not len(combdf):
        raise ValueError('All batches failed')
    xpidsfail1000 = combdf.source_id[combdf.flux.isna()].values.tolist()
    if len(xpidsfail1000) > 1000:
        # checking in steps of nearest 100
        stackdf100 = allcalib(xpidsfail1000, 100, **kwargs)
        stackall = pd.concat([stackall, stackdf100], ignore_index=True)
        combdf = pd.merge(df, stackall, on='source_id', how='left') if stackall is not None else combdf
        xpidsfail100 = combdf.source_id[combdf.flux.isna()].values.tolist()
    else:
        xpidsfail100 = xpidsfail1000
    if len(xpidsfail100) > 100:
        # checking in steps of nearest 10
        stackdf10 = allcalib(xpidsfail100, 10, **kwargs)
        stackall = pd.concat([stackall, stackdf10], ignore_index=True)
        combdf = pd.merge(df, stackall, on='source_id', how='left') if stackall is not None else combdf
        xpidsfail10 = combdf.source_id[combdf.flux.isna()].values.tolist()
    else:
        xpidsfail10 = xpidsfail100
    if len(xpidsfail10):
        stackdf1 = allcalib(xpidsfail10, 1, **kwargs)
        stackall = pd.concat([stackall, stackdf1], ignore_index=True)
    combdf = pd.merge(df, stackall, on='source_id', how='left') if stackall is not None else combdf
    combdf.dropna(subset=['flux', ]).to_csv(fnameout, index=False)
    return combdf


def applytodata(fname: str, combdf: pd.DataFrame):
    # going over full table
    dfull = pd.read_csv(fname)
    dfullm = pd.merge(dfull, combdf[['source_id', 'flux', 'error']], how='left', on='source_id')
    dfullm.to_csv(fname.replace('.csv', '_RP.csv'), index=False)
    return


def main(fname: str, redownload: bool = False, **kwargs):
    rpfname = kwargs.get('rpfname', 'full_xp.csv')
    rpfnameout = rpfname.replace('.csv', '_RP_spectra.csv')
    if redownload or not os.path.exists(rpfnameout):
        combdf = download_rp(rpfname, rpfnameout, **kwargs)
    else:
        combdf = pd.read_csv(rpfnameout)
    applytodata(fname, combdf)
    return


if __name__ == '__main__':
    # main('gaia_ucds.csv')
    main('gucds.csv', redownload=True, rpfname='gucds_xp.csv')
