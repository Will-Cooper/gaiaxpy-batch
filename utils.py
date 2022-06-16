import argparse

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def sysargs():
    """
    These are the system arguments given after calling this python script
    Returns
    -------
    _args
        The different argument parameters, can be grabbed via their long names (e.g. _args.filename)
    """
    _args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    _args.add_argument('-f', dest='filename', required=True, type=str, metavar='Filename',
                       help='Input CSV file containing the Gaia DR3 Source IDs')
    _args.add_argument('-s', dest='sampling', default=None, type=str, nargs=4,
                       metavar=('numpy function', 'start', 'end', 'step'),
                       help='Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120"')
    _args.add_argument('-t', dest='truncate', help='Truncate set of bases?', default=False, action='store_true')
    _args.add_argument('-o', dest='outputstyle', choices=('fits', 'txt'), default=None,
                       help='Output of produced spectra (default: None)')
    _args.add_argument('-i', dest='idcolname', help='Name of the column relating to DR3 Source ID',
                       default='source_id', metavar='Source ID column name')
    _args.add_argument('-n', dest='namecol', default='source_id', metavar='Object name column name',
                       help='Name of the column relating to the object name (for saving spectra)')
    _args.add_argument('-u', dest='uncalibrated', action='store_true', default=False,
                       help='Switch to uncalibrated spectra mode?')
    _args.add_argument('-x', dest='xp', default='RP', type=str, choices=('BP', 'RP'),
                       help='If in uncalibrated mode, RP or BP? (default: RP)')
    _args.add_argument('-v', dest='verbose', default=False, action='store_true', help='Print failure errors?')
    _args = _args.parse_args()
    if _args.sampling is not None:
        exec(f'_args.sampling = '
             f'np.{_args.sampling[0]}({_args.sampling[1]}, {_args.sampling[2]}, {_args.sampling[3]})')
    return _args


def str2float(s: str):
    """
    Convenience function for converting weird format of the flux/ error cells which sometimes occurs in dataframe

    Parameters
    ----------
    s: str
        What should be an array of floats, being interpreted as a string

    Returns
    -------
    _: np.ndarray
        Array of floats as parsed from string
    """
    return np.fromiter(filter(None, s[1:-1].replace('\n', '').split(' ')), float)


def getdispersion(xp: str = 'RP'):
    """
    Creates dispersion spline for converting absolute wavelength to pseudo-wavelength

    Parameters
    ----------
    xp: str
        BP or RP

    Returns
    -------

    """
    xp = xp.upper()
    if xp not in ('BP', 'RP'):
        raise ValueError('Needs to be BP or RP for dispersion')
    df = pd.read_csv('TabulatedDispersionFunction.csv')
    df = df.loc[df[f'{xp}_sample'] > -99]
    worig = df['Wavelength[nm]']
    rporig = df[f'{xp}_sample']
    p = interp1d(worig, rporig, kind=3, fill_value='extrapolate')
    return p
