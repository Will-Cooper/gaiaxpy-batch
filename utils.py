import argparse
import getpass

import numpy as np


class PasswordPromptAction(argparse.Action):
    """
    Taken from https://stackoverflow.com/questions/27921629/python-using-getpass-with-argparse
    For prompting a hidden password
    """
    def __init__(self,
                 option_strings,
                 dest=None,
                 nargs=0,
                 default=None,
                 required=False,
                 type=None,  # if you're looking why this shadows outerscope, blame argparse
                 metavar=None,
                 help=None):  # if you're looking why this shadows outerscope, blame argparse
        super(PasswordPromptAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=nargs,
            default=default,
            required=required,
            metavar=metavar,
            type=type,
            help=help)

    def __call__(self, parser, args, values, option_string=None):
        password = getpass.getpass()
        setattr(args, self.dest, password)


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
    _args.add_argument('-u', dest='username', type=str, metavar='Cosmos username', required=True,
                       help='Username on Gaia Archive')
    _args.add_argument('-p', dest='password', action=PasswordPromptAction, type=str,
                       help='Prompt Password on Gaia Archive', required=True)
    _args.add_argument('-s', dest='sampling', default=None, type=str, nargs=4,
                       metavar=('numpy function', 'start', 'end', 'step'),
                       help='Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120"')
    _args.add_argument('-t', dest='truncate', help='Truncate set of bases?', default=False, action='store_true')
    _args.add_argument('-o', dest='outputstyle', help='Output of produced spectra', choices=('fits', 'txt'),
                       default=None)
    _args.add_argument('-i', dest='idcolname', help='Name of the column relating to DR3 Source ID',
                       default='source_id', metavar='Source ID column name')
    _args.add_argument('-n', dest='namecol', default='source_id', metavar='Object name column name',
                       help='Name of the column relating to the object name (for saving spectra)')
    _args.add_argument('-v', dest='verbose', default=False, action='store_true', help='Print failure errors?')
    _args = _args.parse_args()
    if _args.sampling is not None:
        exec(f'_args.sampling = '
             f'np.{_args.sampling[0]}({_args.sampling[1]}, {_args.sampling[2]}, {_args.sampling[3]})')
    return _args


def str2arr(s: str):
    """
    Convenience function for converting weird format of the flux/ error cells in dataframe

    Parameters
    ----------
    s: str
        What should be an array of floats, being interpreted as a string

    Returns
    -------
    _: np.ndarray
        Array of floats as parsed from string
    """
    return np.fromiter(s[1:-1].replace('\n', '').split('  '), float)
