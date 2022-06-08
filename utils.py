import argparse
import getpass

import numpy as np


class PasswordPromptAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest=None,
                 nargs=0,
                 default=None,
                 required=False,
                 type=None,
                 metavar=None,
                 help=None):
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
        The different argument parameters, can be grabbed via their long names (e.g. _args.host)
    """
    _args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    _args.add_argument('-f', dest='filename', required=True, type=str,
                       help='Input CSV file containing the Gaia DR3 Source IDs')
    _args.add_argument('-s', dest='sampling', default=None, type=str,
                       help='Wavelength [absolute nm] sampling to be retrieved, e.g. "np.linspace(600, 1050, 120)"')
    _args.add_argument('-u', dest='username', default=None, type=str,
                       help='Username on Gaia Archive')
    _args.add_argument('-p', dest='password', action=PasswordPromptAction, default=None, type=str,
                       help='Password on Gaia Archive')
    _args.add_argument('-t', dest='truncate', help='Truncate set of bases?', default=False, action='store_true')
    _args.add_argument('-o', dest='outputstyle', help='Output of produced spectra', choices=('fits', 'txt'),
                       default=None)
    _args.add_argument('-i', dest='idcolname', help='Name of the column relating to DR3 Source ID',
                       default='source_id')
    _args.add_argument('-n', dest='namecol', default='source_id',
                       help='Name of the column relating to the object name (for saving spectra)')
    _args = _args.parse_args()
    if _args.sampling is not None:
        exec('_args.sampling = ' + str(_args.sampling))
    return _args


def str2arr(s):
    return np.fromiter(s[1:-1].replace('\n', '').split('  '), float)

