"""
Plots XP spectra saved by gaiaxpy-batch as fits files, and saves the plot as a png next to the first one.

Usage:
    python examples/plot_spectrum.py examples/outputspectra/example_1.fits [more.fits ...]
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')  # save to a file, no window needed
import matplotlib.pyplot as plt  # noqa: E402
from astropy.table import Table  # noqa: E402

COLOURS = ['black', '#004488', '#DDAA33', '#BB5566']  # colour-blind safe


def main(paths):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_prop_cycle(color=COLOURS)
    for path in paths:
        t = Table.read(path)
        ax.plot(t['wave'], t['flux'], lw=1, label=t.meta.get('OBJECT', os.path.basename(path)))
        ax.fill_between(t['wave'], t['flux'] - t['fluxerror'], t['flux'] + t['fluxerror'], alpha=0.3, lw=0)
    # calibrated spectra have astropy units, uncalibrated ones electron/s (not a FITS unit, so read as a string)
    unit = t['flux'].unit.to_string('latex_inline') if t['flux'].unit is not None else r'e$^{-}$ s$^{-1}$'
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel(f'Flux [{unit}]')
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.minorticks_on()
    ax.legend()
    out = os.path.splitext(paths[0])[0] + '.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f'Saved {out}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
