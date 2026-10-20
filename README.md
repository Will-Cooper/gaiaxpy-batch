# GaiaXPy - Batch
A convenience script for getting XP spectra for lots (>50) of Gaia DR3
source IDs, calibrated with `gaiaxpy.calibrate` or converted with `gaiaxpy.convert`
(but it can handle any number).
It allows one to download the spectra directly to
individual fits or txt files.

## Do you need this?
[GaiaXPy](https://github.com/gaia-dpci/GaiaXPy) itself accepts a list of source IDs
(or an ADQL query), but it sends all of them to the Gaia archive in a single DataLink request,
logging in again each time.
The [Gaia archive user quotas](https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas)
limit DataLink requests to **50 source IDs anonymously and 5000 logged in**,
and allow 5 active jobs per registered user.
A larger request is rejected ("Maximum number of Source Ids reached (max: 50, ...)"),
so GaiaXPy alone cannot retrieve more than 50 sources at once anonymously,
and one failure (the archive is not always reliable) loses the whole list.

This script downloads the XP continuous spectra (the basis function coefficients) itself with `astroquery`,
in the same way GaiaXPy does, and hands them to GaiaXPy as a DataFrame to calibrate or convert,
so the spectra are identical to GaiaXPy's own. It:
* splits the list into chunks (default 50 source IDs anonymously, 1000 logged in) and downloads them in
  parallel (default 4 workers), with one session per worker that logs in once,
* retries a request that fails on the server's side (HTTP 5xx, broken responses, time-outs) with increasing waits
  (30, 60, 120 s by default), and stops if several requests in a row fail every retry, as the archive is then
  probably down,
* when the archive reports its maximum number of source IDs per request, splits all chunks to that size,
* splits a chunk rejected as a bad request in half and tries each half again, down to single source IDs,
* splits spectra that fail to calibrate in half and calibrates each half again (locally, without downloading again),
* records the source IDs without spectra, and why, in a separate csv that can be used as the input to try again,
* writes each spectrum to its own fits or txt file, with the correct wavelength grid and units,
* takes the sampling in absolute wavelength [nm] for uncalibrated (`convert`) spectra too.

For a handful of sources (at most 50 when anonymous), calling GaiaXPy directly is just as good:
```python
from gaiaxpy import calibrate
df, sampling = calibrate([138877256416618112, 142522579844473600], save_file=False)
```
For millions of sources, it is better to download the
[bulk XP continuous files](https://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/)
and pass them to GaiaXPy locally.

## Installation
### Downloading:
```bash
git clone git@github.com:Will-Cooper/gaiaxpy-batch.git
cd gaiaxpy-batch/
```
To keep updated:
```bash
git pull
```
### Set-up:
Requires Python >= 3.11 (as does GaiaXPy >= 2.1.5).
#### Conda
The best way of setting up a clean environment with this
package & `gaiaxpy` is to use `conda`.
```bash
conda env create -f environment.yml
conda activate gaiaxpy-batch
conda deactivate # to leave environment
```
#### Venv (pip)
Alternatively, to create a `pip` environment:
```bash
python3 -m venv env  # creates directory called 'env'
source env/bin/activate
python3 -m pip install .  # or: python3 -m pip install -r requirements.txt
deactivate  # to leave
```
Installing (`pip install .`) adds a `gaiaxpy-batch` command;
`python gaiaxpy_batch.py` works without installing.

## Usage
Pass the `gaiaxpy_batch` a **csv** containing DR3
source IDs as discussed in the usage below:
```
usage: gaiaxpy-batch [-h] -f Filename [-s numpy function start end step] [-t] [-o {fits,txt}]
                     [-i Source ID column name] [-n Object name column name] [-u] [-x {BP,RP}]
                     [-c CHUNKSIZE] [-w WORKERS] [--retries RETRIES] [--retry-delay RETRY_DELAY]
                     [--timeout TIMEOUT] [-a] [-v] [--version]

options:
  -h, --help            show this help message and exit
  -f, --filename Filename
                        Input CSV file containing the Gaia DR3 Source IDs
  -s, --sampling numpy function start end step
                        Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120"
                        (functions: linspace, arange, geomspace)
  -t, --truncate        Truncate set of bases?
  -o, --outputstyle {fits,txt}
                        Output of produced spectra (default: None)
  -i, --idcolname Source ID column name
                        Name of the column relating to DR3 Source ID
  -n, --namecol Object name column name
                        Name of the column relating to the object name (for saving spectra,
                        default: GaiaDR3_<source_id>)
  -u, --uncalibrated    Switch to uncalibrated spectra mode?
  -x, --xp {BP,RP}      If in uncalibrated mode, RP or BP? (default: RP)
  -c, --chunksize CHUNKSIZE
                        Number of source IDs per archive request (default: 50 anonymous, 1000 logged in;
                        the archive allows up to 50 and 5000)
  -w, --workers WORKERS
                        Number of parallel archive requests (default: 4, at most 5 when logged in)
  --retries RETRIES     Retries of a request after a server-side failure (default: 3)
  --retry-delay RETRY_DELAY
                        Seconds to wait before the first retry, doubled for each one after (default: 30)
  --timeout TIMEOUT     Seconds without a response before a connection is given up (default: 900)
  -a, --anonymous       Query the archive without logging in (otherwise uses $GAIA_USERNAME and
                        $GAIA_PASSWORD, or prompts)
  -v, --verbose         Print failure errors?
  --version             show program's version number and exit
```
The filename is the only required input.
The username and password are read from the `GAIA_USERNAME` and `GAIA_PASSWORD`
environment variables, or prompted for, unless `-a` is given.
We recommend logging in, as that allows 5000 rather than 50 source IDs per request, c.f. the
[user quotas](https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas).
Each worker logs in once; if the login fails (astroquery only logs this, the script checks),
it is retried like any other request and the script then stops,
as the username or password is probably wrong, or the archive is down.

It's recommended that you check that your source IDs have XP spectra available.
You can do this by cross-matching with the `gaiadr3.xp_summary` table
(or `has_xp_continuous` in `gaiadr3.gaia_source`).

### Outputs
For an input `mysourceids.csv`:
* `mysourceids_XP_spectra.csv`: the input csv with `flux` and `flux_error` array columns added,
* `mysourceids_XP_spectra_cut.csv`: the same, only for the rows with spectra,
* `mysourceids_XP_spectra_failed.csv` (if any): the source IDs without spectra, with a `reason`:
  `no XP spectrum`, `calibration failed`, `download failed` (failed every retry) or
  `not attempted, archive unavailable` (the script stopped),
* `outputspectra/<name>.fits` or `.txt` (with `-o`): wavelength [nm], flux and flux error per object.

If some sources could not be downloaded, the script says so and exits with status 1 (2 if nothing was downloaded).
The failed csv can be given straight back to the script, later, to try them again:
```bash
gaiaxpy-batch -f mysourceids_XP_spectra_failed.csv -o fits
```

Fluxes are in W m<sup>-2</sup> nm<sup>-1</sup> (calibrated) or electron s<sup>-1</sup> (uncalibrated, `-u`;
`electron/s` is not a FITS standard unit, so astropy warns when reading it).
The default wavelength grid is GaiaXPy's: 336-1020 nm in 2 nm steps when calibrated,
and pseudo-wavelengths 0-60 (converted to nm) when uncalibrated.
The array columns in the csvs can be read back with `gaiaxpy_batch.str2float`.

### Examples
When calling from the command line:
#### Minimum Example
```bash
gaiaxpy-batch -f mysourceids.csv
```

#### Expanded example
```bash
gaiaxpy-batch -f mysourceids.csv -s linspace 600 1050 120 -t -o fits -i source_id -n shortname
```
where `mysourceids.csv` looks like:
```csv
source_id,shortname
1,alice
2,bob
3,charles
...
```
#### Python example
When importing the package within python:
```python
import numpy as np
from gaiaxpy_batch import batch
kwargs = dict(sampling=np.linspace(600, 1050, 120), truncate=False, outputstyle=None)
df = batch('mysourceids.csv', **kwargs)
```

### `calibrate`
For more information on the `gaiaxpy.calibrate` function,
read [the tutorial](https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Calibrator%20tutorial.html).
Additionally, the documentation can be found
[here](https://gaiaxpy.readthedocs.io/en/latest/gaiaxpy.calibrator.html#gaiaxpy.calibrator.calibrator.calibrate).

### `convert`
For more information on the `gaiaxpy.convert` function,
read [the tutorial](https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Converter%20tutorial.html).
Additionally, the documentation can be found
[here](https://gaiaxpy.readthedocs.io/en/latest/gaiaxpy.converter.html#gaiaxpy.converter.converter.convert).

### Gaia DR4
GaiaXPy 2.1.5 retrieves Gaia DR3 spectra only.
This script will be updated for Gaia DR4 (due 2 December 2026) once GaiaXPy supports it.

## Testing
```bash
python -m pip install .[dev]
python -m pytest              # offline, the archive is mocked
python -m pytest -m network   # queries the Gaia archive
```

## Citing
Please cite `gaiaxpy` following the instructions on
the [documentation page](https://gaiaxpy.readthedocs.io/en/latest/cite.html).
This version was tested with GaiaXPy 2.1.5.

If you use the software from this repo, please cite
the latest version through this zenodo badge:
[![DOI](https://zenodo.org/badge/501233453.svg)](https://zenodo.org/badge/latestdoi/501233453)
