# GaiaXPy - Batch
Download the Gaia DR3 BP/RP ("XP") spectra of a list of sources, however long, and save them as
csv tables and, if you like, one fits or txt file per object.
You give it a csv file of Gaia DR3 source IDs; it fetches the spectra from the
[Gaia archive](https://gea.esac.esa.int/archive/) and turns them into flux against wavelength with
[GaiaXPy](https://github.com/gaia-dpci/GaiaXPy), the official Gaia tool for XP spectra.

Contents:
[Quick start](#quick-start) ·
[Your list of sources](#your-list-of-sources) ·
[Logging in](#logging-in) ·
[Examples](#examples) ·
[What you get](#what-you-get) ·
[All the options](#all-the-options) ·
[Troubleshooting](#troubleshooting) ·
[How it works](#how-it-works-and-do-you-need-it) ·
[For developers](#for-developers) ·
[Citing](#citing)

## Quick start
These steps assume no Python set up at all; skip any you have done before.
Type the commands into a terminal (on Windows, the "Miniforge Prompt" installed in step 1).

1. **Install conda** (a manager for Python and its packages), using [Miniforge](https://github.com/conda-forge/miniforge).
   On Linux or macOS:
   ```bash
   curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
   bash Miniforge3-$(uname)-$(uname -m).sh
   ```
   (answer "yes" to the questions, then close and reopen the terminal).
   On Windows, download and run the `Miniforge3-Windows-x86_64.exe` installer from the
   [Miniforge page](https://github.com/conda-forge/miniforge).
2. **Download this repository**, with git:
   ```bash
   git clone https://github.com/Will-Cooper/gaiaxpy-batch.git
   cd gaiaxpy-batch
   ```
   or with the green "Code" button, "Download ZIP" on GitHub, then unzip it and `cd` into the folder.
3. **Create the environment** (a separate Python with exactly the package versions this was tested with;
   this takes a few minutes and is only needed once):
   ```bash
   conda env create -f environment.yml
   ```
4. **Activate it** (needed in every new terminal):
   ```bash
   conda activate gaiaxpy-batch
   ```
5. **Check it works**:
   ```bash
   gaiaxpy-batch --version
   ```
   which prints something like `gaiaxpy-batch 2.0.0 (GaiaXPy 2.1.5)`.
6. **Run the example**, five sources (one of which has no XP spectrum), without logging in:
   ```bash
   gaiaxpy-batch -a -f examples/example_ids.csv -n name -o fits
   ```
   This takes from a few seconds to a few minutes, depending on how busy the archive is, and ends with
   `4 of 5 sources have XP spectra, 1 no XP spectrum`.
   The spectra are in `examples/outputspectra/` (`example_1.fits`, ...) and in the csv files next to the input.
7. **Plot one**:
   ```bash
   python examples/plot_spectrum.py examples/outputspectra/example_1.fits examples/outputspectra/example_2.fits
   ```
   which saves `examples/outputspectra/example_1.png`.

To update later, run `git pull` (or download the ZIP again), then
`conda env update -f environment.yml --prune`.

Alternatively, without conda (Python 3.11 or newer):
```bash
python3 -m venv env            # creates a folder called env
source env/bin/activate        # on Windows: env\Scripts\activate
python -m pip install -r requirements.txt
python -m pip install --no-deps .
```

## Your list of sources
The input is a csv (comma separated values) file with a column of Gaia DR3 source IDs, called `source_id`
(or any name, given with `-i`), and any other columns you like, which are kept in the outputs:
```csv
source_id,shortname
138877256416618112,alice
142522579844473600,bob
Gaia DR3 208760841091097472,charles
```
* The IDs can be plain numbers or `Gaia DR3 <number>`, as SIMBAD and VizieR write them.
* Rows without a valid source ID are listed in the failed csv (see [What you get](#what-you-get)).
  If the IDs look like `1.38877E+17`, a spreadsheet has rounded them: they have lost digits and cannot be used.
  Make the column "text" in the spreadsheet, or keep the file as a csv.
* A column of names (`-n shortname`) names the per-object files; otherwise they are `GaiaDR3_<source_id>`.

Where to get source IDs:
* **The Gaia archive**: [gea.esac.esa.int/archive](https://gea.esac.esa.int/archive/), "Search", "Advanced (ADQL)".
  [`examples/find_ids.adql`](examples/find_ids.adql) is an example query (the 2237 sources within 20 pc with XP spectra);
  download the result as csv.
* **SIMBAD**: the "Gaia DR3 ..." identifier of each object.
* **Cross-matches** of your own catalogue with `gaiadr3.gaia_source`.

Only about 220 million sources (mostly brighter than G = 17.65 mag) have XP spectra:
those with `has_xp_continuous = 'True'` in `gaiadr3.gaia_source`, or in `gaiadr3.xp_summary`.
The others are listed as `no XP spectrum`, so checking first saves time on long lists.

## Logging in
Without logging in (`-a`), the archive allows only 50 source IDs per request.
Logged in, it allows 5000, so long lists download much faster
(see the archive's [user quotas](https://www.cosmos.esa.int/web/gaia-users/archive/collaborate#GaiaArchiveUserQuotas)).
Accounts are free: [register here](https://www.cosmos.esa.int/web/gaia-users/register).

Without `-a`, the script asks for your username and password.
To avoid typing the username (or to run it without a terminal, e.g. overnight with `nohup`), set it first:
```bash
export GAIA_USERNAME=your_username          # Linux, macOS
set GAIA_USERNAME=your_username             # Windows Command Prompt
$env:GAIA_USERNAME = "your_username"        # Windows PowerShell
```
`GAIA_PASSWORD` works the same way, but leaves the password in your shell history, so typing it when asked is safer.

## Examples
All of these are run from the `gaiaxpy-batch` folder, with the environment activated.

Calibrated spectra on the default wavelength grid, only as csv tables:
```bash
gaiaxpy-batch -f mysourceids.csv
```
One fits file per object, named from the `shortname` column:
```bash
gaiaxpy-batch -f mysourceids.csv -n shortname -o fits
```
Your own wavelength grid (120 points from 600 to 1050 nm) with the smoother, truncated spectra, as txt files:
```bash
gaiaxpy-batch -f mysourceids.csv -s linspace 600 1050 120 -t -o txt
```
The uncalibrated (internally calibrated) BP spectrum, in electron/s:
```bash
gaiaxpy-batch -f mysourceids.csv -u -x BP -o fits
```
A long list, logged in, with details of every request:
```bash
export GAIA_USERNAME=your_username
gaiaxpy-batch -f my100000sources.csv -o fits -v
```
Trying again the sources that failed last time (the failed csv can be the input; see [Troubleshooting](#troubleshooting)):
```bash
gaiaxpy-batch -f mysourceids_XP_spectra_failed.csv -o fits
```
From Python (e.g. a Jupyter notebook) instead of the command line, with the same options as keywords:
```python
import numpy as np
from gaiaxpy_batch import batch

df = batch('mysourceids.csv', sampling=np.linspace(600, 1050, 120), outputstyle='fits', anonymous=True)
# df is your table with flux and flux_error arrays added
```

## What you get
For an input `mysourceids.csv`, next to it:

| File | Contents |
|---|---|
| `mysourceids_XP_spectra.csv` | your table, with `flux` and `flux_error` columns added (empty where there is no spectrum) |
| `mysourceids_XP_spectra_cut.csv` | the same, only the rows with spectra |
| `mysourceids_XP_spectra_failed.csv` | the sources without spectra, and why (only if there are any) |
| `outputspectra/<name>.fits` or `.txt` | with `-o`: wavelength, flux and flux error, one file per object |

* **Wavelengths** are in nm. By default, calibrated spectra are sampled from 336 to 1020 nm in 2 nm steps
  (GaiaXPy's default), and uncalibrated ones at pseudo-wavelengths 0 to 60 (converted to nm).
* **Fluxes** are in W m<sup>-2</sup> nm<sup>-1</sup> (calibrated) or electron s<sup>-1</sup> (uncalibrated, `-u`).
  `electron/s` is not a standard FITS unit, so astropy warns when reading it.
* **Partial spectra:** a calibrated spectrum is NaN (not a number) where a band is missing, e.g. below about 640 nm
  (roughly 44 per cent of the default grid) for a source with only an RP spectrum, which is common for very red
  sources such as ultracool dwarfs. These still count as successes; the script warns how many there are.
  **Check for NaNs before integrating or fitting a spectrum** (e.g. `numpy.isfinite(flux)`), or functions such as
  `numpy.trapezoid` will return NaN.
* The wiggles in XP spectra are real: they come from the functions the spectra are made of, and are strongest in
  faint sources and at the blue end. Truncation (`-t`) smooths them.
* The fits files have the header keywords `SOURCEID`, `OBJECT`, `XPMODE` (calibrate, or convert BP/RP) and
  `GAIAXPY` (its version). The txt files have a one line header.
* The csv `flux` and `flux_error` columns are arrays written as text; to read them back:
  ```python
  import pandas as pd
  from gaiaxpy_batch import str2float
  df = pd.read_csv('mysourceids_XP_spectra_cut.csv')
  flux = str2float(df.flux[0])
  ```

The `reason` column of the failed csv is one of:

| Reason | Meaning | What to do |
|---|---|---|
| `no XP spectrum` | Gaia DR3 has no XP spectrum for this source | nothing, most sources fainter than G = 17.65 mag have none |
| `no BP spectrum` / `no RP spectrum` | with `-u`, the source has only the other band | try the other band (`-x`), or calibrated spectra |
| `invalid source ID` | the row has no usable source ID | correct the input |
| `calibration failed` | GaiaXPy could not use the downloaded spectrum | rare, check the source in the archive |
| `download failed` | the archive failed every retry | try again later with the failed csv as input |
| `not attempted` | the script stopped (the archive seemed down, or Ctrl-C) | try again later with the failed csv as input |

Loading a fits spectrum in Python:
```python
from astropy.table import Table
t = Table.read('examples/outputspectra/example_1.fits')
print(t['wave'], t['flux'], t['fluxerror'], t.meta['SOURCEID'])
```

## All the options
```
gaiaxpy-batch -f FILE [options]
```
| Option | What it does |
|---|---|
| `-f FILE`, `--filename` | the input csv (the only one you must give) |
| `-i NAME`, `--idcolname` | the column with the source IDs (default `source_id`) |
| `-n NAME`, `--namecol` | the column with names for the per-object files (default `GaiaDR3_<source_id>`) |
| `-o fits` or `-o txt`, `--outputstyle` | also write one file per object, in `outputspectra/` |
| `-s FUNCTION START END STEP`, `--sampling` | the wavelength grid in nm: `linspace 600 1050 120` (120 points), `arange 400 1000 2` (every 2 nm) or `geomspace ...`; 330 to 1050 nm |
| `-t`, `--truncate` | use only the basis functions GaiaXPy recommends: smoother, less noisy spectra |
| `-u`, `--uncalibrated` | the uncalibrated BP or RP spectrum (`gaiaxpy.convert`) instead of the calibrated one (`gaiaxpy.calibrate`) |
| `-x BP` or `-x RP`, `--xp` | with `-u`, which band (default RP) |
| `-a`, `--anonymous` | don't log in (50 sources per request) |
| `-c N`, `--chunksize` | source IDs per request (default 50 anonymous, 1000 logged in; the archive allows up to 50 and 5000) |
| `-w N`, `--workers` | requests at the same time (default 4, at most 5 logged in) |
| `--retries N` | how often to retry a request the archive failed (default 3) |
| `--retry-delay S` | seconds before the first retry, doubled each time (default 30, so 30, 60, 120) |
| `--timeout S` | seconds without an answer before a request is given up (default 900) |
| `-v`, `--verbose` | print details of every request, login and failure |
| `--version`, `-h` | the version, or this list |

The exit status is 0 if everything worked (or has no XP spectrum), 1 if some sources could not be downloaded,
and 2 if nothing could be downloaded or the input is wrong, for use in scripts.

## Troubleshooting
* **It sits at 0 per cent for a long time.** The progress bar moves when a whole request finishes, and the archive
  can take many minutes per request when it is busy. Use `-v` to see the logins and requests.
* **`Error 500`, `Error 503`, `Could not open JPA EntityManager`, `File is not a zip file`.**
  The archive is struggling. These are retried automatically (after 30, 60 and 120 s);
  if three requests in a row fail every retry, the script stops, saves what it has, and lists the rest as
  `not attempted`. Try again later with the failed csv as the input. The website
  ([gea.esac.esa.int/archive](https://gea.esac.esa.int/archive/)) shows whether the archive is up.
* **`HTTP Status 401 - Unauthorized: Credentials expiration`.** The archive ends logins after about 30 minutes;
  the script logs in again by itself. No need to do anything.
* **`Could not log in to the Gaia archive`.** Check the username and password (on the archive website), or the
  archive is down.
* **Stopping it.** Ctrl-C stops new requests, waits for those in progress, and saves what has finished;
  a second Ctrl-C quits straight away, without saving.
* **`Column 'source_id' not found`.** Give the name of your source ID column with `-i`.
* **`conda: command not found`.** Close and reopen the terminal after installing Miniforge, or on Windows use the
  "Miniforge Prompt".
* **`gaiaxpy-batch: command not found`.** Activate the environment first: `conda activate gaiaxpy-batch`
  (or run `python gaiaxpy_batch.py ...` from this folder).

## How it works, and do you need it?
GaiaXPy itself accepts a list of source IDs, but sends them all to the archive in one request, logging in each time,
so it cannot fetch more than 50 sources at once anonymously (5000 logged in), and one failure loses the whole list.
For a handful of sources, GaiaXPy alone is just as good:
```python
from gaiaxpy import calibrate
df, sampling = calibrate([138877256416618112, 142522579844473600], save_file=False)
```
This script downloads the XP spectra (as the coefficients of basis functions) itself with `astroquery`,
exactly as GaiaXPy does, and hands them to GaiaXPy to calibrate or convert, so the spectra are GaiaXPy's own. It:
* splits the list into requests (50 anonymously, 1000 logged in, smaller if needed to keep every worker busy) and
  sends several at once, each worker logging in once (again when the archive ends the login),
* retries requests the archive fails with increasing waits, and stops if the archive seems to be down,
* splits a request in half when the archive rejects it as too large or as a bad request,
* splits spectra that fail to calibrate in half and calibrates each half again, without downloading them again,
* saves everything it could get, lists the rest and why, and saves what has finished when stopped with Ctrl-C,
* checks the input before asking for a password or downloading anything,
* converts your wavelength grid to the pseudo-wavelengths uncalibrated spectra need.

For millions of sources, downloading the archive's
[bulk XP files](https://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/)
and giving them to GaiaXPy locally is quicker.

More about the GaiaXPy functions used:
[calibrate tutorial](https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Calibrator%20tutorial.html),
[convert tutorial](https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Converter%20tutorial.html).

**Gaia DR4** (due 2 December 2026): GaiaXPy 2.1.5 only retrieves Gaia DR3 spectra;
this script will be updated once GaiaXPy supports DR4.

## For developers
* `requirements.txt` has the exact versions tested (`requirements-dev.txt` adds the test tools), and
  `environment.yml` installs them; `pyproject.toml` has the looser ranges used by `pip install .`.
  To regenerate the pins, in a fresh Python 3.12 environment run `python -m pip install . && python -m pip freeze --exclude-editable`
  and drop gaiaxpy-batch, pip, setuptools and wheel.
  If your conda adds packages to every new environment (`conda config --show create_default_packages`),
  `conda env create` adds them here too; the pinned versions are still the ones installed.
* Tests (the archive and GaiaXPy are mocked, except for the network test):
  ```bash
  python -m pytest --cov=gaiaxpy_batch       # offline
  python -m pytest -m network                # queries the Gaia archive
  ruff check .
  ```
* GitHub Actions runs the tests on the pinned versions and on the newest allowed versions (Python 3.11 to 3.14,
  Linux, macOS and Windows), weekly as well as on every pull request.
* Changes are listed in [CHANGELOG.md](CHANGELOG.md).

## Citing
Please cite GaiaXPy following the instructions on its
[documentation page](https://gaiaxpy.readthedocs.io/en/latest/cite.html); this version was tested with GaiaXPy 2.1.5.

If you use the software from this repository, please cite the latest version through this Zenodo badge:
[![DOI](https://zenodo.org/badge/501233453.svg)](https://zenodo.org/badge/latestdoi/501233453)
