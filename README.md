# GaiaXPy - Batch
A convenience script for submitting lots (>1000) Gaia DR3
source IDs through the `gaiaxpy.calibrate` function
(but can handle any number).
It allows one to download the spectra directly to 
individual fits or txt files.

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
# if venv isn't installed:
python3 -m pip install --user virtualenv
# then / otherwise
python3 -m venv env  # creates directory called 'env'
source env/bin/activate
python3 -m pip install -r requirements.txt
deactivate  # to leave
```
## Usage
Pass the `gaiaxpy_batch` a **csv** containing DR3
source IDs as discussed in the usage below:
```
usage: gaiaxpy_batch.py [-h] -f Filename -u Cosmos username -p [-s numpy function start end step] [-t] [-o {fits,txt}] [-i Source ID column name]
                        [-n Object name column name] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -f Filename           Input CSV file containing the Gaia DR3 Source IDs
  -u Cosmos username    Username on Gaia Archive
  -p                    Prompt Password on Gaia Archive
  -s numpy function start end step
                        Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120"
  -t                    Truncate set of bases?
  -o {fits,txt}         Output of produced spectra
  -i Source ID column name
                        Name of the column relating to DR3 Source ID
  -n Object name column name
                        Name of the column relating to the object name (for saving spectra)
  -v                    Print failure errors?

```
The filename, username and password are the only required inputs
(password will be prompted for, don't type it on the command line). 
We require a username & password as that allows less restrictive
access to the Gaia archive, c.f. 
[FAQ](https://www.cosmos.esa.int/web/gaia-users/archive/faq#account-limits-2020).

### Examples
When calling from the command line:
#### Minimum Example
```bash
python gaiaxpy_batch.py -f mysourceids.csv -u <username> -p
```

#### Expanded example
```bash
python gaiaxpy_batch.py -f mysourceids.csv -u <username> -p -s linspace 600 1050 120 -t -o fits -i source_id -n shortname
```
where `mysourceids.csv` looks like:
```csv
source_id, shortname
1, alice
2, bob
3, charles
...
```
#### Python example
When importing the package within python:
```python
import numpy as np
from gaiaxpy_batch import batch
kwargs = dict(sampling=np.linspace(600, 1050, 120), truncate=False, outputstyle=None)
batch('mysourceids.csv', **kwargs)
```

### `calibrate`
For more information on the `gaiaxpy.calibrate` function,
read [the tutorial](https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Calibrator%20tutorial.html).
Additionally, the documentation can be found
[here](https://gaiaxpy.readthedocs.io/en/latest/gaiaxpy.calibrator.html#gaiaxpy.calibrator.calibrator.calibrate).

## Citing
Please cite `gaiaxpy` following the instructions on
the [documentation page](https://gaiaxpy.readthedocs.io/en/latest/cite.html). The version of that code used here
is:  
Daniela Ruz-Mieres. (2022). gaia-dpci/GaiaXPy: GaiaXPy 1.1.0 (1.1.0). Zenodo. https://doi.org/10.5281/zenodo.6624491  

If you use the software from this repo, please cite
the latest version through this zenodo badge:
[![DOI](https://zenodo.org/badge/501233453.svg)](https://zenodo.org/badge/latestdoi/501233453)
