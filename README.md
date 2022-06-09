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
The best way of setting up a clean environment with this
package & `gaiaxpy` is to use `conda`.
```bash
conda env create -f environment.yml
conda activate gaiaxpy-batch
conda deactivate # to leave environment
```
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
usage: gaiaxpy_batch.py [-h] -f Filename [-s numpy function start end step] [-u Cosmos username] [-p] [-t] [-o {fits,txt}] [-i Source ID column name]
                        [-n Object name column name]

optional arguments:
  -h, --help            show this help message and exit
  -f Filename           Input CSV file containing the Gaia DR3 Source IDs
  -s numpy function start end step
                        Wavelength [absolute nm] sampling to be retrieved, e.g. "linspace 600 1050 120"
  -u Cosmos username    Username on Gaia Archive
  -p                    Prompt Password on Gaia Archive
  -t                    Truncate set of bases?
  -o {fits,txt}         Output of produced spectra
  -i Source ID column name
                        Name of the column relating to DR3 Source ID
  -n Object name column name
                        Name of the column relating to the object name (for saving spectra)

```
The filename is the only required input but most recent
testing showed that Cosmos username & password are also
required.

### Example
When calling from the command line:
```bash
python gaiaxpy_batch.py -f mysourceids.csv -s linspace 600 1050 120 -u wcooper -p -t -o fits -i source_id -n shortname
```
where `mysourceids.csv` looks like:
```csv
source_id, shortname
1, alice
2, bob
3, charles
...
```
### Other Example
When importing the package within python:
```python
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
Daniela Ruz-Mieres. (2022). gaia-dpci/GaiaXPy: GaiaXPy 1.0.2 (1.0.2). Zenodo. https://doi.org/10.5281/zenodo.6569912  

If you use the software from this repo, please use this zenodo badge:  [![DOI](https://zenodo.org/badge/501233453.svg)](https://zenodo.org/badge/latestdoi/501233453)
