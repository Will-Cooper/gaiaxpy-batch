# Changelog

## 2.0.1 (2026-10-23)
For GaiaXPy 2.1.5 and the Gaia archive's October 2026 request limits. Needs Python 3.11 or newer.

### Changed
- Downloads the XP continuous spectra with astroquery, in requests of at most 50 source IDs anonymously or
  1000 logged in (the archive allows 50 and 5000), several at once (`-w`), each worker logging in once,
  and hands them to GaiaXPy to calibrate or convert. GaiaXPy alone sends a whole list in one request.
- Retries requests the archive fails (HTTP 5xx, broken responses, time-outs) with increasing waits
  (`--retries`, `--retry-delay`, `--timeout`), logs in again when the archive ends a login, and stops if the
  archive seems to be down, saving what it has.
- Ctrl-C saves what has finished; a second Ctrl-C quits at once.
- The failed csv lists every source without a spectrum with a reason, and can be used as the input to try again.
  Exit status 1 if some sources could not be downloaded, 2 if none could or the input is wrong.
- Checks the input csv and options before asking for a password or downloading anything.
- Source IDs are read exactly (a blank cell no longer corrupts 19 digit IDs), and `Gaia DR3 <id>` names work.
- `-x` defaults to unset (RP), so giving it without `-u` can be warned about.
- Packaged with `pyproject.toml`, with a `gaiaxpy-batch` command; `python gaiaxpy_batch.py` still works.
  `utils.py` and `TabulatedDispersionFunction.csv` are gone (GaiaXPy has the same dispersion function).
- `requirements.txt` and `environment.yml` install the exact versions tested.

### Fixed
- The wavelength grid of saved spectra is the one GaiaXPy used (the default grid is 336 to 1020 nm,
  not the 336 to 1200 nm previously written).
- Uncalibrated mode without `-s`, and uncalibrated BP spectra, no longer fail.
- Saving without `-n`, with blank, shared or unsafe names, no longer fails or overwrites files.
- Arrays of more than 1000 points are written to the csv in full.
- Output names for paths containing `.csv`, and `.CSV` inputs.
- The flux unit of uncalibrated fits files (`electron/s`) is kept.
- The sampling (`-s`) is no longer passed to `exec`.
- astroquery's temporary files no longer end up in the working directory.

## 1.2.0 (2022-06-16)
The last release for GaiaXPy 1.1.1 and Python 3.9.
