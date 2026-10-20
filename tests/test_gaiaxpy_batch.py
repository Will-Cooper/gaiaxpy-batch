"""
Tests for gaiaxpy_batch. The archive is mocked out except for the tests marked network,
which are deselected by default (run them with pytest -m network).
"""
import os

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table
from gaiaxpy import wl_to_pwl

import gaiaxpy_batch as gb

BAD = 999  # source ID whose requests always fail, like a timeout
NOXP = 888  # source ID with no XP spectrum
DEFAULT_WAVE = np.arange(336., 1021., 2.)  # gaiaxpy.calibrate default sampling


class FakeArchive:
    """Stands in for gaiaxpy.calibrate and gaiaxpy.convert, counting the requests made"""

    def __init__(self, error=ConnectionError):
        self.calls = []
        self.error = error

    def _check(self, idlist):
        self.calls.append(list(idlist))
        if BAD in idlist:
            raise self.error('archive timeout')
        good = [xid for xid in idlist if not 800 <= xid < 900]  # 800-899 have no XP spectra
        if not good:
            raise ValueError('No continuous BP/RP data found for the given sources.')
        return good

    def calibrate(self, idlist, sampling=None, **kwargs):
        good = self._check(idlist)
        sampling = DEFAULT_WAVE if sampling is None else sampling
        return pd.DataFrame({'source_id': good, 'flux': [sampling * xid for xid in good],
                             'flux_error': [np.ones_like(sampling) for _ in good]}), sampling

    def convert(self, idlist, sampling=None, **kwargs):
        good = self._check(idlist)
        rows = [{'source_id': xid, 'xp': xp, 'flux': sampling.copy(), 'flux_error': np.ones_like(sampling)}
                for xid in good for xp in ('BP', 'RP')]
        return pd.DataFrame(rows), sampling


@pytest.fixture
def archive(monkeypatch):
    fake = FakeArchive()
    monkeypatch.setattr(gb, 'calibrate', fake.calibrate)
    monkeypatch.setattr(gb, 'convert', fake.convert)
    return fake


def test_parse_sampling():
    assert np.allclose(gb.parse_sampling('linspace', '600', '1050', '120'), np.linspace(600, 1050, 120))
    assert np.allclose(gb.parse_sampling('arange', 400, 1000, 2), np.arange(400, 1000, 2))
    with pytest.raises(ValueError):
        gb.parse_sampling('__import__("os").system', 1, 2, 3)


def test_sysargs():
    args = gb.sysargs(['-f', 'a.csv', '-s', 'linspace', '600', '1000', '50', '-x', 'bp', '-u'])
    assert len(args.sampling) == 50 and args.xp == 'BP' and args.namecol is None and args.chunksize == 1000


def test_get_sampling():
    wave, pwl = gb.get_sampling(np.linspace(340, 640, 30), uncalibrated=True, xp='BP')
    assert np.all(np.diff(pwl) > 0)  # GaiaXPy needs ascending pseudo-wavelength
    assert np.allclose(pwl, wl_to_pwl('BP', wave))
    assert gb.get_sampling(None) == (None, None)
    with pytest.raises(ValueError):
        gb.get_sampling(np.linspace(300, 1100, 10))
    with pytest.raises(ValueError):
        gb.get_sampling(np.linspace(1000, 600, 10))


def test_bisection(archive):
    ids = list(range(1, 26)) + [BAD]
    stackdf, wave, failed = gb.download_xp(ids, chunksize=10, workers=2)
    assert failed == [BAD]
    assert sorted(stackdf.source_id) == list(range(1, 26))
    assert np.allclose(wave, DEFAULT_WAVE)
    assert all(len(f) == len(DEFAULT_WAVE) for f in stackdf.flux)
    assert [BAD] in archive.calls  # split down to the single bad ID


def test_no_data_not_bisected(archive):
    ids = [1, 2, 801, 802, 802]  # duplicate dropped
    stackdf, wave, failed = gb.download_xp(ids, chunksize=2)
    assert failed == [801, 802] and len(archive.calls) == 2  # the chunk with no XP spectra is not split
    with pytest.raises(ValueError, match='All batches failed'):
        gb.download_xp([NOXP])


def test_fatal_error_not_retried(monkeypatch):
    fake = FakeArchive(error=TypeError)
    monkeypatch.setattr(gb, 'calibrate', fake.calibrate)
    with pytest.raises(TypeError):
        gb.download_xp([1, 2, BAD, 4])
    assert len(fake.calls) == 1


def test_convert_bp_ascending(archive):
    stackdf, wave, failed = gb.download_xp([1, 2], uncalibrated=True, xp='BP', sampling=np.linspace(340, 640, 50))
    assert np.all(np.diff(wave) > 0)
    for flux in stackdf.flux:  # the fake flux is the pseudo-wavelength, so it should follow the wavelengths
        assert np.allclose(flux, wl_to_pwl('BP', wave))
    stackdf, wave, failed = gb.download_xp([1, 2], uncalibrated=True, xp='RP')  # default sampling
    assert len(wave) == len(gb.DEFAULT_PWL) and np.all(np.diff(wave) > 0)


@pytest.mark.parametrize('outputstyle', ['txt', 'fits'])
@pytest.mark.parametrize('namecol', [None, 'shortname'])
def test_batch_saves(archive, tmp_path, outputstyle, namecol):
    fname = str(tmp_path / 'ids.csv')
    pd.DataFrame({'source_id': [1, 2, 3, NOXP], 'shortname': ['alice', 'bob', 'charles', 'dave']}).to_csv(fname,
                                                                                                          index=False)
    sampling = np.linspace(400, 1000, 1500)  # more than numpy's 1000 element summarisation threshold
    combdf = gb.batch(fname, sampling=sampling, outputstyle=outputstyle, namecol=namecol, anonymous=True,
                      chunksize=2)
    assert len(combdf) == 4 and combdf.flux.isna().sum() == 1
    cut = pd.read_csv(fname.replace('.csv', '_XP_spectra_cut.csv'))
    assert len(cut) == 3
    assert np.allclose(gb.str2float(cut.flux.iloc[0]), sampling * cut.source_id.iloc[0])
    assert pd.read_csv(fname.replace('.csv', '_XP_spectra_failed.csv')).source_id.tolist() == [NOXP]
    outdir = tmp_path / 'outputspectra'
    files = sorted(os.listdir(outdir))
    expected = ['alice', 'bob', 'charles'] if namecol else ['GaiaDR3_1', 'GaiaDR3_2', 'GaiaDR3_3']
    assert files == [f'{name}.{outputstyle}' for name in expected]
    path = str(outdir / files[0])
    if outputstyle == 'txt':
        arr = np.loadtxt(path)
    else:
        t = Table.read(path)
        assert t.meta['SOURCEID'] == 1 and str(t['wave'].unit) == 'nm'
        arr = np.array([t['wave'], t['flux'], t['fluxerror']]).T
    assert arr.shape == (len(sampling), 3) and np.allclose(arr[:, 0], sampling)


def test_other_id_column(archive, tmp_path):
    fname = str(tmp_path / 'ids.csv')
    pd.DataFrame({'gaia_id': [1, 2], 'name': ['a', 'b']}).to_csv(fname, index=False)
    combdf = gb.batch(fname, idcolname='gaia_id', anonymous=True)
    assert combdf.flux.notna().all()


@pytest.mark.network
def test_live_archive(tmp_path):
    fname = str(tmp_path / 'ids.csv')
    pd.DataFrame({'source_id': [138877256416618112, 142522579844473600]}).to_csv(fname, index=False)
    combdf = gb.batch(fname, outputstyle='fits', anonymous=True)
    assert combdf.flux.notna().all()
    assert len(os.listdir(tmp_path / 'outputspectra')) == 2
