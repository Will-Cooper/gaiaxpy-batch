"""
Tests for gaiaxpy_batch. The archive (astroquery) and GaiaXPy are mocked out except for the tests marked network,
which are deselected by default (run them with pytest -m network).
"""
import os
from types import SimpleNamespace
from zipfile import BadZipFile

import numpy as np
import pandas as pd
import pytest
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from gaiaxpy import wl_to_pwl

import gaiaxpy_batch as gb

BAD = 999  # source ID that makes the archive reject the request (HTTP 400)
BADCAL = 777  # source ID whose downloaded spectrum fails to calibrate
NOXP = 888  # source IDs 800-899 have no XP spectrum
DEFAULT_WAVE = np.arange(336., 1021., 2.)  # gaiaxpy.calibrate default sampling


class FakeServer:
    """Stands in for the Gaia archive, counting logins and requests"""

    def __init__(self, maxids=None, transient=0, down_after=None, login_ok=True, expire=0):
        self.maxids = maxids  # request size limit
        self.transient = transient  # number of requests that fail with a server error before working
        self.down_after = down_after  # successful requests before the archive goes down
        self.login_ok = login_ok
        self.expire = expire  # number of requests that fail with HTTP 401 (expired session)
        self.logins = 0
        self.calls = []
        self.successes = 0

    def gaiaclass(self, **kwargs):
        return FakeGaia(self)

    def load(self, ids):
        self.calls.append(list(ids))
        if self.down_after is not None and self.successes >= self.down_after:
            raise OSError('Error 503:\nService Unavailable')
        if self.maxids is not None and len(ids) > self.maxids:
            raise OSError(f'Error 500:\nMaximum number of Source Ids reached (max: {self.maxids}, '
                          f'found: {len(ids)})\n\nCannot process request')
        if self.transient:
            self.transient -= 1
            raise BadZipFile('File is not a zip file') if self.transient % 2 else \
                OSError('Error 500:\nCould not open JPA EntityManager for transaction')
        if self.expire:
            self.expire -= 1
            raise OSError('Error 401:\nUnauthorized')
        if BAD in ids:
            raise OSError('Error 400:\nBad request')
        self.successes += 1
        good = [xid for xid in ids if not 800 <= xid < 900]
        if not good:
            return {}
        table = SimpleNamespace(to_pandas=lambda: pd.DataFrame({'source_id': good}))
        return {'XP_CONTINUOUS_DATAMODEL_GAIA.csv': [table]}


class FakeGaia:
    """Stands in for astroquery's GaiaClass"""

    def __init__(self, server):
        self.server = server
        self._TapPlus__isLoggedIn = False
        self._GaiaClass__gaiadata = SimpleNamespace(_TapPlus__isLoggedIn=False)

    def login(self, user=None, password=None):
        self.server.logins += 1
        if self.server.login_ok:  # astroquery does not raise on failure, it just stays logged out
            self._TapPlus__isLoggedIn = True
            self._GaiaClass__gaiadata._TapPlus__isLoggedIn = True

    def load_data(self, ids, **kwargs):
        return self.server.load(ids)


class FakeGaiaXPy:
    """Stands in for gaiaxpy.calibrate and gaiaxpy.convert on downloaded spectra"""

    def __init__(self, error=ValueError):
        self.calls = []
        self.error = error

    def _check(self, coeffs):
        self.calls.append(list(coeffs.source_id))
        if BADCAL in list(coeffs.source_id):
            raise self.error('bad coefficients')
        return list(coeffs.source_id)

    def calibrate(self, coeffs, sampling=None, **kwargs):
        ids = self._check(coeffs)
        sampling = DEFAULT_WAVE if sampling is None else sampling
        return pd.DataFrame({'source_id': ids, 'flux': [sampling * xid for xid in ids],
                             'flux_error': [np.ones_like(sampling) for _ in ids]}), sampling

    def convert(self, coeffs, sampling=None, **kwargs):
        ids = self._check(coeffs)
        rows = [{'source_id': xid, 'xp': xp, 'flux': sampling.copy(), 'flux_error': np.ones_like(sampling)}
                for xid in ids for xp in ('BP', 'RP')]
        return pd.DataFrame(rows), sampling


def install(monkeypatch, server=None, gaiaxpy=None):
    server = server or FakeServer()
    gaiaxpy = gaiaxpy or FakeGaiaXPy()
    monkeypatch.setattr(gb, 'GaiaClass', server.gaiaclass)
    monkeypatch.setattr(gb, 'calibrate', gaiaxpy.calibrate)
    monkeypatch.setattr(gb, 'convert', gaiaxpy.convert)
    monkeypatch.setattr(gb, 'RETRY_DELAY', 0.)
    return server, gaiaxpy


@pytest.fixture
def archive(monkeypatch):
    server, _ = install(monkeypatch)
    return server


def test_parse_sampling():
    assert np.allclose(gb.parse_sampling('linspace', '600', '1050', '120'), np.linspace(600, 1050, 120))
    assert np.allclose(gb.parse_sampling('arange', 400, 1000, 2), np.arange(400, 1000, 2))
    with pytest.raises(ValueError):
        gb.parse_sampling('__import__("os").system', 1, 2, 3)


def test_sysargs():
    args = gb.sysargs(['-f', 'a.csv', '-s', 'linspace', '600', '1000', '50', '-x', 'bp', '-u', '--retries', '5'])
    assert len(args.sampling) == 50 and args.xp == 'BP' and args.namecol is None
    assert args.chunksize is None and args.retries == 5 and args.retry_delay == gb.RETRY_DELAY


def test_get_sampling():
    wave, pwl = gb.get_sampling(np.linspace(340, 640, 30), uncalibrated=True, xp='BP')
    assert np.all(np.diff(pwl) > 0)  # GaiaXPy needs ascending pseudo-wavelength
    assert np.allclose(pwl, wl_to_pwl('BP', wave))
    assert gb.get_sampling(None) == (None, None)
    with pytest.raises(ValueError):
        gb.get_sampling(np.linspace(300, 1100, 10))
    with pytest.raises(ValueError):
        gb.get_sampling(np.linspace(1000, 600, 10))


def test_is_transient():
    assert gb.is_transient(BadZipFile('File is not a zip file'))
    assert gb.is_transient(OSError('Error 500:\nCould not open JPA EntityManager for transaction'))
    assert gb.is_transient(OSError('Error 503:\nService Unavailable'))
    assert gb.is_transient(TimeoutError('timed out'))
    assert not gb.is_transient(OSError('Error 400:\nBad request'))
    assert gb.http_status(OSError('Error 401:\nUnauthorized')) == 401


def test_bad_request_bisection(archive):
    ids = list(range(1, 26)) + [BAD]
    stackdf, wave, failed = gb.download_xp(ids, chunksize=10, workers=2)
    assert failed == {BAD: gb.DOWNLOAD_FAILED}
    assert sorted(stackdf.source_id) == list(range(1, 26))
    assert np.allclose(wave, DEFAULT_WAVE)
    assert all(len(f) == len(DEFAULT_WAVE) for f in stackdf.flux)
    assert [BAD] in archive.calls  # split down to the single bad ID


def test_no_xp_not_bisected(archive):
    ids = [1, 2, 801, 802, 802]  # duplicate dropped
    stackdf, wave, failed = gb.download_xp(ids, chunksize=2)
    assert failed == {801: gb.NO_XP, 802: gb.NO_XP} and len(archive.calls) == 2
    with pytest.raises(ValueError, match='No XP spectra found'):
        gb.download_xp([NOXP])


def test_calibration_failure_bisected_locally(monkeypatch):
    server, gaiaxpy = install(monkeypatch)
    stackdf, wave, failed = gb.download_xp([1, 2, BADCAL, 4], chunksize=10)
    assert failed == {BADCAL: gb.CALIBRATION_FAILED} and sorted(stackdf.source_id) == [1, 2, 4]
    assert len(server.calls) == 1  # one download, the splitting happens on the downloaded rows
    assert len(gaiaxpy.calls) > 1


def test_transient_failures_retried(monkeypatch):
    server, _ = install(monkeypatch, FakeServer(transient=2))
    stackdf, wave, failed = gb.download_xp([1, 2, 3], chunksize=10, retries=3)
    assert failed == {} and len(server.calls) == 3  # the same chunk three times, not split
    assert all(call == [1, 2, 3] for call in server.calls)


def test_retries_exhausted(monkeypatch):
    server, _ = install(monkeypatch, FakeServer(transient=2))
    stackdf, wave, failed = gb.download_xp([1, 2, 3, 4], chunksize=2, workers=1, retries=1)
    # the first chunk fails twice and is given up, the second works
    assert failed == {1: gb.DOWNLOAD_FAILED, 2: gb.DOWNLOAD_FAILED} and len(server.calls) == 3
    server, _ = install(monkeypatch, FakeServer(transient=4))
    with pytest.raises(gb.ArchiveError, match='No XP spectra downloaded'):
        gb.download_xp([1, 2, 3, 4], chunksize=2, workers=1, retries=1)
    assert len(server.calls) == 4


def test_archive_down_stops(monkeypatch):
    server, _ = install(monkeypatch, FakeServer(down_after=2))
    stackdf, wave, failed = gb.download_xp(list(range(1, 11)), chunksize=1, workers=1, retries=1)
    assert sorted(stackdf.source_id) == [1, 2]
    reasons = list(failed.values())
    assert reasons.count(gb.DOWNLOAD_FAILED) == gb.ABORT_AFTER
    assert reasons.count(gb.NOT_ATTEMPTED) == 10 - 2 - gb.ABORT_AFTER
    assert len(server.calls) == 2 + gb.ABORT_AFTER * 2  # no requests after stopping
    server, _ = install(monkeypatch, FakeServer(down_after=0))
    with pytest.raises(gb.ArchiveUnavailable):
        gb.download_xp(list(range(1, 11)), chunksize=1, workers=1, retries=0)


def test_login_once_per_worker(monkeypatch):
    server, _ = install(monkeypatch)
    stackdf, wave, failed = gb.download_xp(list(range(1, 13)), chunksize=2, workers=2, username='me', password='pw')
    assert failed == {} and len(server.calls) == 6 and 1 <= server.logins <= 2


def test_login_failure_stops(monkeypatch):
    server, _ = install(monkeypatch, FakeServer(login_ok=False))
    with pytest.raises(gb.ArchiveUnavailable, match='Could not log in'):
        gb.download_xp(list(range(1, 7)), chunksize=2, workers=1, username='me', password='pw', retries=2)
    assert server.logins == 3 and server.calls == []


def test_expired_session_logs_in_again(monkeypatch):
    server, _ = install(monkeypatch, FakeServer(expire=1))
    stackdf, wave, failed = gb.download_xp([1, 2], username='me', password='pw')
    assert failed == {} and server.logins == 2


def test_default_chunksize_and_workers(monkeypatch):
    server, _ = install(monkeypatch)
    gb.download_xp(list(range(1000, 2000)))
    assert max(len(call) for call in server.calls) == gb.ANONYMOUS_CHUNKSIZE
    server, _ = install(monkeypatch)
    gb.download_xp(list(range(1000, 11000)), username='me', password='pw', workers=10)
    assert max(len(call) for call in server.calls) == gb.LOGGED_IN_CHUNKSIZE
    assert server.logins <= gb.LOGGED_IN_MAX_WORKERS
    server, _ = install(monkeypatch)  # fewer sources than workers times the chunk size: spread over the workers
    gb.download_xp(list(range(1000, 2000)), username='me', password='pw', workers=4)
    assert [len(call) for call in server.calls] == [250] * 4
    server, _ = install(monkeypatch)  # unless the chunk size is given
    gb.download_xp(list(range(1000, 2000)), username='me', password='pw', workers=4, chunksize=1000)
    assert [len(call) for call in server.calls] == [1000]


def test_interrupt_keeps_finished(monkeypatch):
    install(monkeypatch)
    real = gb.as_completed

    def interrupted(futures):  # Ctrl-C after the first chunk is collected
        iterator = real(futures)
        yield next(iterator)
        raise KeyboardInterrupt

    monkeypatch.setattr(gb, 'as_completed', interrupted)
    stackdf, wave, failed = gb.download_xp(list(range(1, 21)), chunksize=1, workers=1)
    assert len(stackdf) >= 1 and len(stackdf) + len(failed) == 20
    assert set(failed.values()) == {gb.NOT_ATTEMPTED}


def test_archive_limit(monkeypatch):
    server, _ = install(monkeypatch, FakeServer(maxids=3))
    ids = list(range(1, 21)) + [BAD]
    stackdf, wave, failed = gb.download_xp(ids, chunksize=10, workers=1)
    assert failed == {BAD: gb.DOWNLOAD_FAILED} and sorted(stackdf.source_id) == list(range(1, 21))
    assert all(len(call) <= 3 for call in server.calls[1:])  # only the first request is too large
    assert len(server.calls) <= 1 + 7 + 2  # one rejection, then chunks of 3, plus the bad chunk


def test_fatal_error_not_retried(monkeypatch):
    server, gaiaxpy = install(monkeypatch, gaiaxpy=FakeGaiaXPy(error=TypeError))
    with pytest.raises(TypeError):
        gb.download_xp([1, 2, BADCAL, 4], chunksize=10)
    assert len(gaiaxpy.calls) == 1


def test_fits_units(archive, tmp_path):
    for uncalibrated, unit in ((False, 'W m-2 nm-1'), (True, 'electron/s')):
        fname = str(tmp_path / f'ids{uncalibrated}.csv')
        pd.DataFrame({'source_id': [1]}).to_csv(fname, index=False)
        gb.batch(fname, outputstyle='fits', anonymous=True, uncalibrated=uncalibrated,
                 sampling=np.linspace(650, 1000, 20))
        header = fits.getheader(str(tmp_path / 'outputspectra' / 'GaiaDR3_1.fits'), 1)
        assert header['TUNIT1'] == 'nm' and header['SOURCEID'] == 1
        if uncalibrated:
            assert header['TUNIT2'] == header['TUNIT3'] == unit
        else:
            assert u.Unit(header['TUNIT2'], format='fits') == u.Unit(unit, format='fits')


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
    assert len(combdf) == 4 and combdf.flux.isna().sum() == 1 and combdf.attrs['incomplete'] == 0
    cut = pd.read_csv(fname.replace('.csv', '_XP_spectra_cut.csv'))
    assert len(cut) == 3
    assert np.allclose(gb.str2float(cut.flux.iloc[0]), sampling * cut.source_id.iloc[0])
    failedcsv = pd.read_csv(fname.replace('.csv', '_XP_spectra_failed.csv'))
    assert failedcsv.source_id.tolist() == [NOXP] and failedcsv.reason.tolist() == [gb.NO_XP]
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


def test_batch_incomplete(monkeypatch, tmp_path):
    install(monkeypatch, FakeServer(down_after=1))
    fname = str(tmp_path / 'ids.csv')
    pd.DataFrame({'source_id': list(range(1, 9))}).to_csv(fname, index=False)
    combdf = gb.batch(fname, anonymous=True, chunksize=1, workers=1, retries=0)
    failedcsv = pd.read_csv(fname.replace('.csv', '_XP_spectra_failed.csv'))
    assert combdf.attrs['incomplete'] == 7 == len(failedcsv)
    install(monkeypatch)  # the archive is back
    combdf = gb.batch(fname.replace('.csv', '_XP_spectra_failed.csv'), anonymous=True)  # the failed csv is an input
    assert combdf.flux.notna().sum() == 7 and combdf.attrs['incomplete'] == 0


def test_other_id_column(archive, tmp_path):
    fname = str(tmp_path / 'ids.csv')
    pd.DataFrame({'gaia_id': [1, 2], 'name': ['a', 'b']}).to_csv(fname, index=False)
    combdf = gb.batch(fname, idcolname='gaia_id', anonymous=True)
    assert combdf.flux.notna().all()


@pytest.mark.network
def test_live_archive(tmp_path):
    fname = str(tmp_path / 'ids.csv')
    pd.DataFrame({'source_id': [138877256416618112, 142522579844473600, 923972363019345536]}).to_csv(fname,
                                                                                                     index=False)
    combdf = gb.batch(fname, outputstyle='fits', anonymous=True)
    assert combdf.flux.notna().sum() == 2 and combdf.attrs['incomplete'] == 0
    assert len(os.listdir(tmp_path / 'outputspectra')) == 2
