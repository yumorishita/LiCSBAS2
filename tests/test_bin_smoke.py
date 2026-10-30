"""End-to-end smoke tests of the LiCSBAS time series analysis chain
(steps 11-16) on a small synthetic dataset.

The step fixtures in conftest.py run each script once per session via
subprocess; the tests here check their outputs. All data is synthetic
(tests/synth.py): 10x10 pixels, 5 epochs, 7 exactly loop-consistent
interferograms of a known linear velocity field, so the inverted
velocity can be compared against the truth.
"""
import re
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest

import LiCSBAS_io_lib as io_lib
import synth

pytestmark = pytest.mark.smoke


#%% Step 11: check unw
def test_step11_outputs(ts11, geocml):
    infodir = ts11 / 'info'
    assert (infodir / '11bad_ifg.txt').read_text() == ''
    stats = (infodir / '11ifg_stats.txt').read_text().splitlines()
    data_lines = [l for l in stats if re.match(r'^2\d{7}_2\d{7}', l)]
    assert len(data_lines) == len(geocml.ifgdates)
    assert (infodir / 'slc.mli.par').exists()
    assert (infodir / 'EQA.dem_par').exists()
    network_pngs = list((ts11 / 'network').glob('network11*.png'))
    assert network_pngs and all(p.stat().st_size > 0 for p in network_pngs)


#%% Step 12: loop closure
def test_step12_outputs(ts12, geocml):
    infodir = ts12 / 'info'
    assert (infodir / '12bad_ifg.txt').read_text() == ''

    ref = (infodir / '12ref.txt').read_text().strip()
    x1, x2, y1, y2 = [int(s) for s in re.split('[:/]', ref)]
    assert 0 <= x1 < x2 <= geocml.width
    assert 0 <= y1 < y2 <= geocml.length

    resultsdir = ts12 / 'results'
    n_unw = io_lib.read_img(str(resultsdir / 'n_unw'),
                            geocml.length, geocml.width)
    np.testing.assert_array_equal(n_unw, len(geocml.ifgdates))

    coh_avg = io_lib.read_img(str(resultsdir / 'coh_avg'),
                              geocml.length, geocml.width)
    np.testing.assert_allclose(coh_avg, 180 / 255, atol=0.01)

    n_loop_err = io_lib.read_img(str(resultsdir / 'n_loop_err'),
                                 geocml.length, geocml.width)
    np.testing.assert_array_equal(n_loop_err, 0)


#%% Step 13: small baseline inversion
def test_step13_cum_h5(ts13, geocml):
    cumfile = ts13 / 'cum.h5'
    assert cumfile.exists()
    with h5py.File(str(cumfile), 'r') as f:
        cum = f['cum'][()]
        vel = f['vel'][()]
        imdates = [str(d) for d in f['imdates'][()]]
        refarea = f['refarea'][()]

    n_im = len(geocml.imdates)
    assert cum.shape == (n_im, geocml.length, geocml.width)
    assert vel.shape == (geocml.length, geocml.width)
    assert imdates == geocml.imdates
    assert np.all(np.isfinite(vel))

    # Velocity must match the synthetic truth after subtracting the
    # reference area (mm/yr)
    refarea = refarea.decode() if isinstance(refarea, bytes) else str(refarea)
    x1, x2, y1, y2 = [int(s) for s in re.split('[:/]', refarea)]
    vel_ref = np.nanmean(vel[y1:y2, x1:x2])
    truth_ref = np.nanmean(geocml.vel_mm[y1:y2, x1:x2])
    np.testing.assert_allclose(vel - vel_ref, geocml.vel_mm - truth_ref,
                               atol=0.1)

    params = (ts13 / 'info' / '13parameters.txt').read_text()
    assert 'wavelength' in params
    assert 'pixel_spacing_r' in params


def test_step13_ref(ts13):
    ref = (ts13 / 'info' / '13ref.txt').read_text().strip()
    assert re.match(r'^\d+:\d+/\d+:\d+$', ref)


#%% Step 14: velocity std (vstd is calculated only with --vstd, which is
#   tested in test_bin_defect.py)
def test_step14_outputs(ts14, geocml):
    resultsdir = ts14 / 'results'
    assert not (resultsdir / 'vstd').exists()
    assert not (resultsdir / 'vstd.png').exists()

    stc = io_lib.read_img(str(resultsdir / 'stc'),
                          geocml.length, geocml.width)
    assert np.all(np.isfinite(stc))
    assert np.all(stc < 1)  # exact linear data -> tiny stc


#%% Step 15: mask
def test_step15_outputs(ts15, geocml):
    resultsdir = ts15 / 'results'
    mask = io_lib.read_img(str(resultsdir / 'mask'),
                           geocml.length, geocml.width)
    assert set(np.unique(mask)) <= {0.0, 1.0}
    assert np.any(mask == 1)

    vel_mskd = io_lib.read_img(str(resultsdir / 'vel.mskd'),
                               geocml.length, geocml.width)
    assert np.any(np.isfinite(vel_mskd))

    params = (ts15 / 'info' / '15parameters.txt').read_text()
    assert re.search(r'^- vstd +: not used', params, re.M)
    # vstd keeps its place in the list, as in the panels of mask_ts.png
    names = re.findall(r'^- (\S+)', params, re.M)
    assert names[:4] == ['coh_avg', 'n_unw', 'vstd', 'maxTlen']


def test_step15_ignores_vstd_thre_without_vstd(ts15, geocml, run_script):
    mask_file = str(ts15 / 'results' / 'mask')
    mask = io_lib.read_img(mask_file, geocml.length, geocml.width)
    res = run_script('LiCSBAS15_mask_ts.py', '-t', 'TS_GEOCml1', '-v', '100',
                     cwd=geocml.workdir)
    assert '-v is ignored' in res.stderr
    params = (ts15 / 'info' / '15parameters.txt').read_text()
    assert re.search(r'^- vstd +: not used', params, re.M)
    np.testing.assert_array_equal(
        io_lib.read_img(mask_file, geocml.length, geocml.width), mask)


#%% Step 14 re-run (after the chain, as it rewrites the results dir)
def test_step14_removes_old_vstd(ts16, geocml, run_script):
    """Without --vstd, a vstd left by a previous run must be removed
    because step15 would use it otherwise."""
    resultsdir = ts16 / 'results'
    np.full((geocml.length, geocml.width), 1000,
            dtype=np.float32).tofile(str(resultsdir / 'vstd'))
    (resultsdir / 'vstd.png').write_bytes(b'')
    run_script('LiCSBAS14_vel_std.py', '-t', 'TS_GEOCml1', cwd=geocml.workdir)
    assert not (resultsdir / 'vstd').exists()
    assert not (resultsdir / 'vstd.png').exists()


#%% Step 16: filter
def test_step16_outputs(ts16, geocml):
    cumffile = ts16 / 'cum_filt.h5'
    assert cumffile.exists()
    with h5py.File(str(cumffile), 'r') as f:
        cum_filt = f['cum'][()]
        # Default is the linear fit in time, n_lsq_t 9 cut to 5 epochs
        assert f['n_lsq_t'][()] == len(geocml.imdates)
        assert np.isnan(f['filtwidth_yr'][()])
    n_im = len(geocml.imdates)
    assert cum_filt.shape == (n_im, geocml.length, geocml.width)

    vel_filt = io_lib.read_img(str(ts16 / 'results' / 'vel.filt'),
                               geocml.length, geocml.width)
    assert np.all(np.isfinite(vel_filt))
    assert (ts16 / 'info' / '16parameters.txt').exists()


def nonlinear_signal(geocml):
    """Signal g (mm) nonlinear in any 3 epochs, of amplitude w (1-2) in
    each pixel, to be added to the linear synthetic data."""
    g = np.array([0., 4., -3., 5., 1.])
    assert len(g) == len(geocml.imdates)
    x, y = np.meshgrid(np.arange(geocml.width), np.arange(geocml.length))
    w = 1 + (x / (geocml.width - 1) + y / (geocml.length - 1)) / 2
    return g, w


def hp_in_time(g, n):
    """HP in time of g with --n_lsq_t n: the residual from the linear fit
    to n epochs around each epoch, shifted inward near both ends."""
    t = synth.dt_cum_years()
    hp = np.zeros_like(g)
    for i in range(len(t)):
        i1 = min(max(i - n // 2, 0), len(t) - n)
        fit = np.polyfit(t[i1:i1 + n], g[i1:i1 + n], 1)
        hp[i] = g[i] - np.polyval(fit, t[i])
    return hp


def run16_copy(ts16, geocml, run_script, tsname, *args, add_cum=None):
    """Run step 16 on a copy of the TS dir so that ts16 is not overwritten.
    add_cum (n_im x length x width) is added to cum.h5 of the copy."""
    tsdir = geocml.workdir / tsname
    shutil.copytree(str(ts16), str(tsdir))
    if add_cum is not None:
        with h5py.File(str(tsdir / 'cum.h5'), 'r+') as f:
            f['cum'][...] = f['cum'][()] + add_cum.astype(np.float32)
    run_script('LiCSBAS16_filt_ts.py', '-t', tsname, '--n_para', '1', *args,
               cwd=geocml.workdir)
    return tsdir


def test_step16_lsq_hp_in_time(ts16, geocml, run_script):
    """LP in time is the value at each epoch of the linear fit to n_lsq_t
    epochs around it, shifted inward near both ends. Add a nonlinear signal
    to the linear synthetic data and make the spatial LP identity (tiny -s):
    the filter must remove from each pixel exactly the residual of the signal
    from that fit (HP in time), and leave the linear part unchanged."""
    n = 3
    g, w = nonlinear_signal(geocml)
    hp = hp_in_time(g, n)
    assert np.all(np.abs(hp) > 1)  # every epoch is filtered

    with h5py.File(str(ts16 / 'cum.h5'), 'r') as f:
        cum = f['cum'][()]
    tsdir = run16_copy(ts16, geocml, run_script, 'TS_GEOCml1_lsq3',
                       '--n_lsq_t', str(n), '-s', '1e-5', '--nomask',
                       add_cum=g[:, None, None] * w)

    with h5py.File(str(tsdir / 'cum_filt.h5'), 'r') as f:
        cum_filt = f['cum'][()]
        assert f['n_lsq_t'][()] == n

    # Compare up to a constant in each epoch, which depends on the ref
    # point. cum of step 13 is linear to ~2e-3 mm (float32).
    expected = cum + (g - hp)[:, None, None] * w
    np.testing.assert_allclose(
        cum_filt - cum_filt.mean(axis=(1, 2), keepdims=True),
        expected - expected.mean(axis=(1, 2), keepdims=True), atol=3e-3)

    params = (tsdir / 'info' / '16parameters.txt').read_text()
    assert re.search(r'^n_lsq_t: +3$', params, re.M)
    assert 'filtwidth_yr' not in params


def test_step16_gauss_auto(ts16, geocml, run_script):
    """-y auto uses the gaussian filter in time of avg_interval*3."""
    tsdir = run16_copy(ts16, geocml, run_script, 'TS_GEOCml1_gauss',
                       '-y', 'auto')

    dt_yr = np.array(synth.dt_cum_years())
    filtwidth_yr = dt_yr[-1] / (len(dt_yr) - 1) * 3
    with h5py.File(str(tsdir / 'cum_filt.h5'), 'r') as f:
        assert 'n_lsq_t' not in f
        np.testing.assert_allclose(f['filtwidth_yr'][()], filtwidth_yr,
                                   rtol=1e-6)

    params = (tsdir / 'info' / '16parameters.txt').read_text()
    assert re.search(r'^filtwidth_yr: +{}$'.format(
        re.escape('{:.3f}'.format(filtwidth_yr))), params, re.M)
    assert 'n_lsq_t' not in params


def test_step16_n_lsq_t_and_y(ts15, geocml, bin_env, repo_root):
    """-y and --n_lsq_t together is an error before processing."""
    res = subprocess.run(
        [sys.executable, str(repo_root / 'bin' / 'LiCSBAS16_filt_ts.py'),
         '-t', 'TS_GEOCml1', '-y', '0.1', '--n_lsq_t', '9'],
        env=bin_env, cwd=str(geocml.workdir), capture_output=True, text=True,
        timeout=300)
    assert res.returncode == 2
    assert 'Both -y and --n_lsq_t' in res.stderr


@pytest.fixture(scope='module')
def ts16_nl(ts16, geocml, run_script):
    """Step 16 with --n_lsq_t 5 on the data plus nonlinear_signal.

    HP in time of the signal is >0.5 mm in every epoch and pixel (w >= 1),
    so the filter phase is far from 0 also at the ref point; the tests below
    rely on it to catch re-referencing with the unfiltered cum.
    """
    g, w = nonlinear_signal(geocml)
    assert np.all(np.abs(hp_in_time(g, 5)) > 0.5)
    return run16_copy(ts16, geocml, run_script, 'TS_GEOCml1_nl',
                      '--n_lsq_t', '5', add_cum=g[:, None, None] * w)


def read_ref16(tsdir):
    ref = (tsdir / 'info' / '16ref.txt').read_text().strip()
    return [int(s) for s in re.split('[:/]', ref)]


def test_step16_zero_at_ref(ts16_nl, geocml):
    """The filtered time series, velocity and intercept are referenced to the
    stable ref point of step 16 (16ref.txt), so they are 0 there."""
    x1, x2, y1, y2 = read_ref16(ts16_nl)
    with h5py.File(str(ts16_nl / 'cum_filt.h5'), 'r') as f:
        cum_filt = f['cum'][()]
    np.testing.assert_array_equal(cum_filt[:, y1:y2, x1:x2], 0)

    for name in ('vel.filt', 'vintercept.filt'):
        data = io_lib.read_img(str(ts16_nl / 'results' / name),
                               geocml.length, geocml.width)
        np.testing.assert_allclose(data[y1:y2, x1:x2], 0, atol=1e-6,
                                   err_msg=name)


def test_step16_independent_of_cum_ref(ts16, ts16_nl, geocml, run_script):
    """The results must not depend on the reference of cum.h5 from step 13.
    Re-reference cum.h5 to another pixel (an offset in each epoch, nonlinear
    in time) and run step 16 again: the same ref point is selected and the
    results are the same."""
    g, w = nonlinear_signal(geocml)
    with h5py.File(str(ts16 / 'cum.h5'), 'r') as f:
        cum_nl = f['cum'][()] + g[:, None, None] * w  # input of ts16_nl
    py, px = geocml.length - 1, 0
    tsdir = run16_copy(ts16, geocml, run_script, 'TS_GEOCml1_nl_reref',
                       '--n_lsq_t', '5', add_cum=g[:, None, None] * w
                       - cum_nl[:, py:py+1, px:px+1])

    assert read_ref16(tsdir) == read_ref16(ts16_nl)
    with h5py.File(str(ts16_nl / 'cum_filt.h5'), 'r') as f:
        cum_filt = f['cum'][()]
    with h5py.File(str(tsdir / 'cum_filt.h5'), 'r') as f:
        cum_filt_reref = f['cum'][()]
    np.testing.assert_allclose(cum_filt_reref, cum_filt, atol=1e-3)


#%% Downstream tool: cum2vel
def test_cum2vel(ts13, geocml, run_script):
    run_script('LiCSBAS_cum2vel.py', '-i', 'TS_GEOCml1/cum.h5',
               '-o', 'vel.cum2vel', cwd=geocml.workdir)
    with h5py.File(str(ts13 / 'cum.h5'), 'r') as f:
        vel_h5 = f['vel'][()]
    vel_out = io_lib.read_img(str(geocml.workdir / 'vel.cum2vel'),
                              geocml.length, geocml.width)
    np.testing.assert_allclose(vel_out, vel_h5, atol=0.01)
