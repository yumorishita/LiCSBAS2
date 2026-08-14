"""End-to-end smoke tests of the LiCSBAS time series analysis chain
(steps 11-16) on a small synthetic dataset.

The step fixtures in conftest.py run each script once per session via
subprocess; the tests here check their outputs. All data is synthetic
(tests/synth.py): 10x10 pixels, 5 epochs, 7 exactly loop-consistent
interferograms of a known linear velocity field, so the inverted
velocity can be compared against the truth.
"""
import re

import h5py
import numpy as np
import pytest

import LiCSBAS_io_lib as io_lib

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


#%% Step 14: velocity std
def test_step14_outputs(ts14, geocml):
    resultsdir = ts14 / 'results'
    vstd = io_lib.read_img(str(resultsdir / 'vstd'),
                           geocml.length, geocml.width)
    assert np.all(np.isfinite(vstd))
    assert np.all(vstd < 1)  # exact linear data -> tiny bootstrap std

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


#%% Step 16: filter
def test_step16_outputs(ts16, geocml):
    cumffile = ts16 / 'cum_filt.h5'
    assert cumffile.exists()
    with h5py.File(str(cumffile), 'r') as f:
        cum_filt = f['cum'][()]
    n_im = len(geocml.imdates)
    assert cum_filt.shape == (n_im, geocml.length, geocml.width)

    vel_filt = io_lib.read_img(str(ts16 / 'results' / 'vel.filt'),
                               geocml.length, geocml.width)
    assert np.all(np.isfinite(vel_filt))
    assert (ts16 / 'info' / '16parameters.txt').exists()


#%% Downstream tool: cum2vel
def test_cum2vel(ts13, geocml, run_script):
    run_script('LiCSBAS_cum2vel.py', '-i', 'TS_GEOCml1/cum.h5',
               '-o', 'vel.cum2vel', cwd=geocml.workdir)
    with h5py.File(str(ts13 / 'cum.h5'), 'r') as f:
        vel_h5 = f['vel'][()]
    vel_out = io_lib.read_img(str(geocml.workdir / 'vel.cum2vel'),
                              geocml.length, geocml.width)
    np.testing.assert_allclose(vel_out, vel_h5, atol=0.01)
