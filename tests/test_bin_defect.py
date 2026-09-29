"""End-to-end tests of steps 11-16 on a synthetic dataset with defects.

test_bin_smoke.py runs the same chain on flawless data, so every rejection
path stays unexecuted there: no ifg is ever discarded, no pixel is ever
nodata, and no epoch ever drops out. This module runs the chain on
tests/synth.py's defective dataset instead (see build_geocml_defect for the
design of the network and of each defect) and checks that

  - step 11 discards the low coverage and the low coherence ifgs,
  - step 12 pins the 2pi unwrapping error on the single ifg carrying it,
  - the epoch whose only ifg was discarded drops out of the time series,
  - nodata stays nan all the way through the noise indices, the mask and
    the spatio-temporal filter, leaking in neither direction, and
  - the inverted velocity still matches the truth once the bad ifgs are
    gone, which it would not if the rejection had missed any of them.
"""
import re

import h5py
import numpy as np
import pytest

import LiCSBAS_io_lib as io_lib

pytestmark = pytest.mark.smoke


#%% Step 11: reject by coverage and by coherence
def test_step11_rejects_defective_ifgs(ts11d, geocml_defect):
    bad = (ts11d / 'info' / '11bad_ifg.txt').read_text().split()
    assert sorted(bad) == geocml_defect.bad_ifg11


def test_step11_rejection_reasons(ts11d, geocml_defect):
    """The stats file must show which threshold each bad ifg failed."""
    stats = {}
    for line in (ts11d / 'info' / '11ifg_stats.txt').read_text().splitlines():
        if re.match(r'^2\d{7}_2\d{7}', line):
            f = line.split()
            stats[f[0]] = {'cov': float(f[3]), 'coh': float(f[4])}

    assert len(stats) == len(geocml_defect.ifgdates)

    from synth import LOW_COV_IFG, LOW_COH_IFG, ISOLATED_IFG
    # Low coverage ifgs fail on coverage alone, not on coherence
    for ifgd in (LOW_COV_IFG, ISOLATED_IFG):
        assert stats[ifgd]['cov'] < 0.3
        assert stats[ifgd]['coh'] >= 0.05
    # and the low coherence ifg the other way round
    assert stats[LOW_COH_IFG]['coh'] < 0.05
    assert stats[LOW_COH_IFG]['cov'] >= 0.3

    # Every other ifg must survive both thresholds
    good = set(stats) - set(geocml_defect.bad_ifg11)
    for ifgd in good:
        assert stats[ifgd]['cov'] >= 0.3 and stats[ifgd]['coh'] >= 0.05


#%% Step 12: pin the unwrapping error on the ifg that carries it
def test_step12_identifies_loop_error_ifg(ts12d, geocml_defect):
    bad = (ts12d / 'info' / '12bad_ifg.txt').read_text().split()
    assert bad == geocml_defect.bad_ifg12


def test_step12_network_stays_connected(ts12d, geocml_defect):
    """Removing the bad ifg must leave a single connected network."""
    gap_info = (ts12d / 'info' / '12network_gap_info.txt').read_text()
    assert 'Gaps' not in gap_info

    # One "yyyymmdd-yyyymmdd (years, n_image)" line per connected component
    comps = re.findall(r'^(\d{8})-(\d{8}) \([\d.]+, (\d+)\)$',
                       gap_info, re.MULTILINE)
    assert len(comps) == 1
    first, last, n_im = comps[0]
    kept = geocml_defect.imdates_kept
    assert (first, last, int(n_im)) == (kept[0], kept[-1], len(kept))

    # No image is reported as removed: the isolated epoch is already gone
    # after step 11, so step 12 never sees it.
    assert (ts12d / 'info' / '12removed_image.txt').read_text().strip() == ''


def test_step12_n_loop_err_zero_after_removal(ts12d, geocml_defect):
    """No loop error may remain once the bad ifg is discarded."""
    n_loop_err = io_lib.read_img(str(ts12d / 'results' / 'n_loop_err'),
                                 geocml_defect.length, geocml_defect.width)
    np.testing.assert_array_equal(n_loop_err, 0)


@pytest.mark.xfail(strict=True, reason='calc_n_unw in LiCSBAS12_loop_closure '
                                       'ignores its _ifgdates argument and '
                                       'loops over the global ifgdates, so '
                                       'results/n_unw counts the ifgs that '
                                       'step 12 discarded')
def test_step12_n_unw_counts_good_ifgs_only(ts12d, geocml_defect):
    n_unw = io_lib.read_img(str(ts12d / 'results' / 'n_unw'),
                            geocml_defect.length, geocml_defect.width)
    n_good = (len(geocml_defect.ifgdates) - len(geocml_defect.bad_ifg11)
              - len(geocml_defect.bad_ifg12))
    land = geocml_defect.land
    # The block where the low coverage ifgs do have data is excluded: there
    # the count is higher by construction.
    from synth import COV_BLOCK
    plain = land.copy()
    plain[:COV_BLOCK, :COV_BLOCK] = False
    np.testing.assert_array_equal(n_unw[plain], n_good)


def test_step12_nodata_stays_nodata(ts12d, geocml_defect):
    """Pixels with no valid unw in any ifg must not be counted."""
    n_unw = io_lib.read_img(str(ts12d / 'results' / 'n_unw'),
                            geocml_defect.length, geocml_defect.width)
    np.testing.assert_array_equal(n_unw[~geocml_defect.land], 0)

    coh_avg = io_lib.read_img(str(ts12d / 'results' / 'coh_avg'),
                              geocml_defect.length, geocml_defect.width)
    assert np.all(np.isnan(coh_avg[~geocml_defect.land]))


#%% Step 13: the isolated epoch drops out, nodata stays nan, truth holds
def test_step13_drops_isolated_epoch(ts13d, geocml_defect):
    """The epoch whose only ifg was rejected must leave the time series."""
    with h5py.File(str(ts13d / 'cum.h5'), 'r') as f:
        imdates = [str(d) for d in f['imdates'][()]]
        cum = f['cum'][()]

    assert imdates == geocml_defect.imdates_kept
    assert geocml_defect.isolated_imd not in imdates
    assert cum.shape == (len(geocml_defect.imdates_kept),
                         geocml_defect.length, geocml_defect.width)


def test_step13_nodata_is_nan(ts13d, geocml_defect):
    with h5py.File(str(ts13d / 'cum.h5'), 'r') as f:
        cum = f['cum'][()]
        vel = f['vel'][()]

    sea = ~geocml_defect.land
    assert np.all(np.isnan(vel[sea]))
    assert np.all(np.isnan(cum[:, sea]))
    assert np.all(np.isfinite(vel[geocml_defect.land]))


def test_step13_velocity_matches_truth_after_rejection(ts13d, geocml_defect):
    """The 2pi error must be gone from the solution.

    Had step 12 kept the corrupted ifg, the velocity over the columns it
    damages would be wrong by far more than this tolerance.
    """
    with h5py.File(str(ts13d / 'cum.h5'), 'r') as f:
        vel = f['vel'][()]
        refarea = f['refarea'][()]

    refarea = refarea.decode() if isinstance(refarea, bytes) else str(refarea)
    x1, x2, y1, y2 = [int(s) for s in re.split('[:/]', refarea)]

    land = geocml_defect.land
    truth = geocml_defect.vel_mm
    vel_ref = np.nanmean(vel[y1:y2, x1:x2])
    truth_ref = np.nanmean(truth[y1:y2, x1:x2])

    np.testing.assert_allclose((vel - vel_ref)[land],
                               (truth - truth_ref)[land], atol=0.1)


def test_step13_no_gap_reported(ts13d, geocml_defect):
    with h5py.File(str(ts13d / 'cum.h5'), 'r') as f:
        n_gap = f['n_gap'][()]
        n_ifg_noloop = f['n_ifg_noloop'][()]

    land = geocml_defect.land
    np.testing.assert_array_equal(n_gap[land], 0)
    np.testing.assert_array_equal(n_ifg_noloop[land], 0)


#%% Steps 14-16: nodata must survive the noise indices, the mask and the
#   spatio-temporal filter without leaking either way
def test_step14_noise_indices_keep_nodata(ts14d, geocml_defect):
    land, sea = geocml_defect.land, ~geocml_defect.land
    for name in ('vstd', 'stc'):
        a = io_lib.read_img(str(ts14d / 'results' / name),
                            geocml_defect.length, geocml_defect.width)
        assert np.all(np.isnan(a[sea])), '{} leaked into nodata'.format(name)
        assert np.all(np.isfinite(a[land])), '{} nan on valid data'.format(name)
        assert np.all(a[land] < 1)  # near-exact linear data -> tiny


def test_step15_mask_keeps_nodata_nan(ts15d, geocml_defect):
    land, sea = geocml_defect.land, ~geocml_defect.land

    mask = io_lib.read_img(str(ts15d / 'results' / 'mask'),
                           geocml_defect.length, geocml_defect.width)
    assert np.all(np.isnan(mask[sea]))
    assert set(np.unique(mask[land])) <= {0.0, 1.0}
    assert np.any(mask[land] == 1)

    vel_mskd = io_lib.read_img(str(ts15d / 'results' / 'vel.mskd'),
                               geocml_defect.length, geocml_defect.width)
    assert np.all(np.isnan(vel_mskd[sea]))
    assert np.any(np.isfinite(vel_mskd[land]))


def test_step16_filter_does_not_leak_across_nodata(ts16d, geocml_defect):
    """The gaussian spatial filter must neither spread into the nodata area
    nor let nodata poison the valid pixels next to it."""
    land, sea = geocml_defect.land, ~geocml_defect.land

    with h5py.File(str(ts16d / 'cum_filt.h5'), 'r') as f:
        cum_filt = f['cum'][()]
        imdates = [str(d) for d in f['imdates'][()]]

    assert imdates == geocml_defect.imdates_kept
    assert np.all(np.isnan(cum_filt[:, sea]))
    assert np.all(np.isfinite(cum_filt[:, land]))

    vel_filt = io_lib.read_img(str(ts16d / 'results' / 'vel.filt'),
                               geocml_defect.length, geocml_defect.width)
    assert np.all(np.isnan(vel_filt[sea]))
    assert np.all(np.isfinite(vel_filt[land]))


def test_step16_velocity_still_matches_truth(ts16d, geocml_defect):
    """End of the chain: the 2pi error must not have survived into vel.filt.

    The tolerance is looser than in step 13 because the spatial filter
    smooths the linear ramp near the edges of the valid area. It is still
    an order of magnitude tighter than the error the rejected ifg would
    introduce (its 2pi jump is ~28 mm over 24 days).
    """
    vel_filt = io_lib.read_img(str(ts16d / 'results' / 'vel.filt'),
                               geocml_defect.length, geocml_defect.width)
    ref = (ts16d / 'info' / '16ref.txt').read_text().strip()
    x1, x2, y1, y2 = [int(s) for s in re.split('[:/]', ref)]

    land = geocml_defect.land
    truth = geocml_defect.vel_mm
    d = ((vel_filt - np.nanmean(vel_filt[y1:y2, x1:x2]))
         - (truth - np.nanmean(truth[y1:y2, x1:x2])))
    assert np.nanmax(np.abs(d[land])) < 0.5
