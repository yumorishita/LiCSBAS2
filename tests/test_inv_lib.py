"""Unit tests for LiCSBAS_inv_lib (SBAS design matrices and inversion)."""
import numpy as np

import LiCSBAS_inv_lib as inv_lib
import synth


IFGDATES = synth.IFGDATES  # 5 epochs, 7 ifgs (4 consecutive + 3 skip-1)
N_IM = len(synth.IMDATES)
N_IFG = len(IFGDATES)
DT_CUM = synth.dt_cum_years()


def forward_model_unw(vel, n_pt):
    """unw (n_pt, n_ifg) for point velocities vel (n_pt,) [unit/yr]."""
    G = inv_lib.make_sb_matrix(IFGDATES)
    inc_true = vel[:, None] * np.diff(DT_CUM)[None, :]  # (n_pt, n_im-1)
    return inc_true @ G.T.astype(np.float64), inc_true


#%% Design matrices
def test_make_sb_matrix():
    G = inv_lib.make_sb_matrix(IFGDATES)
    assert G.shape == (N_IFG, N_IM - 1)
    for i, ifgd in enumerate(IFGDATES):
        ix_p = synth.IMDATES.index(ifgd[:8])
        ix_s = synth.IMDATES.index(ifgd[-8:])
        expected = np.zeros(N_IM - 1)
        expected[ix_p:ix_s] = 1
        np.testing.assert_array_equal(G[i], expected)


def test_make_sb_matrix2():
    A = inv_lib.make_sb_matrix2(IFGDATES)
    assert A.shape == (N_IFG, N_IM)
    assert np.all(A.sum(axis=1) == 0)  # one -1 and one +1 per row
    assert np.all((A == -1).sum(axis=1) == 1)
    assert np.all((A == 1).sum(axis=1) == 1)


#%% Velocity estimation
def test_calc_vel_exact_linear():
    vel_true = np.array([0.0, 5.0, -3.0, 12.5])
    vconst_true = np.array([1.0, 0.0, -2.0, 4.0])
    cum = (vconst_true[:, None] + vel_true[:, None] * DT_CUM[None, :]
           ).astype(np.float32)
    vel, vconst = inv_lib.calc_vel(cum, DT_CUM)
    np.testing.assert_allclose(vel, vel_true, atol=1e-4)
    np.testing.assert_allclose(vconst, vconst_true, atol=1e-4)


def test_calc_vel_with_nan():
    vel_true = np.array([5.0, -3.0])
    cum = (vel_true[:, None] * DT_CUM[None, :]).astype(np.float32)
    cum[0, 2] = np.nan
    vel, vconst = inv_lib.calc_vel(cum, DT_CUM)
    np.testing.assert_allclose(vel, vel_true, atol=1e-4)


def test_calc_velsin_recovers_sine():
    # 3 years of 12-day sampling for a stable sin fit
    import datetime as dtmod
    imd0 = '20200101'
    days = np.arange(0, 3 * 365, 12)
    dt_cum = days / 365.25
    vel_true, amp_true = 4.0, 6.0
    delta_t_true = 100.0  # day of year of the sine peak offset
    phi = 2 * np.pi * (-delta_t_true) / 365.25
    cum = (vel_true * dt_cum
           + amp_true * np.sin(2 * np.pi * dt_cum + phi)
           ).astype(np.float32)[None, :]
    vel, vconst, amp, delta_t = inv_lib.calc_velsin(cum, dt_cum, imd0)
    np.testing.assert_allclose(vel[0], vel_true, atol=0.05)
    np.testing.assert_allclose(amp[0], amp_true, atol=0.05)
    np.testing.assert_allclose(delta_t[0], delta_t_true, atol=2)


#%% STC
def test_calc_stc_finite_for_smooth_field():
    rng = np.random.default_rng(1)
    cum = rng.normal(size=(4, 5, 5)).astype(np.float32)
    stc = inv_lib.calc_stc(cum)
    assert stc.shape == (5, 5)
    assert np.all(np.isfinite(stc))
    assert np.all(stc > 0)


def test_calc_stc_isolated_pixel_nan():
    rng = np.random.default_rng(1)
    cum = np.full((4, 5, 5), np.nan, dtype=np.float32)
    cum[:, 0, 0] = rng.normal(size=4)  # isolated pixel, no valid neighbor
    stc = inv_lib.calc_stc(cum)
    assert np.isnan(stc[0, 0])


#%% Censored least squares
def test_censored_lstsq_slow_recovers():
    rng = np.random.default_rng(2)
    A = rng.normal(size=(20, 3))
    X_true = rng.normal(size=(3, 5))
    B = A @ X_true
    M = rng.random((20, 5)) > 0.2  # ~20% censored
    X = inv_lib.censored_lstsq_slow(A, B, M)
    np.testing.assert_allclose(X, X_true, atol=1e-8)


def test_censored_lstsq2_multicolumn():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(20, 2))
    X_true = rng.normal(size=(2, 6))
    B = A @ X_true
    M = rng.random((20, 6)) > 0.2
    inv_lib.bootcount, inv_lib.bootnum = 0, 1  # globals used by a print
    X = inv_lib.censored_lstsq2(A, B, M)
    np.testing.assert_allclose(X, X_true, atol=1e-8)


def test_censored_lstsq2_single_column():
    inv_lib.bootcount, inv_lib.bootnum = 0, 1
    A = np.eye(4)
    B = np.arange(4.0)[:, None]
    M = np.ones((4, 1), dtype=bool)
    X = inv_lib.censored_lstsq2(A, B, M)
    np.testing.assert_allclose(np.ravel(X), np.arange(4.0))


def test_censored_lstsq2_1d():
    inv_lib.bootcount, inv_lib.bootnum = 0, 1
    A = np.eye(4)
    b = np.arange(4.0)
    m = np.ones(4, dtype=bool)
    X = inv_lib.censored_lstsq2(A, b, m)
    np.testing.assert_allclose(np.ravel(X), np.arange(4.0))


#%% NSBAS inversion
def test_invert_nsbas_forward_model():
    rng = np.random.default_rng(4)
    n_pt = 50
    vel_true = rng.uniform(-20, 20, n_pt)
    unw, inc_true = forward_model_unw(vel_true, n_pt)
    unw = unw.astype(np.float32)
    # Knock out some observations (keep enough for inversion)
    unw[rng.random(unw.shape) < 0.1] = np.nan
    unw[0, :] = np.nan  # one point with all nan
    G = inv_lib.make_sb_matrix(IFGDATES)
    inc, vel, vconst = inv_lib.invert_nsbas(unw, G, DT_CUM, 0.0001,
                                            n_core=1, gpu=False)
    assert inc.shape == (N_IM - 1, n_pt)
    np.testing.assert_allclose(vel[1:], vel_true[1:], atol=0.01)
    np.testing.assert_allclose(inc[:, 1:], inc_true[1:, :].T, atol=0.01)


def test_invert_nsbas_ncore2_matches_ncore1():
    rng = np.random.default_rng(5)
    n_pt = 30
    vel_true = rng.uniform(-20, 20, n_pt)
    unw, _ = forward_model_unw(vel_true, n_pt)
    unw = unw.astype(np.float32)
    unw[rng.random(unw.shape) < 0.15] = np.nan  # ensure point-by-point path
    G = inv_lib.make_sb_matrix(IFGDATES)
    res1 = inv_lib.invert_nsbas(unw.copy(), G, DT_CUM, 0.0001, 1, False)
    res2 = inv_lib.invert_nsbas(unw.copy(), G, DT_CUM, 0.0001, 2, False)
    for a, b in zip(res1, res2):
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-6)


def test_invert_nsbas_wls():
    rng = np.random.default_rng(6)
    n_pt = 20
    vel_true = rng.uniform(-20, 20, n_pt)
    unw, _ = forward_model_unw(vel_true, n_pt)
    unw = unw.astype(np.float32)
    var = np.full_like(unw, 2.0)
    G = inv_lib.make_sb_matrix(IFGDATES)
    inc, vel, vconst = inv_lib.invert_nsbas_wls(unw, var, G, DT_CUM,
                                                0.0001, n_core=1)
    np.testing.assert_allclose(vel, vel_true, atol=0.01)


def test_invert_nsbas_wls_ncore2_matches_ncore1():
    rng = np.random.default_rng(8)
    n_pt = 30
    vel_true = rng.uniform(-20, 20, n_pt)
    unw, _ = forward_model_unw(vel_true, n_pt)
    unw = unw.astype(np.float32)
    unw[rng.random(unw.shape) < 0.15] = np.nan
    var = np.full_like(unw, 2.0)
    G = inv_lib.make_sb_matrix(IFGDATES)
    res1 = inv_lib.invert_nsbas_wls(unw.copy(), var, G, DT_CUM, 0.0001, 1)
    res2 = inv_lib.invert_nsbas_wls(unw.copy(), var, G, DT_CUM, 0.0001, 2)
    for a, b in zip(res1, res2):
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-6)


#%% Bootstrap velocity std (seeded -> deterministic)
def test_calc_velstd_withnan():
    rng = np.random.default_rng(7)
    n_pt, n_im = 4, 20
    dt_cum = np.arange(n_im) * 12 / 365.25
    vel_true = np.array([0.0, 10.0, -5.0, 3.0])
    noise = rng.normal(0, 0.5, (n_pt, n_im))
    cum = (vel_true[:, None] * dt_cum[None, :] + noise).astype(np.float32)
    vstd = inv_lib.calc_velstd_withnan(cum, dt_cum)
    assert vstd.shape == (n_pt,)
    assert np.all(np.isfinite(vstd))
    assert np.all(vstd > 0)
    assert np.all(vstd < 5)  # noise 0.5 -> vstd of order 1 mm/yr
