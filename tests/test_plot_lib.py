"""Smoke tests for LiCSBAS_plot_lib (and loop_lib.make_loop_png).

These only assert that a non-empty output file is produced; pixel
content is not compared.
"""
import numpy as np
import pytest

import LiCSBAS_plot_lib as plot_lib
import LiCSBAS_loop_lib as loop_lib
import synth


def _assert_png(path):
    assert path.exists() and path.stat().st_size > 0


def test_make_im_png(tmp_path):
    data = np.arange(100, dtype=np.float32).reshape(10, 10)
    png = tmp_path / 'im.png'
    plot_lib.make_im_png(data, str(png), 'viridis', 'title')
    _assert_png(png)

    # Documents current behavior: float32 input is mutated in place
    # (0 replaced by nan in the caller's array)
    assert np.isnan(data[0, 0])
    assert data[0, 1] == 1


def test_make_im_png_insar_cmap(tmp_path):
    data = np.random.default_rng(0).uniform(
        -np.pi, np.pi, (10, 10)).astype(np.float32)
    png = tmp_path / 'im_insar.png'
    plot_lib.make_im_png(data, str(png), 'insar', 'wrapped',
                         vmin=-np.pi, vmax=np.pi)
    _assert_png(png)


def test_make_3im_png(tmp_path):
    data3 = [np.ones((10, 10), dtype=np.float32) * i for i in range(1, 4)]
    png = tmp_path / '3im.png'
    plot_lib.make_3im_png(data3, str(png), 'viridis', ['a', 'b', 'c'])
    _assert_png(png)

    # Documents current behavior: the input list entries are emptied
    # to free memory, so callers cannot reuse them
    assert all(d == [] for d in data3)


def test_plot_network(tmp_path):
    png = tmp_path / 'network.png'
    plot_lib.plot_network(synth.IFGDATES, synth.BPERP, [], str(png))
    _assert_png(png)


def test_plot_network_with_removed(tmp_path):
    png = tmp_path / 'network_rm.png'
    plot_lib.plot_network(synth.IFGDATES, synth.BPERP,
                          [synth.IFGDATES[0]], str(png))
    _assert_png(png)


def test_plot_network_passes_locator_to_date_formatter(tmp_path):
    # plot_network used to pass the return value of set_major_locator()
    # (which is always None) into ConciseDateFormatter, leaving the
    # formatter without a locator. matplotlib <=3.10 tolerated it, 3.11
    # raises AttributeError while drawing. Capture what plot_network
    # actually hands to the formatter so the bug cannot come back.
    import matplotlib.dates as mdates

    captured = []
    orig = mdates.ConciseDateFormatter

    class Spy(orig):
        def __init__(self, locator, *args, **kwargs):
            captured.append(locator)
            super().__init__(locator, *args, **kwargs)

    mdates.ConciseDateFormatter = Spy
    try:
        plot_lib.plot_network(synth.IFGDATES, synth.BPERP, [],
                              str(tmp_path / 'network_loc.png'))
    finally:
        mdates.ConciseDateFormatter = orig

    assert captured, 'ConciseDateFormatter was never constructed'
    assert isinstance(captured[0], mdates.DateLocator), \
        'plot_network passed {!r} instead of a locator'.format(captured[0])


def test_plot_hgt_corr(tmp_path):
    rng = np.random.default_rng(1)
    n = 200
    hgt = rng.uniform(0, 1000, n).astype(np.float32)
    fit_hgt = (0.01 * hgt).astype(np.float32)
    data_bf = fit_hgt + rng.normal(0, 1, n).astype(np.float32)
    data_bf[5] = np.nan
    png = tmp_path / 'hgt_corr.png'
    plot_lib.plot_hgt_corr(data_bf, fit_hgt, hgt, 'test', str(png))
    _assert_png(png)


def test_plot_gacos_info(tmp_path):
    info = tmp_path / 'gacos_info.txt'
    info.write_text(
        'date       std_bf std_af rate\n'
        '20200101   2.5    1.2    52.0%\n'
        '20200113   3.1    1.5    51.6%\n'
        '20200125   0.0    0.0    0.0%\n'   # skipped
        '20200206   nan    nan    nan%\n'   # skipped
        '20200218   1.8    2.0    -11.1%\n')
    png = tmp_path / 'gacos_info.png'
    plot_lib.plot_gacos_info(str(info), str(png))
    _assert_png(png)


def test_make_loop_png(tmp_path):
    unw12 = synth.unw_phase('20200101', '20200113')
    unw23 = synth.unw_phase('20200113', '20200125')
    unw13 = synth.unw_phase('20200101', '20200125')
    loop_ph = unw12 + unw23 - unw13
    png = tmp_path / 'loop.png'
    loop_lib.make_loop_png(unw12, unw23, unw13, loop_ph, str(png),
                           ['12', '23', '13', 'loop'], 3)
    _assert_png(png)
