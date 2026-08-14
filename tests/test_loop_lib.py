"""Unit tests for LiCSBAS_loop_lib (loop closure)."""
import numpy as np

import LiCSBAS_loop_lib as loop_lib
import synth


def test_make_loop_matrix_triangle():
    ifgdates = ['20200101_20200113', '20200113_20200125',
                '20200101_20200125']
    Aloop = loop_lib.make_loop_matrix(ifgdates)
    np.testing.assert_array_equal(Aloop, [[1, 1, -1]])


def test_make_loop_matrix_no_loop():
    ifgdates = ['20200101_20200113', '20200113_20200125']
    Aloop = loop_lib.make_loop_matrix(ifgdates)
    assert len(Aloop) == 0


def test_make_loop_matrix_network():
    Aloop = loop_lib.make_loop_matrix(synth.IFGDATES)
    # 5 epochs, 4 consecutive + 3 skip-1 ifgs -> 3+3=6 triplet loops:
    # (i,i+1,i+2) via skip ifg for i=0..2, each counted once, plus
    # loops with the skip ifgs as legs
    assert Aloop.shape[1] == len(synth.IFGDATES)
    assert len(Aloop) >= 3
    assert np.all((Aloop == 1).sum(axis=1) == 2)
    assert np.all((Aloop == -1).sum(axis=1) == 1)

    # A consistent phase field must close every loop exactly
    phase = np.array([synth.unw_phase(ifgd[:8], ifgd[-8:])[5, 5]
                      for ifgd in synth.IFGDATES])
    np.testing.assert_allclose(Aloop @ phase, 0, atol=1e-5)


def test_identify_bad_ifg():
    bad_cand = ['b', 'a', 'c', 'a']
    good = ['b']
    assert loop_lib.identify_bad_ifg(bad_cand, good) == ['a', 'c']


def test_read_unw_loop_ph(tmp_path):
    ifgdates = ['20200101_20200113', '20200113_20200125',
                '20200101_20200125']
    length = width = 10
    for ifgd in ifgdates:
        d = tmp_path / ifgd
        d.mkdir()
        synth.write_img(d / (ifgd + '.unw'),
                        synth.unw_phase(ifgd[:8], ifgd[-8:]))

    Aloop = loop_lib.make_loop_matrix(ifgdates)
    unw12, unw23, unw13, ifgd12, ifgd23, ifgd13 = loop_lib.read_unw_loop_ph(
        Aloop[0], ifgdates, str(tmp_path), length, width)

    assert (ifgd12, ifgd23, ifgd13) == tuple(ifgdates)
    loop_ph = unw12 + unw23 - unw13
    np.testing.assert_allclose(loop_ph, 0, atol=1e-5)


def test_read_unw_loop_ph_with_error(tmp_path):
    ifgdates = ['20200101_20200113', '20200113_20200125',
                '20200101_20200125']
    length = width = 10
    for ifgd in ifgdates:
        d = tmp_path / ifgd
        d.mkdir()
        unw = synth.unw_phase(ifgd[:8], ifgd[-8:])
        if ifgd == '20200113_20200125':
            unw = unw + 2 * np.pi  # simulated unwrapping error
        synth.write_img(d / (ifgd + '.unw'), unw)

    Aloop = loop_lib.make_loop_matrix(ifgdates)
    unw12, unw23, unw13, *_ = loop_lib.read_unw_loop_ph(
        Aloop[0], ifgdates, str(tmp_path), length, width)
    loop_ph = unw12 + unw23 - unw13
    np.testing.assert_allclose(loop_ph, 2 * np.pi, atol=1e-5)


def test_read_unw_loop_ph_zero_to_nan(tmp_path):
    # 0 in unw files is interpreted as nodata -> nan
    ifgdates = ['20200101_20200113', '20200113_20200125',
                '20200101_20200125']
    length = width = 5
    for ifgd in ifgdates:
        d = tmp_path / ifgd
        d.mkdir()
        unw = np.ones((length, width), dtype=np.float32)
        unw[0, 0] = 0
        synth.write_img(d / (ifgd + '.unw'), unw)

    Aloop = loop_lib.make_loop_matrix(ifgdates)
    unw12, *_ = loop_lib.read_unw_loop_ph(
        Aloop[0], ifgdates, str(tmp_path), length, width)
    assert np.isnan(unw12[0, 0])
    assert unw12[1, 1] == 1
