"""End-to-end tests of --keep_ifg_list of step 12 (#119).

When a loop does not close, step 12 removes every ifg of the loop that is in
no closed loop, which can include good ifgs. tests/synth.py's
build_geocml_keep makes the simplest such case: a 2pi error in KEEP_ERR_IFG
also takes the good KEEP_GOOD_IFG and the first epoch with it. Each test
runs step 12 on a copy of the TS dir made by step 11, so that the runs with
different lists do not interfere.
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
from synth import KEEP_GOOD_IFG, KEEP_ERR_IFG, IMDATES, LOOP_ERR_COLS

pytestmark = pytest.mark.smoke


@pytest.fixture(scope='module')
def geocml_keep(tmp_path_factory, run_script):
    workdir = tmp_path_factory.mktemp('licsbas_keep')
    truth = synth.build_geocml_keep(workdir)
    run_script('LiCSBAS11_check_unw.py', '-d', 'GEOCml1', cwd=workdir)
    return truth


def run12(geocml_keep, run_script, tsname, keep_ifg=None):
    """Run step 12 into a copy of TS_GEOCml1 with keep_ifg as the list."""
    workdir = geocml_keep.workdir
    shutil.copytree(str(workdir / 'TS_GEOCml1'), str(workdir / tsname))
    args = []
    if keep_ifg is not None:
        listfile = workdir / (tsname + '_keep.txt')
        listfile.write_text(''.join(i + '\n' for i in keep_ifg))
        args = ['--keep_ifg_list', listfile.name]
    res = run_script('LiCSBAS12_loop_closure.py', '-d', 'GEOCml1',
                     '-t', tsname, '--n_para', '1', *args, cwd=workdir)
    return workdir / tsname, res


def read_list(tsdir, name):
    return (tsdir / name).read_text().split()


def first_loop_flag(tsdir):
    """Flag of the failing loop (0, 1, 2) in loop_info.txt."""
    for line in (tsdir / '12loop' / 'loop_info.txt').read_text().splitlines():
        f = line.split()
        if f[:3] == IMDATES[:3]:
            return f[4]
    raise AssertionError('No loop {} in loop_info.txt'.format(IMDATES[:3]))


def read_n_loop_err(tsdir, geocml_keep):
    return io_lib.read_img(str(tsdir / 'results' / 'n_loop_err'),
                           geocml_keep.length, geocml_keep.width)


def test_step12_removes_good_ifg_without_keep(geocml_keep, run_script):
    """The case of #119: the good ifg goes with the bad one."""
    tsdir, _ = run12(geocml_keep, run_script, 'TS_nokeep')
    assert read_list(tsdir, 'info/12bad_ifg.txt') == sorted(
        [KEEP_GOOD_IFG, KEEP_ERR_IFG])
    assert read_list(tsdir, 'info/12removed_image.txt') == [IMDATES[0]]


def test_step12_keep_good_ifg(geocml_keep, run_script):
    """Keeping the good ifg removes only the bad one and saves the epoch,
    and the time series still matches the truth."""
    absent = '20190101_20190113'  # not in GEOCml1
    tsdir, res = run12(geocml_keep, run_script, 'TS_keep_good',
                       [KEEP_GOOD_IFG, absent])

    assert read_list(tsdir, 'info/12bad_ifg.txt') == [KEEP_ERR_IFG]
    assert read_list(tsdir, 'info/12removed_image.txt') == []
    assert read_list(tsdir, '12loop/keep_ifg_man.txt') == [KEEP_GOOD_IFG,
                                                           absent]
    assert re.search(absent + r' \(not used', res.stdout)
    # The loop closure itself is recorded as is
    assert read_list(tsdir, '12loop/bad_ifg_loop.txt') == sorted(
        [KEEP_GOOD_IFG, KEEP_ERR_IFG])
    # The loop is still bad because of the removed ifg
    assert first_loop_flag(tsdir) == '*'
    np.testing.assert_array_equal(read_n_loop_err(tsdir, geocml_keep), 0)

    run_script('LiCSBAS13_sb_inv.py', '-d', 'GEOCml1', '-t', tsdir.name,
               '--n_para', '1', cwd=geocml_keep.workdir)
    with h5py.File(str(tsdir / 'cum.h5'), 'r') as f:
        imdates = [str(d) for d in f['imdates'][()]]
        vel = f['vel'][()]
        refarea = f['refarea'][()]
    assert imdates == IMDATES

    refarea = refarea.decode() if isinstance(refarea, bytes) else str(refarea)
    x1, x2, y1, y2 = [int(s) for s in re.split('[:/]', refarea)]
    truth = geocml_keep.vel_mm
    np.testing.assert_allclose(
        vel - np.nanmean(vel[y1:y2, x1:x2]),
        truth - np.nanmean(truth[y1:y2, x1:x2]), atol=0.1)


def test_step12_keep_all_ifgs_of_bad_loop(geocml_keep, run_script):
    """Kept ifgs stay even when the check with the ref point finds them bad
    again; the unclosed loop is then left as a candidate and counted in
    n_loop_err."""
    tsdir, _ = run12(geocml_keep, run_script, 'TS_keep_all',
                     [KEEP_GOOD_IFG, KEEP_ERR_IFG])

    assert read_list(tsdir, 'info/12bad_ifg.txt') == []
    assert read_list(tsdir, '12loop/bad_ifg_loopref.txt') == sorted(
        [KEEP_GOOD_IFG, KEEP_ERR_IFG])
    ifg12 = '{}_{}'.format(IMDATES[1], IMDATES[2])
    assert read_list(tsdir, 'info/12bad_ifg_cand.txt') == sorted(
        [KEEP_GOOD_IFG, KEEP_ERR_IFG, ifg12])
    assert first_loop_flag(tsdir) == '/'

    n_loop_err = read_n_loop_err(tsdir, geocml_keep)
    np.testing.assert_array_equal(n_loop_err[:, :LOOP_ERR_COLS], 1)
    np.testing.assert_array_equal(n_loop_err[:, LOOP_ERR_COLS:], 0)


def test_step12_keep_and_rm_same_ifg(geocml_keep, bin_env, repo_root):
    """An ifg both to be removed and kept is an error before processing."""
    workdir = geocml_keep.workdir
    (workdir / 'both.txt').write_text(KEEP_GOOD_IFG + '\n')
    res = subprocess.run(
        [sys.executable, str(repo_root / 'bin' / 'LiCSBAS12_loop_closure.py'),
         '-d', 'GEOCml1', '-t', 'TS_GEOCml1', '--rm_ifg_list', 'both.txt',
         '--keep_ifg_list', 'both.txt'],
        env=bin_env, cwd=str(workdir), capture_output=True, text=True,
        timeout=300)
    assert res.returncode == 2
    assert KEEP_GOOD_IFG + ' in both' in res.stderr
    assert not (workdir / 'TS_GEOCml1' / '12loop').exists()
