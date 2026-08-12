"""
Shared pytest configuration and fixtures for the LiCSBAS test suite.

LiCSBAS is not an installable package; LiCSBAS_lib must be on sys.path
(normally done by sourcing bashrc_LiCSBAS.sh). pytest.ini sets
pythonpath = LiCSBAS_lib, and the sys.path insertion below makes the
tests work even when invoked in ways that bypass pytest.ini.
"""
import os
import sys
import subprocess
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

REPO = Path(__file__).resolve().parents[1]
LIB = REPO / 'LiCSBAS_lib'
if str(LIB) not in sys.path:
    sys.path.insert(0, str(LIB))

import pytest  # noqa: E402


@pytest.fixture(scope='session')
def repo_root():
    return REPO


@pytest.fixture(scope='session')
def bin_env():
    """Environment for running bin/ scripts as subprocesses."""
    env = os.environ.copy()
    env['PYTHONPATH'] = str(LIB) + os.pathsep + env.get('PYTHONPATH', '')
    env['MPLBACKEND'] = 'Agg'
    env['QT_QPA_PLATFORM'] = 'offscreen'
    return env


@pytest.fixture(scope='session')
def run_script(bin_env):
    """Run bin/<script> with args in cwd and assert exit code 0."""
    def run(script, *args, cwd):
        cmd = [sys.executable, str(REPO / 'bin' / script)] + [str(a) for a in args]
        res = subprocess.run(cmd, env=bin_env, cwd=str(cwd),
                             capture_output=True, text=True, timeout=300)
        assert res.returncode == 0, (
            '{} failed with code {}\n--- stdout ---\n{}\n--- stderr ---\n{}'
            .format(script, res.returncode, res.stdout, res.stderr))
        return res
    return run


# The synthetic dataset builder lives in tests/synth.py. Import it here so
# test modules and fixtures can share it regardless of invocation directory.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import synth  # noqa: E402


@pytest.fixture(scope='session')
def geocml(tmp_path_factory):
    """Build the synthetic GEOCml1 dataset once per session.

    Returns a SimpleNamespace with workdir, geocdir and the truth model
    (vel_mm, imdates, ifgdates, coef_r2m, ...).
    """
    workdir = tmp_path_factory.mktemp('licsbas_smoke')
    truth = synth.build_geocml(workdir)
    return truth


@pytest.fixture(scope='session')
def ts11(geocml, run_script):
    run_script('LiCSBAS11_check_unw.py', '-d', 'GEOCml1', cwd=geocml.workdir)
    return geocml.workdir / 'TS_GEOCml1'


@pytest.fixture(scope='session')
def ts12(ts11, geocml, run_script):
    run_script('LiCSBAS12_loop_closure.py', '-d', 'GEOCml1', '--n_para', '1',
               cwd=geocml.workdir)
    return ts11


@pytest.fixture(scope='session')
def ts13(ts12, geocml, run_script):
    run_script('LiCSBAS13_sb_inv.py', '-d', 'GEOCml1', '--n_para', '1',
               cwd=geocml.workdir)
    return ts12


@pytest.fixture(scope='session')
def ts14(ts13, geocml, run_script):
    run_script('LiCSBAS14_vel_std.py', '-t', 'TS_GEOCml1', cwd=geocml.workdir)
    return ts13


@pytest.fixture(scope='session')
def ts15(ts14, geocml, run_script):
    run_script('LiCSBAS15_mask_ts.py', '-t', 'TS_GEOCml1', cwd=geocml.workdir)
    return ts14


@pytest.fixture(scope='session')
def ts16(ts15, geocml, run_script):
    run_script('LiCSBAS16_filt_ts.py', '-t', 'TS_GEOCml1', '--n_para', '1',
               cwd=geocml.workdir)
    return ts15
