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


# --- Defective dataset: exercises the rejection paths of steps 11-13 ---

@pytest.fixture(scope='session')
def geocml_defect(tmp_path_factory):
    """Build the synthetic GEOCml1 dataset with defective ifgs.

    Separate workdir from the geocml fixture so the two chains do not
    interfere. See synth.build_geocml_defect for what each defect is.
    """
    workdir = tmp_path_factory.mktemp('licsbas_defect')
    return synth.build_geocml_defect(workdir)


@pytest.fixture(scope='session')
def ts11d(geocml_defect, run_script):
    run_script('LiCSBAS11_check_unw.py', '-d', 'GEOCml1',
               cwd=geocml_defect.workdir)
    return geocml_defect.workdir / 'TS_GEOCml1'


@pytest.fixture(scope='session')
def ts12d(ts11d, geocml_defect, run_script):
    run_script('LiCSBAS12_loop_closure.py', '-d', 'GEOCml1', '--n_para', '1',
               cwd=geocml_defect.workdir)
    return ts11d


@pytest.fixture(scope='session')
def ts13d(ts12d, geocml_defect, run_script):
    run_script('LiCSBAS13_sb_inv.py', '-d', 'GEOCml1', '--n_para', '1',
               cwd=geocml_defect.workdir)
    return ts12d


@pytest.fixture(scope='session')
def ts14d(ts13d, geocml_defect, run_script):
    run_script('LiCSBAS14_vel_std.py', '-t', 'TS_GEOCml1',
               cwd=geocml_defect.workdir)
    return ts13d


@pytest.fixture(scope='session')
def ts15d(ts14d, geocml_defect, run_script):
    run_script('LiCSBAS15_mask_ts.py', '-t', 'TS_GEOCml1',
               cwd=geocml_defect.workdir)
    return ts14d


@pytest.fixture(scope='session')
def ts16d(ts15d, geocml_defect, run_script):
    run_script('LiCSBAS16_filt_ts.py', '-t', 'TS_GEOCml1', '--n_para', '1',
               cwd=geocml_defect.workdir)
    return ts15d


# --- Steps 02, 04 and 05 (preparation of GEOCml) ---

@pytest.fixture(scope='session')
def geoc(tmp_path_factory):
    """Synthetic GEOC (GeoTIFF) dataset. See synth.build_geoc."""
    return synth.build_geoc(tmp_path_factory.mktemp('licsbas_geoc'))


@pytest.fixture(scope='session')
def geocml1_02(geoc, run_script):
    run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
               '--n_para', '1', cwd=geoc.workdir)
    return geoc.workdir / 'GEOCml1'


@pytest.fixture(scope='session')
def geocml2_02(geoc, run_script):
    run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', str(synth.NLOOK),
               '--n_para', '1', cwd=geoc.workdir)
    return geoc.workdir / 'GEOCml{}'.format(synth.NLOOK)


@pytest.fixture(scope='session')
def geoc_nometa(tmp_path_factory):
    """GEOC without metadata.txt nor baselines."""
    return synth.build_geoc(tmp_path_factory.mktemp('licsbas_geoc_nometa'),
                            metadata=False, baselines=False)


@pytest.fixture(scope='session')
def geocml1_02_nometa(geoc_nometa, run_script):
    """Step 02 without metadata.txt and without --freq: all defaults."""
    run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
               '--n_para', '1', cwd=geoc_nometa.workdir)
    return geoc_nometa.workdir / 'GEOCml1'


@pytest.fixture(scope='session')
def geocml1_02_freq(geoc_nometa, run_script):
    """Step 02 without metadata.txt, with --freq."""
    run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
               '-o', 'GEOCml1freq', '--freq', '1.27e9', '--n_para', '1',
               cwd=geoc_nometa.workdir)
    return geoc_nometa.workdir / 'GEOCml1freq'


@pytest.fixture(scope='session')
def geocml1_02_freq_meta(geoc, run_script):
    """Step 02 with both metadata.txt (with radar_freq) and --freq.

    Returns the output dir and the CompletedProcess (for the warning)."""
    res = run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
                     '-o', 'GEOCml1freq', '--freq', '1.27e9', '--n_para', '1',
                     cwd=geoc.workdir)
    return geoc.workdir / 'GEOCml1freq', res


@pytest.fixture(scope='session')
def geoc_meta_nofreq(tmp_path_factory):
    """GEOC whose metadata.txt has center_time but no radar_freq."""
    return synth.build_geoc(tmp_path_factory.mktemp('licsbas_geoc_nofreq'),
                            metadata_freq=False)


@pytest.fixture(scope='session')
def geocml1_02_meta_nofreq(geoc_meta_nofreq, run_script):
    run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
               '--n_para', '1', cwd=geoc_meta_nofreq.workdir)
    return geoc_meta_nofreq.workdir / 'GEOCml1'


@pytest.fixture(scope='session')
def geocml1_02_meta_nofreq_freq(geoc_meta_nofreq, run_script):
    run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
               '-o', 'GEOCml1freq', '--freq', '1.27e9', '--n_para', '1',
               cwd=geoc_meta_nofreq.workdir)
    return geoc_meta_nofreq.workdir / 'GEOCml1freq'


@pytest.fixture(scope='session')
def geocml_prep(tmp_path_factory):
    """Input GEOCml1 of steps 04 and 05. See synth.build_geocml_prep."""
    return synth.build_geocml_prep(tmp_path_factory.mktemp('licsbas_prep'))


@pytest.fixture(scope='session')
def mask04(geocml_prep, run_script):
    run_script('LiCSBAS04op_mask_unw.py', '-i', 'GEOCml1', '-o', 'GEOCml1mask',
               '-c', '0.2', '-r', '0:2/0:3', '--n_para', '1',
               cwd=geocml_prep.workdir)
    return geocml_prep.workdir / 'GEOCml1mask'


@pytest.fixture(scope='session')
def mask04_file(geocml_prep, run_script):
    rangefile = geocml_prep.workdir / 'mask_ranges.txt'
    rangefile.write_text('3:5/4:6\n6:8/0:2\n')
    run_script('LiCSBAS04op_mask_unw.py', '-i', 'GEOCml1',
               '-o', 'GEOCml1maskf', '-f', rangefile.name, '--n_para', '1',
               cwd=geocml_prep.workdir)
    return geocml_prep.workdir / 'GEOCml1maskf'


@pytest.fixture(scope='session')
def clip05(geocml_prep, run_script):
    run_script('LiCSBAS05op_clip_unw.py', '-i', 'GEOCml1', '-o', 'GEOCml1clip',
               '-r', '1:7/1:10', '--n_para', '1', cwd=geocml_prep.workdir)
    return geocml_prep.workdir / 'GEOCml1clip'


@pytest.fixture(scope='session')
def clip05_geo(geocml_prep, run_script):
    # Same area as clip05 in lon/lat (grid registration). lat_s is on the
    # frame edge, so the clamping branch of tools_lib.read_range_geo is
    # exercised as well as the unclamped ones.
    run_script('LiCSBAS05op_clip_unw.py', '-i', 'GEOCml1',
               '-o', 'GEOCml1clipg', '-g', '132.001/132.006/33.991/33.999',
               '--n_para', '1', cwd=geocml_prep.workdir)
    return geocml_prep.workdir / 'GEOCml1clipg'


# --- Step 03 (GACOS) ---

@pytest.fixture(scope='session')
def geocml_gacos(tmp_path_factory):
    """GEOCml1 with U.geo and GACOS sltd. See synth.build_geocml_gacos."""
    return synth.build_geocml_gacos(tmp_path_factory.mktemp('licsbas_gacos'))


@pytest.fixture(scope='session')
def gacos03(geocml_gacos, run_script):
    run_script('LiCSBAS03op_GACOS.py', '-i', 'GEOCml1', '-o', 'GEOCml1GACOS',
               '-g', 'GACOS', '--n_para', '1', cwd=geocml_gacos.workdir)
    return geocml_gacos.workdir / 'GEOCml1GACOS'


@pytest.fixture(scope='session')
def gacos03_ztd(geocml_gacos, run_script):
    """Step 03 from ztd.tif (m) instead of sltd.geo.tif (rad)."""
    run_script('LiCSBAS03op_GACOS.py', '-i', 'GEOCml1', '-o', 'GEOCml1GACOSztd',
               '-g', 'GACOS_ztd', '--n_para', '1', cwd=geocml_gacos.workdir)
    return geocml_gacos.workdir / 'GEOCml1GACOSztd'
