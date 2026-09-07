"""Unit tests for LiCSBAS_tools_lib."""
import os
import pathlib
import subprocess
import sys

import numpy as np
import pytest

import LiCSBAS_tools_lib as tools_lib


# Geometry used in several tests (grid registration, north-up)
LAT1, POSTLAT = 34.0, -0.001
LON1, POSTLON = 132.0, 0.001
WIDTH, LENGTH = 10, 10


#%% Coordinate conversions
def test_bl2xy_corner():
    assert tools_lib.bl2xy(LON1, LAT1, WIDTH, LENGTH,
                           LAT1, POSTLAT, LON1, POSTLON) == [0, 0]


def test_bl2xy_xy2bl_roundtrip():
    x, y = 7, 3
    lat, lon = tools_lib.xy2bl(x, y, LAT1, POSTLAT, LON1, POSTLON)
    assert tools_lib.bl2xy(lon, lat, WIDTH, LENGTH,
                           LAT1, POSTLAT, LON1, POSTLON) == [x, y]


def test_read_range_geo():
    # Full-extent bbox covers the whole image
    lat2 = LAT1 + POSTLAT * (LENGTH - 1)
    lon2 = LON1 + POSTLON * (WIDTH - 1)
    range_str = '{}/{}/{}/{}'.format(LON1, lon2, lat2, LAT1)
    assert tools_lib.read_range_geo(range_str, WIDTH, LENGTH,
                                    LAT1, POSTLAT, LON1, POSTLON) \
        == [0, WIDTH, 0, LENGTH]

    # Interior bbox: x2/y2 are exclusive-style (+1)
    lon_w, lon_e = tools_lib.xy2bl(2, 0, LAT1, POSTLAT, LON1, POSTLON)[1], \
        tools_lib.xy2bl(5, 0, LAT1, POSTLAT, LON1, POSTLON)[1]
    lat_n = tools_lib.xy2bl(0, 3, LAT1, POSTLAT, LON1, POSTLON)[0]
    lat_s = tools_lib.xy2bl(0, 6, LAT1, POSTLAT, LON1, POSTLON)[0]
    range_str = '{}/{}/{}/{}'.format(lon_w, lon_e, lat_s, lat_n)
    assert tools_lib.read_range_geo(range_str, WIDTH, LENGTH,
                                    LAT1, POSTLAT, LON1, POSTLON) \
        == [2, 6, 3, 7]


#%% String parsers
def test_read_point():
    assert tools_lib.read_point('3/5', WIDTH, LENGTH) == [3, 5]
    assert tools_lib.read_point('10/5', WIDTH, LENGTH) is False  # out of range
    assert tools_lib.read_point('bad', WIDTH, LENGTH) is False


def test_read_range():
    assert tools_lib.read_range('2:5/3:8', WIDTH, LENGTH) == [2, 5, 3, 8]
    # 0 for x2/y2 means to the edge
    assert tools_lib.read_range('2:0/3:0', WIDTH, LENGTH) == [2, WIDTH, 3, LENGTH]
    assert tools_lib.read_range('5:2/3:8', WIDTH, LENGTH) is False  # x1 >= x2
    assert tools_lib.read_range('2:11/3:8', WIDTH, LENGTH) is False  # exceeds
    assert tools_lib.read_range('bad', WIDTH, LENGTH) is False


def test_read_range_line():
    # Input is x1,y1/x2,y2 but return is [x1, x2, y1, y2]
    assert tools_lib.read_range_line('1,2/8,9', WIDTH, LENGTH) == [1, 8, 2, 9]
    assert tools_lib.read_range_line('1,2/8,10', WIDTH, LENGTH) is False


def test_ifgdates2imdates():
    ifgdates = ['20200113_20200125', '20200101_20200113', '20200101_20200125']
    assert tools_lib.ifgdates2imdates(ifgdates) \
        == ['20200101', '20200113', '20200125']


def test_get_ifgdates(tmp_path):
    for d in ['20200101_20200113', '20200113_20200125']:
        (tmp_path / d).mkdir()
    (tmp_path / '20200101_20200113.txt').touch()  # file, not dir
    (tmp_path / 'not_a_pair_dir').mkdir()
    (tmp_path / '19990101_19990113').mkdir()  # not starting with 2
    assert tools_lib.get_ifgdates(str(tmp_path)) \
        == ['20200101_20200113', '20200113_20200125']


def test_get_pair_folders(tmp_path):
    for d in ['20200113_20200125', '20200101_20200113']:
        (tmp_path / d).mkdir()
    (tmp_path / 'network').mkdir()
    assert tools_lib.get_pair_folders(str(tmp_path)) \
        == ['20200101_20200113', '20200113_20200125']


#%% Numeric utilities
def test_multilook_mean():
    array = np.arange(16, dtype=np.float32).reshape(4, 4)
    ml = tools_lib.multilook(array, 2, 2)
    expected = np.array([[2.5, 4.5], [10.5, 12.5]])
    np.testing.assert_allclose(ml, expected)


def test_multilook_nan_threshold():
    array = np.arange(16, dtype=np.float32).reshape(4, 4)
    array[0, 0] = array[0, 1] = array[1, 0] = np.nan  # 3/4 invalid in block
    ml = tools_lib.multilook(array, 2, 2, n_valid_thre=0.5)
    assert np.isnan(ml[0, 0])
    assert not np.isnan(ml[0, 1])
    # With a lower threshold the single valid pixel is kept
    ml2 = tools_lib.multilook(array, 2, 2, n_valid_thre=0.25)
    assert ml2[0, 0] == 5.0


def test_get_patchrow():
    width, length = 100, 1000
    n_data, memory_size = 100, 10  # 100*1000*100*4B = ~38 MiB -> 4 patches
    n_patch, patchrow = tools_lib.get_patchrow(width, length, n_data,
                                               memory_size)
    assert n_patch == len(patchrow)
    assert patchrow[0][0] == 0
    assert patchrow[-1][-1] == length
    for (a, b), (c, d) in zip(patchrow[:-1], patchrow[1:]):
        assert b == c  # contiguous


def test_convert_size():
    assert tools_lib.convert_size(0) == '0B'
    assert tools_lib.convert_size(1024) == '1.0KB'
    assert tools_lib.convert_size(1536) == '1.5KB'
    assert tools_lib.convert_size(5 * 1024**3) == '5.0GB'


#%% run_pool
def _square(i):
    return i * i


def _row(i):
    return np.full(4, i, dtype=np.float32)


def _pair(i):
    return i, i * i


@pytest.mark.parametrize('n_para', [1, 3])  # serial path and parallel path
def test_run_pool_preserves_order(n_para):
    assert tools_lib.run_pool(_square, range(50), n_para) \
        == [i * i for i in range(50)]


def test_run_pool_empty_args():
    assert tools_lib.run_pool(_square, [], 4) == []


def test_run_pool_out_streaming():
    out = np.empty((10, 4), dtype=np.float32)
    ret = tools_lib.run_pool(_row, range(10), 3, chunksize=1, out=out)
    assert ret is out
    expected = np.repeat(np.arange(10, dtype=np.float32)[:, None], 4, axis=1)
    np.testing.assert_array_equal(out, expected)


def test_run_pool_store():
    firsts = np.empty(10, dtype=np.int64)
    seconds = np.empty(10, dtype=np.int64)

    def store(i, result):
        firsts[i], seconds[i] = result

    tools_lib.run_pool(_pair, range(10), 3, store=store)
    np.testing.assert_array_equal(firsts, np.arange(10))
    np.testing.assert_array_equal(seconds, np.arange(10) ** 2)


def test_run_pool_mem_cap(capsys):
    # Absurdly large per-worker memory must cap n_para to 1 but still work
    result = tools_lib.run_pool(_square, range(10), 4,
                                mem_per_worker_mb=10**12)
    assert result == [i * i for i in range(10)]
    assert 'Reduce n_para' in capsys.readouterr().out


@pytest.mark.skipif(sys.platform != 'linux', reason='needs SIGKILL')
def test_run_pool_dead_worker_raises_instead_of_hanging():
    # A worker killed mid-map (like by the OOM killer) must raise a clear
    # error. Run in a subprocess with a hard timeout so that a regression
    # (hang, as plain Pool.map does) can never hang the test session.
    lib_dir = str(pathlib.Path(__file__).resolve().parent.parent
                  / 'LiCSBAS_lib')
    code = '\n'.join([
        'import os, signal, sys, time',
        'sys.path.insert(0, {!r})'.format(lib_dir),
        'import LiCSBAS_tools_lib as tools_lib',
        'def f(i):',
        '    if i == 3: os.kill(os.getpid(), signal.SIGKILL)',
        '    time.sleep(0.01)',
        '    return i',
        'try:',
        '    tools_lib.run_pool(f, range(8), 2, chunksize=1)',
        'except RuntimeError as e:',
        '    assert "n_para" in str(e)',
        '    print("OK")',
    ])
    res = subprocess.run([sys.executable, '-c', code], timeout=60,
                         capture_output=True, text=True)
    assert 'OK' in res.stdout, res.stderr


#%% cgroup awareness
def _make_cgroup_v2(path, memory_max=None, memory_current=None,
                    inactive_file=None, cpu_max=None):
    path.mkdir(parents=True, exist_ok=True)
    (path / 'memory.max').write_text('{}\n'.format(
        'max' if memory_max is None else memory_max))
    if memory_current is not None:
        (path / 'memory.current').write_text('{}\n'.format(memory_current))
    if inactive_file is not None:
        (path / 'memory.stat').write_text(
            'anon 12345\ninactive_file {}\n'.format(inactive_file))
    (path / 'cpu.max').write_text('{}\n'.format(
        'max 100000' if cpu_max is None else cpu_max))
    return path


def test_read_cgroup_int(tmp_path):
    (tmp_path / 'unlimited').write_text('max\n')
    (tmp_path / 'unlimited_v1').write_text('-1\n')
    (tmp_path / 'huge_v1').write_text('9223372036854771712\n')
    (tmp_path / 'limited').write_text('1234\n')
    assert tools_lib._read_cgroup_int(tmp_path / 'unlimited') is None
    assert tools_lib._read_cgroup_int(tmp_path / 'unlimited_v1') is None
    assert tools_lib._read_cgroup_int(tmp_path / 'huge_v1') is None
    assert tools_lib._read_cgroup_int(tmp_path / 'nonexistent') is None
    assert tools_lib._read_cgroup_int(tmp_path / 'limited') == 1234


def test_get_cgroup_dirs_no_crash():
    # Whatever (or no) cgroup this test runs in, it must return a list of
    # existing directories without raising
    dirs = tools_lib._get_cgroup_dirs()
    assert isinstance(dirs, list)
    assert all(os.path.isdir(d) for d in dirs)
    assert len(dirs) == len(set(dirs))


def test_get_mem_avail_mb_no_cgroup(monkeypatch):
    import psutil
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs', lambda: [])
    mem_avail = tools_lib.get_mem_avail_mb()
    assert mem_avail == pytest.approx(
        psutil.virtual_memory().available / 2**20, rel=0.1)


def test_get_mem_avail_mb_cgroup_v2(monkeypatch, tmp_path):
    # 100 MB limit, 40 MB used of which 10 MB is reclaimable page cache
    cg = _make_cgroup_v2(tmp_path / 'cg', memory_max=100 * 2**20,
                         memory_current=40 * 2**20,
                         inactive_file=10 * 2**20)
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs', lambda: [str(cg)])
    assert tools_lib.get_mem_avail_mb() == pytest.approx(70, abs=1)


def test_get_mem_avail_mb_cgroup_ancestor_limit(monkeypatch, tmp_path):
    # The tightest limit in the hierarchy wins, even if it is an ancestor
    parent = _make_cgroup_v2(tmp_path / 'parent', memory_max=50 * 2**20,
                             memory_current=10 * 2**20)
    child = _make_cgroup_v2(parent / 'child', memory_current=10 * 2**20)
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs',
                        lambda: [str(child), str(parent)])
    assert tools_lib.get_mem_avail_mb() == pytest.approx(40, abs=1)


def test_get_mem_avail_mb_cgroup_v1(monkeypatch, tmp_path):
    cg = tmp_path / 'cg'
    cg.mkdir()
    (cg / 'memory.limit_in_bytes').write_text('{}\n'.format(200 * 2**20))
    (cg / 'memory.usage_in_bytes').write_text('{}\n'.format(80 * 2**20))
    (cg / 'memory.stat').write_text(
        'cache 0\ntotal_inactive_file {}\n'.format(30 * 2**20))
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs', lambda: [str(cg)])
    assert tools_lib.get_mem_avail_mb() == pytest.approx(150, abs=1)


@pytest.mark.skipif(not hasattr(os, 'sched_getaffinity'),
                    reason='needs sched_getaffinity')
def test_get_n_cpu_avail_no_quota(monkeypatch):
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs', lambda: [])
    assert tools_lib.get_n_cpu_avail() == len(os.sched_getaffinity(0))


def test_get_n_cpu_avail_cgroup_v2(monkeypatch, tmp_path):
    cg = _make_cgroup_v2(tmp_path / 'cg', cpu_max='150000 100000')
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs', lambda: [str(cg)])
    assert tools_lib.get_n_cpu_avail() == 2  # 1.5 cores, rounded up


def test_get_n_cpu_avail_cgroup_v1(monkeypatch, tmp_path):
    cg = tmp_path / 'cg'
    cg.mkdir()
    (cg / 'cpu.cfs_quota_us').write_text('200000\n')
    (cg / 'cpu.cfs_period_us').write_text('100000\n')
    monkeypatch.setattr(tools_lib, '_get_cgroup_dirs', lambda: [str(cg)])
    assert tools_lib.get_n_cpu_avail() == 2


def test_run_pool_cpu_quota_cap(monkeypatch, capsys):
    monkeypatch.setattr(tools_lib, 'get_n_cpu_avail', lambda: 1)
    assert tools_lib.run_pool(_square, range(10), 4) \
        == [i * i for i in range(10)]
    assert 'CPU quota' in capsys.readouterr().out


def test_run_pool_mem_cap_uses_cgroup(monkeypatch, capsys):
    # 100 MB available in the cgroup allows only 100/2/40 = 1 worker,
    # regardless of the memory of the host
    monkeypatch.setattr(tools_lib, 'get_mem_avail_mb', lambda: 100)
    result = tools_lib.run_pool(_square, range(10), 4, mem_per_worker_mb=40)
    assert result == [i * i for i in range(10)]
    assert 'Reduce n_para from 4 to 1' in capsys.readouterr().out


#%% calculate_common_geometry
def test_calculate_common_geometry_overlap():
    gt1 = (132.0, 0.001, 0, 34.0, 0, -0.001)
    gt2 = (132.005, 0.001, 0, 33.995, 0, -0.001)
    new_gt, (ny, nx) = tools_lib.calculate_common_geometry(
        gt1, (10, 10), gt2, (10, 10))
    assert new_gt == (132.005, 0.001, 0, 33.995, 0, -0.001)
    assert (ny, nx) == (5, 5)


def test_calculate_common_geometry_no_overlap():
    gt1 = (132.0, 0.001, 0, 34.0, 0, -0.001)
    gt2 = (133.0, 0.001, 0, 34.0, 0, -0.001)
    with pytest.raises(ValueError):
        tools_lib.calculate_common_geometry(gt1, (10, 10), gt2, (10, 10))


#%% 2D fitting
def test_fit2d_exact_plane():
    a, b, c = 3.0, 0.5, -0.25
    X, Y = np.meshgrid(np.arange(20), np.arange(15))
    A = (a + b * X + c * Y).astype(np.float32)
    Afit, m = tools_lib.fit2d(A, deg='1')
    np.testing.assert_allclose(Afit, A, atol=1e-3)
    np.testing.assert_allclose(m, [a, b, c], atol=1e-5)


def test_fit2d_with_nan():
    a, b, c = 3.0, 0.5, -0.25
    X, Y = np.meshgrid(np.arange(20), np.arange(15))
    A = (a + b * X + c * Y).astype(np.float32)
    A[3, 4] = A[10, 10] = np.nan
    Afit, m = tools_lib.fit2d(A, deg='1')
    np.testing.assert_allclose(m, [a, b, c], atol=1e-4)
    assert not np.any(np.isnan(Afit))  # fit fills the nan holes


@pytest.mark.parametrize('deg, n_param', [('1', 3), ('bl', 4), ('2', 6)])
def test_fit2d_degrees(deg, n_param):
    X, Y = np.meshgrid(np.arange(20), np.arange(15))
    A = (1.0 + 0.1 * X - 0.2 * Y + 0.01 * X * Y).astype(np.float32)
    Afit, m = tools_lib.fit2d(A, deg=deg)
    assert len(m) == n_param
    if deg != '1':  # bilinear term only representable for bl and 2
        np.testing.assert_allclose(Afit, A, atol=1e-3)


def test_fit2d_bad_deg_returns_false():
    # Documents current behavior: a bare False is returned (not a tuple),
    # so callers that unpack two values raise TypeError. A future fix
    # should probably raise instead.
    A = np.ones((5, 5), dtype=np.float32)
    assert tools_lib.fit2d(A, deg='7') is False


def test_fit2dh_recovers_hgt_coefficient():
    k = 0.05
    X, Y = np.meshgrid(np.arange(20), np.arange(15))
    hgt = (np.random.default_rng(0).integers(0, 1000, (15, 20))
           ).astype(np.float32)
    A = (2.0 + 0.1 * X - 0.2 * Y + k * hgt).astype(np.float32)
    Afit, m = tools_lib.fit2dh(A, '1', hgt, -100, 10000)
    np.testing.assert_allclose(m[-1], k, atol=1e-6)
    np.testing.assert_allclose(Afit, A, atol=1e-3)


#%% Colormaps
@pytest.mark.parametrize('name', ['viridis', 'SCM.roma', 'GMT.polar',
                                  'cm_insar', 'cm_insar_r', 'cm_isce'])
def test_get_cmap(name):
    cmap = tools_lib.get_cmap(name, 256)
    assert cmap.N == 256
    rgba = cmap(np.linspace(0, 1, 256))
    assert rgba.shape == (256, 4)
    assert np.all((rgba >= 0) & (rgba <= 1))


def test_cmap_insar_cdict():
    cdict = tools_lib.cmap_insar()
    assert set(cdict.keys()) == {'red', 'green', 'blue'}
    for tup in cdict.values():
        assert tup[0][0] == 0.0
        assert tup[-1][0] == 1.0
