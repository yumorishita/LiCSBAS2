"""Unit tests for LiCSBAS_tools_lib."""
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
