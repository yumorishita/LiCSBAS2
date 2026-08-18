"""Unit tests for LiCSBAS_io_lib (file I/O, GAMMA par files, GeoTIFF)."""
import numpy as np
import pytest

import LiCSBAS_io_lib as io_lib
import synth


#%% Raw binary images
def test_read_img_roundtrip_float32(tmp_path):
    data = np.arange(20, dtype=np.float32).reshape(4, 5)
    file = tmp_path / 'img.raw'
    data.tofile(str(file))
    np.testing.assert_array_equal(io_lib.read_img(str(file), 4, 5), data)


def test_read_img_uint8(tmp_path):
    data = np.arange(20, dtype=np.uint8).reshape(4, 5)
    file = tmp_path / 'img.raw'
    data.tofile(str(file))
    out = io_lib.read_img(str(file), 4, 5, dtype=np.uint8)
    np.testing.assert_array_equal(out, data)


def test_read_img_big_endian(tmp_path):
    data = np.arange(20, dtype=np.float32).reshape(4, 5)
    file = tmp_path / 'img.raw'
    data.astype('>f4').tofile(str(file))
    out = io_lib.read_img(str(file), 4, 5, endian='big')
    np.testing.assert_array_equal(out, data)


def test_read_img_wrong_size_raises(tmp_path):
    data = np.arange(20, dtype=np.float32)
    file = tmp_path / 'img.raw'
    data.tofile(str(file))
    with pytest.raises(ValueError):
        io_lib.read_img(str(file), 4, 4)


#%% ifg list
def test_read_ifg_list(tmp_path):
    file = tmp_path / 'ifg.txt'
    file.write_text('# comment line\n'
                    '20200101_20200113\n'
                    '20200113_20200125  trailing tokens ignored\n')
    assert io_lib.read_ifg_list(str(file)) \
        == ['20200101_20200113', '20200113_20200125']


#%% GAMMA par files
def test_get_param_par(tmp_path):
    par = tmp_path / 'slc.mli.par'
    synth.write_mli_par(par, width=123, length=456)
    assert io_lib.get_param_par(str(par), 'range_samples') == '123'
    assert io_lib.get_param_par(str(par), 'azimuth_lines') == '456'


def test_get_param_par_dem_par(tmp_path):
    par = tmp_path / 'EQA.dem_par'
    synth.write_dem_par(par)
    assert io_lib.get_param_par(str(par), 'corner_lat') == '34.0'
    assert float(io_lib.get_param_par(str(par), 'post_lat')) == -0.001


def test_get_param_par_first_match_wins(tmp_path):
    # Documents current behavior: grep substring matching returns the
    # first matching line (relevant e.g. for n_im vs n_im_all).
    par = tmp_path / 'params.txt'
    par.write_text('n_im_all: 10\nn_im: 7\n')
    assert io_lib.get_param_par(str(par), 'n_im') == '10'


#%% baselines files
def test_read_bperp_file_new_format(tmp_path):
    file = tmp_path / 'baselines'
    synth.write_baselines(file)
    bperp = io_lib.read_bperp_file(str(file), synth.IMDATES)
    np.testing.assert_allclose(bperp, synth.BPERP, atol=0.1)


def test_read_bperp_file_old_format(tmp_path):
    file = tmp_path / 'baselines'
    imdates = synth.IMDATES
    with open(file, 'w') as f:
        # num mdate sdate bp dt dt_m_sm dt_s_sm bp_m_sm bp_s_sm
        for i, imd in enumerate(imdates[1:], 1):
            print('{} {} {} {} {} {} {} {} {}'.format(
                i, imdates[0], imd, synth.BPERP[i], 12*i,
                0.0, 12.0*i, 0.0, synth.BPERP[i]), file=f)
    bperp = io_lib.read_bperp_file(str(file), imdates)
    np.testing.assert_allclose(bperp, synth.BPERP, atol=0.1)


def test_read_bperp_file_missing_date(tmp_path, capsys):
    # Documents current behavior: prints ERROR and returns False
    file = tmp_path / 'baselines'
    synth.write_baselines(file)
    result = io_lib.read_bperp_file(str(file), synth.IMDATES + ['20990101'])
    assert result is False
    assert 'ERROR' in capsys.readouterr().err


def test_make_bperp_file_roundtrip(tmp_path):
    file = tmp_path / 'baselines'
    # dict referenced to an arbitrary common reference; re-referenced to
    # the first imdate on write
    bperp_dict = {imd: bp + 100.0
                  for imd, bp in zip(synth.IMDATES, synth.BPERP)}
    io_lib.make_bperp_file(str(file), synth.IMDATES, bperp_dict)
    bperp = io_lib.read_bperp_file(str(file), synth.IMDATES)
    np.testing.assert_allclose(bperp, synth.BPERP, atol=0.1)


def test_make_bperp_file_missing_date_raises(tmp_path):
    bperp_dict = {imd: 0.0 for imd in synth.IMDATES[:-1]}
    with pytest.raises(ValueError, match=synth.IMDATES[-1]):
        io_lib.make_bperp_file(str(tmp_path / 'baselines'),
                               synth.IMDATES, bperp_dict)


def test_make_dummy_bperp(tmp_path):
    file = tmp_path / 'baselines'
    io_lib.make_dummy_bperp(str(file), synth.IMDATES)
    bperp = io_lib.read_bperp_file(str(file), synth.IMDATES)
    assert len(bperp) == len(synth.IMDATES)
    assert bperp[0] == 0.0
    assert np.all(np.abs(bperp) <= 1)  # dummy values are within -1..1


#%% Time series text
def test_make_tstxt(tmp_path):
    file = tmp_path / 'ts.txt'
    imdates = synth.IMDATES
    dt_cum = synth.dt_cum_years()
    vel_true = 25.0
    ts = (vel_true * dt_cum).astype(np.float32)
    gap = np.zeros(len(imdates) - 1)
    io_lib.make_tstxt(3, 5, imdates, ts, str(file), 1, 2, 3, 4, gap,
                      lat=34.0, lon=132.0)
    text = file.read_text()
    assert '# x, y    : 3, 5' in text
    assert '# ref     : 1:2/3:4' in text
    for imd in imdates:
        assert imd in text
    # linear model line contains the recovered slope
    model_line = [l for l in text.splitlines()
                  if l.startswith('# linear model:')][0]
    slope = float(model_line.split(':')[1].split('*')[0])
    assert slope == pytest.approx(vel_true, abs=0.01)


def test_make_point_kml(tmp_path):
    file = tmp_path / 'point.kml'
    io_lib.make_point_kml(34.5, 132.5, str(file))
    assert '<coordinates>132.5,34.5</coordinates>' in file.read_text()


#%% GeoTIFF (GDAL)
def test_geotiff_roundtrip(tmp_path):
    data = np.arange(80, dtype=np.float32).reshape(8, 10)
    file = tmp_path / 'test.tif'
    io_lib.make_geotiff(data, 34.0, 132.0, -0.001, 0.001, str(file),
                        ['COMPRESS=DEFLATE'], nodata=np.nan)
    out = io_lib.read_geotiff(str(file))
    np.testing.assert_array_equal(out, data)

    gt, proj, shape = io_lib.get_geotiff_info(str(file))
    assert shape == (8, 10)
    np.testing.assert_allclose(gt, (132.0, 0.001, 0, 34.0, 0, -0.001))
    assert '4326' in proj


def test_geotiff_uint8(tmp_path):
    data = np.arange(80, dtype=np.uint8).reshape(8, 10)
    file = tmp_path / 'test.tif'
    io_lib.make_geotiff(data, 34.0, 132.0, -0.001, 0.001, str(file), [])
    np.testing.assert_array_equal(io_lib.read_geotiff(str(file)), data)


def test_make_geotiff_unsupported_dtype(tmp_path):
    data = np.arange(80, dtype=np.int32).reshape(8, 10)
    with pytest.raises(ValueError, match='not supported'):
        io_lib.make_geotiff(data, 34.0, 132.0, -0.001, 0.001,
                            str(tmp_path / 'test.tif'), [])


def test_read_geotiff_with_ref(tmp_path):
    data = np.arange(80, dtype=np.float32).reshape(8, 10)
    file1 = tmp_path / 'a.tif'
    file2 = tmp_path / 'b.tif'
    io_lib.make_geotiff(data, 34.0, 132.0, -0.001, 0.001, str(file1), [])
    io_lib.make_geotiff(data * 2, 34.0, 132.0, -0.001, 0.001, str(file2), [])
    out = io_lib.read_geotiff(str(file1), file_ref=str(file2))
    np.testing.assert_array_equal(out, data)

    # Different geometry -> exception
    file3 = tmp_path / 'c.tif'
    io_lib.make_geotiff(data, 35.0, 132.0, -0.001, 0.001, str(file3), [])
    with pytest.raises(Exception, match='not identical'):
        io_lib.read_geotiff(str(file1), file_ref=str(file3))


def test_resample_geotiff(tmp_path):
    data = np.ones((8, 10), dtype=np.float32) * 7
    src = tmp_path / 'src.tif'
    dst = tmp_path / 'dst.tif'
    io_lib.make_geotiff(data, 34.0, 132.0, -0.001, 0.001, str(src), [])
    gt, proj, _ = io_lib.get_geotiff_info(str(src))
    new_gt = (132.0, 0.002, 0, 34.0, 0, -0.002)
    ok = io_lib.resample_geotiff(str(src), str(dst), new_gt, (4, 5), proj,
                                 resample_alg='near')
    assert ok is True
    out_gt, _, out_shape = io_lib.get_geotiff_info(str(dst))
    assert out_shape == (4, 5)
    np.testing.assert_allclose(out_gt, new_gt)
    np.testing.assert_array_equal(io_lib.read_geotiff(str(dst)),
                                  np.ones((4, 5), dtype=np.float32) * 7)
