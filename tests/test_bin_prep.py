"""End-to-end tests of the preparation steps 02, 04 and 05 on synthetic data.

Step 02 converts LiCSAR-like GeoTIFFs (GEOC) to the raw GEOCml format,
optionally multilooking them; steps 04 (mask) and 05 (clip) rewrite
GEOCml. All three change the data every later step works on, so their
outputs are checked value by value against a model of the input
(tests/synth.py: build_geoc, build_geocml_prep), not only for existence.
"""
import filecmp
import os
import re

import numpy as np
import pytest

import LiCSBAS_io_lib as io_lib
import synth

pytestmark = pytest.mark.smoke

L, W = synth.LENGTH, synth.WIDTH                 # GEOCml2 and steps 04/05
GL, GW = synth.GEOC_LENGTH, synth.GEOC_WIDTH     # GEOCml1 from step 02


def read_unw(d, ifgd, length, width):
    return io_lib.read_img(str(d / ifgd / (ifgd + '.unw')), length, width)


def read_cc(d, ifgd, length, width):
    return io_lib.read_img(str(d / ifgd / (ifgd + '.cc')), length, width,
                           np.uint8)


def par(file, field):
    return io_lib.get_param_par(str(file), field)


#%% Step 02: GeoTIFF -> GEOCml
def test_step02_nlook1_unw_cc(geocml1_02):
    """Without multilooking, unw and cc are the tifs as they are, with 0
    (nodata) in unw turned into nan and a float cc of 0-1 into uint8."""
    for ifgd in synth.IFGDATES:
        unw_exp = synth.geoc_unw(ifgd)
        unw_exp[unw_exp == 0] = np.nan
        np.testing.assert_array_equal(read_unw(geocml1_02, ifgd, GL, GW),
                                      unw_exp)

        cc_exp = synth.geoc_cc(ifgd)
        if cc_exp.dtype == np.float32:
            cc_exp = np.floor(cc_exp * 255).astype(np.uint8)
        np.testing.assert_array_equal(read_cc(geocml1_02, ifgd, GL, GW),
                                      cc_exp)


def test_step02_nlook2_multilook(geocml2_02):
    """Multilooking must average valid pixels only, give nan where fewer
    than half of a block is valid, and floor the averaged cc."""
    for ifgd in synth.IFGDATES:
        np.testing.assert_allclose(read_unw(geocml2_02, ifgd, L, W),
                                   synth.geoc_unw_ml_expected(ifgd),
                                   atol=1e-6)  # nan positions must match too
        np.testing.assert_array_equal(read_cc(geocml2_02, ifgd, L, W),
                                      synth.geoc_cc_ml_expected(ifgd))


@pytest.mark.parametrize('fixture, length, width', [
    ('geocml1_02', GL, GW), ('geocml2_02', L, W)])
def test_step02_auxiliary_files(fixture, length, width, request):
    d = request.getfixturevalue(fixture)
    expected = {'hgt': synth.GEOC_HGT, 'slc.mli': synth.GEOC_MLI}
    expected.update({enu + '.geo': v for enu, v in synth.GEOC_ENU.items()})
    for name, value in expected.items():
        a = io_lib.read_img(str(d / name), length, width)
        np.testing.assert_allclose(a, value, rtol=1e-6, err_msg=name)
    assert (d / 'hgt.png').stat().st_size > 0
    assert (d / 'slc.mli.png').stat().st_size > 0


def test_step02_par_files_nlook1(geocml1_02):
    mlipar = geocml1_02 / 'slc.mli.par'
    assert int(par(mlipar, 'range_samples')) == GW
    assert int(par(mlipar, 'azimuth_lines')) == GL
    assert float(par(mlipar, 'radar_frequency')) == pytest.approx(
        synth.METADATA_FREQ)  # from metadata.txt
    assert par(mlipar, 'center_time') == synth.METADATA_CENTER_TIME

    dempar = geocml1_02 / 'EQA.dem_par'
    assert int(par(dempar, 'width')) == GW
    assert int(par(dempar, 'nlines')) == GL
    assert float(par(dempar, 'post_lat')) == pytest.approx(synth.GEOC_DLAT)
    assert float(par(dempar, 'post_lon')) == pytest.approx(synth.GEOC_DLON)
    # EQA.dem_par is in grid registration: center of the first pixel
    assert float(par(dempar, 'corner_lat')) == pytest.approx(
        synth.GEOC_LAT_N + synth.GEOC_DLAT / 2)
    assert float(par(dempar, 'corner_lon')) == pytest.approx(
        synth.GEOC_LON_W + synth.GEOC_DLON / 2)


def test_step02_par_files_nlook2_size(geocml2_02):
    assert int(par(geocml2_02 / 'slc.mli.par', 'range_samples')) == W
    assert int(par(geocml2_02 / 'slc.mli.par', 'azimuth_lines')) == L
    dempar = geocml2_02 / 'EQA.dem_par'
    assert int(par(dempar, 'width')) == W
    assert int(par(dempar, 'nlines')) == L
    assert float(par(dempar, 'post_lat')) == pytest.approx(
        synth.GEOC_DLAT * synth.NLOOK)
    assert float(par(dempar, 'post_lon')) == pytest.approx(
        synth.GEOC_DLON * synth.NLOOK)


@pytest.mark.xfail(strict=True, reason='#174: step 02 shifts the pixel-'
                                       'registration corner by half an '
                                       'ORIGINAL pixel even when multilooking, '
                                       'instead of half a multilooked pixel, '
                                       'so the grid is off by (nlook-1)/2 '
                                       'original pixels')
def test_step02_corner_nlook2(geocml2_02):
    """The first multilooked pixel averages the first NLOOK x NLOOK
    original pixels, so its center is NLOOK/2 original pixels inside the
    frame edge."""
    dempar = geocml2_02 / 'EQA.dem_par'
    assert float(par(dempar, 'corner_lat')) == pytest.approx(
        synth.GEOC_LAT_N + synth.GEOC_DLAT * synth.NLOOK / 2)
    assert float(par(dempar, 'corner_lon')) == pytest.approx(
        synth.GEOC_LON_W + synth.GEOC_DLON * synth.NLOOK / 2)


@pytest.mark.parametrize('fixture', ['geocml1_02', 'geocml2_02'])
def test_step02_missing_cc(fixture, request):
    """An ifg without cc tif is skipped and listed, not half converted."""
    d = request.getfixturevalue(fixture)
    assert (d / 'no_unw_list.txt').read_text().split() == [synth.GEOC_NOCC_IFG]
    assert not (d / synth.GEOC_NOCC_IFG).exists()
    for ifgd in synth.IFGDATES:
        assert (d / ifgd / (ifgd + '.unw')).exists()
        assert (d / ifgd / (ifgd + '.unw.png')).stat().st_size > 0


def test_step02_baselines_copied(geoc, geocml1_02):
    assert filecmp.cmp(str(geoc.geocdir / 'baselines'),
                       str(geocml1_02 / 'baselines'), shallow=False)


def test_step02_dummy_baselines(geocml1_02_nometa):
    """Without baselines, a dummy one covering all epochs is made."""
    bperp = io_lib.read_bperp_file(str(geocml1_02_nometa / 'baselines'),
                                   synth.IMDATES)
    assert len(bperp) == len(synth.IMDATES)


@pytest.mark.xfail(strict=True, reason='#173: the default radar_freq '
                                       '(5.405e9) is assigned after the '
                                       'options are parsed, so --freq is '
                                       'always overwritten')
def test_step02_freq_option(geocml1_02_nometa):
    radar_freq = float(par(geocml1_02_nometa / 'slc.mli.par',
                           'radar_frequency'))
    assert radar_freq == pytest.approx(1.27e9)


def test_step02_rerun_skips_existing(geoc, geocml1_02, run_script):
    """Existing outputs are not recreated; the ifg without cc is retried
    and listed again."""
    unwfiles = [geocml1_02 / d / (d + '.unw') for d in synth.IFGDATES]
    mtimes = [os.stat(f).st_mtime_ns for f in unwfiles]
    parfiles = [geocml1_02 / 'slc.mli.par', geocml1_02 / 'EQA.dem_par']
    pars = [f.read_text() for f in parfiles]

    res = run_script('LiCSBAS02_ml_prep.py', '-i', 'GEOC', '-n', '1',
                     '--n_para', '1', cwd=geoc.workdir)

    assert re.search(r'\b{}/\s*{} unw and cc already exist'.format(
        len(synth.IFGDATES), len(geoc.ifgdates)), res.stdout)
    assert [os.stat(f).st_mtime_ns for f in unwfiles] == mtimes
    assert [f.read_text() for f in parfiles] == pars
    assert (geocml1_02 / 'no_unw_list.txt').read_text().split() == \
        [synth.GEOC_NOCC_IFG]


#%% Step 04: mask
def prep_unw(geocml_prep, ifgd):
    """Input unw of steps 04/05 with nodata as nan, as they must read it."""
    unw = read_unw(geocml_prep.geocdir, ifgd, L, W)
    unw[unw == 0] = np.nan
    return unw


def test_step04_mask_by_coherence_and_range(geocml_prep, mask04):
    """-c masks the low coherence column and -r the given range; nodata
    stays nan and everything else is untouched."""
    bool_mask = np.zeros((L, W), dtype=bool)
    bool_mask[:, synth.PREP_LOWCC_COL] = True   # -c 0.2 (30/255 = 0.12)
    bool_mask[0:3, 0:2] = True                  # -r 0:2/0:3

    mask = io_lib.read_img(str(mask04 / 'mask'), L, W)
    np.testing.assert_array_equal(mask, np.float32(~bool_mask))

    for ifgd in synth.IFGDATES:
        unw_exp = prep_unw(geocml_prep, ifgd)
        unw_exp[bool_mask] = np.nan
        np.testing.assert_array_equal(read_unw(mask04, ifgd, L, W), unw_exp)

    # The nodata pixel is outside the mask but must still be nan
    assert not bool_mask[synth.PREP_NODATA_YX]


def test_step04_coh_avg(mask04):
    coh_avg = io_lib.read_img(str(mask04 / 'coh_avg'), L, W)
    exp = np.full((L, W), 180 / 255, dtype=np.float32)
    exp[:, synth.PREP_LOWCC_COL] = synth.PREP_LOWCC / 255
    np.testing.assert_allclose(coh_avg, exp, atol=1e-6)


def test_step04_cc_and_other_files(geocml_prep, mask04):
    """cc is linked (or copied) unchanged, and other files are copied."""
    for ifgd in synth.IFGDATES:
        assert filecmp.cmp(str(geocml_prep.geocdir / ifgd / (ifgd + '.cc')),
                           str(mask04 / ifgd / (ifgd + '.cc')), shallow=False)
    for name in ('slc.mli.par', 'EQA.dem_par', 'baselines', 'hgt', 'slc.mli'):
        assert filecmp.cmp(str(geocml_prep.geocdir / name),
                           str(mask04 / name), shallow=False), name


def test_step04_mask_by_range_file(geocml_prep, mask04_file):
    bool_mask = np.zeros((L, W), dtype=bool)
    bool_mask[4:6, 3:5] = True   # 3:5/4:6
    bool_mask[0:2, 6:8] = True   # 6:8/0:2

    mask = io_lib.read_img(str(mask04_file / 'mask'), L, W)
    np.testing.assert_array_equal(mask, np.float32(~bool_mask))

    for ifgd in synth.IFGDATES:
        unw_exp = prep_unw(geocml_prep, ifgd)
        unw_exp[bool_mask] = np.nan
        np.testing.assert_array_equal(read_unw(mask04_file, ifgd, L, W),
                                      unw_exp)
    assert not (mask04_file / 'coh_avg').exists()  # -c not given


#%% Step 05: clip
X1, X2, Y1, Y2 = 0, 7, 1, 10   # clip05: -r 0:7/1:10


def test_step05_clip_unw_cc(geocml_prep, clip05):
    lc, wc = Y2 - Y1, X2 - X1
    for ifgd in synth.IFGDATES:
        unw_exp = prep_unw(geocml_prep, ifgd)[Y1:Y2, X1:X2]
        np.testing.assert_array_equal(read_unw(clip05, ifgd, lc, wc), unw_exp)

        cc_exp = read_cc(geocml_prep.geocdir, ifgd, L, W)[Y1:Y2, X1:X2]
        np.testing.assert_array_equal(read_cc(clip05, ifgd, lc, wc), cc_exp)
        assert (clip05 / ifgd / (ifgd + '.unw.png')).stat().st_size > 0

    # The nodata pixel is inside the clipped area and must be nan there
    y, x = synth.PREP_NODATA_YX
    assert np.isnan(read_unw(clip05, synth.IFGDATES[0], lc, wc)[y - Y1, x - X1])


def test_step05_clip_float_files(clip05):
    lc, wc = Y2 - Y1, X2 - X1
    hgt = io_lib.read_img(str(clip05 / 'hgt'), lc, wc)
    np.testing.assert_array_equal(hgt, synth.prep_hgt()[Y1:Y2, X1:X2])
    mli = io_lib.read_img(str(clip05 / 'slc.mli'), lc, wc)
    np.testing.assert_array_equal(mli, synth.prep_hgt()[Y1:Y2, X1:X2] + 1000)
    # pngs are recreated from the clipped data, not copied (inputs are empty)
    assert (clip05 / 'hgt.png').stat().st_size > 0
    assert (clip05 / 'slc.mli.png').stat().st_size > 0


def test_step05_par_files(geocml_prep, clip05):
    lc, wc = Y2 - Y1, X2 - X1
    mlipar_in = geocml_prep.geocdir / 'slc.mli.par'
    mlipar = clip05 / 'slc.mli.par'
    assert int(par(mlipar, 'range_samples')) == wc
    assert int(par(mlipar, 'azimuth_lines')) == lc
    assert par(mlipar, 'radar_frequency') == par(mlipar_in, 'radar_frequency')

    dempar_in = geocml_prep.geocdir / 'EQA.dem_par'
    dempar = clip05 / 'EQA.dem_par'
    assert int(par(dempar, 'width')) == wc
    assert int(par(dempar, 'nlines')) == lc
    post_lat = float(par(dempar_in, 'post_lat'))
    post_lon = float(par(dempar_in, 'post_lon'))
    assert float(par(dempar, 'post_lat')) == post_lat
    assert float(par(dempar, 'post_lon')) == post_lon
    assert float(par(dempar, 'corner_lat')) == pytest.approx(
        float(par(dempar_in, 'corner_lat')) + post_lat * Y1)
    assert float(par(dempar, 'corner_lon')) == pytest.approx(
        float(par(dempar_in, 'corner_lon')) + post_lon * X1)

    assert (clip05 / 'cliparea.txt').read_text() == \
        '{}:{}/{}:{}'.format(X1, X2, Y1, Y2)
    assert filecmp.cmp(str(geocml_prep.geocdir / 'baselines'),
                       str(clip05 / 'baselines'), shallow=False)


def test_step05_geo_range_equals_index_range(clip05, clip05_geo):
    """-g with the lon/lat of the -r area must give identical output."""
    cmp = filecmp.dircmp(str(clip05), str(clip05_geo))
    pngs = lambda files: [f for f in files if not f.endswith('.png')]
    assert not cmp.left_only and not cmp.right_only
    assert pngs(cmp.diff_files) == []
    for ifgd in synth.IFGDATES:
        for ext in ('.unw', '.cc'):
            f = os.path.join(ifgd, ifgd + ext)
            assert filecmp.cmp(str(clip05 / f), str(clip05_geo / f),
                               shallow=False), f
