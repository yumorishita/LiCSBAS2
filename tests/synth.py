"""
Builder of a small synthetic GEOCml-format dataset for testing.

The dataset is 10x10 pixels, geocoded (mli size == dem size), with 5 epochs
and 7 interferograms whose phases are exactly consistent (zero loop phase)
and follow a linear-in-time deformation with a known velocity field.
"""
import datetime as dt
from types import SimpleNamespace

import numpy as np

WIDTH = 10
LENGTH = 10
RADAR_FREQUENCY = 5405000000.0  # Hz (Sentinel-1 C-band)
SPEED_OF_LIGHT = 299792458.0  # m/s
WAVELENGTH = SPEED_OF_LIGHT / RADAR_FREQUENCY  # ~0.0555 m
COEF_R2M = -WAVELENGTH / 4 / np.pi * 1000  # rad -> mm (LOS)

IMDATES = ['20200101', '20200113', '20200125', '20200206', '20200218']
IFG_PAIRS = [(0, 1), (1, 2), (2, 3), (3, 4), (0, 2), (1, 3), (2, 4)]
IFGDATES = sorted('{}_{}'.format(IMDATES[i], IMDATES[j]) for i, j in IFG_PAIRS)

# bperp (m) of each epoch relative to the first
BPERP = [0.0, 30.0, -55.0, 12.0, 80.0]

GEOCML_CC = 180  # uint8 coherence of every pixel of the clean dataset


def vel_truth_mm():
    """True velocity field (mm/yr), strictly positive everywhere."""
    x, y = np.meshgrid(np.arange(WIDTH), np.arange(LENGTH))
    return (5.0 + 20.0 * (x + y) / 18).astype(np.float32)


def dt_cum_years(imdates=IMDATES):
    """Years since the first epoch for each epoch."""
    d0 = dt.datetime.strptime(imdates[0], '%Y%m%d')
    return np.array([(dt.datetime.strptime(imd, '%Y%m%d') - d0).days / 365.25
                     for imd in imdates])


def unw_phase(imd1, imd2, imdates=IMDATES):
    """Unwrapped phase (rad, float32) of the pair imd1_imd2."""
    dt_cum = dt_cum_years(imdates)
    t1 = dt_cum[imdates.index(imd1)]
    t2 = dt_cum[imdates.index(imd2)]
    dmm = vel_truth_mm() * (t2 - t1)
    return (dmm / COEF_R2M).astype(np.float32)


def write_mli_par(path, width=WIDTH, length=LENGTH,
                  radar_frequency=RADAR_FREQUENCY):
    with open(path, 'w') as f:
        print('range_samples:    {}'.format(width), file=f)
        print('azimuth_lines:    {}'.format(length), file=f)
        print('radar_frequency:  {} Hz'.format(radar_frequency), file=f)


def write_dem_par(path, width=WIDTH, length=LENGTH, corner_lat=34.0,
                  corner_lon=132.0, post_lat=-0.001, post_lon=0.001):
    with open(path, 'w') as f:
        print('width:          {}'.format(width), file=f)
        print('nlines:         {}'.format(length), file=f)
        print('corner_lat:     {}  decimal degrees'.format(corner_lat), file=f)
        print('corner_lon:     {}  decimal degrees'.format(corner_lon), file=f)
        print('post_lat:       {}  decimal degrees'.format(post_lat), file=f)
        print('post_lon:       {}  decimal degrees'.format(post_lon), file=f)
        print('ellipsoid_ra:   6378137.000  m', file=f)
        print('ellipsoid_reciprocal_flattening:  298.2572236', file=f)


def write_img(path, array):
    """Write an array as little-endian raw binary (GEOCml format)."""
    array.tofile(str(path))


def write_baselines(path, imdates=IMDATES, bperp=BPERP):
    """Write a baselines file in the new 4-column format."""
    d0 = dt.datetime.strptime(imdates[0], '%Y%m%d')
    with open(path, 'w') as f:
        for imd, bp in zip(imdates[1:], bperp[1:]):
            days = (dt.datetime.strptime(imd, '%Y%m%d') - d0).days
            print('{} {} {:6.1f} {:5.1f}'.format(imdates[0], imd, bp,
                                                 float(days)), file=f)


def build_geocml(workdir):
    """Create workdir/GEOCml1 and return the truth model."""
    geocdir = workdir / 'GEOCml1'
    geocdir.mkdir()

    write_mli_par(geocdir / 'slc.mli.par')
    write_dem_par(geocdir / 'EQA.dem_par')
    write_baselines(geocdir / 'baselines')

    # Tiny deterministic noise: with exactly zero loop phase, the ref
    # selection in step 12 masks rms==0 pixels as nodata and fails.
    # sigma = 0.0003 rad (~0.001 mm) keeps the truth check meaningful.
    rng = np.random.default_rng(42)

    cc = np.full((LENGTH, WIDTH), GEOCML_CC, dtype=np.uint8)
    for ifgd in IFGDATES:
        d = geocdir / ifgd
        d.mkdir()
        unw = unw_phase(ifgd[:8], ifgd[-8:]) \
            + rng.normal(0, 0.0003, (LENGTH, WIDTH)).astype(np.float32)
        write_img(d / (ifgd + '.unw'), unw)
        write_img(d / (ifgd + '.cc'), cc)
        (d / (ifgd + '.unw.png')).touch()  # existence only is checked

    return SimpleNamespace(
        workdir=workdir, geocdir=geocdir, width=WIDTH, length=LENGTH,
        imdates=IMDATES, ifgdates=IFGDATES, bperp=BPERP,
        vel_mm=vel_truth_mm(), dt_cum=dt_cum_years(),
        wavelength=WAVELENGTH, coef_r2m=COEF_R2M)


#%% Defective dataset
# A second dataset whose defects exercise the rejection paths of steps
# 11-13, which the clean dataset above never reaches.
#
#   - bottom SEA_ROWS rows are nodata (0) in every ifg -> nan handling
#   - LOW_COV_IFG    : almost no valid pixels     -> rejected by step 11 (-u)
#   - LOW_COH_IFG    : coherence below coh_thre   -> rejected by step 11 (-c)
#   - ISOLATED_IFG   : almost no valid pixels, and the only ifg of the last
#                      epoch, so that epoch drops out of the time series
#   - LOOP_ERR_IFG   : carries a 2pi unwrapping error -> rejected by step 12
#
# The nine good pairs form exactly four triangular loops (A-D) over the
# first six epochs:
#   A: 0_1, 1_2, 0_2    B: 1_2, 2_3, 1_3
#   C: 2_3, 3_4, 2_4    D: 3_4, 4_5, 3_5
# LOOP_ERR_IFG (1_3) belongs to loop B only, so loop B is the only one that
# fails and identify_bad_ifg() can pin the error on that single ifg. The
# network stays fully connected after it is removed, so the inverted
# velocity must still match the truth.

IMDATES_DEF = ['20200101', '20200113', '20200125', '20200206', '20200218',
               '20200301', '20200313']
BPERP_DEF = [0.0, 30.0, -55.0, 12.0, 80.0, -20.0, 45.0]

IFG_PAIRS_DEF_GOOD = [(0, 1), (1, 2), (0, 2), (2, 3), (1, 3),
                      (3, 4), (2, 4), (4, 5), (3, 5)]
IFG_PAIRS_DEF_BAD = [(0, 4), (1, 4), (5, 6)]


def _pair_def(i, j):
    return '{}_{}'.format(IMDATES_DEF[i], IMDATES_DEF[j])


IFGDATES_DEF = sorted(_pair_def(i, j) for i, j
                      in IFG_PAIRS_DEF_GOOD + IFG_PAIRS_DEF_BAD)

LOW_COV_IFG = _pair_def(0, 4)
LOW_COH_IFG = _pair_def(1, 4)
ISOLATED_IFG = _pair_def(5, 6)
ISOLATED_IMD = IMDATES_DEF[6]
LOOP_ERR_IFG = _pair_def(1, 3)

SEA_ROWS = 2          # bottom rows, nodata in every ifg
COV_BLOCK = 2         # side of the only valid block in the low coverage ifgs
LOOP_ERR_COLS = 4     # columns of land carrying the 2pi error
LOW_COH_CC = 5        # 5/255 = 0.02 < coh_thre (0.05)
GOOD_CC = 180


def build_geocml_defect(workdir):
    """Create workdir/GEOCml1 with defective ifgs; return the truth model."""
    geocdir = workdir / 'GEOCml1'
    geocdir.mkdir()

    write_mli_par(geocdir / 'slc.mli.par')
    write_dem_par(geocdir / 'EQA.dem_par')
    write_baselines(geocdir / 'baselines', IMDATES_DEF, BPERP_DEF)

    # Same tiny noise as build_geocml: exactly zero loop phase makes the
    # step 12 reference selection fail.
    rng = np.random.default_rng(42)

    land = np.ones((LENGTH, WIDTH), dtype=bool)
    land[LENGTH - SEA_ROWS:, :] = False

    for ifgd in IFGDATES_DEF:
        d = geocdir / ifgd
        d.mkdir()

        unw = unw_phase(ifgd[:8], ifgd[-8:], IMDATES_DEF) \
            + rng.normal(0, 0.0003, (LENGTH, WIDTH)).astype(np.float32)

        if ifgd == LOOP_ERR_IFG:
            unw[:LENGTH - SEA_ROWS, :LOOP_ERR_COLS] += np.float32(2 * np.pi)

        unw[~land] = 0  # nodata

        if ifgd in (LOW_COV_IFG, ISOLATED_IFG):
            keep = np.zeros_like(land)
            keep[:COV_BLOCK, :COV_BLOCK] = True
            unw[~keep] = 0

        # cc is nodata (0) wherever unw is, as in real data
        cc_val = LOW_COH_CC if ifgd == LOW_COH_IFG else GOOD_CC
        cc = np.where(unw == 0, 0, cc_val).astype(np.uint8)

        write_img(d / (ifgd + '.unw'), unw)
        write_img(d / (ifgd + '.cc'), cc)
        (d / (ifgd + '.unw.png')).touch()

    # Epochs surviving step 11: the last one loses its only ifg
    imdates_kept = IMDATES_DEF[:6]

    return SimpleNamespace(
        workdir=workdir, geocdir=geocdir, width=WIDTH, length=LENGTH,
        imdates=IMDATES_DEF, imdates_kept=imdates_kept,
        ifgdates=IFGDATES_DEF, bperp=BPERP_DEF,
        land=land, sea_rows=SEA_ROWS,
        bad_ifg11=sorted([LOW_COV_IFG, LOW_COH_IFG, ISOLATED_IFG]),
        bad_ifg12=[LOOP_ERR_IFG], isolated_imd=ISOLATED_IMD,
        vel_mm=vel_truth_mm(), dt_cum=dt_cum_years(imdates_kept),
        wavelength=WAVELENGTH, coef_r2m=COEF_R2M)


#%% GEOC dataset: GeoTIFF input of step 02
# 2x the size of the GEOCml dataset, so that multilooking by NLOOK must give
# back exactly unw_phase(). Each 2x2 block holds its GEOCml value plus a
# checkerboard of +/-CHECKER, which averages out over a full block but not
# over a partial one. That separates a correct nanmean from nodata being
# treated as 0 or a block being taken by its first pixel.
#
# Defects, all in GEOC_DEFECT_IFG (cc is nodata, 0, where unw is, as in
# real data):
#   - block (0, 0): 1 of 4 pixels valid         -> nan (< n_valid_thre 0.5)
#                                                  cc 0
#   - block (0, 1): the 2 diagonal +CHECKER pixels valid -> value + CHECKER
#                                                  cc GEOC_CC (not halved)
#   - block (0, 2): no pixel valid              -> nan, cc 0
#   - cc block (1, 0): 100, 101, 102, 103       -> 101 (mean 101.5, floored)
# and
#   - GEOC_FLOATCC_IFG has a float32 cc of 0-1 instead of uint8
#   - GEOC_NOCC_IFG has no cc tif at all       -> no_unw_list.txt

NLOOK = 2
GEOC_WIDTH = WIDTH * NLOOK
GEOC_LENGTH = LENGTH * NLOOK
# Pixel registration (GeoTIFF convention): outer edges of the frame
GEOC_LAT_N = 34.0
GEOC_LON_W = 132.0
GEOC_DLAT = -0.0005
GEOC_DLON = 0.0005

CHECKER = 0.1  # rad
GEOC_CC = 200
GEOC_FLOATCC = 0.5  # -> 127 after *255 and flooring
GEOC_HGT = 100.0
GEOC_MLI = 1000.0
GEOC_ENU = {'E': 0.6, 'N': -0.1, 'U': 0.8}
METADATA_FREQ = 5.40500045433e9
METADATA_CENTER_TIME = '09:26:43.500000'

GEOC_DEFECT_IFG = IFGDATES[0]
GEOC_FLOATCC_IFG = IFGDATES[1]
GEOC_NOCC_IFG = '{}_{}'.format(IMDATES[0], IMDATES[3])  # dates already used


def _checker(length, width):
    y, x = np.mgrid[0:length, 0:width]
    return np.where((x + y) % 2 == 0, CHECKER, -CHECKER).astype(np.float32)


def _upsample(a):
    return np.kron(a, np.ones((NLOOK, NLOOK), dtype=a.dtype))


def geoc_unw(ifgd):
    """Full resolution unw (rad, float32, 0 as nodata) of ifgd in GEOC."""
    unw = _upsample(unw_phase(ifgd[:8], ifgd[-8:])) \
        + _checker(GEOC_LENGTH, GEOC_WIDTH)
    if ifgd == GEOC_DEFECT_IFG:
        # block (0,0): keep pixel (0,0) only
        unw[0, 1] = unw[1, 0] = unw[1, 1] = 0
        # block (0,1): keep the +CHECKER diagonal, pixels (0,2) and (1,3)
        unw[0, 3] = unw[1, 2] = 0
        # block (0,2): nothing valid
        unw[0:2, 4:6] = 0
    return unw.astype(np.float32)


def geoc_unw_ml_expected(ifgd):
    """What step 02 -n NLOOK must write for ifgd (nan as nodata)."""
    unw = unw_phase(ifgd[:8], ifgd[-8:]).copy()
    if ifgd == GEOC_DEFECT_IFG:
        unw[0, 0] = np.nan
        unw[0, 1] += CHECKER
        unw[0, 2] = np.nan
    return unw


def geoc_cc(ifgd):
    """Full resolution cc as written to the GEOC tif."""
    if ifgd == GEOC_FLOATCC_IFG:
        return np.full((GEOC_LENGTH, GEOC_WIDTH), GEOC_FLOATCC, np.float32)
    cc = np.full((GEOC_LENGTH, GEOC_WIDTH), GEOC_CC, np.uint8)
    if ifgd == GEOC_DEFECT_IFG:
        cc[geoc_unw(ifgd) == 0] = 0
        cc[2:4, 0:2] = [[100, 101], [102, 103]]
    return cc


def geoc_cc_ml_expected(ifgd):
    """What step 02 -n NLOOK must write for the cc of ifgd (uint8)."""
    if ifgd == GEOC_FLOATCC_IFG:
        return np.full((LENGTH, WIDTH), int(GEOC_FLOATCC * 255), np.uint8)
    cc = np.full((LENGTH, WIDTH), GEOC_CC, np.uint8)
    if ifgd == GEOC_DEFECT_IFG:
        cc[0, 0] = cc[0, 2] = 0  # below n_valid_thre -> nan -> 0
        cc[1, 0] = 101
    return cc


def _write_geotiff(path, data):
    import LiCSBAS_io_lib as io_lib
    io_lib.make_geotiff(data, GEOC_LAT_N, GEOC_LON_W, GEOC_DLAT, GEOC_DLON,
                        str(path), [])


def build_geoc(workdir, metadata=True, baselines=True, metadata_freq=True):
    """Create workdir/GEOC (LiCSAR-like GeoTIFFs) and return its model.

    metadata_freq=False writes metadata.txt without radar_freq (as for the
    LiCSAR frames that lack it)."""
    geocdir = workdir / 'GEOC'
    geocdir.mkdir()

    for ifgd in IFGDATES + [GEOC_NOCC_IFG]:
        d = geocdir / ifgd
        d.mkdir()
        _write_geotiff(d / (ifgd + '.geo.unw.tif'), geoc_unw(ifgd))
        if ifgd != GEOC_NOCC_IFG:
            _write_geotiff(d / (ifgd + '.geo.cc.tif'), geoc_cc(ifgd))

    frame = 'test_frame'
    full = np.ones((GEOC_LENGTH, GEOC_WIDTH), np.float32)
    _write_geotiff(geocdir / (frame + '.geo.hgt.tif'), full * GEOC_HGT)
    _write_geotiff(geocdir / (frame + '.geo.mli.tif'), full * GEOC_MLI)
    for enu, v in GEOC_ENU.items():
        _write_geotiff(geocdir / '{}.geo.{}.tif'.format(frame, enu), full * v)

    if metadata:
        with open(geocdir / 'metadata.txt', 'w') as f:
            print('center_time={}'.format(METADATA_CENTER_TIME), file=f)
            if metadata_freq:
                print('radar_freq={}'.format(METADATA_FREQ), file=f)
    if baselines:
        write_baselines(geocdir / 'baselines')

    return SimpleNamespace(workdir=workdir, geocdir=geocdir,
                           ifgdates=IFGDATES + [GEOC_NOCC_IFG])


#%% GEOCml dataset for steps 04 and 05
# The clean GEOCml dataset plus what steps 04 and 05 act on:
#   - one column of low coherence in every ifg      -> step 04 -c
#   - one pixel of nodata (0) in every unw          -> must come out nan
#   - float files (hgt, slc.mli) with a ramp        -> step 05 clips them
#   - slc.mli.png and hgt.png                        -> step 05 recreates them

PREP_LOWCC_COL = WIDTH - 1
PREP_LOWCC = 30          # 30/255 = 0.12
PREP_NODATA_YX = (LENGTH - 1, 1)


def prep_hgt():
    """hgt (m) with a distinct value in every pixel, to check clipping."""
    return np.arange(LENGTH * WIDTH, dtype=np.float32).reshape(LENGTH, WIDTH)


def prep_mli():
    """slc.mli, distinct in every pixel and different from hgt."""
    return prep_hgt() + 1000


def build_geocml_prep(workdir):
    """build_geocml plus the features steps 04 and 05 act on.

    The returned model also holds the unw and cc written (unw_in, cc_in),
    so that tests take their expectations from memory and would notice a
    step modifying its input in place.
    """
    truth = build_geocml(workdir)
    geocdir = truth.geocdir
    truth.unw_in, truth.cc_in = {}, {}

    for ifgd in IFGDATES:
        unwfile = geocdir / ifgd / (ifgd + '.unw')
        unw = np.fromfile(str(unwfile), dtype=np.float32).reshape(LENGTH, WIDTH)
        unw[PREP_NODATA_YX] = 0
        write_img(unwfile, unw)
        truth.unw_in[ifgd] = unw

        ccfile = geocdir / ifgd / (ifgd + '.cc')
        cc = np.fromfile(str(ccfile), dtype=np.uint8).reshape(LENGTH, WIDTH)
        cc[:, PREP_LOWCC_COL] = PREP_LOWCC
        write_img(ccfile, cc)
        truth.cc_in[ifgd] = cc

    write_img(geocdir / 'hgt', prep_hgt())
    write_img(geocdir / 'slc.mli', prep_mli())
    (geocdir / 'hgt.png').touch()
    (geocdir / 'slc.mli.png').touch()

    return truth


#%% GACOS dataset for step 03
# The clean GEOCml dataset plus U.geo and GACOS/yyyymmdd.sltd.geo.tif.
# Each sltd (rad) is linear in lon and lat with a slope that differs per
# epoch, on a finer grid (GACOS_SUB per GEOCml pixel) extending GACOS_MARGIN
# GEOCml pixels beyond the frame. Resampling reproduces a linear field
# exactly, so the sltd step 03 writes must equal gacos_sltd() at the
# GEOCml pixel centers; a grid shifted or shrunk by a fraction of a pixel
# shows up as an error of slope x shift.

GACOS_SUB = 2
GACOS_MARGIN = 5
LOS_U = 0.8


def gacos_sltd(imd, lon, lat):
    """sltd (rad) of epoch imd at lon/lat, linear; never 0 (nodata)."""
    e = IMDATES.index(imd)
    x = (lon - 132.0) / 0.001      # GEOCml pixel coordinates (synth dem_par)
    y = (34.0 - lat) / 0.001
    return 100.0 + e + 0.2 * e * x + 0.1 * e * y


def gacos_sltd_expected(imd):
    """sltd at the GEOCml pixel centers, as step 03 must write it."""
    lon, lat = np.meshgrid(132.0 + 0.001 * np.arange(WIDTH),
                           34.0 - 0.001 * np.arange(LENGTH))
    return gacos_sltd(imd, lon, lat).astype(np.float32)


def build_geocml_gacos(workdir):
    """build_geocml plus U.geo and GACOS/*.sltd.geo.tif."""
    import LiCSBAS_io_lib as io_lib
    truth = build_geocml(workdir)
    write_img(truth.geocdir / 'U.geo',
              np.full((LENGTH, WIDTH), LOS_U, dtype=np.float32))

    gacosdir = workdir / 'GACOS'
    gacosdir.mkdir()
    d = 0.001 / GACOS_SUB
    n = GACOS_SUB * (WIDTH + 2 * GACOS_MARGIN)
    m = GACOS_SUB * (LENGTH + 2 * GACOS_MARGIN)
    lon_w = 132.0 - 0.001 / 2 - GACOS_MARGIN * 0.001   # outer edges
    lat_n = 34.0 + 0.001 / 2 + GACOS_MARGIN * 0.001
    lon, lat = np.meshgrid(lon_w + (np.arange(n) + 0.5) * d,
                           lat_n - (np.arange(m) + 0.5) * d)
    for imd in IMDATES:
        io_lib.make_geotiff(gacos_sltd(imd, lon, lat).astype(np.float32),
                            lat_n, lon_w, -d, d,
                            str(gacosdir / (imd + '.sltd.geo.tif')), [])

    truth.gacosdir = gacosdir
    return truth
