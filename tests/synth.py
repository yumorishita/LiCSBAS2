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

    cc = np.full((LENGTH, WIDTH), 180, dtype=np.uint8)
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
