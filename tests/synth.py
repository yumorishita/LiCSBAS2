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


def dt_cum_years():
    """Years since the first epoch for each epoch."""
    d0 = dt.datetime.strptime(IMDATES[0], '%Y%m%d')
    return np.array([(dt.datetime.strptime(imd, '%Y%m%d') - d0).days / 365.25
                     for imd in IMDATES])


def unw_phase(imd1, imd2):
    """Unwrapped phase (rad, float32) of the pair imd1_imd2."""
    dt_cum = dt_cum_years()
    t1 = dt_cum[IMDATES.index(imd1)]
    t2 = dt_cum[IMDATES.index(imd2)]
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
