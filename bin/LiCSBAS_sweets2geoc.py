#!/usr/bin/env python3
"""
Convert a sweets/dolphin output directory into a LiCSBAS-style `GEOC/` layout.

Inputs (sweets work directory, or the dolphin directory itself):
- dolphin/
    - unwrapped/{yyyymmdd_yyyymmdd}.unw.tif          -> unw (radians)
    - interferograms/{yyyymmdd_yyyymmdd}.int.cor.tif -> cc (0-1)
- dolphin_config.yaml (wavelength; cslc_file_list for center_time [NISAR])
- geometry/height.tif (or dem.tif)                   -> hgt
- geometry/los_east.tif, los_north.tif (S1)          -> ENU

Outputs:
- GEOC/{yyyymmdd_yyyymmdd}/{yyyymmdd_yyyymmdd}.geo.unw.tif
- GEOC/{yyyymmdd_yyyymmdd}/{yyyymmdd_yyyymmdd}.geo.cc.tif
- GEOC/X.geo.hgt.tif (if geometry/height.tif or dem.tif found)
- GEOC/X.geo.E.tif, X.geo.N.tif, X.geo.U.tif (if geometry files found)
- GEOC/metadata.txt (center_time if known, and radar_freq)
- GEOC/baselines (bperp fetched from the ASF baseline API; Sentinel-1 only)

Notes:
- All outputs are reprojected from the dolphin projection (usually UTM) to
  EPSG:4326 on a single common grid.
- The dolphin unw phase sign convention is the same as LiCSBAS (positive
  phase = motion away from satellite); no sign flip, unlike ARIA products.
- baselines are fetched from the ASF baseline API (no authentication needed);
  on any failure (NISAR, no network, incomplete stack) no file is created
  and LiCSBAS02 generates a dummy bperp instead.
- Not converted (candidates for future options): wrapped int, conncomp,
  temporal_coherence/ps_mask masking, timeseries/velocity products.

"""

import argparse
import datetime
import glob
import multiprocessing as mp
import os
import re
import shutil
import sys
import tempfile
import time

from osgeo import gdal
import numpy as np


def write_geotiff(path, arr, gt, proj, dtype=gdal.GDT_Float32, nodata=0.0, options=None):
    """Write a GeoTIFF via a local temporary file.

    libtiff creation involves random-access writes which can silently corrupt
    files on shared/network file systems (e.g. vboxsf); moving a finished file
    is a sequential copy and safe.
    """
    if options is None:
        options = ['COMPRESS=DEFLATE', 'TILED=YES', 'PREDICTOR=2']
    driver = gdal.GetDriverByName('GTiff')
    ny, nx = arr.shape
    with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmpf:
        tmp_path = tmpf.name
    try:
        outds = driver.Create(tmp_path, nx, ny, 1, dtype, options=options)
        outds.SetGeoTransform(gt)
        outds.SetProjection(proj)
        band = outds.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        band.WriteArray(arr)
        outds.FlushCache()
        outds = None
        shutil.move(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    # make failures loud rather than leaving a silently broken file
    ds_check = gdal.Open(path)
    if ds_check is None:
        raise IOError(f'Failed to write a valid GeoTIFF: {path}')
    ds_check = None


def find_dirs(indir):
    """Return (workdir, dolphin_dir) from indir.

    Accepts either a sweets work directory containing dolphin/ or a dolphin
    directory itself containing unwrapped/.
    """
    indir = os.path.abspath(indir)
    if os.path.isdir(os.path.join(indir, 'unwrapped')):
        return os.path.dirname(indir), indir
    if os.path.isdir(os.path.join(indir, 'dolphin', 'unwrapped')):
        return indir, os.path.join(indir, 'dolphin')
    raise FileNotFoundError(
        f"{indir} contains neither dolphin/unwrapped/ (sweets work dir) "
        f"nor unwrapped/ (dolphin dir)")


def read_dolphin_config(workdir, dolphin_dir):
    """Return (wavelength [m] or None, list of cslc file names) from dolphin_config.yaml."""
    config = None
    for d in [workdir, dolphin_dir]:
        path = os.path.join(d, 'dolphin_config.yaml')
        if os.path.exists(path):
            config = path
            break
    if config is None:
        print('No dolphin_config.yaml found; wavelength/center_time unavailable')
        return None, []
    print(f'Config: {config}')

    with open(config) as f:
        txt = f.read()

    try:
        import yaml
        cfg = yaml.safe_load(txt)
        wavelength = (cfg.get('input_options') or {}).get('wavelength')
        cslc_names = cfg.get('cslc_file_list') or []
    except ImportError:
        print('pyyaml not available; parsing config with regex')
        m = re.search(r'^\s*wavelength:\s*([\d.eE+-]+)', txt, re.M)
        wavelength = m.group(1) if m else None
        cslc_names = []
        in_list = False
        for line in txt.splitlines():
            if re.match(r'^cslc_file_list:', line):
                in_list = True
            elif in_list:
                m2 = re.match(r'^\s*-\s+(\S+)', line)
                if m2:
                    cslc_names.append(m2.group(1))
                elif line and line[0] not in ' \t#':  # next top-level key
                    break

    if wavelength is not None:
        wavelength = float(wavelength)

    return wavelength, [os.path.basename(str(name)) for name in cslc_names]


def get_center_time(cslc_names):
    """Return 'HH:MM:SS' averaged over the acquisition timestamps in the first
    cslc file name, or None if unavailable (e.g. no cslc_file_list for S1)."""
    if not cslc_names:
        return None
    matches = re.findall(r'(\d{8})T(\d{6})', cslc_names[0])
    if not matches:
        return None
    # NISAR names contain start and end time; average those on the acquisition
    # date and ignore other timestamps (e.g. production time on another date)
    date0 = matches[0][0]
    secs = [int(t[0:2]) * 3600 + int(t[2:4]) * 60 + int(t[4:6])
            for d, t in matches if d == date0]
    sec = round(sum(secs) / len(secs))
    return f'{sec // 3600:02d}:{sec % 3600 // 60:02d}:{sec % 60:02d}'


def find_pairs(dolphin_dir):
    """Return sorted list of (pair, unw_path, cor_path or None)."""
    pairs = []
    for unw in sorted(glob.glob(os.path.join(dolphin_dir, 'unwrapped', '*.unw.tif'))):
        m = re.match(r'^(\d{8})_(\d{8})\.unw\.tif$', os.path.basename(unw))
        if not m:
            continue
        pair = '_'.join(sorted(m.groups()))  # earlier date first
        cor = os.path.join(dolphin_dir, 'interferograms',
                           f'{m.group(1)}_{m.group(2)}.int.cor.tif')
        pairs.append((pair, unw, cor if os.path.exists(cor) else None))
    return pairs


def get_target_grid(unw_path):
    """Determine the common output grid in EPSG:4326 from the first unw."""
    # Note: gdal<3.7 return xres=yres
    # TODO upgrade gdal to >=3.7 so that the resolution of the source file(s) is preserved
    ds = gdal.Warp('', unw_path, format='MEM', dstSRS='EPSG:4326')
    gt = ds.GetGeoTransform()
    width = ds.RasterXSize
    length = ds.RasterYSize
    proj = ds.GetProjection()
    ds = None
    bounds = (gt[0], gt[3] + length * gt[5], gt[0] + width * gt[1], gt[3])
    return {'gt': gt, 'width': width, 'length': length, 'proj': proj,
            'bounds': bounds}


def warp_to_grid(src_path, grid, resample_alg):
    """Reproject a raster onto the common EPSG:4326 grid; return float32 array.

    outputBounds/width/height are set explicitly so that every output shares
    the identical grid even if the sources are gridded differently (e.g. dem).
    """
    ds = gdal.Warp('', src_path, format='MEM', dstSRS='EPSG:4326',
                   outputBounds=grid['bounds'], width=grid['width'],
                   height=grid['length'], resampleAlg=resample_alg)
    arr = ds.ReadAsArray().astype(np.float32)
    ds = None
    return arr


def process_pair(args):
    """Convert one pair to GEOC unw and cc tifs (parallel worker)."""
    pair, unw_path, cor_path, out_geoc, grid, cc_thresh, resample_alg, \
        overwrite = args

    tgt_dir = os.path.join(out_geoc, pair)
    tgt_unw = os.path.join(tgt_dir, f'{pair}.geo.unw.tif')
    tgt_cc = os.path.join(tgt_dir, f'{pair}.geo.cc.tif')
    if not overwrite and os.path.exists(tgt_unw) and os.path.exists(tgt_cc):
        print(f'  {pair} already exists. Skipping', flush=True)
        return 'skip'
    if cor_path is None:
        print(f'  WARN: no int.cor.tif for {pair}. Skipping', flush=True)
        return 'no_cc'

    try:
        print(f'Processing {pair}...', flush=True)
        unw = warp_to_grid(unw_path, grid, resample_alg)
        cc = warp_to_grid(cor_path, grid, gdal.GRA_NearestNeighbour)

        # scale cc if 0-255
        if np.nanmax(cc) > 1.5:
            cc = cc / 255.0

        unw[~np.isfinite(unw)] = 0.0
        cc[~np.isfinite(cc)] = 0.0

        # mask unw by coherence (nodata cc=0 is always < cc_thresh)
        unw[cc < cc_thresh] = 0.0

        os.makedirs(tgt_dir, exist_ok=True)
        write_geotiff(tgt_unw, unw, grid['gt'], grid['proj'])
        write_geotiff(tgt_cc, cc, grid['gt'], grid['proj'])
    except Exception as e:
        # do not let one broken pair kill the whole run
        print(f'  ERROR: {pair} failed: {e}', flush=True)
        return 'fail'

    return 'ok'


def find_geometry_dir(workdir, dolphin_dir):
    """Return the directory containing los_east.tif/los_north.tif, or None.

    sweets creates geometry/ in the work directory for Sentinel-1 (from OPERA
    CSLC static layers); no geometry is available for NISAR at the moment.
    """
    for d in [workdir, dolphin_dir]:
        for sub in ['geometry', '.']:
            geom_dir = os.path.normpath(os.path.join(d, sub))
            if (os.path.exists(os.path.join(geom_dir, 'los_east.tif')) and
                    os.path.exists(os.path.join(geom_dir, 'los_north.tif'))):
                return geom_dir
    return None


def make_hgt(dem_path, grid, out_geoc, frameid='X'):
    hgt = warp_to_grid(dem_path, grid, gdal.GRA_Bilinear)
    hgt[~np.isfinite(hgt)] = 0.0
    write_geotiff(os.path.join(out_geoc, f'{frameid}.geo.hgt.tif'), hgt,
                  grid['gt'], grid['proj'])
    print(f'Wrote {frameid}.geo.hgt.tif from {dem_path}')


def make_enu(geom_dir, grid, out_geoc, frameid='X'):
    """Create ENU LOS unit vector tifs from sweets geometry files.

    los_east/los_north are the E/N components of the LOS unit vector from
    target to satellite, the same convention as LiCSBAS *.geo.E/N/U tifs
    (verified against LiCSAR frame 046D over the same area); U is recovered
    as sqrt(1 - E^2 - N^2) since the vector always points above the horizon.
    """
    if geom_dir is None:
        print('No geometry files (los_east/los_north) found; skip ENU creation')
        return

    print(f'Geometry files found in: {geom_dir}')
    E = warp_to_grid(os.path.join(geom_dir, 'los_east.tif'), grid,
                     gdal.GRA_Bilinear)
    N = warp_to_grid(os.path.join(geom_dir, 'los_north.tif'), grid,
                     gdal.GRA_Bilinear)
    invalid = ~np.isfinite(E) | ~np.isfinite(N) | ((E == 0) & (N == 0))
    E[invalid] = 0.0
    N[invalid] = 0.0
    U = np.sqrt(np.clip(1.0 - E**2 - N**2, 0.0, 1.0))
    U[invalid] = 0.0

    for name, arr in zip(['E', 'N', 'U'], [E, N, U]):
        write_geotiff(os.path.join(out_geoc, f'{frameid}.geo.{name}.tif'),
                      arr, grid['gt'], grid['proj'])
    print(f'Wrote {frameid}.geo.E/N/U.tif')


def write_metadata(out_geoc, wavelength, center_time):
    if wavelength and wavelength > 0:
        radar_freq = 2.99792458e8 / wavelength
    else:
        radar_freq = 5.405e9
        print('WARN: wavelength unknown; using default Sentinel-1 radar_freq')

    metadata_path = os.path.join(out_geoc, 'metadata.txt')
    with open(metadata_path, 'w') as f:
        if center_time:
            f.write(f'center_time={center_time}\n')
        f.write(f'radar_freq={radar_freq}\n')
    print(f'Wrote metadata: {metadata_path}')


def make_baselines(out_geoc, pair_names, bounds, center_time, wavelength,
                   cslc_names, overwrite=False):
    """Create GEOC/baselines with bperp from the ASF baseline API (S1 only).

    Never raises: any failure prints a warning and returns without a file,
    in which case LiCSBAS02 generates a dummy bperp.
    """
    bperp_file = os.path.join(out_geoc, 'baselines')
    if os.path.exists(bperp_file) and not overwrite:
        print('baselines already exists. Skipping')
        return

    if any('NISAR' in name for name in cslc_names):
        print('NISAR data: ASF baseline API not supported; '
              'LiCSBAS02 will generate dummy bperp')
        return
    if wavelength and not 0.05 < wavelength < 0.06:
        print(f'Non-S1 wavelength ({wavelength} m): ASF baseline API not '
              'supported; LiCSBAS02 will generate dummy bperp')
        return

    try:
        import LiCSBAS_tools_lib as tools_lib
        import LiCSBAS_io_lib as io_lib
    except ImportError:
        print('WARN: LiCSBAS_lib not found in PYTHONPATH; skip baselines; '
              'LiCSBAS02 will generate dummy bperp')
        return

    try:
        imdates = tools_lib.ifgdates2imdates(pair_names)
        clon = (bounds[0] + bounds[2]) / 2
        clat = (bounds[1] + bounds[3]) / 2

        # track (path) number from OPERA CSLC burst names if available
        m = re.search(r'_T(\d{3})-\d{6}-IW\d', cslc_names[0]) if cslc_names else None
        path = int(m.group(1)) if m else None

        print(f'\nGetting baseline info from ASF for POINT({clon:.3f} {clat:.3f})...')
        bperp_dict = tools_lib.get_bperp_asf(clon, clat, imdates,
                                             center_time=center_time, path=path)
        if not bperp_dict:
            print('WARN: Could not get baseline info from ASF; '
                  'LiCSBAS02 will generate dummy bperp')
            return

        missing = [imd for imd in imdates if imd not in bperp_dict]
        if missing:
            print(f'WARN: {len(missing)} of {len(imdates)} epochs missing from '
                  f'ASF baseline stack: {" ".join(missing)}')
            print('No baselines file created; LiCSBAS02 will generate dummy bperp')
            return

        io_lib.make_bperp_file(bperp_file, imdates, bperp_dict)
        print(f'Wrote baselines with bperp of {len(imdates)} epochs from ASF')

    except Exception as e:
        print(f'WARN: baselines creation failed ({e}); '
              'LiCSBAS02 will generate dummy bperp')


def main(indir='.', outdir='GEOC', cc_thresh=0.3, n_workers=None,
         resampling='nearest', dem=None, wavelength=None, center_time=None,
         overwrite=False, no_baselines=False):

    # %% Setting
    workdir, dolphin_dir = find_dirs(indir)
    print(f'Work directory: {workdir}')
    print(f'dolphin directory: {dolphin_dir}')

    out_geoc = os.path.abspath(outdir)
    os.makedirs(out_geoc, exist_ok=True)

    # metadata from dolphin_config.yaml (CLI options take precedence)
    wl_config, cslc_names = read_dolphin_config(workdir, dolphin_dir)
    if wavelength is None:
        wavelength = wl_config
    if center_time is None:
        center_time = get_center_time(cslc_names)
        if center_time is None:
            print('No acquisition time found (no cslc_file_list?); skip center_time')

    # map resampling choice to GDAL constant
    resample_map = {
        'nearest': gdal.GRA_NearestNeighbour,
        'bilinear': gdal.GRA_Bilinear,
        'cubic': gdal.GRA_Cubic,
        'lanczos': gdal.GRA_Lanczos,
        'average': gdal.GRA_Average,
        'mode': gdal.GRA_Mode,
    }
    resample_alg = resample_map.get(resampling, gdal.GRA_NearestNeighbour)
    print(f'Resampling method for unw: {resampling}')
    print(f'cc_thresh: {cc_thresh}')

    pairs = find_pairs(dolphin_dir)
    if not pairs:
        raise FileNotFoundError(f'No *.unw.tif found in {dolphin_dir}/unwrapped')
    n_no_cc = sum(1 for pair in pairs if pair[2] is None)
    print(f'Found {len(pairs)} IFGs to process ({n_no_cc} without int.cor.tif)')

    # common output grid for all products, from the first unw
    grid = get_target_grid(pairs[0][1])
    print(f"Output grid: {grid['width']} x {grid['length']}, EPSG:4326")

    n_cpu = max(1, mp.cpu_count() - 1)
    if n_workers is not None and n_workers > 0:
        n_workers = min(len(pairs), int(n_workers))
    else:
        n_workers = min(len(pairs), n_cpu)
    print(f'Number of parallel workers: {n_workers}')


    # %% Process unw and cc
    print(f'\nStart {n_workers} parallel processing of {len(pairs)} unw and cc...\n')
    args_list = [(pair, unw, cor, out_geoc, grid, cc_thresh, resample_alg,
                  overwrite) for pair, unw, cor in pairs]
    ctx = mp.get_context('fork')
    with ctx.Pool(n_workers) as pool:
        results = pool.map(process_pair, args_list)

    n_ok = results.count('ok')
    n_skip = results.count('skip')
    n_fail = results.count('fail')
    print(f'\n{n_ok} IFGs converted, {n_skip} skipped (existing), '
          f'{results.count("no_cc")} skipped (no cc), {n_fail} failed')

    if n_fail > 0:
        failed = [pairs[i][0] for i, r in enumerate(results) if r == 'fail']
        raise RuntimeError(f'{n_fail} IFGs failed ({" ".join(failed)}); '
                           're-run to retry them (existing IFGs are skipped)')
    if n_ok + n_skip == 0:
        print('No IFGs converted; aborting top-level file creation')
        return


    # %% Top-level files: hgt, ENU, metadata
    geom_dir = find_geometry_dir(workdir, dolphin_dir)

    if dem is None:
        # prefer geometry/height.tif: already geocoded on the processing grid
        candidates = [os.path.join(geom_dir, 'height.tif')] if geom_dir else []
        candidates += [os.path.join(d, 'dem.tif') for d in [workdir, dolphin_dir]]
        dem = next((path for path in candidates if os.path.exists(path)), None)
    if dem and os.path.exists(dem):
        make_hgt(dem, grid, out_geoc)
    else:
        print('No height.tif or dem.tif found; skip hgt creation (optional for LiCSBAS)')

    make_enu(geom_dir, grid, out_geoc)

    write_metadata(out_geoc, wavelength, center_time)

    print('No mli in dolphin output; skip (optional for LiCSBAS)')

    if no_baselines:
        print('--no_baselines specified; LiCSBAS02 will generate dummy bperp')
    else:
        make_baselines(out_geoc, [pair[0] for pair in pairs], grid['bounds'],
                       center_time, wavelength, cslc_names, overwrite)


if __name__ == '__main__':
    start = time.time()
    prog = os.path.basename(sys.argv[0])
    print(f"\n{prog} ver1.1.0 20260811 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    # Show argparse defaults in help messages
    p = argparse.ArgumentParser(description='Convert sweets/dolphin output directory to LiCSBAS GEOC layout',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--indir', default='.', help='sweets work directory (containing dolphin/) or dolphin directory itself')
    p.add_argument('-o', '--outdir', default='GEOC', help='Target GEOC directory')
    p.add_argument('-t', '--cc_thresh', type=float, default=0.3, help='Coherence threshold for masking unw')
    p.add_argument('-n', '--n_workers', type=int, default=None, help='Number of parallel workers (default: num CPU cores - 1)')
    p.add_argument('-r', '--resampling', default='nearest',
                   choices=['nearest', 'bilinear', 'cubic', 'lanczos', 'average', 'mode'],
                   help='Resampling method for unw reprojection')
    p.add_argument('--dem', default=None, help='DEM GeoTIFF for hgt (default: geometry/height.tif or dem.tif in the work directory)')
    p.add_argument('--wavelength', type=float, default=None, help='Radar wavelength in m (default: from dolphin_config.yaml)')
    p.add_argument('--center_time', default=None, help='Center time HH:MM:SS (default: from cslc file names in dolphin_config.yaml)')
    p.add_argument('--overwrite', action='store_true', help='Overwrite existing outputs')
    p.add_argument('--no_baselines', action='store_true', help='Do not query ASF for bperp (no network access); LiCSBAS02 will generate dummy bperp')
    args = p.parse_args()

    main(args.indir, args.outdir, args.cc_thresh, args.n_workers,
         args.resampling, args.dem, args.wavelength, args.center_time,
         args.overwrite, args.no_baselines)

    # Finish
    elapsed_time = datetime.timedelta(seconds=(time.time()-start))
    print(f"\nElapsed time: {elapsed_time}")
    print(f'\n{prog} Successfully finished!!\n')

    print(f"Output directory: {args.outdir}\n")
