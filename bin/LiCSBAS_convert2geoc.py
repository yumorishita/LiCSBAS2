#!/usr/bin/env python3
"""
bin/LiCSBAS_asf2geoc.py

Convert a parent directory of ASF product folders into a LiCSBAS-style `GEOC/` layout.

This script scans a parent directory for ASF product subdirectories and, for each
product that contains an interferogram (unwrapped phase) and coherence maps, it
produces the following under the output `GEOC/` directory:

- Per-IFG folder: GEOC/YYYYMMDD_YYYYMMDD/
  - <yyyymmdd>_<yyyymmdd>.geo.unw.tif  (masked by coherence threshold)
  - <yyyymmdd>_<yyyymmdd>.geo.cc.tif   (coherence map)

- Top-level files (from the first successful product):
  - X.geo.E.tif, X.geo.N.tif, X.geo.U.tif  (look-vector-derived ENU unit maps, or copied if available)
  - X.geo.hgt.tif  (DEM)
  - X.geo.mli.tif  (amplitude/mli)
  - metadata.txt  (contains center_time and radar_freq)

Example: expected layout of an input parent directory
- a parent directory (one example product folder name pattern): *YYYYMMDD*YYYYMMDD* (e.g., S1AA_20220111T204320_20220216T204318_VVP036_INT80_G_ueF_263B)
  - *_unw_phase.tif
  - *_corr.tif (filtered coherence, not common coherence?)
  - *_lv_phi.tif
  - *_lv_theta.tif
  - *_dem.tif
  - *_amp.tif

Key behaviors and implementation notes:
- The script reprojects inputs to EPSG:4326 (WGS84) and writes all outputs
  as GeoTIFFs in EPSG:4326.
- Unwrapped phase (UNW) pixels are masked where coherence (CC) is below a
  configurable threshold (default 0.5). Masked UNW pixels are set to 0.
- ENU (E/N/U) maps are computed from lv_phi and lv_theta when available; the
  script normalizes to unit vectors and writes the three bands as separate
  GeoTIFFs named `X.geo.E.tif`, `X.geo.N.tif`, `X.geo.U.tif`.
- Outputs are written with DEFLATE compression, tiling and predictor enabled.
- The script supports an optional geographic crop bbox to restrict outputs
  (format: lonmin/lonmax/latmin/latmax).
- Processing of individual products is parallelized (multiprocessing Pool).
- Resampling algorithm for reprojection is selectable (nearest, bilinear, etc.).
- metadata.txt is created with `center_time` parsed from the product TXT file
  (looks for lines like "UTC time: 74601.153494") and a fixed `radar_freq=5.405e9`.

CLI (helper usage is printed by running the script with -h): typical options
include: -i/--indir, -o/--outdir, -c/--crop_geo, -t/--cc_thresh,
-n/--n_workers, -r/--resampling.

Behavior on errors:
- The script currently prefers to let exceptions propagate (fail-fast) so that
  errors are visible during processing. This can be changed to more tolerant
  behavior if desired.

"""

# %% Import
import argparse
import datetime
import glob
import multiprocessing as mp
import os
import re
import shutil
import sys
import time

from osgeo import gdal
import numpy as np

# %% Functions
def parse_ifg_dates_from_name(name):
    """Return (yyyymmdd, yyyymmdd) or None if not found."""
    # try pattern with time: 20220111T204320_20220216T204318 or with other separators
    m = re.search(r"(\d{8})T\d{6}.*?(\d{8})T\d{6}", name)
    if m:
        return m.group(1), m.group(2)
    # try simple YYYYMMDD_YYYYMMDD
    m = re.search(r"(\d{8})[_-](\d{8})", name)
    if m:
        return m.group(1), m.group(2)
    return None


def find_file_by_suffix(directory, suffixes):
    """Search for files under directory (non-recursive) whose name contains any of suffixes.
    Return the first matching path or None.
    """
    for suf in suffixes:
        pattern = os.path.join(directory, f"*{suf}")
        hits = glob.glob(pattern)
        if hits:
            return hits[0]
    return None


def lonlat_to_xy(lon, lat, gt):
    """Convert lon,lat (degrees) to pixel x,y using GDAL geoTransform gt.
    gt: (originX, pixelWidth, 0, originY, 0, pixelHeight)
    Returns integer (x,y) indices (0-based).
    Uses pixel-center convention by shifting origin to the center of the pixel.
    """
    originX, pixelWidth, _, originY, _, pixelHeight = gt
    # pixel center reference
    lon0 = originX + pixelWidth / 2.0
    lat0 = originY + pixelHeight / 2.0
    x = int(round((lon - lon0) / pixelWidth))
    y = int(round((lat - lat0) / pixelHeight))
    return x, y


def crop_array_and_gt(arr, gt, crop_geo):
    """Crop a 2D or band-first 3D array to the geographic bbox using gt.
    Returns (cropped_array, new_gt). If bbox is outside or invalid, returns (arr, gt).
    crop_geo is (lon_min, lon_max, lat_min, lat_max)
    """
    if crop_geo is None:
        return arr, gt
    lon_min, lon_max, lat_min, lat_max = crop_geo
    # determine array shape
    if arr.ndim == 3:
        _, ny, nx = arr.shape
    else:
        ny, nx = arr.shape
    # compute pixel indices
    x1, y1 = lonlat_to_xy(lon_min, lat_max, gt)  # top-left
    x2, y2 = lonlat_to_xy(lon_max, lat_min, gt)  # bottom-right
    # normalize order
    if x1 > x2: x1, x2 = x2, x1
    if y1 > y2: y1, y2 = y2, y1
    # clamp
    x1 = max(0, x1); y1 = max(0, y1)
    x2 = min(nx-1, x2); y2 = min(ny-1, y2)
    # if empty region, return original
    if x2 < x1 or y2 < y1:
        return arr, gt
    # slice
    if arr.ndim == 3:
        sub = arr[:, y1:y2+1, x1:x2+1]
    else:
        sub = arr[y1:y2+1, x1:x2+1]
    # new origin
    new_originX = gt[0] + x1 * gt[1]
    new_originY = gt[3] + y1 * gt[5]
    new_gt = (new_originX, gt[1], 0.0, new_originY, 0.0, gt[5])
    return sub, new_gt


def write_geotiff(path, arr, gt, proj, dtype=gdal.GDT_Float32, nodata=0.0, options=None):
    """Create a single-band GeoTIFF at path from 2D array (or band-first 3D by taking first band).
    arr: 2D or 3D (bands, y, x)
    gt: geotransform tuple
    proj: WKT projection string
    """
    if options is None:
        options = ['COMPRESS=DEFLATE', 'TILED=YES', 'PREDICTOR=2']
    driver = gdal.GetDriverByName('GTiff')
    # if 3D (bands, y, x) use first band
    if getattr(arr, 'ndim', 2) == 3:
        arr_write = arr[0]
    else:
        arr_write = arr
    ny, nx = arr_write.shape
    outds = driver.Create(path, nx, ny, 1, dtype, options=options)
    outds.SetGeoTransform(gt)
    outds.SetProjection(proj)
    outds.GetRasterBand(1).WriteArray(arr_write)
    outds.GetRasterBand(1).SetNoDataValue(nodata)
    outds.FlushCache()
    outds = None


# Worker for parallel processing of individual products
def process_prod(args):
    prod, out_geoc, cc_thresh, crop_geo, resample_alg = args
    basename = os.path.basename(prod.rstrip(os.sep))
    dates = parse_ifg_dates_from_name(basename)
    unw = find_file_by_suffix(prod, ['_unw_phase.tif'])
    cc = find_file_by_suffix(prod, ['_corr.tif'])
    if not (dates and unw and cc):
        print(f'{basename} is not a product directory. Skip.')
        return {'prod': prod, 'basename': basename, 'has_ifg': False}

    ifgdir = f"{dates[0]}_{dates[1]}"
    target_ifg_dir = os.path.join(out_geoc, ifgdir)
    os.makedirs(target_ifg_dir, exist_ok=True)

    target_unw = os.path.join(target_ifg_dir, f"{ifgdir}.geo.unw.tif")
    target_cc = os.path.join(target_ifg_dir, f"{ifgdir}.geo.cc.tif")

    print('Processing:', ifgdir)

    # Open and reproject to WGS84 (EPSG:4326) in-memory so all outputs are in that CRS
    # Note: gdal<3.7 return xres=yres
    # TODO upgrade gdal to >=3.7 so that the resolution of the source file(s) is preserved
    src_ds_unw = gdal.Open(unw)
    ds_unw = gdal.Warp('', src_ds_unw, format='MEM', dstSRS='EPSG:4326')
    unw_arr = ds_unw.ReadAsArray().astype(np.float32)
    src_ds_cc = gdal.Open(cc)
    ds_cc = gdal.Warp('', src_ds_cc, format='MEM', dstSRS='EPSG:4326')
    cc_arr = ds_cc.ReadAsArray().astype(np.float32)

    # scale cc if in 0-255
    if np.nanmax(cc_arr) > 1.5:
        cc_arr = cc_arr / 255.0

    # prepare output GeoTransform (may be adjusted if cropping)
    out_gt = ds_unw.GetGeoTransform()
    cropped = False

    # Apply crop_geo if provided: crop arrays and update out_gt
    if crop_geo is not None:
        gt = ds_unw.GetGeoTransform()
        unw_arr, out_gt = crop_array_and_gt(unw_arr, gt, crop_geo)
        cc_arr, _ = crop_array_and_gt(cc_arr, gt, crop_geo)
        cropped = True

    # compute mask after cropping so shapes match
    mask = (cc_arr < cc_thresh) | np.isnan(cc_arr)

    # mask unw
    unw_arr[mask] = 0.0

    # write masked unw (ds_unw is already reprojected to EPSG:4326)
    proj = ds_unw.GetProjection()
    if unw_arr.ndim == 3:
        unw_write = unw_arr[0]
    else:
        unw_write = unw_arr
    write_geotiff(target_unw, unw_write, out_gt, proj)

    # write cc (cropped if crop_geo applied) or copy original if not warped
    if 'cc_arr' in locals():
        if cc_arr.ndim == 3:
            cc_write = cc_arr[0]
        else:
            cc_write = cc_arr
        if cropped:
            gt_cc = out_gt
        else:
            gt_cc = ds_cc.GetGeoTransform()
        proj_cc = ds_cc.GetProjection()
        write_geotiff(target_cc, cc_write, gt_cc, proj_cc)

    return {'prod': prod, 'basename': basename, 'has_ifg': True, 'ifgdir': ifgdir}


# %% main
def main(indir='.', outdir='GEOC', crop_geo=None, cc_thresh=0.5,
         n_workers=None, resampling='nearest'):

    # %% Setting
    indir = os.path.abspath(indir)
    if not os.path.isdir(indir):
        raise FileNotFoundError(f"ERROR: {indir} is not a directory")

    out_geoc = os.path.abspath(outdir)
    os.makedirs(out_geoc, exist_ok=True)

    # collect product subdirectories (non-recursive)
    # Exclude the output GEOC directory if it resides inside indir to avoid processing it
    prod_dirs = sorted([os.path.join(indir, name) for name in os.listdir(indir)
                        if os.path.isdir(os.path.join(indir, name)) and
                        os.path.abspath(os.path.join(indir, name)) != out_geoc])
    if len(prod_dirs) == 0:
        raise FileNotFoundError(f"ERROR: No subdirectories found in {indir}")
    print(f'Number of product subdirectories: {len(prod_dirs)}')

    # Parallel processing of products for IFG (unw/cc) handling
    n_cpu = max(1, mp.cpu_count() - 1)

    # allow user override via CLI
    if n_workers is not None and n_workers > 0:
        n_workers = min(len(prod_dirs), int(n_workers))
    else:
        n_workers = min(len(prod_dirs), n_cpu)
    print(f'Number of parallel workers: {n_workers}')

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
    print(f'Resampling method: {resampling}')

    # parse crop_geo arg if provided (format: lonmin/lonmax/latmin/latmax)
    if crop_geo is not None:
        parts = crop_geo.split('/')
        if len(parts) != 4:
            raise ValueError('crop_geo must be lonmin/lonmax/latmin/latmax')
        lon_min, lon_max, lat_min, lat_max = [float(v) for v in parts]
        crop_geo = (lon_min, lon_max, lat_min, lat_max)
        print(f'Crop area: {crop_geo}')


    # %% Process unw and cc
    print(f'\nStart {n_workers} parallel processing of {len(prod_dirs)} unw and cc...\n')
    args_list = [(prod, out_geoc, cc_thresh, crop_geo, resample_alg) for prod in prod_dirs]
    ctx = mp.get_context('fork')
    with ctx.Pool(n_workers) as pool:
        results = pool.map(process_prod, args_list)

    print('unw and cc processing finished.\n')


    # %% Find first product that had IFG processed successfully to create top-level files
    first = next((r for r in results if r.get('has_ifg')), None)
    if first:
        prod_for_top = first['prod']
        basename = os.path.basename(prod_for_top.rstrip(os.sep))
        # frameID is fixed to 'X'
        frameid = 'X'

        # Compute ENU from lv_phi and lv_theta (look vector maps) if available
        src_phi = find_file_by_suffix(prod_for_top, ['_lv_phi.tif'])
        src_theta = find_file_by_suffix(prod_for_top, ['_lv_theta.tif'])
        if src_phi and src_theta:
            # Reproject look vector maps to EPSG:4326 before computation
            src_ds_phi = gdal.Open(src_phi)
            src_ds_theta = gdal.Open(src_theta)
            ds_phi = gdal.Warp('', src_ds_phi, format='MEM', dstSRS='EPSG:4326')
            ds_theta = gdal.Warp('', src_ds_theta, format='MEM', dstSRS='EPSG:4326')
            phi = ds_phi.ReadAsArray().astype(np.float32)
            theta = ds_theta.ReadAsArray().astype(np.float32)

            mask = np.isnan(phi) | np.isnan(theta) | (phi == 0) | (theta == 0)
            E = np.cos(theta) * np.cos(phi)
            N = np.cos(theta) * np.sin(phi)
            U = np.sin(theta)
            norm = np.sqrt(E**2 + N**2 + U**2)
            nz = norm > 0
            E[nz] = E[nz] / norm[nz]
            N[nz] = N[nz] / norm[nz]
            U[nz] = U[nz] / norm[nz]
            E[mask] = 0.0
            N[mask] = 0.0
            U[mask] = 0.0

            gt = ds_phi.GetGeoTransform()
            proj = ds_phi.GetProjection()
            ny, nx = phi.shape

            dst_E = os.path.join(out_geoc, f"{frameid}.geo.E.tif")
            dst_N = os.path.join(out_geoc, f"{frameid}.geo.N.tif")
            dst_U = os.path.join(out_geoc, f"{frameid}.geo.U.tif")

            # Apply crop for ENU if requested using helper to avoid duplication
            if crop_geo is not None:
                E_c, new_gt = crop_array_and_gt(E, gt, crop_geo)
                N_c, _ = crop_array_and_gt(N, gt, crop_geo)
                U_c, _ = crop_array_and_gt(U, gt, crop_geo)
                new_ny, new_nx = E_c.shape
                arrs_to_write = [(E_c, dst_E, new_nx, new_ny, new_gt),
                                 (N_c, dst_N, new_nx, new_ny, new_gt),
                                 (U_c, dst_U, new_nx, new_ny, new_gt)]
            else:
                arrs_to_write = [(E, dst_E, nx, ny, gt), (N, dst_N, nx, ny, gt), (U, dst_U, nx, ny, gt)]

            for arr, dst, out_nx, out_ny, out_gt in arrs_to_write:
                write_geotiff(dst, arr, out_gt, proj)

            ds_phi = None
            ds_theta = None

            print(f"Computed ENU and wrote {os.path.basename(dst_E)}, {os.path.basename(dst_N)}, {os.path.basename(dst_U)}")
        else:
            for ENU in ['E', 'N', 'U']:
                src = find_file_by_suffix(prod_for_top, [f".geo.{ENU}.tif", f".{ENU}.geo.tif"])  # try common variants
                if src:
                    dst = os.path.join(out_geoc, f"{frameid}.geo.{ENU}.tif")
                    print(f"Copying/reprojecting {os.path.basename(src)} -> {os.path.basename(dst)}")
                    # ensure output is EPSG:4326
                    gdal.Warp(dst, src, format='GTiff', dstSRS='EPSG:4326',
                              resampleAlg=resample_alg,
                              options=['COMPRESS=DEFLATE','TILED=YES','PREDICTOR=2'])

        # hgt: use suffix _dem.tif
        src_hgt = find_file_by_suffix(prod_for_top, ['_dem.tif'])
        if src_hgt:
            dst_hgt = os.path.join(out_geoc, f"{frameid}.geo.hgt.tif")
            print(f"Copying/reprojecting {os.path.basename(src_hgt)} -> {os.path.basename(dst_hgt)}")
            if crop_geo is not None:
                # reproject to EPSG:4326 first and crop using helper
                src_ds_hgt = gdal.Open(src_hgt)
                ds_hgt = gdal.Warp('', src_ds_hgt, format='MEM', dstSRS='EPSG:4326')
                hgt_arr = ds_hgt.ReadAsArray().astype(np.float32)
                gt_h = ds_hgt.GetGeoTransform()
                proj_h = ds_hgt.GetProjection()
                sub, new_gt = crop_array_and_gt(hgt_arr, gt_h, crop_geo)
                write_geotiff(dst_hgt, sub, new_gt, proj_h)
                ds_hgt = None

        # mli: use suffix _amp.tif
        src_mli = find_file_by_suffix(prod_for_top, ['_amp.tif'])
        if src_mli:
            dst_mli = os.path.join(out_geoc, f"{frameid}.geo.mli.tif")
            print(f"Copying/reprojecting {os.path.basename(src_mli)} -> {os.path.basename(dst_mli)}")
            if crop_geo is not None:
                # reproject to EPSG:4326 first and crop using helper
                src_ds_m = gdal.Open(src_mli)
                ds_m = gdal.Warp('', src_ds_m, format='MEM', dstSRS='EPSG:4326')
                mli_arr = ds_m.ReadAsArray().astype(np.float32)
                gt_m = ds_m.GetGeoTransform()
                proj_m = ds_m.GetProjection()
                subm, new_gt = crop_array_and_gt(mli_arr, gt_m, crop_geo)
                write_geotiff(dst_mli, subm, new_gt, proj_m)
                ds_m = None

        # Create metadata.txt in GEOC using UTC time found in the product's fixed TXT file
        txt_file = os.path.join(prod_for_top, f"{basename}.txt")
        utc_secs = []
        if os.path.exists(txt_file):
            try:
                with open(txt_file, 'r') as fh:
                    for line in fh:
                        # match either 'UTC time: 74601.153494' or 'UTC: 74601.153494' (case-insensitive)
                        m = re.search(r'UTC(?:\s*time)?\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)', line, flags=re.IGNORECASE)
                        if m:
                            try:
                                utc_secs.append(float(m.group(1)))
                            except Exception:
                                pass
            except Exception:
                pass
        else:
            print(f"WARN: Expected product TXT {txt_file} not found.")

        if utc_secs:
            center_sec = float(sum(utc_secs)) / len(utc_secs)
            hours = int(center_sec // 3600)
            rem = center_sec - hours * 3600
            minutes = int(rem // 60)
            seconds = rem - minutes * 60
            center_time_str = f"{hours}:{minutes:02d}:{seconds:09.6f}"
            print(f"Computed center_time from {os.path.basename(txt_file)}: {center_time_str}")
        else:
            center_time_str = None
            print('WARN: No UTC time found in product TXT file. metadata.txt will contain default radar_freq only.')

        metadata_path = os.path.join(out_geoc, 'metadata.txt')
        try:
            with open(metadata_path, 'w') as mf:
                if center_time_str is not None:
                    mf.write(f"center_time={center_time_str}\n")
                mf.write("radar_freq=5.405e9\n")
            print(f"Wrote metadata: {metadata_path}")
        except Exception as e:
            print(f"ERROR writing metadata.txt: {e}")

        # Note: baseline file generation is intentionally not implemented here.
        # TODO: implement baselines file creation if required.

    else:
        print('No IFG products were found/processed. No top-level files created.')

    print('\nDone.')


# %% main
if __name__ == '__main__':
    start = time.time()
    prog = os.path.basename(sys.argv[0])
    print(f"\n{prog} ver1.0.0 20251016 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    # Read args
    p = argparse.ArgumentParser(description='Convert ASF product parent dir to LiCSBAS GEOC layout')
    p.add_argument('-i', '--indir', default='.',
                   help='Parent directory that contains ASF product subdirectories (default: .)')
    p.add_argument('-o', '--outdir', default='GEOC', help='Target GEOC directory (default: ./GEOC)')
    p.add_argument('-c', '--crop_geo', default=None,
                   help="Crop outputs to geographic bbox in degrees: 'lonmin/lonmax/latmin/latmax' (e.g. 138/140/34/36)")
    p.add_argument('-t', '--cc_thresh', type=float, default=0.5,
                   help='Coherence threshold (0-1). Pixels with cc < threshold will be masked in unw (default: 0.5)')
    p.add_argument('-n', '--n_workers', type=int, default=None,
                   help='Number of parallel workers to use for IFG processing (default: cpu_count()-1)')
    p.add_argument('-r', '--resampling', type=str, default='nearest',
                   choices=['nearest', 'bilinear', 'cubic', 'lanczos', 'average', 'mode'],
                   help='Resampling algorithm for reprojection/cropping (default: nearest)')
    args = p.parse_args()

    main(args.indir,
         args.outdir,
         args.crop_geo,
         args.cc_thresh,
         args.n_workers,
         args.resampling)

    # Finish
    elapsed_time = datetime.timedelta(seconds=(time.time()-start))
    print(f"\nElapsed time: {elapsed_time}")
    print(f'\n{prog} Successfully finished!!\n')

    print(f"Output directory: {args.outdir}\n")