#!/usr/bin/env python3
"""
Convert an ARIA-style parent directory into a LiCSBAS-style `GEOC/` layout.

This is a simplified converter compared to the ASF converter: it assumes ARIA
VRT inputs are already in EPSG:4326 and therefore does not perform reprojection
or resampling. It also does not implement geographic cropping.

Inputs (expected under ARIA parent dir):
- amplitude/{yyyymmdd_yyyymmdd}.vrt    -> X.geo.mli.tif (per-frame mli)
- azimuthAngle/{yyyymmdd_yyyymmdd}.vrt -> lv_phi equivalent
- incidenceAngle/{yyyymmdd_yyyymmdd}.vrt -> lv_theta equivalent
- coherence/{yyyymmdd_yyyymmdd}.vrt    -> cc
- unwrappedPhase/{yyyymmdd_yyyymmdd}.vrt -> unw
- DEM/glo_90.dem.vrt                   -> hgt
- products/*.nc                         -> used to parse center_time (hhmmss in filename)

Outputs:
- GEOC/{yyyymmdd_yyyymmdd}/{yyyymmdd_yyyymmdd}.geo.unw.tif
- GEOC/{yyyymmdd_yyyymmdd}/{yyyymmdd_yyyymmdd}.geo.cc.tif
- GEOC/X.geo.E.tif, X.geo.N.tif, X.geo.U.tif (computed from azimuth/incidence)
- GEOC/X.geo.hgt.tif
- GEOC/X.geo.mli.tif
- GEOC/metadata.txt (center_time from products filenames, radar_freq fixed)

"""

import argparse
import datetime
import glob
import os
import re
import sys
import time

from osgeo import gdal
import numpy as np


def find_vrt(folder, key):
    """Return path to a VRT for key (tries exact, then _uncropped)."""
    if not os.path.isdir(folder):
        return None
    p1 = os.path.join(folder, f"{key}.vrt")
    if os.path.exists(p1):
        return p1
    # sometimes files don't have .vrt extension listed exactly; try glob
    hits = glob.glob(os.path.join(folder, f"*{key}*.vrt"))
    return hits[0] if hits else None


def write_geotiff(path, arr, gt, proj, dtype=gdal.GDT_Float32, nodata=0.0, options=None):
    if options is None:
        options = ['COMPRESS=DEFLATE', 'TILED=YES', 'PREDICTOR=2']
    driver = gdal.GetDriverByName('GTiff')
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


def parse_center_time_from_products(products_dir, date_key=None):
    """Search products for a filename that contains date_key and extract hhmmss.
    If date_key is None, use the first product file found.
    Returns time string 'HH:MM:SS' or None.
    """
    if not os.path.isdir(products_dir):
        return None
    files = sorted([os.path.basename(p) for p in glob.glob(os.path.join(products_dir, '*'))])
    for fn in files:
        if date_key and date_key not in fn:
            continue
        # look for a 6-digit time like -204350- or _204350_
        m = re.search(r'[-_](\d{6})[-_]', fn)
        if m:
            hhmmss = m.group(1)
            hh = hhmmss[0:2]; mm = hhmmss[2:4]; ss = hhmmss[4:6]
            return f"{hh}:{mm}:{ss}"
    return None


def main(indir='.', outdir='GEOC', cc_thresh=0.5, overwrite=False):
    indir = os.path.abspath(indir)
    if not os.path.isdir(indir):
        raise FileNotFoundError(f"{indir} not found")

    out_geoc = os.path.abspath(outdir)
    os.makedirs(out_geoc, exist_ok=True)

    # expected ARIA subfolders
    amp_dir = os.path.join(indir, 'amplitude')
    az_dir = os.path.join(indir, 'azimuthAngle')
    inc_dir = os.path.join(indir, 'incidenceAngle')
    coh_dir = os.path.join(indir, 'coherence')
    unw_dir = os.path.join(indir, 'unwrappedPhase')
    dem_dir = os.path.join(indir, 'DEM')
    products_dir = os.path.join(indir, 'products')

    # collect available date keys from unwrappedPhase VRTs
    date_keys = []
    if os.path.isdir(unw_dir):
        for p in glob.glob(os.path.join(unw_dir, '*.vrt')):
            bn = os.path.basename(p)
            key = bn.replace('.vrt', '')
            date_keys.append(key)
    date_keys = sorted(list(dict.fromkeys(date_keys)))
    if not date_keys:
        raise FileNotFoundError('No unwrappedPhase VRTs found in ARIA directory')

    print(f'Found {len(date_keys)} IFGs to process')

    processed_any = False
    for key in date_keys:
        print(f'Processing {key}...')
        src_unw = find_vrt(unw_dir, key)
        src_cc = find_vrt(coh_dir, key)
        if not (src_unw and src_cc):
            print(f'  WARN: missing unw or cc for {key}. Skipping')
            continue

        # open unw and cc (inputs are assumed to be in EPSG:4326 already)
        ds_unw = gdal.Open(src_unw)
        ds_cc = gdal.Open(src_cc)
        if ds_unw is None or ds_cc is None:
            print(f'  ERROR opening datasets for {key}. Skipping')
            continue

        unw_arr = ds_unw.ReadAsArray().astype(np.float32)
        cc_arr = ds_cc.ReadAsArray().astype(np.float32)
        # scale cc if 0-255
        if np.nanmax(cc_arr) > 1.5:
            cc_arr = cc_arr / 255.0

        # ARIA unwrappedPhase uses opposite sign compared to GEOC/ASF: flip sign
        unw_arr = -1.0 * unw_arr

        mask = (cc_arr < cc_thresh) | np.isnan(cc_arr)
        unw_arr[mask] = 0.0

        gt = ds_unw.GetGeoTransform()
        proj = ds_unw.GetProjection()

        # ARIA date key is reversed relative to GEOC/ASF (ARIA: y2_y1). Swap
        # order for GEOC output names so they become y1_y2.
        parts = key.split('_', 1)
        if len(parts) == 2:
            out_key = f"{parts[1]}_{parts[0]}"
        else:
            out_key = key

        # create IFG folder with swapped date order
        tgt_dir = os.path.join(out_geoc, out_key)
        os.makedirs(tgt_dir, exist_ok=True)
        tgt_unw = os.path.join(tgt_dir, f"{out_key}.geo.unw.tif")
        tgt_cc = os.path.join(tgt_dir, f"{out_key}.geo.cc.tif")

        write_geotiff(tgt_unw, unw_arr, gt, proj)
        write_geotiff(tgt_cc, cc_arr, gt, proj)
        print(f'  Wrote {os.path.basename(tgt_unw)} and {os.path.basename(tgt_cc)}')
        processed_any = True

    # create top-level files using the first successful key
    if not processed_any:
        print('No IFGs processed; aborting top-level file creation')
        return

    first_key = date_keys[0]

    # ENU from azimuthAngle (phi) and incidenceAngle (theta)
    src_phi = find_vrt(az_dir, first_key)
    src_theta = find_vrt(inc_dir, first_key)
    frameid = 'X'
    if src_phi and src_theta:
        # ARIA: azimuthAngle in degrees, 0 = North, positive clockwise
        # ARIA: incidenceAngle in degrees, 0 = zenith (up), 90 = horizontal
        ds_phi = gdal.Open(src_phi)
        ds_theta = gdal.Open(src_theta)
        az_deg = ds_phi.ReadAsArray().astype(np.float32)
        inc_deg = ds_theta.ReadAsArray().astype(np.float32)

        # Convert ARIA az/inc into angles for ENU computation.
        # Use phi = deg2rad(az_deg) and elevation = 90 - inc_deg (degrees).
        # Then theta = deg2rad(elevation).
        phi = np.deg2rad(az_deg)
        theta = np.deg2rad(90.0 - inc_deg)

        # mask invalid where source arrays are nan
        mask = np.isnan(az_deg) | np.isnan(inc_deg)

        E = np.cos(theta) * np.cos(phi)
        N = np.cos(theta) * np.sin(phi)
        U = np.sin(theta)
        norm = np.sqrt(E**2 + N**2 + U**2)
        nz = norm > 0
        E[nz] = E[nz] / norm[nz]
        N[nz] = N[nz] / norm[nz]
        U[nz] = U[nz] / norm[nz]
        E[mask] = 0.0; N[mask] = 0.0; U[mask] = 0.0

        gt_phi = ds_phi.GetGeoTransform()
        proj_phi = ds_phi.GetProjection()

        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.E.tif"), E, gt_phi, proj_phi)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.N.tif"), N, gt_phi, proj_phi)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.U.tif"), U, gt_phi, proj_phi)
        print('Wrote ENU unit maps')
    else:
        print('WARN: azimuthAngle or incidenceAngle missing for top-level ENU creation')

    # mli: amplitude first_key
    src_mli = find_vrt(amp_dir, first_key)
    if src_mli:
        ds_m = gdal.Open(src_mli)
        mli_arr = ds_m.ReadAsArray().astype(np.float32)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.mli.tif"), mli_arr, ds_m.GetGeoTransform(), ds_m.GetProjection())
        print('Wrote mli')
    else:
        print('WARN: amplitude mli not found for top-level mli')

    # hgt: DEM/glo_90.dem.vrt
    src_hgt = find_vrt(dem_dir, 'glo_90.dem') or find_vrt(dem_dir, 'glo_90') or find_vrt(dem_dir, 'glo_90.dem.vrt')
    # also try explicit file
    if not src_hgt:
        possible = glob.glob(os.path.join(dem_dir, '*.vrt'))
        src_hgt = possible[0] if possible else None
    if src_hgt:
        ds_h = gdal.Open(src_hgt)
        hgt_arr = ds_h.ReadAsArray().astype(np.float32)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.hgt.tif"), hgt_arr, ds_h.GetGeoTransform(), ds_h.GetProjection())
        print('Wrote hgt')
    else:
        print('WARN: DEM not found for top-level hgt')

    # metadata: parse center_time from products
    center_time = parse_center_time_from_products(products_dir, first_key)
    metadata_path = os.path.join(out_geoc, 'metadata.txt')
    with open(metadata_path, 'w') as mf:
        if center_time:
            mf.write(f"center_time={center_time}\n")
        mf.write('radar_freq=5.405e9\n')
    print(f'Wrote metadata: {metadata_path}')


if __name__ == '__main__':
    start = time.time()
    prog = os.path.basename(sys.argv[0])
    print(f"\n{prog} ver1.0.0 20251204 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    # Show argparse defaults in help messages
    p = argparse.ArgumentParser(description='Convert ARIA-style folder to LiCSBAS GEOC layout',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--indir', default='.', help='ARIA parent directory')
    p.add_argument('-o', '--outdir', default='GEOC', help='Target GEOC directory')
    p.add_argument('-t', '--cc_thresh', type=float, default=0.5, help='Coherence threshold for masking unw')
    p.add_argument('--overwrite', action='store_true', help='Overwrite existing outputs')
    args = p.parse_args()

    main(args.indir, args.outdir, args.cc_thresh, args.overwrite)

    # Finish
    elapsed_time = datetime.timedelta(seconds=(time.time()-start))
    print(f"\nElapsed time: {elapsed_time}")
    print(f'\n{prog} Successfully finished!!\n')

    print(f"Output directory: {args.outdir}\n")