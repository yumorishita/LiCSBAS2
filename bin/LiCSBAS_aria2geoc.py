#!/usr/bin/env python3
"""
Convert an ARIA-style parent directory into a LiCSBAS-style `GEOC/` layout.

Inputs (expected under ARIA parent dir):
- stack/
    - ampStack.vrt (amplitude/{yyyymmdd_yyyymmdd}.vrt) -> X.geo.mli.tif (per-frame mli)
    - cohStack.vrt (coherence/{yyyymmdd_yyyymmdd}.vrt) -> cc
    - unwrapStack.vrt (unwrappedPhase/{yyyymmdd_yyyymmdd}.vrt) -> unw
      (contains metadata including `Wavelength (m)` used to derive radar frequency)
- azimuthAngle/{yyyymmdd_yyyymmdd}.vrt -> lv_phi equivalent
- incidenceAngle/{yyyymmdd_yyyymmdd}.vrt -> lv_theta equivalent
- DEM/glo_90.dem.vrt                   -> hgt

Outputs:
- GEOC/{yyyymmdd_yyyymmdd}/{yyyymmdd_yyyymmdd}.geo.unw.tif
- GEOC/{yyyymmdd_yyyymmdd}/{yyyymmdd_yyyymmdd}.geo.cc.tif
- GEOC/X.geo.E.tif, X.geo.N.tif, X.geo.U.tif (computed from azimuth/incidence)
- GEOC/X.geo.hgt.tif
- GEOC/X.geo.mli.tif
- GEOC/metadata.txt (center_time and radar_freq from stack metadata)

"""

import argparse
import datetime
import glob
import os
import re
import sys
import time

from osgeo import gdal

### Raise exceptions on GDAL errors instead of returning None.
### Required from GDAL 3.7 (FutureWarning) and the default in GDAL 4.0.
gdal.UseExceptions()
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


def _keys_from_stack(vrt_path):
    """Return a list of date keys found in the band metadata of *vrt_path*.

    The ARIA stack VRTs store per-band metadata with a key named ``Dates``
    containing the IFG identifier.  We use this to enumerate all available
    interferograms without touching the individual subdirectories.
    """
    ks = []
    if os.path.exists(vrt_path):
        ds = gdal.Open(vrt_path)
        if ds is not None:
            for ib in range(1, ds.RasterCount + 1):
                band = ds.GetRasterBand(ib)
                # look through all domains until we find a Dates entry
                for dom in band.GetMetadataDomainList() or []:
                    md = band.GetMetadata(dom)
                    if 'Dates' in md:
                        ks.append(md['Dates'])
                        break
            ds = None
    return ks


def _canonical_ifg_key(key):
    """Return a normalized key with the earlier date first.

    ARIA-style keys are sometimes stored with the later date first (``y2_y1``)
    whereas many workflows (including the LiCSBAS GEOC layout) conventionally
    use the earlier date first (``y1_y2``).  We cannot rely on the input
    order, so this helper takes an IFG key
    and, when it looks like two eight-digit dates separated by an underscore,
    it returns a string with the two dates sorted in chronological order.  If
    the key does not match that simple pattern it is returned unchanged.

    This function is used both for sorting the list of keys we process and for
    constructing the output directory/file names.
    """

    parts = key.split('_', 1)
    if len(parts) == 2:
        try:
            d0 = datetime.datetime.strptime(parts[0], '%Y%m%d')
            d1 = datetime.datetime.strptime(parts[1], '%Y%m%d')
        except ValueError:
            # not two parseable dates, just return the original key
            return key
        # assemble with earlier date first
        if d0 <= d1:
            return f"{parts[0]}_{parts[1]}"
        else:
            return f"{parts[1]}_{parts[0]}"
    return key


def _open_band(stack_path, key):
    """Open *stack_path* and return the dataset plus band index for *key*.

    If the stack VRT exists and contains a band whose ``Dates`` metadata
    matches *key* the dataset object and the 1-based band number are
    returned.  Otherwise ``(None, None)`` is returned.
    """
    if os.path.exists(stack_path):
        ds = gdal.Open(stack_path)
        if ds is not None:
            for ib in range(1, ds.RasterCount + 1):
                band = ds.GetRasterBand(ib)
                for dom in band.GetMetadataDomainList() or []:
                    md = band.GetMetadata(dom)
                    if md.get('Dates') == key:
                        return ds, ib
            ds = None
    return None, None


def _wavelength_from_stack(vrt_path):
    """Return the first wavelength value (in metres) found in *vrt_path*.

    ARIA stack VRTs include a metadata field ``Wavelength (m)`` in the
    interferogram-domain metadata for each band.  We inspect the first band
    and return the float value, or ``None`` if it cannot be determined.
    """
    if os.path.exists(vrt_path):
        ds = gdal.Open(vrt_path)
        if ds is not None and ds.RasterCount >= 1:
            band = ds.GetRasterBand(1)
            for dom in band.GetMetadataDomainList() or []:
                md = band.GetMetadata(dom)
                if 'Wavelength (m)' in md:
                    try:
                        return float(md['Wavelength (m)'])
                    except ValueError:
                        pass
            ds = None
    return None


def _utc_time_from_stack(vrt_path):
    """Return the UTC time string found in the first band metadata of *vrt_path*.

    ARIA stack VRTs include a metadata field ``UTCTime (HH:MM:SS.ss)`` in the
    interferogram-domain metadata for each band.  We inspect the first band
    and return the time (truncated to HH:MM:SS), or ``None`` if it cannot be determined.
    """
    if os.path.exists(vrt_path):
        ds = gdal.Open(vrt_path)
        if ds is not None and ds.RasterCount >= 1:
            band = ds.GetRasterBand(1)
            for dom in band.GetMetadataDomainList() or []:
                md = band.GetMetadata(dom)
                if 'UTCTime (HH:MM:SS.ss)' in md:
                    utc_str = md['UTCTime (HH:MM:SS.ss)']
                    # extract HH:MM:SS part (first 8 chars)
                    if len(utc_str) >= 8:
                        return utc_str[:8]
                    return utc_str
            ds = None
    return None


def main(indir='.', outdir='GEOC', cc_thresh=None, overwrite=False):
    indir = os.path.abspath(indir)
    if not os.path.isdir(indir):
        raise FileNotFoundError(f"{indir} not found")

    out_geoc = os.path.abspath(outdir)
    os.makedirs(out_geoc, exist_ok=True)

    # expected ARIA subfolders.  amplitude/coherence/unwrappedPhase data are
    # read via the stack VRT files; azimuth/incidence/DEM are read directly.
    stack_dir = os.path.join(indir, 'stack')
    az_dir = os.path.join(indir, 'azimuthAngle')
    inc_dir = os.path.join(indir, 'incidenceAngle')
    dem_dir = os.path.join(indir, 'DEM')

    # convenience paths for stack-based VRTs; may not all exist but we will
    # try to open them when needed.
    amp_stack = os.path.join(stack_dir, 'ampStack.vrt')
    coh_stack = os.path.join(stack_dir, 'cohStack.vrt')
    unw_stack = os.path.join(stack_dir, 'unwrapStack.vrt')

    # determine the list of IFG keys exclusively from the unwrappedPhase
    # stack VRT.  absence of a valid stack is considered an error.
    date_keys = _keys_from_stack(unw_stack)
    if not date_keys:
        raise FileNotFoundError('No unwrappedPhase stack VRTs found in ARIA directory')
    # sort according to canonical ordering so earliest IFG is first regardless
    # of whether the input key was ``y1_y2`` or ``y2_y1``
    date_keys = sorted(date_keys, key=_canonical_ifg_key)

    # Determine coherence threshold based on Satellite.
    wl = _wavelength_from_stack(unw_stack)
    # Estimate satellite from wavelength: C-band <= 0.1m, else L-band
    is_cband = (wl is None) or (wl > 0 and wl <= 0.1)

    if cc_thresh is None:
        if is_cband:
            cc_thresh = 0.5  # Sentinel-1 default because of filtered coherence
        else:
            cc_thresh = 0.2  # NISAR default because of unfiltered coherence
        print(f'Auto-detected cc_thresh={cc_thresh} (wavelength={wl}m if known)')
    else:
        print(f'Using user-specified cc_thresh={cc_thresh}')

    print(f'Found {len(date_keys)} IFGs to process')

    processed_any = False
    for key in date_keys:
        print(f'Processing {key}...')
        # open the required IFG bands from the stacks; missing bands are
        # reported and the IFG is skipped.
        ds_unw, band_unw = _open_band(unw_stack, key)
        ds_cc, band_cc = _open_band(coh_stack, key)
        if not (ds_unw and ds_cc):
            print(f'  WARN: missing unw or cc for {key}. Skipping')
            if ds_unw:
                ds_unw = None
            if ds_cc:
                ds_cc = None
            continue

        unw_arr = ds_unw.GetRasterBand(band_unw).ReadAsArray().astype(np.float32)
        cc_arr = ds_cc.GetRasterBand(band_cc).ReadAsArray().astype(np.float32)
        # scale cc if 0-255
        if np.nanmax(cc_arr) > 1.5:
            cc_arr = cc_arr / 255.0

        # ARIA unwrappedPhase uses opposite sign compared to GEOC/ASF
        unw_arr = -1.0 * unw_arr

        mask = (cc_arr < cc_thresh) | np.isnan(cc_arr)
        unw_arr[mask] = 0.0

        gt = ds_unw.GetGeoTransform()
        proj = ds_unw.GetProjection()

        # construct output key with earlier date first; this handles both the
        # traditional ``y2_y1`` ordering and the reversed ``y1_y2`` form.
        # ``_canonical_ifg_key`` will leave other
        # non-date-like keys untouched.
        out_key = _canonical_ifg_key(key)

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

        if not is_cband:
            # For left-lokking NISAR/L-band
            E = -1.0 * E
            N = -1.0 * N

        gt_phi = ds_phi.GetGeoTransform()
        proj_phi = ds_phi.GetProjection()

        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.E.tif"), E, gt_phi, proj_phi)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.N.tif"), N, gt_phi, proj_phi)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.U.tif"), U, gt_phi, proj_phi)
        print('Wrote ENU unit maps')
    else:
        print('WARN: azimuthAngle or incidenceAngle missing for top-level ENU creation')

    # mli: read the first-band of the amplitude stack.  if it's missing we
    # just warn since we already processed IFGs above.
    src_mli_ds, band_mli = _open_band(amp_stack, first_key)
    if src_mli_ds:
        mli_arr = src_mli_ds.GetRasterBand(band_mli).ReadAsArray().astype(np.float32)
        write_geotiff(os.path.join(out_geoc, f"{frameid}.geo.mli.tif"),
                      mli_arr,
                      src_mli_ds.GetGeoTransform(),
                      src_mli_ds.GetProjection())
        print('Wrote mli')
        src_mli_ds = None
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

    # metadata: read center_time and radar_freq from stack metadata
    center_time = _utc_time_from_stack(unw_stack)
    # compute radar frequency from wavelength (already read above)
    radar_freq = None
    if wl and wl > 0:
        c = 2.99792458e8  # speed of light m/s
        radar_freq = c / wl
    if radar_freq is None:
        # fall back to default sentinel-1 frequency
        radar_freq = 5.405e9
        print('WARN: could not determine wavelength from stack; using default radar_freq')

    metadata_path = os.path.join(out_geoc, 'metadata.txt')
    with open(metadata_path, 'w') as mf:
        if center_time:
            mf.write(f"center_time={center_time}\n")
        mf.write(f"radar_freq={radar_freq}\n")
    print(f'Wrote metadata: {metadata_path}')


if __name__ == '__main__':
    start = time.time()
    prog = os.path.basename(sys.argv[0])
    print(f"\n{prog} ver2.0.1 20260829 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    # Show argparse defaults in help messages
    p = argparse.ArgumentParser(description='Convert ARIA-style folder to LiCSBAS GEOC layout',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-i', '--indir', default='.', help='ARIA parent directory')
    p.add_argument('-o', '--outdir', default='GEOC', help='Target GEOC directory')
    p.add_argument('-t', '--cc_thresh', type=float, default=None, help='Coherence threshold for masking unw; 0.5 for S1 filtered coherence, 0.2 for NISAR unfiltered coherence')
    p.add_argument('--overwrite', action='store_true', help='Overwrite existing outputs')
    args = p.parse_args()

    main(args.indir, args.outdir, args.cc_thresh, args.overwrite)

    # Finish
    elapsed_time = datetime.timedelta(seconds=(time.time()-start))
    print(f"\nElapsed time: {elapsed_time}")
    print(f'\n{prog} Successfully finished!!\n')

    print(f"Output directory: {args.outdir}\n")