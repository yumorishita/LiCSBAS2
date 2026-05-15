#!/usr/bin/env python3
"""
Align two GEOC folders by resampling the second folder's geotiffs to match the geometry of the first folder.

Inputs:
- Two GEOC folders with different resolutions and date pairs
- Reference geotiff: unw.geo.unw.tif from the first pair folder

Outputs:
- Aligned GEOC folder with unified geometry and integrated date pairs
- GEOC-level geotiffs copied from first folder
- metadata.txt copied from first folder
- baselines: TODO (not implemented yet)

"""

import argparse
import glob
import math
import os
import re
import shutil
import sys
import time

from osgeo import gdal
import numpy as np

import LiCSBAS_io_lib as io_lib
import LiCSBAS_tools_lib as tools_lib


def main(geoc1_dir, geoc2_dir, out_geoc_dir, overwrite=False):
    geoc1_dir = os.path.abspath(geoc1_dir)
    geoc2_dir = os.path.abspath(geoc2_dir)
    out_geoc_dir = os.path.abspath(out_geoc_dir)

    if not os.path.isdir(geoc1_dir):
        raise FileNotFoundError(f"{geoc1_dir} not found")
    if not os.path.isdir(geoc2_dir):
        raise FileNotFoundError(f"{geoc2_dir} not found")

    os.makedirs(out_geoc_dir, exist_ok=True)

    # Get reference geotiff from first pair in geoc1
    pairs1 = tools_lib.get_pair_folders(geoc1_dir)
    if not pairs1:
        raise ValueError(f"No pair folders in {geoc1_dir}")
    ref_pair = pairs1[0]
    ref_tif = os.path.join(geoc1_dir, ref_pair, f"{ref_pair}.geo.unw.tif")
    if not os.path.exists(ref_tif):
        raise FileNotFoundError(f"Reference tif {ref_tif} not found")

    gt1, proj1, shape1 = io_lib.get_geotiff_info(ref_tif)

    # Get pairs from geoc2
    pairs2 = tools_lib.get_pair_folders(geoc2_dir)

    # Determine common geometry using geoc2's first pair if exists
    if pairs2:
        ref_tif2 = os.path.join(geoc2_dir, pairs2[0], f"{pairs2[0]}.geo.unw.tif")
        if os.path.exists(ref_tif2):
            gt2, _, shape2 = io_lib.get_geotiff_info(ref_tif2)
            common_gt, common_shape = tools_lib.calculate_common_geometry(gt1, shape1, gt2, shape2)
        else:
            common_gt, common_shape = gt1, shape1
    else:
        common_gt, common_shape = gt1, shape1

    # Collect all unique pairs
    all_pairs = set(pairs1 + pairs2)

    for pair in sorted(all_pairs):
        out_pair_dir = os.path.join(out_geoc_dir, pair)
        os.makedirs(out_pair_dir, exist_ok=True)

        # Determine source: prefer geoc1
        if pair in pairs1:
            src_dir = geoc1_dir
        elif pair in pairs2:
            src_dir = geoc2_dir
        else:
            continue

        src_tif = os.path.join(src_dir, pair, f"{pair}.geo.unw.tif")
        if not os.path.exists(src_tif):
            continue

        out_tif = os.path.join(out_pair_dir, f"{pair}.geo.unw.tif")
        if os.path.exists(out_tif) and not overwrite:
            print(f"Skipping {out_tif} (exists)")
            continue

        success = io_lib.resample_geotiff(src_tif, out_tif, common_gt, common_shape, proj1)
        if success:
            print(f"Created {out_tif}")
        else:
            print(f"Failed to create {out_tif}")

        # Handle other tifs: cc, etc. if exist
        for ext in ['cc', 'mli']:  # add more if needed
            src_other = os.path.join(src_dir, pair, f"{pair}.geo.{ext}.tif")
            if os.path.exists(src_other):
                out_other = os.path.join(out_pair_dir, f"{pair}.geo.{ext}.tif")
                io_lib.resample_geotiff(src_other, out_other, common_gt, common_shape, proj1)

    # Resample GEOC-level geotiffs from geoc1 to the common output geometry
    for tif in glob.glob(os.path.join(geoc1_dir, "*.geo.*.tif")):
        basename = os.path.basename(tif)
        out_tif = os.path.join(out_geoc_dir, basename)
        if os.path.exists(out_tif) and not overwrite:
            print(f"Skipping {out_tif} (exists)")
            continue

        if io_lib.resample_geotiff(tif, out_tif, common_gt, common_shape, proj1):
            print(f"Resampled {basename}")
        else:
            print(f"Failed to resample {basename}")

    # Copy metadata.txt from geoc1
    metadata_src = os.path.join(geoc1_dir, "metadata.txt")
    metadata_dst = os.path.join(out_geoc_dir, "metadata.txt")
    if os.path.exists(metadata_src):
        shutil.copy2(metadata_src, metadata_dst)
        print("Copied metadata.txt")

    # baselines: TODO


if __name__ == '__main__':
    start = time.time()
    prog = os.path.basename(sys.argv[0])
    print(f"\n{prog} ver1.0.0 20260515 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    p = argparse.ArgumentParser(description='Align two GEOC folders to unified geometry',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-g1', '--geoc1', required=True, help='First GEOC directory (reference)')
    p.add_argument('-g2', '--geoc2', required=True, help='Second GEOC directory')
    p.add_argument('-o', '--outdir', required=True, help='Output GEOC directory')
    p.add_argument('--overwrite', action='store_true', help='Overwrite existing outputs')
    args = p.parse_args()

    main(args.geoc1, args.geoc2, args.outdir, args.overwrite)

    # Finish
    elapsed = time.time() - start
    print(f"\n{prog} finished in {elapsed:.2f} sec")