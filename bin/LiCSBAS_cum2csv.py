#!/usr/bin/env python3
"""
This script outputs csv file containing all time series data from cum*.h5.
The output csv can be opened in QGIS (as a delimited text layer)
and the time series can be interactively plotted by InSAR Explorer plugin
(https://insar-explorer.readthedocs.io/en/latest/).
"""

#%% Import
import argparse
import datetime
import os
import re
import sys
import time

import h5py as h5
import numpy as np

import LiCSBAS_tools_lib as tools_lib


# %% Main
def main(cumfile='cum_filt.h5',
        outcsv='cum_filt.csv',
        refarea=None,
        refarea_geo=None,
        nomask=False,
        crop_str=None,
        crop_geo_str=None):


    # %% Read cumfile
    if not os.path.exists(cumfile):
        raise FileNotFoundError(f'ERROR: {cumfile} not found!')
    cumh5 = h5.File(cumfile, 'r')
    cum = cumh5['cum']
    vel = cumh5['vel']
    imdates = cumh5['imdates'][()].astype(str).tolist()
    n_im, length, width = cum.shape

    lat1 = float(cumh5['corner_lat'][()])
    lon1 = float(cumh5['corner_lon'][()])
    dlat = float(cumh5['post_lat'][()])
    dlon = float(cumh5['post_lon'][()])


    # %% Set ref area
    if refarea:
        if not tools_lib.read_range(refarea, width, length):
            raise ValueError(f'ERROR with -r {refarea}')
        refx1, refx2, refy1, refy2 = tools_lib.read_range(refarea, width, length)
    elif refarea_geo:
        if not tools_lib.read_range_geo(refarea_geo, width, length, lat1, dlat, lon1, dlon):
            raise ValueError(f'ERROR with -ref_geo {refarea_geo}')
        refx1, refx2, refy1, refy2 = tools_lib.read_range_geo(refarea_geo, width, length, lat1, dlat, lon1, dlon)
    else:
        refarea_val = cumh5['refarea'][()]
        if type(refarea_val) is bytes:
            refarea_val = refarea_val.decode('utf-8')
        refx1, refx2, refy1, refy2 = [int(s) for s in re.split('[:/]', refarea_val)]


    # %% Mask
    if nomask or 'mask' not in cumh5.keys():
        mask = np.ones((length, width), dtype=np.float32)
    else:
        mask = cumh5['mask'][()]
        mask[mask==0] = np.nan


    # %% Read noise indices
    noise_names = ['mask', 'coh_avg', 'n_unw', 'vstd', 'maxTlen', 'n_gap', 'stc',
                   'n_ifg_noloop', 'n_loop_err', 'resid_rms', 'slc.mli', 'hgt']
    noise_names_exist = []
    noise_dict = {}
    for name in noise_names:
        if name in cumh5.keys():
            noise_names_exist.append(name)
            noise_dict[name] = cumh5[name][()]


    # %% Determine crop area if specified
    if crop_str:
        if not tools_lib.read_range(crop_str, width, length):
            raise ValueError(f'ERROR with --crop {crop_str}')
        crop_x1, crop_x2, crop_y1, crop_y2 = tools_lib.read_range(crop_str, width, length)
    elif crop_geo_str:
        if not tools_lib.read_range_geo(crop_geo_str, width, length, lat1, dlat, lon1, dlon):
            raise ValueError(f'ERROR with --crop_geo {crop_geo_str}')
        crop_x1, crop_x2, crop_y1, crop_y2 = tools_lib.read_range_geo(crop_geo_str, width, length, lat1, dlat, lon1, dlon)
    else:
        crop_x1, crop_x2, crop_y1, crop_y2 = 0, width, 0, length


    # %% Compute reference values
    ts_ref = np.nanmean(cum[:, refy1:refy2, refx1:refx2]*mask[refy1:refy2, refx1:refx2], axis=(1,2))


    # %% Output csv for each point
    with open(outcsv, 'w') as f:
        # Write header
        header = ['latitude', 'longitude', 'velocity', 'line', 'pixel']
        header += noise_names
        header += imdates
        f.write(','.join(header) + '\n')

        # Loop over cropped area
        for y in range(crop_y1, crop_y2):
            for x in range(crop_x1, crop_x2):
                if np.isnan(mask[y, x]):
                    continue

                lat, lon = tools_lib.xy2bl(x, y, lat1, dlat, lon1, dlon)
                velocity = vel[y, x]
                if np.isnan(velocity):
                    continue

                # Noise indices
                noise_vals = []
                for name in noise_names_exist:
                    arr = noise_dict[name]
                    val = arr[y, x]
                    if np.isnan(val):
                        noise_vals.append('')
                    else:
                        if name in ['mask', 'n_unw', 'n_gap', 'n_ifg_noloop', 'n_loop_err']:
                            val = str(int(val))
                        else:
                            val = f'{val:.3f}'
                        noise_vals.append(val)

                # Time series (ref area subtraction)
                ts = cum[:, y, x]
                ts_dif = ts - ts_ref
                ts_dif = ts_dif - ts_dif[0]
                ts_vals = [f'{ts_dif[i]:.1f}' if not np.isnan(ts_dif[i]) else '' for i in range(n_im)]

                # Write row
                row = [f'{lat:.6f}', f'{lon:.6f}', f'{velocity:.3f}',
                       str(y), str(x)] + noise_vals + ts_vals
                f.write(','.join(row) + '\n')


# %% if main
if __name__ == "__main__":
    # Read arg
    description = 'Output csv from cum*h5.'
    parser = argparse.ArgumentParser(description=description)
    addarg = parser.add_argument
    addarg('-i', '--input', default='cum_filt.h5',
           help='Input cum*.h5 file (default: %(default)s)')
    addarg('-o', '--output', default='cum_filt.csv',
           help='Output CSV file (default: %(default)s)')
    addarg('-r', '--refarea', default=None, help='Reference area x1:x2/y1:y2')
    addarg('--ref_geo', default=None,
           help='Reference area in geographical coordinates lon1/lon2/lat1/lat2')
    addarg('--nomask', action='store_true', help='Do not apply mask')
    addarg('--crop', default=None, help='Crop area x1:x2/y1:y2')
    addarg('--crop_geo', default=None,
           help='Crop area in geographical coordinates lon1/lon2/lat1/lat2')

    args = parser.parse_args()

    start = time.time()
    prog = os.path.basename(sys.argv[0])
    print(f"\n{prog} ver1.0.0 20250910 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    main(args.input,
        args.output,
        args.refarea,
        args.ref_geo,
        args.nomask,
        args.crop,
        args.crop_geo)

    # Finish
    elapsed_time = datetime.timedelta(seconds=(time.time()-start))
    print(f"\nElapsed time: {elapsed_time}")
    print(f'\n{prog} Successfully finished!!\n')

    print(f"Output CSV: {args.output}\n")
