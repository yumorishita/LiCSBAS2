#!/usr/bin/env python3
"""
This script creates a colored GeoTIFF from a data GeoTIFF.

=====
Usage
=====
LiCSBAS_color_geotiff.py -i infile [-c cmap] [-o outfile] [--cmin float] [--cmax float] [--n_color int] [--no_colorbar]

 -i  Input data GeoTIFF file
 -c  Colormap name (Default: SCM.roma_r, reverse of SCM.roma)
     Available colormaps (all cmap can be reversed with "_r"):
     - Matplotlib predefined name (e.g. viridis)
       https://matplotlib.org/tutorials/colors/colormaps.html
     - Scientific colour maps (e.g. SCM.roma)
       http://www.fabiocrameri.ch/colourmaps.php
     - Generic Mapping Tools (e.g. GMT.polar)
       https://docs.generic-mapping-tools.org/dev/cookbook/cpts.html
     - cmocean (e.g. cmocean.phase)
       https://matplotlib.org/cmocean/
     - colorcet (e.g. colorcet.CET_C1)
       https://colorcet.holoviz.org/
     - cm_insar (GAMMA standard rainbow color for wrapped phase)
     - cm_isce (ISCE standard rainbow color for wrapped phase)
 -o  Output colored GeoTIFF file (Default: [infile%.tif].cmap_cmin_cmax.tif)
 --cmin|cmax  Min|max values of color (Default: None (auto))
 --n_color    Number of rgb quantization levels (Default: 256)
 --no_colorbar   Do not create colorbar image (name: cmap_cmin_cmax.pdf)

"""

#%% Import
import getopt
import os
import sys
#import gdal  # may cause segmentation fault in savefig?
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import LiCSBAS_tools_lib as tools_lib
from osgeo import gdal

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


#%% Main
## Not use def main to use global valuables
if __name__ == "__main__":
    argv = sys.argv

    start = time.time()
    ver="1.4.1"; date=20250715; author="Y. Morishita"
    print("\n{} ver{} {} {}".format(os.path.basename(argv[0]), ver, date, author), flush=True)
    print("{} {}".format(os.path.basename(argv[0]), ' '.join(argv[1:])), flush=True)

    #%% Set default
    infile = []
    cmap_name = "SCM.roma_r"
    outfile = []
    nodata = np.nan
    cmin = None
    cmax = None
    n_color = 256
    cbar_flag = True

    gdal_option = ['COMPRESS=DEFLATE']

    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hi:c:o:", ["help", "cmin=", "cmax=", "n_color=", "no_colorbar"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                sys.exit(0)
            elif o == '-i':
                infile = a
            elif o == '-c':
                cmap_name = a
            elif o == '-o':
                outfile = a
            elif o == '--nodata':
                nodata = float(a)
            elif o == '--cmin':
                cmin = float(a)
            elif o == '--cmax':
                cmax = float(a)
            elif o == '--n_color':
                n_color = int(a)
            elif o == '--no_colorbar':
                cbar_flag = False

        if not infile:
            raise Usage('No input file given, -i is not optional!')
        elif not os.path.exists(infile):
            raise Usage('No {} exists!'.format(infile))

    except Usage as err:
        print("\nERROR:", file=sys.stderr, end='')
        print("  "+str(err.msg), file=sys.stderr)
        print("\nFor help, use -h or --help.\n", file=sys.stderr)
        sys.exit(2)


    #%% Set cmap
    cmap = tools_lib.get_cmap(cmap_name, n_color)


    #%% Set color range
    ### Auto
    if cmin is None:
        cmin = gdal.Info(infile, computeMinMax=True, format="json")["bands"][0]['computedMin']
    if cmax is None:
        cmax = gdal.Info(infile, computeMinMax=True, format="json")["bands"][0]['computedMax']

    ### float -> int if int
    if np.mod(cmin, 1) == 0: cmin = int(cmin)
    if np.mod(cmax, 1) == 0: cmax = int(cmax)


    #%% Set file name
    if not outfile:
        outfile = infile.replace('.tif', '.{}_{}_{}.tif'.format(cmap_name, cmin, cmax))


    #%% Create color table
    ### Format: value R G B alpha
    colorfile = '{}_{}_{}.txt'.format(cmap_name, cmin, cmax)
    cmap_RGB = np.int16(np.round(cmap(np.linspace(0, 1, n_color))*255))
    with open(colorfile, "w") as f:
        for i in range(n_color):
            print("{} {} {} {} 255".format(cmin+i*(cmax-cmin)/(n_color-1),
                   cmap_RGB[i, 0], cmap_RGB[i, 1], cmap_RGB[i, 2]), file=f)
        print("nv 0 0 0 0", file=f)


    #%% gdal dem
    gdal.DEMProcessing(outfile, infile, "color-relief", colorFilename=colorfile,
                       format="GTiff", creationOptions=gdal_option, addAlpha=4)


    #%% colorbar. Move before import gdal to avoid segmentation fault
    if cbar_flag:
        fig, ax = plt.subplots(figsize=(3, 1))
        norm = mpl.colors.Normalize(cmin, cmax)
        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
        cbarfile = "{}_{}_{}.pdf".format(cmap_name, cmin, cmax)
        plt.tight_layout()
        plt.savefig(cbarfile, transparent=True)
        plt.close(fig)


    #%% Remove intermediate files
    os.remove(colorfile)


    #%% Finish
    elapsed_time = time.time()-start
    hour = int(elapsed_time/3600)
    minite = int(np.mod((elapsed_time/60),60))
    sec = int(np.mod(elapsed_time,60))
    print("\nElapsed time: {0:02}h {1:02}m {2:02}s".format(hour,minite,sec))

    print('\n{} Successfully finished!!\n'.format(os.path.basename(argv[0])))
    print('Output: {}'.format(outfile), flush=True)
    if cbar_flag:
        print('        {}'.format(cbarfile), flush=True)
    print('')

    sys.exit(0)
