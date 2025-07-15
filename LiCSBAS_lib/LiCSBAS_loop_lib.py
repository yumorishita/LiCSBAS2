#!/usr/bin/env python3
"""
Python3 library of loop closure check functions for LiCSBAS.

v1.5.2 20210303 Yu Morishita
"""

import os
import numpy as np
import LiCSBAS_io_lib as io_lib
import LiCSBAS_tools_lib as tools_lib

os.environ['QT_QPA_PLATFORM']='offscreen'
import warnings
import matplotlib as mpl
with warnings.catch_warnings(): ## To silence user warning
    warnings.simplefilter('ignore', UserWarning)
    mpl.use('Agg')
from matplotlib import pyplot as plt


#%%
def make_loop_matrix(ifgdates):
    """
    Make loop matrix (containing 1, -1, 0) from ifgdates.

    Inputs:
      ifgdates : Unwrapped phase vector at a point without nan (n_ifg)

    Returns:
      Aloop : Loop matrix with 1 for ifg12/ifg23 and -1 for ifg13
              (n_loop, n_ifg)

    """
    n_ifg = len(ifgdates)
    Aloop = []

    for ix_ifg12, ifgd12 in enumerate(ifgdates):
        primary12 = ifgd12[0:8]
        secondary12 = ifgd12[9:17]
        ifgdates23 = [ ifgd for ifgd in ifgdates if
                      ifgd.startswith(secondary12)] # all candidates of ifg23

        for ifgd23 in ifgdates23: # for each candidate of ifg23
            secondary23 = ifgd23[9:17]
            try:
                ## Search ifg13
                ix_ifg13 = ifgdates.index(primary12+'_'+secondary23)
            except: # no loop for this ifg23. Next.
                continue

            ## Loop found
            ix_ifg23 = ifgdates.index(ifgd23)

            Aline = [0]*n_ifg
            Aline[ix_ifg12] = 1
            Aline[ix_ifg23] = 1
            Aline[ix_ifg13] = -1
            Aloop.append(Aline)

    Aloop = np.array(Aloop)

    return Aloop


#%%
def read_unw_loop_ph(Aloop1, ifgdates, ifgdir, length, width, bad_ifg=[]):
    ### Find index of ifg
    ix_ifg12, ix_ifg23 = np.where(Aloop1 == 1)[0]
    ix_ifg13 = np.where(Aloop1 == -1)[0][0]
    ifgd12 = ifgdates[ix_ifg12]
    ifgd23 = ifgdates[ix_ifg23]
    ifgd13 = ifgdates[ix_ifg13]

    ### Read unw data
    unw12file = os.path.join(ifgdir, ifgd12, ifgd12+'.unw')
    unw12 = io_lib.read_img(unw12file, length, width)
    unw12[unw12 == 0] = np.nan # Fill 0 with nan
    unw23file = os.path.join(ifgdir, ifgd23, ifgd23+'.unw')
    unw23 = io_lib.read_img(unw23file, length, width)
    unw23[unw23 == 0] = np.nan # Fill 0 with nan
    unw13file = os.path.join(ifgdir, ifgd13, ifgd13+'.unw')
    unw13 = io_lib.read_img(unw13file, length, width)
    unw13[unw13 == 0] = np.nan # Fill 0 with n

    return unw12, unw23, unw13, ifgd12, ifgd23, ifgd13


#%%
def identify_bad_ifg(bad_ifg_cand, good_ifg):
    ### Identify bad ifgs and output text
    good_ifg = list(set(good_ifg))
    good_ifg.sort()
    bad_ifg_cand = list(set(bad_ifg_cand))
    bad_ifg_cand.sort()

    bad_ifg = list(set(bad_ifg_cand)-set(good_ifg)) # difference
    bad_ifg.sort()

    return bad_ifg


#%%
def make_loop_png(unw12, unw23, unw13, loop_ph, png, titles4, cycle):
    cmap_wrap = tools_lib.get_cmap('cm_insar')
    cmap_loop = tools_lib.get_cmap('SCM.vik')

    ### Settings
    plt.rcParams['axes.titlesize'] = 10
    data = [unw12, unw23, unw13]

    length, width = unw12.shape
    if length > width:
        figsize_y = 10
        figsize_x = int((figsize_y-1)*width/length)
        if figsize_x < 5: figsize_x = 5
    else:
        figsize_x = 10
        figsize_y = int(figsize_x*length/width+1)
        if figsize_y < 3: figsize_y = 3

    ### Plot
    fig = plt.figure(figsize = (figsize_x, figsize_y))

    ## 3 ifgs
    for i in range(3):
        data_wrapped = np.angle(np.exp(1j*(data[i]/cycle))*cycle)
        ax = fig.add_subplot(2, 2, i+1) #index start from 1
        im = ax.imshow(data_wrapped, vmin=-np.pi, vmax=+np.pi, cmap=cmap_wrap,
                  interpolation='nearest')
        ax.set_title('{}'.format(titles4[i]))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        cax = plt.colorbar(im)
        cax.set_ticks([])

    ## loop phase
    ax = fig.add_subplot(2, 2, 4) #index start from 1
    im = ax.imshow(loop_ph, vmin=-np.pi, vmax=+np.pi, cmap=cmap_loop,
              interpolation='nearest')
    ax.set_title('{}'.format(titles4[3]))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    cax = plt.colorbar(im)

    plt.tight_layout()
    plt.savefig(png)
    plt.close(fig)

