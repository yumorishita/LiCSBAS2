#!/usr/bin/env python3
"""
Python3 library of plot functions for LiCSBAS.

v1.3.3 20250918 Yu Morishita

"""
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import numpy as np
import datetime as dt

import warnings
import matplotlib as mpl
with warnings.catch_warnings(): ## To silence user warning
    warnings.simplefilter('ignore', UserWarning)
    mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import dates as mdates

import LiCSBAS_tools_lib as tools_lib
import LiCSBAS_inv_lib as inv_lib


#%%
def make_im_png(data, pngfile, cmap, title, vmin=None, vmax=None, cbar=True):
    """
    Make png image.
    cmap can be 'insar'. To wrap data, np.angle(np.exp(1j*x/cycle)*cycle)
    """

    if cmap=='insar':
        cdict = tools_lib.cmap_insar()
        plt.register_cmap(cmap=mpl.colors.LinearSegmentedColormap('insar', cdict))
        interp = 'nearest'
    else:
        interp = 'nearest' #'antialiased'

    length, width = data.shape
    figsizex = 8
    xmergin = 2 if cbar else 0
    figsizey = int((figsizex-xmergin)*(length/width))+1

    if data.dtype==np.float32:
        data[data==0] = np.nan

    ### Plot
    fig, ax = plt.subplots(1, 1, figsize=(figsizex, figsizey))
    plt.tight_layout()

    im = ax.imshow(data, vmin=vmin, vmax=vmax, cmap=cmap, interpolation=interp)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(title)
    if cbar: fig.colorbar(im)

    plt.savefig(pngfile)
    plt.close(fig)

    return


#%%
def make_3im_png(data3, pngfile, cmap, title3, vmin=None, vmax=None, cbar=True):
    """
    Make png with 3 images for comparison.
    data3 and title3 must be list with 3 elements.
    cmap can be 'insar'. To wrap data, np.angle(np.exp(1j*x/cycle)*cycle)
    """
    ### Plot setting
    if cmap=='insar':
        cdict = tools_lib.cmap_insar()
        plt.register_cmap(cmap=mpl.colors.LinearSegmentedColormap('insar', cdict))
        interp = 'nearest'
    else:
        interp = 'nearest' #'antialiased'

    length, width = data3[0].shape
    figsizex = 12
    xmergin = 4 if cbar else 0
    figsizey = int((figsizex-xmergin)/3*length/width)+2

    fig = plt.figure(figsize = (figsizex, figsizey))

    for i in range(3):
        if data3[i].dtype==np.float32:
            data3[i][data3[i]==0] = np.nan # as no data
        ax = fig.add_subplot(1, 3, i+1) #index start from 1
        im = ax.imshow(data3[i], vmin=vmin, vmax=vmax, cmap=cmap, interpolation=interp)
        data3[i] = [] # delete to save memory
        ax.set_title(title3[i])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if cbar: fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig(pngfile)
    plt.close(fig)

    return


#%%
def plot_gacos_info(gacos_infofile, pngfile):
    figsize = (7, 3) #3x7
    sizec, colorc, markerc, alphac = 2, 'k', 'o', 0.8

    ### Figure
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 2, 1) #index start from 1
    ax2 = fig.add_subplot(1, 2, 2) #index start from 1

    ### Read data
    with open(gacos_infofile, "r") as f:
        info = f.readlines()[1:]

    std_bf, std_af, rate = [], [], [];

    for line in info:
        date, std_bf1, std_af1, rate1 = line.split()
        if std_bf1=='0.0' or std_bf1=='nan' or std_af1=='0.0' or std_af1=='nan':
            continue
        std_bf.append(float(std_bf1))
        std_af.append(float(std_af1))
        rate.append(float(rate1[:-1]))

    std_bf = np.array(std_bf)
    std_af = np.array(std_af)
    rate = np.array(rate)
    rate[rate>99] = 99
    rate[rate< -99] = -99

    ### Plot
    xylim1 = np.max(np.concatenate((std_bf, std_af)))+1
    ax1.scatter(std_bf, std_af, s=sizec, c=colorc, marker=markerc, alpha=alphac, zorder=4)
    ax1.set_xlim(0, xylim1)
    ax1.set_ylim(0, xylim1)
    ax1.plot([0, xylim1], [0, xylim1], linewidth=2, color='grey', alpha=0.5, zorder=2)
    ax1.grid(zorder=0)
    ax1.set_xlabel('STD before GACOS (rad)')
    ax1.set_ylabel('STD after GACOS (rad)')

    ### Plot
    ax2.scatter(std_bf, rate, s=sizec, c=colorc, marker=markerc, alpha=alphac, zorder=4)
    ax2.plot([0, xylim1], [0, 0], linewidth=2, color='grey', alpha=0.5, zorder=2)
    ax2.grid(zorder=0)
    ax2.set_xlim(0, xylim1)
    ax2.set_xlabel('STD before GACOS (rad)')
    ax2.set_ylabel('STD reduction rate (%)')

    fig.tight_layout()
    fig.savefig(pngfile)


#%%
def plot_hgt_corr(data_bf, fit_hgt, hgt, title, pngfile):
    # Select non-NaN data
    bool_nan = np.isnan(data_bf)
    valid_hgt = hgt[~bool_nan]
    valid_data_bf = data_bf[~bool_nan]
    valid_fit_hgt = fit_hgt[~bool_nan]

    # Calculate the data after correction
    data_af = valid_data_bf - valid_fit_hgt

    # Calculate the endpoints for the correction line
    ix_hgt0 = np.nanargmin(valid_hgt)
    ix_hgt1 = np.nanargmax(valid_hgt)
    hgt0 = valid_hgt[ix_hgt0]
    hgt1 = valid_hgt[ix_hgt1]
    fit_hgt0 = valid_fit_hgt[ix_hgt0]
    fit_hgt1 = valid_fit_hgt[ix_hgt1]

    # --- Plotting setup ---
    # Create a 1-row, 2-column subplot to display before and after correction side-by-side
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4), sharey=True)
    ax1, ax2 = axes

    # Set Y-axis limits to ensure all data is visible
    # Find the min and max across BOTH 'before' and 'after' datasets
    ymin = min(np.min(valid_data_bf), np.min(data_af))
    ymax = max(np.max(valid_data_bf), np.max(data_af))
    if ymax == ymin:
        ymax = ymax + 1
        ymin = ymin - 1

    # Add a 5% padding to the top and bottom for better visualization
    padding = (ymax - ymin) * 0.05

    # --- Plot before correction (left side) ---
    # Create a 2D histogram (heatmap)
    # 'bins' determines the resolution of the histogram.
    # cmin=1 ensures that bins with no points are not colored.
    h1 = ax1.hist2d(valid_hgt, valid_data_bf, bins=100, cmap='viridis', cmin=1)
    fig.colorbar(h1[3], ax=ax1, label='Point Density') # Add a colorbar
    ax1.plot([hgt0, hgt1], [fit_hgt0, fit_hgt1], linewidth=2, color='r', alpha=0.8, zorder=8, label='Correction')
    ax1.grid(zorder=0)
    ax1.set_title('Before Correction', fontsize=10)
    ax1.set_xlabel('Height (m)')
    ax1.set_ylabel('Displacement (mm)')
    ax1.legend(loc='upper right')

    # --- Plot after correction (right side) ---
    h2 = ax2.hist2d(valid_hgt, data_af, bins=100, cmap='viridis', cmin=1)
    fig.colorbar(h2[3], ax=ax2, label='Point Density') # Add a colorbar
    ax2.grid(zorder=0)
    ax2.set_title('After Correction', fontsize=10)
    ax2.set_xlabel('Height (m)')
    ax2.set_ylim(ymin - padding, ymax + padding)

    # Overall title and layout adjustment
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to prevent suptitle overlap
    fig.savefig(pngfile)
    plt.close(fig)

    return


#%%
def plot_network(ifgdates, bperp, rm_ifgdates, pngfile, plot_bad=True):
    """
    Plot network of interferometric pairs.

    bperp can be dummy (-1~1).
    Suffix of pngfile can be png, ps, pdf, or svg.
    plot_bad
        True  : Plot bad ifgs by red lines
        False : Do not plot bad ifgs
    """

    imdates_all = tools_lib.ifgdates2imdates(ifgdates)
    n_im_all = len(imdates_all)
    imdates_dt_all = np.array(([dt.datetime.strptime(imd, '%Y%m%d') for imd in imdates_all])) ##datetime

    ifgdates = list(set(ifgdates)-set(rm_ifgdates))
    ifgdates.sort()
    imdates = tools_lib.ifgdates2imdates(ifgdates)
    n_im = len(imdates)
    imdates_dt = np.array(([dt.datetime.strptime(imd, '%Y%m%d') for imd in imdates])) ##datetime

    ### Identify gaps
    G = inv_lib.make_sb_matrix(ifgdates)
    ixs_inc_gap = np.where(G.sum(axis=0)==0)[0]

    ### Plot fig
    figsize_x = np.round(((imdates_dt_all[-1]-imdates_dt_all[0]).days)/80)+2
    fig = plt.figure(figsize=(figsize_x, 6))
    ax = fig.add_axes([0.06, 0.12, 0.92,0.85])

    ### IFG blue lines
    for i, ifgd in enumerate(ifgdates):
        ix_m = imdates_all.index(ifgd[:8])
        ix_s = imdates_all.index(ifgd[-8:])
        label = 'IFG' if i==0 else '' #label only first
        plt.plot([imdates_dt_all[ix_m], imdates_dt_all[ix_s]], [bperp[ix_m],
                bperp[ix_s]], color='b', alpha=0.6, zorder=2, label=label)

    ### IFG bad red lines
    if plot_bad:
        for i, ifgd in enumerate(rm_ifgdates):
            ix_m = imdates_all.index(ifgd[:8])
            ix_s = imdates_all.index(ifgd[-8:])
            label = 'Removed IFG' if i==0 else '' #label only first
            plt.plot([imdates_dt_all[ix_m], imdates_dt_all[ix_s]], [bperp[ix_m],
                    bperp[ix_s]], color='r', alpha=0.6, zorder=6, label=label)

    ### Image points and dates
    ax.scatter(imdates_dt_all, bperp, alpha=0.6, zorder=4)
    for i in range(n_im_all):
        if bperp[i] > np.median(bperp): va='bottom'
        else: va = 'top'
        ax.annotate(imdates_all[i][4:6]+'/'+imdates_all[i][6:],
                    (imdates_dt_all[i], bperp[i]), ha='center', va=va, zorder=8)

    ### gaps
    if len(ixs_inc_gap)!=0:
        gap_dates_dt = []
        for ix_gap in ixs_inc_gap:
            ddays_td = imdates_dt[ix_gap+1]-imdates_dt[ix_gap]
            gap_dates_dt.append(imdates_dt[ix_gap]+ddays_td/2)
        plt.vlines(gap_dates_dt, 0, 1, transform=ax.get_xaxis_transform(),
                   zorder=1, label='Gap', alpha=0.6, colors='k', linewidth=3)

    ### Locater
    loc = ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    try:  # Only support from Matplotlib 3.1
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    except:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
        for label in ax.get_xticklabels():
            label.set_rotation(20)
            label.set_horizontalalignment('right')
    ax.grid(which='major')

    ### Add bold line every 1yr
    ax.xaxis.set_minor_locator(mdates.YearLocator())
    ax.grid(which='minor', linewidth=2)

    ax.set_xlim((imdates_dt_all[0]-dt.timedelta(days=10),
                 imdates_dt_all[-1]+dt.timedelta(days=10)))

    ### Labels and legend
    plt.xlabel('Time')
    if np.all(np.abs(np.array(bperp))<=1): ## dummy
        plt.ylabel('dummy')
    else:
        plt.ylabel('Bperp [m]')

    plt.legend()

    ### Save
    plt.savefig(pngfile)
    plt.close(fig)
