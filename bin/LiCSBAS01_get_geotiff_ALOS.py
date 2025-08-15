#!/usr/bin/env python3
"""
============
Output files
============
 - GEOC/
   - yyyymmdd_yyyymmdd/
     - yyyymmdd_yyyymmdd.geo.unw.tif
     - yyyymmdd_yyyymmdd.geo.cc.tif
   - *.geo.mli.tif (using just one first epoch)
   - *.geo.E.tif
   - *.geo.N.tif
   - *.geo.U.tif
   - *.geo.hgt.tif
   - baselines
   - metadata.txt
   - network.png

"""

#%% download_wrapper
def download_wrapper(args):
    i, ifgd, n_dl, url_data, path_data = args
    dir_data = os.path.dirname(path_data)
    print('  Donwnloading {} ({}/{})...'.format(ifgd, i+1, n_dl), flush=True)
    if not os.path.exists(dir_data): os.mkdir(dir_data)
    tools_lib.download_data(url_data, path_data)
    return


#%% check_exist_wrapper
def check_exist_wrapper(args):
    """
    Returns :
        0 : Local exist, complete, and new (no need to donwload)
        1 : Local incomplete (need to re-donwload)
        2 : Local old (no need to re-donwload)
        3 : Remote not exist  (can not compare, no download)
        4 : Local not exist (need to download)
    """

    i, n_data, url_data, path_data = args
    bname_data = os.path.basename(path_data)

    if os.path.exists(path_data):
        rc = tools_lib.comp_size_time(url_data, path_data)
        if rc == 1:
            print("Size of {} is not identical.".format(bname_data), flush=True)
        elif rc == 2:
            print("Newer {} available.".format(bname_data), flush=True)
        return rc
    else:
        return 4

#%% Import
import argparse
import os
import re
import sys
import time
import datetime
import subprocess
from osgeo import gdal
import multiprocessing as multi
import LiCSBAS_tools_lib as tools_lib


#%% Main
def main(argv=None):

    # %% Read arg
    start = time.time()
    prog = os.path.basename(sys.argv[0])
    description = 'Download AIST ALOS InSAR products for LiCSBAS.'
    print(f"\n{prog} ver1.0.1 20250815 Y. Morishita")
    print(f"{prog} {' '.join(sys.argv[1:])}\n")

    parser = argparse.ArgumentParser(description=description)
    addarg = parser.add_argument
    addarg('-f', '--frameID', type=str, help='Frame id (e.g., 045_2700_343)')
    addarg('-s', '--startdate', type=int, default=20060201,
           help='Start date (default: %(default)s)')
    addarg('-e', '--enddate', type=int, default=20110801,
           help='End date (default: %(default)s)')
    addarg('-p', '--n_para', type=int, default=4,
           help='Number of parallel downloading (default: %(default)s)')
    addarg('-u', '--unwrate_min', type=float, default=0.3,
           help='Minimum unwrap rate to be downloaded (default: %(default)s)')
    addarg('-b', '--bperp_max', type=int, default=10000,
           help='Maximum bperp (m) to be downloaded (default: %(default)s)')
    addarg('-t', '--dt_max', type=int, default=2000,
           help='Maximum time span (day) to be downloaded (default: %(default)s)')

    args = parser.parse_args()
    frameID = args.frameID
    startdate = args.startdate
    enddate = args.enddate
    n_para = args.n_para
    unwrate_min = args.unwrate_min
    bperp_max = args.bperp_max
    dt_max = args.dt_max

    # %% Setting
    q = multi.get_context('fork')
    gunw_url = 'https://gsvrg.ipri.aist.go.jp/pds/palsar-insar-pds/P1INSAR/GUNW/'
    doc_url = 'https://gsvrg.ipri.aist.go.jp/insarbrowser/doc/'

    #%% Determine frameID
    wd = os.getcwd()
    if not frameID: ## if frameID not indicated
        _tmp = re.findall(r'\d{3}_\d{4}_\d{3}', wd)
        ##e.g., 406_0710_343
        if len(_tmp)==0:
            raise Exception('Frame ID cannot be identified from dir name! '
                            'Use -f option')
        else:
            frameID = _tmp[0]
            print('\nFrame ID is {}\n'.format(frameID), flush=True)
    else:
        print('\nFrame ID is {}\n'.format(frameID), flush=True)


    #%% Directory and file setting
    outdir = os.path.join(wd, 'GEOC')
    if not os.path.exists(outdir): os.mkdir(outdir)
    os.chdir(outdir)


    # %% baselines, network.png, products_list.txt and unwrap_rates_list.txt
    url = os.path.join(gunw_url, frameID, f'{frameID}_GUNW.baselines')
    print(f'Download {url} as baselines')
    tools_lib.download_data(url, 'baselines')

    url = os.path.join(doc_url, frameID, 'network.png')
    print(f'Download {url} as network.png')
    tools_lib.download_data(url, 'network.png')

    url = os.path.join(doc_url, frameID, 'products_list.txt')
    print(f'Download {url} as products_list.txt')
    tools_lib.download_data(url, 'products_list.txt')

    url = os.path.join(doc_url, frameID, 'unwrap_rates_list.txt')
    print(f'Download {url} as unwrap_rates_list.txt')
    tools_lib.download_data(url, 'unwrap_rates_list.txt')

    print('', flush=True)

    # %% Create metadata.txt containing radar_freq
    print(f'Download metadata.txt')
    with open('metadata.txt', 'w') as f:
        print(f'radar_freq=1.27e9', file=f)

    # %% Read baselines, unwrap_rate, and products_list
    pdates = []
    sdates = []
    bperps = []
    dts = []
    with open('baselines', 'r') as f:
        for line in f.readlines():
            line = line.split()[1:5]
            pdates.append(int(line[0]))
            sdates.append(int(line[1]))
            bperps.append(float(line[2]))
            dts.append(float(line[3]))

    ifgdates_all = []
    unw_rates = []
    with open('unwrap_rates_list.txt', 'r') as f:
        for line in f.readlines():
            line = line.split(',')
            ifgdates_all.append(line[0])
            unw_rates.append(float(line[1]))

    unw_rate_dict = dict(zip(ifgdates_all, unw_rates))
    bperp_dict = dict(zip(ifgdates_all, bperps))
    dt_dict = dict(zip(ifgdates_all, dts))

    print(f'\nNumber of all available ifgs: {len(ifgdates_all)}')

    with open('products_list.txt', 'r') as f:
        products_urls = [ url.rstrip() for url in f.readlines() ]

    for _url in products_urls:
        if '_GUNW.txt' in _url:
            sceneID = os.path.basename(_url)[:17]
            print(f'\nScene ID is {sceneID}')
            break
    print('')


    # %% ENU and hgt for only 1st ifg
    ifgd1 = ifgdates_all[0]

    for ENU in ['E', 'N', 'U', 'hgt']:
        if ENU in ['E', 'N', 'U']:
            enutif_local = f'{frameID}.geo.{ENU}inv.tif'
            ENU = 'los'+ENU # for remote file name
        else:
            enutif_local = f'{frameID}.geo.{ENU}.tif'
        enutif = f'{sceneID}_{ifgd1}_GUNW_{ENU}.tif'
        url = os.path.join(gunw_url, frameID, ifgd1, enutif)
        if os.path.exists(enutif_local):
            rc = tools_lib.comp_size_time(url, enutif_local)
            if rc == 0:
                print('{} already exist. Skip download.'.format(enutif_local), flush=True)
                continue
            elif rc == 3:
                print('{} not available. Skip download.'.format(url), flush=True)
                continue
            else:
                if rc == 1:
                    print("Size of {} is not identical.".format(enutif))
                elif rc == 2:
                    print("Newer {} available.".format(enutif))

        print(f'Download {url} as {enutif_local}')
        tools_lib.download_data(url, enutif_local)

    print('')

    # %% Invert sign of ENU because the definition is opposite
    for ENU in ['E', 'N', 'U']:
        enutif_inv = f'{frameID}.geo.{ENU}inv.tif'
        enutif = f'{frameID}.geo.{ENU}.tif'
        call = ['gdal_calc.py', '--calc="-A"', '-A', enutif_inv,
                f'--outfile={enutif}', '--overwrite', '--co="COMPRESS=DEFLATE"',
                '--NoDataValue=0']
        subprocess.run(' '.join(call), shell=True, check=True)

    #%% amp as mli only for 1st epoch
    mlitif = frameID+'.geo.mli.tif'
    if os.path.exists(mlitif):
        print('{} already exist. Skip.'.format(mlitif), flush=True)
    else:
        ### Download
        for _url in products_urls:
            if '_GUNW_amp.tif' in _url:
                url_mli = _url
                break
        print(f'Donwnloading {url_mli}...')
        tools_lib.download_data(url_mli, mlitif)

    print('')

    #%% unw and cc
    ### Get available dates
    print('\nDownload geotiff of unw and cc', flush=True)

    ### Extract during start_date to end_date
    ifgdates = []
    for ifgd in ifgdates_all:
        mimd = int(ifgd[:8])
        simd = int(ifgd[-8:])
        if mimd >= startdate and simd <= enddate:
            ifgdates.append(ifgd)

    n_ifg = len(ifgdates)
    imdates = tools_lib.ifgdates2imdates(ifgdates)
    print('{} IFGs available from {} to {}'.format(n_ifg, imdates[0], imdates[-1]), flush=True)

    ### Check if both unw and cc already donwloaded, new, and same size
    print('Checking existing unw and cc ({} parallel, may take time)...'.format(n_para), flush=True)

    ## unw
    args = [(i, n_ifg,
             os.path.join(gunw_url, frameID, ifgd, f'{sceneID}_{ifgd}_GUNW_unw.tif'),
             os.path.join(ifgd, '{}.geo.unw.tif'.format(ifgd))
             ) for i, ifgd in enumerate(ifgdates)]

    p = q.Pool(n_para)
    rc = p.map(check_exist_wrapper, args)
    p.close()

    n_unw_existing = 0
    n_unw_bad = 0
    unwdates_dl = []
    for i, rc1 in enumerate(rc):
        ifgd = ifgdates[i]

        if unw_rate_dict[ifgd] < unwrate_min*100:
            print(f'  Not download {ifgd} because unw_rate'
                  f' {unw_rate_dict[ifgd]}% < {unwrate_min*100}%')
            n_unw_bad = n_unw_bad + 1
        elif abs(bperp_dict[ifgd]) > bperp_max:
            print(f'  Not download {ifgd} because bperp'
                  f' {bperp_dict[ifgd]} > {bperp_max}')
            n_unw_bad = n_unw_bad + 1
        elif dt_dict[ifgd] > dt_max:
            print(f'  Not download {ifgd} because dt'
                  f' {dt_dict[ifgd]} > {dt_max}')
            n_unw_bad = n_unw_bad + 1
        elif rc1 == 0:  ## No need to download
            n_unw_existing = n_unw_existing + 1
        elif rc1 == 3 or rc1 == 5:  ## Can not download
            print('  {}.geo.unw.tif not available.'.format(ifgd), flush=True)
        elif rc1 == 1 or rc1 == 2  or rc1 == 4:  ## Need download
            unwdates_dl.append(ifgd)

    ## cc
    args = [(i, n_ifg,
             os.path.join(gunw_url, frameID, ifgd, f'{sceneID}_{ifgd}_GUNW_coh.tif'),
             os.path.join(ifgd, '{}.geo.cc.tif'.format(ifgd))
             ) for i, ifgd in enumerate(ifgdates)]

    p = q.Pool(n_para)
    rc = p.map(check_exist_wrapper, args)
    p.close()

    n_cc_existing = 0
    ccdates_dl = []
    for i, rc1 in enumerate(rc):
        if unw_rate_dict[ifgdates[i]] < unwrate_min*100 or \
            abs(bperp_dict[ifgd]) > bperp_max or dt_dict[ifgd] > dt_max:
            pass
        elif rc1 == 0:  ## No need to download
            n_cc_existing = n_cc_existing + 1
        elif rc1 == 3 or rc1 == 5:  ## Can not download
            print('  {}.geo.cc.tif not available.'.format(ifgdates[i]), flush=True)
        elif rc1 == 1 or rc1 == 2  or rc1 == 4:  ## Need download
            ccdates_dl.append(ifgdates[i])

    n_unw_dl = len(unwdates_dl)
    n_cc_dl = len(ccdates_dl)
    print(f'{n_unw_existing} unw already downloaded')
    print(f'{n_unw_bad} unw rate < {unwrate_min*100}')
    print(f'{n_unw_dl} unw will be downloaded')
    print(f'{n_cc_existing} cc already downloaded')
    print(f'{n_cc_dl} cc will be downloaded')

    ### Download unw with parallel
    if n_unw_dl != 0:
        print('\nDownload unw ({} parallel)...'.format(n_para), flush=True)
        args = [(i, ifgd, n_unw_dl,
                 os.path.join(gunw_url, frameID, ifgd, f'{sceneID}_{ifgd}_GUNW_unw.tif'),
                 os.path.join(ifgd, '{}.geo.unw.tif'.format(ifgd))
                 ) for i, ifgd in enumerate(unwdates_dl)]

        p = q.Pool(n_para)
        p.map(download_wrapper, args)
        p.close()

    ### Download cc with parallel
    if n_cc_dl != 0:
        print('\nDownload cc ({} parallel)...'.format(n_para), flush=True)
        args = [(i, ifgd, n_cc_dl,
                 os.path.join(gunw_url, frameID, ifgd, f'{sceneID}_{ifgd}_GUNW_coh.tif'),
                 os.path.join(ifgd, '{}.geo.cc.tif'.format(ifgd))
                 ) for i, ifgd in enumerate(ccdates_dl)]

        p = q.Pool(n_para)
        p.map(download_wrapper, args)
        p.close()


    #%% Finish
    elapsed_time = datetime.timedelta(seconds=(time.time()-start))
    print(f"\nElapsed time: {elapsed_time}")
    print(f'\n{prog} Successfully finished!!\n')

    print(f"Output directory: {outdir}\n")


#%% main
if __name__ == "__main__":
    sys.exit(main())

