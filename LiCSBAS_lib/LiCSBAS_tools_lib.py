#!/usr/bin/env python3
"""
Python3 library of time series analysis tools for LiCSBAS.

v1.12.0 20260818 Yu Morishita
"""
import os
import sys
import re
import time
import math
import requests
import dateutil
import datetime as dt
import numpy as np
import warnings
from bs4 import BeautifulSoup
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap as LSC
from matplotlib import pyplot as plt


#%%
def bl2xy(lon, lat, width, length, lat1, postlat, lon1, postlon):
    """
    lat1 is north edge and postlat is negative value.
    lat lon values are in grid registration
    x/y index start from 0, end with width-1
    """
    x = int(np.round((lon - lon1)/postlon))
    y = int(np.round((lat - lat1)/postlat))

    return [x, y]


#%%
def comp_size_time(file_remote, file_local):
    """
    Compare size and time of remote and local files.
    Returns:
        0 : Size is identical and local time is new
        1 : Size is not identical
        2 : Size is identical but remote time is new
        3 : Remote not exist
    """

    response = requests.head(file_remote, allow_redirects=True)

    if response.status_code != 200:
        return 3
    elif response.headers.get("Content-Length") is None:
        ## Remote file exist but empty
        return 3

    size_remote = int(response.headers.get("Content-Length"))
    size_local = os.path.getsize(file_local)
    if size_remote != size_local: ### Different size
        return 1

    time_remote = dateutil.parser.parse(response.headers.get("Last-Modified"))
    time_local = dt.datetime.fromtimestamp(os.path.getmtime(file_local), dt.timezone.utc)
    if time_remote > time_local: ### New file
        return 2

    return 0


#%%
def cm_insar():
    """
    Perceptually uniform InSAR colormap
    """
    rgbs = np.zeros((256,3), dtype=np.uint8)

    for kk in range(106):
        rgbs[kk,0] = int(kk/106*255)
        rgbs[kk,1] = 255-rgbs[kk,0]
        rgbs[kk,2] = 255

    rgbs[106:212,0] = rgbs[0:106,2]
    rgbs[106:212,1] = rgbs[0:106,0]
    rgbs[106:212,2] = rgbs[0:106,1]

    for kk in range(43):
        rgbs[kk+212,2] = int(kk/43*255)
        rgbs[kk+212,0] = 255-rgbs[kk+212,2]
        rgbs[kk+212,1] = 255

    rgbs[255,0] = 0
    rgbs[255,1] = 255
    rgbs[255,2] = 255

    rgbs = np.roll(rgbs, int(256/2), axis=0)  #shift green to the center
    rgbs = (rgbs/255.0*200+55)/255.0

    ma_size = 51 # moving average, should be odd number
    ma_size_half = int(ma_size/2)
    rgb_ex = np.zeros((256+ma_size_half*2, 3), dtype=np.float32)
    rgb_ex[0:ma_size_half, :] = rgbs[-ma_size_half:, :]
    rgb_ex[ma_size_half:ma_size_half+256, :] = rgbs
    rgb_ex[-ma_size_half:, :] = rgbs[:ma_size_half, :]
    rgb_ma = np.zeros_like(rgbs)
    for i in range(3):
        rgb_ma[:, i] = np.convolve(rgb_ex[:, i], np.ones(ma_size)/ma_size, 'valid')

    rgbs = rgb_ma
    rgbs[rgbs>1] = 1 # values should be within 0-1 range

    return LSC.from_list('cm_insar', rgbs)


#%%
def cm_isce():
    # https://github.com/insarlab/MintPy/blob/main/mintpy/objects/colors.py#L199
    # Almost same as GAMMA except for bottom RGB (0vs0.2), shift and flip
    rgbs = np.zeros((256,3), dtype=np.uint8)

    for kk in range(85):
        rgbs[kk,0] = kk*3
        rgbs[kk,1] = 255-kk*3
        rgbs[kk,2] = 255

    rgbs[85:170,0] = rgbs[0:85,2]
    rgbs[85:170,1] = rgbs[0:85,0]
    rgbs[85:170,2] = rgbs[0:85,1]

    rgbs[170:255,0] = rgbs[0:85,1]
    rgbs[170:255,1] = rgbs[0:85,2]
    rgbs[170:255,2] = rgbs[0:85,0]

    rgbs[255,0] = 0
    rgbs[255,1] = 255
    rgbs[255,2] = 255

    rgbs = np.roll(rgbs, int(256/2-214), axis=0)  #shift green to the center
    rgbs = np.flipud(rgbs)   #flip up-down so that orange is in the later half (positive)
    rgbs = rgbs/255.0

    return LSC.from_list('cm_isce', rgbs)


#%%
def cmap_insar():
    """
    How to use cmap_insar (GAMMA standard rainbow cmap):
        import matplotlib as mpl
        from matplotlib import pyplot as plt
        cdict = cmap_insar()
        plt.register_cmap(cmap=mpl.colors.LinearSegmentedColormap('insar', cdict))
        plt.register_cmap(name='insar', data=cdict)
        plt.imshow(array, cmap='insar', vmin=-np.pi, vmax=np.pi, interpolation='nearest')

    Note:
        - Input array should be wrapped and in radian
        - To wrap unwrapped phase, np.angle(np.exp(1j*unw/cycle)*cycle)
    """

#    These are for 0-2pi, not -pi to pi
    red = [255,255,255,255,240,203,165,127,90,55,92,130,167,205,243,255,255,255]
    green = [118,156,193,231,255,255,255,255,255,255,217,179,142,104,66,80,118,156]
    blue = [191,153,116,78,69,106,144,182,219,255,255,255,255,255,255,229,191,153]
    phase = [k/32 for k in range(1,33,2)]
    phase = [0] + phase + [1]

    red_norm = [ k/255 for k in red ] + [ red[0]/255 ]
    green_norm = [ k/255 for k in green ] + [ green[0]/255 ]
    blue_norm = [ k/255 for k in blue ] + [ blue[0]/255 ]

    redtuple=[]
    greentuple=[]
    bluetuple=[]
    for j in range(18):
        redtuple.append((phase[j],red_norm[j],red_norm[j+1]))
        greentuple.append((phase[j],green_norm[j],green_norm[j+1]))
        bluetuple.append((phase[j],blue_norm[j],blue_norm[j+1]))

    redtuple=tuple(redtuple)
    greentuple=tuple(greentuple)
    bluetuple=tuple(bluetuple)

    cdict = { 'red': redtuple, 'green': greentuple, 'blue': bluetuple }

    return cdict


#%%
def convert_size(size_bytes):
   if size_bytes == 0:
       return "0B"
   size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
   i = int(np.floor(np.log(size_bytes)/np.log(1024)))
   p = np.power(1024, i)
   s = round(size_bytes / p, 2)
   return "%s%s" % (s, size_name[i])


def extract_url_licsar(url):
    ''' new since Aug 2025+: URL gets different from direct to HTML file with a link..
    '''
    # first try if the direct url actually works (old version):
    response = requests.head(url, allow_redirects=True)
    if response.status_code == 200:
        return url
    else:
        fname = os.path.basename(url)
        url = os.path.dirname(url)
        # transition period - both old and new version exist - try the temporary url:
        url = url.replace('LiCSAR_products/', 'LiCSAR_products.future/')
        response = requests.get(url)
        if response.status_code != 200:
            # perhaps already after the transition? Getting it back:
            url = url.replace('LiCSAR_products.future/', 'LiCSAR_products/')
            response = requests.get(url)
            if response.status_code != 200:
                return None # this also does not exist
        response.encoding = response.apparent_encoding  # avoid garble
        html_doc = response.text
        soup = BeautifulSoup(html_doc, "html.parser")
        tag = soup.find(href=re.compile(fname))
        if tag:
            return tag.get('href')
        else:
            return None


#%%
def download_data(url, file, n_retry=3):
    for i in range(n_retry):
        try:
            start = time.time()
            with requests.get(url) as res:
                res.raise_for_status()
                size_remote = int(res.headers.get("Content-Length"))
                with open(file, "wb") as output:
                     output.write(res.content)
                size_local = os.path.getsize(file)
                if size_remote != size_local: ### Different size
                    raise IOError('Downloaded file not complete')

            elapsed = int(time.time()-start)
            if elapsed==0: elapsed=elapsed+1
            fsize = convert_size(os.path.getsize(file))
            speed = convert_size(int(os.path.getsize(file)/elapsed))
            print('    {}, {}, {}s, {}/s'.format(os.path.basename(file),
                   fsize, elapsed, speed), flush=True)
            return # success

        except Exception as e:
            print(    '    {} for {} in {}th try'.format(
                e.__class__.__name__, os.path.basename(url), i+1), flush=True)
            pass # try again

    print("    Error while downloading from {}".format(url),
          file=sys.stderr, flush=True)
    return # fail


#%%
def get_asf_reference_scene(lon, lat, imdates, center_time=None, path=None,
                            platform='SENTINEL-1', timeout=60):
    """
    Find a granule (scene) name at ASF usable as the reference of a baseline
    stack query, by searching scenes containing (lon, lat) on the epoch dates.
    No authentication is required.

    Inputs:
      lon, lat    : Point within the area of interest (deg, EPSG:4326)
      imdates     : List of epoch dates ['yyyymmdd',...]; tried in order
                    until a scene is found
      center_time : 'HH:MM:SS' approx acquisition time (UTC) to disambiguate
                    multiple paths covering the point (optional)
      path        : Path (track, relative orbit) number to disambiguate
                    (optional)
      platform    : ASF platform keyword (only SENTINEL-1 supported)

    Returns:
      Scene name (str), or None if not found or network error
      (a warning is printed).
    """
    url = 'https://api.daac.asf.alaska.edu/services/search/param'

    if center_time is not None:
        h, m, s = center_time.split(':')[0:3]
        center_sec = int(h)*3600 + int(m)*60 + int(float(s))

    for imd in imdates:
        params = {'platform': platform,
                  'processingLevel': 'SLC',
                  'beamMode': 'IW',
                  'intersectsWith': 'POINT({} {})'.format(lon, lat),
                  'start': '{}-{}-{}T00:00:00Z'.format(imd[0:4], imd[4:6], imd[6:8]),
                  'end': '{}-{}-{}T23:59:59Z'.format(imd[0:4], imd[4:6], imd[6:8]),
                  'output': 'geojson',
                  'maxResults': 20}
        try:
            with requests.get(url, params=params, timeout=timeout) as res:
                res.raise_for_status()
                features = res.json().get('features', [])
        except Exception as e:
            print('WARN: ASF scene search failed ({}: {})'.format(
                e.__class__.__name__, e), file=sys.stderr, flush=True)
            return None

        cands = [f['properties'] for f in features if
                 f['properties'].get('startTime', '')[0:10].replace('-', '') == imd]

        if path is not None:
            cands = [c for c in cands if int(c['pathNumber']) == int(path)]

        if center_time is not None:
            for c in cands:
                h, m, s = c['startTime'][11:19].split(':')
                scene_sec = int(h)*3600 + int(m)*60 + int(s)
                d = abs(scene_sec - center_sec)
                c['_dt_sec'] = min(d, 86400 - d)
            cands = sorted([c for c in cands if c['_dt_sec'] <= 600],
                           key=lambda c: c['_dt_sec'])

        if not cands:
            continue

        paths = sorted(set(int(c['pathNumber']) for c in cands))
        if len(paths) > 1:
            print('WARN: Multiple paths {} cover the point on {}; '
                  'use path {}'.format(paths, imd, cands[0]['pathNumber']),
                  flush=True)

        return cands[0]['sceneName']

    print('WARN: No SLC scene found at ASF for any epoch date',
          file=sys.stderr, flush=True)
    return None


#%%
def get_bperp_asf(lon, lat, imdates, center_time=None, path=None,
                  platform='SENTINEL-1', timeout=60, n_retry=3):
    """
    Get perpendicular baselines (bperp, m) from the ASF baseline API for the
    stack covering (lon, lat) on imdates. No authentication is required.

    Inputs: same as get_asf_reference_scene.

    Returns:
      bperp_dict : dict {'yyyymmdd': bperp (float)} of ALL dates available
                   in the ASF stack (can be more than imdates), relative to
                   the automatically selected reference scene, or
      None       : if the stack could not be retrieved (warning printed).
                   The caller should check that all needed imdates are in
                   the returned dict.
    """
    refscene = get_asf_reference_scene(lon, lat, imdates, center_time, path,
                                       platform, timeout)
    if refscene is None:
        return None

    print('Reference scene for ASF baseline stack: {}'.format(refscene),
          flush=True)

    url = 'https://api.daac.asf.alaska.edu/services/search/baseline'
    params = {'reference': refscene,
              'processingLevel': 'SLC',
              'output': 'geojson'}
    features = None
    for i in range(n_retry):
        try:
            with requests.get(url, params=params, timeout=timeout) as res:
                res.raise_for_status()
                features = res.json().get('features', [])
            break # success
        except Exception as e:
            print('    {} for baseline stack in {}th try'.format(
                e.__class__.__name__, i+1), flush=True)
            pass # try again

    if features is None:
        print('WARN: Error while getting baseline stack from ASF',
              file=sys.stderr, flush=True)
        return None

    if center_time is not None:
        h, m, s = center_time.split(':')[0:3]
        center_sec = int(h)*3600 + int(m)*60 + int(float(s))

    bperp_dict = {}
    dt_sec_dict = {} ## to keep the scene closest in time per date
    for f in features:
        props = f['properties']
        if props.get('perpendicularBaseline') is None:
            continue ## missing state vectors
        imd = props['startTime'][0:10].replace('-', '')

        if center_time is not None:
            h, m, s = props['startTime'][11:19].split(':')
            scene_sec = int(h)*3600 + int(m)*60 + int(s)
            d = abs(scene_sec - center_sec)
            dt_sec = min(d, 86400 - d)
        else:
            dt_sec = 0

        ## A date can appear more than once (neighboring frames); dedup
        if imd not in bperp_dict or dt_sec < dt_sec_dict[imd]:
            bperp_dict[imd] = float(props['perpendicularBaseline'])
            dt_sec_dict[imd] = dt_sec

    return bperp_dict


#%%
def fit2d(A,w=None,deg="1"):
    """
    Estimate best fit plain with indicated degree of polynomial.

    Inputs:
        A : Input ndarray (can include nan)
        w : Wieight (1/std**2) for each element of A (with the same dimention as A)
        deg : degree of polynomial of fitting plain
         - 1   -> a+bx+cy (ramp)
         - bl -> a+bx+cy+dxy (biliner)
         - 2   -> a+bx+cy+dxy+ex**2_fy**2 (2d polynomial)

    Returns:
        Afit : Best fit plain with the same demention as A
        m    : set of parameters of best fit plain (a,b,c...)

    """
    import statsmodels.api as sm

    ### Make design matrix G
    length,width = A.shape #read dimension
    Xgrid,Ygrid = np.meshgrid(np.arange(width),np.arange(length)) #mesh grid

    if str(deg) == "1":
        G = np.stack((np.ones((length*width)), Xgrid.flatten(), Ygrid.flatten())).T
    elif str(deg) == "bl":
        G = np.stack((np.ones((length*width)), Xgrid.flatten(), Ygrid.flatten(), Xgrid.flatten()*Ygrid.flatten())).T
    elif str(deg) == "2":
        G = np.stack((np.ones((length*width)), Xgrid.flatten(), Ygrid.flatten(), Xgrid.flatten()*Ygrid.flatten(), Xgrid.flatten()**2, Ygrid.flatten()**2)).T
    else:
        print('\nERROR: Not proper deg ({}) is used\n'.format(deg), file=sys.stderr)
        return False

    ### Handle nan by 0 padding and 0 weight
    # Not drop in sm because cannot return Afit
    if np.any(np.isnan(A)):
        bool_nan = np.isnan(A)
        A = A.copy() # to avoid change original value in main
        A[bool_nan] = 0
        if w is None:
            w = np.ones_like(A)
        w = w.copy() # to avoid change original value in main
        w[bool_nan] = 0

    ### Invert
    if w is None: ## Ordinary LS
        results = sm.OLS(A.ravel(), G).fit()
    else: ## Weighted LS
        results = sm.WLS(A.ravel(), G, weights=w.ravel()).fit()

    m = results.params
    Afit = np.float32(results.predict().reshape((length,width)))

    return Afit,m


#%%
def fit2dh(A, deg, hgt, hgt_min, hgt_max, gpu=False):
    """
    Estimate best fit 2d ramp and topography-correlated component simultaneously.

    Inputs:
        A   : Input ndarray (can include nan)
        deg : degree of polynomial for fitting ramp
          - 1  -> a+bx+cy (ramp, default)
          - bl -> a+bx+cy+dxy (biliner)
          - 2  -> a+bx+cy+dxy+ex**2_fy**2 (2d polynomial)
          - []  -> a (must be used with hgt)
        hgt : Input hgt to estimate coefficient of topo-corr component
              If blank, don*t estimate topo-corr component.
        hgt_min : Minimum hgt to take into account in hgt-linear
        hgt_max : Maximum hgt to take into account in hgt-linear
        gpu     : GPU flag

    Returns:
        Afit : Best fit solution with the same demention as A
        m    : Set of parameters of best fit plain (a,b,c...)

    Note: GPU option seems slow and may return error when the size is large.
          Not recommended.

    """
    if gpu:
        import cupy as xp
        A = xp.asarray(A)
        hgt = xp.asarray(hgt)
        hgt_min = xp.asarray(hgt_min)
        hgt_max = xp.asarray(hgt_max)
    else:
        xp = np

    ### Make design matrix G
    length, width = A.shape

    if not deg:
        G = xp.ones((length*width))
    else:
        Xgrid, Ygrid = xp.meshgrid(xp.arange(width), xp.arange(length))
        Xgrid1 = Xgrid.ravel()
        Ygrid1 = Ygrid.ravel()

        if str(deg) == "1":
            G = xp.stack((Xgrid1, Ygrid1)).T
        elif str(deg) == "bl":
            G = xp.stack((Xgrid1, Ygrid1, Xgrid1*Ygrid1)).T
        elif str(deg) == "2":
            G = xp.stack((Xgrid1, Ygrid1, Xgrid1*Ygrid1,
                          Xgrid1**2, Ygrid1**2)).T
        else:
            print('\nERROR: Not proper deg ({}) is used\n'.format(deg), file=sys.stderr)
            return False
        del Xgrid, Ygrid, Xgrid1, Ygrid1

        G = xp.hstack([xp.ones((length*width, 1)), G])

    if len(hgt) > 0:
        _hgt = hgt.copy()  ## Not to overwrite hgt in main
        _hgt[xp.isnan(_hgt)] = 0
        _hgt[_hgt<hgt_min] = 0
        _hgt[_hgt>hgt_max] = 0
        G2 = xp.vstack((G.T, hgt.ravel())).T ## for Afit
        G = xp.vstack((G.T, _hgt.ravel())).T
        del _hgt
    else:
        G2 = G

    G = G.astype(xp.int32)

    ### Invert
    mask = xp.isnan(A.ravel())
    m = xp.linalg.lstsq(G[~mask, :], A.ravel()[~mask], rcond=None)[0]

    Afit = ((xp.matmul(G2, m)).reshape((length, width))).astype(xp.float32)

    if gpu:
        Afit = xp.asnumpy(Afit)
        m = xp.asnumpy(m)
        del A, hgt, hgt_min, hgt_max, length, width, G, G2, mask

    return Afit, m


#%%
def get_cmap(cmap_name, cmapN=256):
    """
    Return cmap (LinearSegmentedColormap).
    cmap_name can be:
        - Matplotlib predefined name (e.g. viridis)
        - Scientific colour maps (e.g. SCM.roma)
        - Generic Mapping Tools (e.g. GMT.polar)
        - cmocean (e.g. cmocean.phase)
        - colorcet (e.g. colorcet.CET_C1)
        - cm_insar, cm_isce
    All cmap can be reversed with "_r"
    """
    flag = 0 # mpl
    if cmap_name.startswith('cm_insar'):
        _cmap = cm_insar()
        flag = 1
    elif cmap_name.startswith('cm_isce'):
        _cmap = cm_isce()
        flag = 1
    elif cmap_name.startswith('SCM'):
        import SCM as CMAP
        flag = 2
    elif cmap_name.startswith('GMT'):
        import GMT as CMAP
        flag = 2
    elif cmap_name.startswith('cmocean'):
        import cmocean as CMAP
        flag = 2
    elif cmap_name.startswith('colorcet'):
        import colorcet as CMAP
        flag = 3

    if flag >= 2:
        cmdir = os.path.dirname(CMAP.__file__)
        name = cmap_name.split('.')[1]
        if cmap_name.endswith('_r'):
            name = name[:-2]
        if flag == 2:
            file = os.path.join(cmdir, name, name+'.txt')
            cm_data = np.loadtxt(file)
        elif flag == 3:
            file = os.path.join(cmdir, name+'.csv')
            cm_data = np.loadtxt(file, delimiter=',')
        _cmap = LSC.from_list(name, cm_data)

    if flag >= 1:
        if cmap_name.endswith('_r'):
            _cmap = _cmap.reversed()
        if cmap_name not in plt.colormaps():
            mpl.colormaps.register(_cmap, name=cmap_name)

    cmap = plt.get_cmap(cmap_name, cmapN)

    return cmap


#%%
def get_ifgdates(ifgdir):
    """
    Get ifgdates and imdates in ifgdir.

    Returns:
        ifgdates : List of dates of ifgs
        imdates  : List of dates of images
    """
    ifgdates = [str(k) for k in sorted(os.listdir(ifgdir))
                if len(k) == 17
                and k[0] =='2'
                and k[8] =='_'
                and k[9] =='2'
                and os.path.isdir(os.path.join(ifgdir, k))]

    return ifgdates


#%%
def get_patchrow(width, length, n_data, memory_size):
    """
    Get patch number of rows for memory size (in MB).

    Returns:
        n_patch : Number of patches (int)
        patchrow : List of the number of rows for each patch.
                ex) [[0, 1234], [1235, 2469],... ]
    """
    data_size = width*length*n_data*4/2**20  #in MiB, 4byte float
    n_patch = int(np.ceil(data_size/memory_size))
    ### Devide only row direction

    patchrow= []
    for i in range(n_patch):
        rowspacing = int(np.ceil(length/n_patch))
        patchrow.append([i*rowspacing,(i+1)*rowspacing])
        if i == n_patch-1:
            patchrow[-1][-1] = length

    return n_patch, patchrow


#%%
def run_pool(func, args, n_para, mem_per_worker_mb=None, chunksize=None,
             out=None, store=None):
    """
    Robust replacement of multiprocessing fork Pool(n_para).map(func, args).

    - Forks workers so that func can use global variables set in the caller
      (same semantics as the conventional fork Pool in LiCSBAS).
    - If mem_per_worker_mb is given, n_para is capped just before execution
      based on the currently available memory to avoid OOM.
    - If a worker dies unexpectedly (e.g., killed by the OS out-of-memory
      killer), raise RuntimeError with advice instead of hanging forever.
    - If out is given, results are stored one by one as out[i] = result
      (streaming; avoids materializing all results at a time).
    - If store is given, store(i, result) is called instead (e.g., for
      functions returning tuples). store overrides out.
    - If n_para == 1, run serially without forking (also should be used
      after GPU use because forking after CUDA initialization is unsafe).

    Returns:
        list of results in the order of args (or out if out is given)
    """
    from concurrent.futures import ProcessPoolExecutor
    from concurrent.futures.process import BrokenProcessPool
    import multiprocessing as multi
    import psutil

    args = list(args)
    n_task = len(args)
    if n_task == 0:
        return [] if (out is None and store is None) else out
    n_para = max(1, min(int(n_para), n_task))

    ### Cap n_para by memory available now
    if mem_per_worker_mb:
        mem_avail = psutil.virtual_memory().available / 2**20  # MB
        n_para_mem = max(1, int(mem_avail / 2 / mem_per_worker_mb))
        if n_para_mem < n_para:
            print('  Reduce n_para from {} to {} due to available memory '
                  '({:.0f} MB, ~{:.0f} MB/worker)'.format(
                      n_para, n_para_mem, mem_avail, mem_per_worker_mb),
                  flush=True)
            n_para = n_para_mem

    if store is None and out is not None:
        def store(i, r): out[i] = r

    def _consume(iterator):
        if store is None:
            return list(iterator)
        for i, r in enumerate(iterator):
            store(i, r)
        return out

    if n_para == 1:  ## Serial, no fork
        return _consume(map(func, args))

    if chunksize is None:  ## Mimic default of Pool.map
        chunksize = -(-n_task // (n_para*4))

    ctx = multi.get_context('fork')
    try:
        with ProcessPoolExecutor(max_workers=n_para, mp_context=ctx) as p:
            return _consume(p.map(func, args, chunksize=chunksize))
    except BrokenProcessPool as e:
        raise RuntimeError(
            'A worker process died unexpectedly, most likely killed by the '
            'OS out-of-memory (OOM) killer because memory ran out\n'
            '(check with e.g. `dmesg -T | grep -i "killed process"`).\n'
            'Try one of the following and rerun:\n'
            '  - Reduce the number of parallel workers (--n_para)\n'
            '  - Reduce memory use (e.g., --mem_size, smaller data size)\n'
            '  - Close other programs or use a machine with more memory'
        ) from e


#%%
def ifgdates2imdates(ifgdates):
    primarylist = []
    secondarylist = []
    for ifgd in ifgdates:
        primarylist.append(ifgd[:8])
        secondarylist.append(ifgd[-8:])

    imdates = list(set(primarylist+secondarylist)) # set is a unique operator
    imdates.sort()

    return imdates


#%%
def multilook(array, nlook_r, nlook_c, n_valid_thre=0.5):
    """
    Nodata in input array must be filled with nan beforehand.
    if the number of valid data is less than n_valid_thre*nlook_r*nlook_c, return nan.
    """
    length, width = array.shape
    length_ml = int(np.floor(length/nlook_r))
    width_ml = int(np.floor(width/nlook_c))

    array_reshape = array[:length_ml*nlook_r,:width_ml*nlook_c].reshape(length_ml, nlook_r, width_ml, nlook_c)

    with warnings.catch_warnings(): ## To silence RuntimeWarning: Mean of empty slice
        warnings.simplefilter('ignore', RuntimeWarning)
        array_ml = np.nanmean(array_reshape, axis=(1, 3))

    n_valid = np.sum(~np.isnan(array_reshape), axis=(1, 3))
    bool_invalid = n_valid < n_valid_thre*nlook_r*nlook_c
    array_ml[bool_invalid] = np.nan

    return array_ml


#%%
def read_point(point_str, width, length):
    if re.match('[0-9]*/[0-9]*', point_str):
        x, y = [int(s) for s in re.split('[/]', point_str)]
        if x > width-1 or y > length-1:
            print("\nERROR:", file=sys.stderr)
            print("Index exceed input dimension ({0},{1})!".format(width,length), file=sys.stderr)
            return False
    else:
        print("\nERROR:", file=sys.stderr)
        print("Point format seems to be wrong (should be x/y)", file=sys.stderr)
        return False

    return [x, y]


#%%
def read_range(range_str, width, length):
    if re.match('[0-9]*:[0-9]*/[0-9]*:[0-9]', range_str):
        x1, x2, y1, y2 = [int(s) for s in re.split('[:/]', range_str)]
        if x2 == 0:
            x2 = width
        if y2 == 0:
            y2 = length
        if x1 > width-1 or x2 > width or y1 > length-1 or y2 > length:
            print("\nERROR:", file=sys.stderr)
            print("Index exceed input dimension ({0},{1})!".format(width,length), file=sys.stderr)
            return False
        if x1 >= x2 or y1 >= y2:
            print("\nERROR: x2/y2 must be larger than x1/y1", file=sys.stderr)
            return False
    else:
        print("\nERROR:", file=sys.stderr)
        print("Range format seems to be wrong (should be x1:x2/y1:y2)", file=sys.stderr)
        return False

    return [x1, x2, y1, y2]


#%%
def read_range_line(range_str, width, length):
    if re.match('[0-9]*,[0-9]*/[0-9]*,[0-9]', range_str):
        x1, y1, x2, y2 = [int(s) for s in re.split('[,/]', range_str)]
        if x1 > width-1 or x2 > width-1 or y1 > length-1 or y2 > length-1:
            print("\nERROR:", file=sys.stderr)
            print("Index exceed input dimension ({0},{1})!".format(width ,length), file=sys.stderr)
            return False
    else:
        print("\nERROR:", file=sys.stderr)
        print("Range format seems to be wrong (must be x1,y1/x2,y2)", file=sys.stderr)
        return False

    return [x1, x2, y1, y2]



#%%
def read_range_geo(range_str, width, length, lat1, postlat, lon1, postlon):
    """
    lat1 is north edge and postlat is negative value.
    lat lon values are in grid registration
    Note: x1/y1 range 0 to width-1, while x2/y2 range 1 to width
    """
    lat2 = lat1+postlat*(length-1)
    lon2 = lon1+postlon*(width-1)

    if re.match('[+-]?\d+(?:\.\d+)?/[+-]?\d+(?:\.\d+)?/[+-]?\d+(?:\.\d+)?/[+-]?\d+(?:\.\d+)?', range_str):
        lon_w, lon_e, lat_s, lat_n = [float(s) for s in range_str.split('/')]
        x1 = int(np.round((lon_w - lon1)/postlon)) if lon_w > lon1 else 0
        x2 = int(np.round((lon_e - lon1)/postlon))+1 if lon_e < lon2 else width
        y1 = int(np.round((lat_n - lat1)/postlat)) if lat_n < lat1 else 0
        y2 = int(np.round((lat_s - lat1)/postlat))+1 if lat_s > lat2 else length
    else:
        print("\nERROR:", file=sys.stderr)
        print("Range format seems to be wrong (should be lon1/lon2/lat1/lat2)", file=sys.stderr)
        return False

    return [x1, x2, y1, y2]


#%%
def read_range_line_geo(range_str, width, length, lat_n, postlat, lon_w, postlon):
    """
    lat_n is north edge and postlat is negative value.
    lat lon values are in grid registration
    """

    if re.match('[+-]?\d+(?:\.\d+)?,[+-]?\d+(?:\.\d+)?/[+-]?\d+(?:\.\d+)?,[+-]?\d+(?:\.\d+)?', range_str):
        lon1, lat1, lon2, lat2 = [float(s) for s in re.split('[,/]', range_str)]
        x1 = int(np.round((lon1 - lon_w)/postlon))
        x2 = int(np.round((lon2 - lon_w)/postlon))
        y1 = int(np.round((lat1 - lat_n)/postlat))
        y2 = int(np.round((lat2 - lat_n)/postlat))
    else:
        print("\nERROR:", file=sys.stderr)
        print("Range format seems to be wrong (should be lon1,lat1/lon2,lat2)", file=sys.stderr)
        return False

    return [x1, x2, y1, y2]


#%%
def xy2bl(x, y, lat1, dlat, lon1, dlon):
    """
    xy index starts from 0, end with width/length-1
    """
    lat = lat1+dlat*y
    lon = lon1+dlon*x

    return lat, lon


#%%
def get_pair_folders(geoc_dir):
    """Get list of pair folders (yyyymmdd_yyyymmdd) in GEOC dir."""
    pairs = []
    for item in os.listdir(geoc_dir):
        path = os.path.join(geoc_dir, item)
        if os.path.isdir(path) and re.match(r'\d{8}_\d{8}', item):
            pairs.append(item)
    return sorted(pairs)


#%%
def calculate_common_geometry(gt1, shape1, gt2, shape2):
    """Calculate common geometry by resampling to gt1 resolution and overlapping extent."""
    # gt1 is reference: (x0, dx, 0, y0, 0, dy)
    x0_1, dx1, _, y0_1, _, dy1 = gt1
    x0_2, dx2, _, y0_2, _, dy2 = gt2

    # Calculate overlapping extent
    x1_end = x0_1 + shape1[1] * dx1
    x2_end = x0_2 + shape2[1] * dx2
    y1_end = y0_1 + shape1[0] * dy1
    y2_end = y0_2 + shape2[0] * dy2

    x_min = max(x0_1, x0_2)
    x_max = min(x1_end, x2_end)
    y_top = min(y0_1, y0_2)
    y_bottom = max(y1_end, y2_end)

    if x_min >= x_max or y_bottom >= y_top:
        raise ValueError("No overlapping area")

    # Align overlap to the reference grid
    eps = 1e-9
    x_offset = int(math.ceil((x_min - x0_1) / dx1 - eps))
    x_end_index = int(math.floor((x_max - x0_1) / dx1 + eps))
    y_offset = int(math.ceil((y0_1 - y_top) / abs(dy1) - eps))
    y_end_index = int(math.floor((y0_1 - y_bottom) / abs(dy1) + eps))

    nx = x_end_index - x_offset
    ny = y_end_index - y_offset
    if nx <= 0 or ny <= 0:
        raise ValueError("No overlapping area after alignment to reference grid")

    new_gt = (x0_1 + x_offset * dx1, dx1, 0, y0_1 + y_offset * dy1, 0, dy1)
    return new_gt, (ny, nx)