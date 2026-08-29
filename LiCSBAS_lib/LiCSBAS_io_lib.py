#!/usr/bin/env python3
"""
Python3 library of input/output functions for LiCSBAS.

v1.6 20260829 Yu Morishita

"""
import sys
import numpy as np
import subprocess as subp
import datetime as dt
import statsmodels.api as sm
from osgeo import gdal, osr

### Raise exceptions on GDAL errors instead of returning None.
### Required from GDAL 3.7 (FutureWarning) and the default in GDAL 4.0.
gdal.UseExceptions()


#%%
def make_dummy_bperp(bperp_file, imdates):
    with open(bperp_file, 'w') as f:
        for i, imd in enumerate(imdates):
            if i==0: bp = 0
            elif np.mod(i, 4)==1: bp = np.random.rand()/2+0.5 #0.5~1
            elif np.mod(i, 4)==2: bp = -np.random.rand()/2-0.5 #-1~-0.5
            elif np.mod(i, 4)==3: bp = np.random.rand()/2 #0~0.5
            elif np.mod(i, 4)==0: bp = -np.random.rand()/2 #-0.5~0

            ifg_dt = dt.datetime.strptime(imd, '%Y%m%d').toordinal() - dt.datetime.strptime(imdates[0], '%Y%m%d').toordinal()

            print('{:3d} {} {} {:5.2f} {:4d} {} {:4d} {} {:5.2f}'.format(i, imdates[0], imd, bp, ifg_dt, 0, ifg_dt, 0, bp), file=f)


#%%
def make_bperp_file(bperp_file, imdates, bperp_dict):
    """
    Make a baselines file in the new 4-column format readable by
    read_bperp_file:
          smdate    sdate    bp    dt
        20170302 20170326 130.9  24.0
    smdate is imdates[0]; bp is bperp (m) of sdate relative to smdate;
    dt is temporal baseline (days).

    bperp_dict is {'yyyymmdd': bperp} with values relative to any common
    reference (e.g. from tools_lib.get_bperp_asf); they are re-referenced
    to imdates[0] here. Raise ValueError if any imdate is missing from
    bperp_dict.
    """
    missing = [imd for imd in imdates if imd not in bperp_dict]
    if missing:
        raise ValueError('bperp not available for: {}'.format(' '.join(missing)))

    smdate = imdates[0]
    bp0 = float(bperp_dict[smdate])
    with open(bperp_file, 'w') as f:
        for imd in imdates[1:]:
            bp = float(bperp_dict[imd]) - bp0
            ifg_dt = dt.datetime.strptime(imd, '%Y%m%d').toordinal() - dt.datetime.strptime(smdate, '%Y%m%d').toordinal()
            print('{} {} {:6.1f} {:5.1f}'.format(smdate, imd, bp, float(ifg_dt)), file=f)


#%%
def make_geotiff(data, latn_p, lonw_p, dlat, dlon, outfile, compress_option, nodata=None):
    length, width = data.shape
    if data.dtype == np.float32:
        dtype = gdal.GDT_Float32
    elif data.dtype == np.uint8:
        dtype = gdal.GDT_Byte
    else:
        raise ValueError('dtype {} is not supported (must be float32 or uint8)'.format(data.dtype))

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(outfile, width, length, 1, dtype, options=compress_option)
    outRaster.SetGeoTransform((lonw_p, dlon, 0, latn_p, 0, dlat))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(data)
    if nodata is not None: outband.SetNoDataValue(nodata)
    outRaster.SetMetadataItem('AREA_OR_POINT', 'Area')
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


#%%
def make_point_kml(lat, lon, kmlfile):
    with open(kmlfile, "w") as f:
        print('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document><Placemark><Point>\n<coordinates>{},{}</coordinates>\n</Point></Placemark></Document>\n</kml>'.format(lon, lat), file=f)


#%%
def make_tstxt(x, y, imdates, ts, tsfile, refx1, refx2, refy1, refy2, gap, lat=None, lon=None, reflat1=None, reflat2=None, reflon1=None, reflon2=None, deramp_flag=None, hgt_linear_flag=None, filtwidth_km=None, filtwidth_yr=None):
    """
    Make txt of time series.
    Format example:
    # x, y    : 432, 532
    # lat, lon: 34.65466, 136.65432
    # ref     : 21:22/54:55
    # refgeo  : 136.98767/136.98767/34.95364/34.95364
    # deramp, filtwidth_km, filtwidth_yr: 1, 2, 0.653
    # hgt_linear_flag: 1
    # gap     : 20160104_20160116, 20170204_20170216
    # linear model: -3.643*t+4.254
    20141030    0.00
    20150216   -3.50
    20160716   -3.5
    """
    ### Calc model
    imdates_ordinal = np.array(([dt.datetime.strptime(imd, '%Y%m%d').toordinal() for imd in imdates])) ##73????
    imdates_yr = (imdates_ordinal-imdates_ordinal[0])/365.25
    A = sm.add_constant(imdates_yr) #[1, t]
    vconst, vel = sm.OLS(ts, A, missing='drop').fit().params

    ### Identify gaps
    ixs_gap = np.where(gap==1)[0] # n_im-1, bool
    gap_str = ''
    for ix_gap in ixs_gap:
        gap_str = gap_str+imdates[ix_gap]+'_'+imdates[ix_gap+1]+' '

    ### Output
    with open(tsfile, 'w') as f:
        print('# x, y    : {}, {}'.format(x, y), file=f)
        if all(v is not None for v in [lat, lon]):
            print('# lat, lon: {:.5f}, {:.5f}'.format(lat, lon), file=f)
        print('# ref     : {}:{}/{}:{}'.format(refx1, refx2, refy1, refy2), file=f)
        if all(v is not None for v in [reflon1, reflon2, reflat1, reflat2]):
            print('# refgeo  : {:.5f}/{:.5f}/{:.5f}/{:.5f}'.format(reflon1, reflon2, reflat1, reflat2), file=f)
        if filtwidth_yr is not None:
            print('# deramp, filtwidth_km, filtwidth_yr : {}, {}, {:.3f}'.format(deramp_flag, filtwidth_km, filtwidth_yr), file=f)
        if hgt_linear_flag is not None:
            print('# hgt_linear_flag : {}'.format(hgt_linear_flag), file=f)
        print('# gap     : {}'.format(gap_str), file=f)
        print('# linear model: {:.3f}*t{:+.3f}'.format(vel, vconst), file=f)

        for i, imd in enumerate(imdates):
            print('{} {:7.2f}'.format(imd, ts[i]), file=f)


#%%
def read_bperp_file(bperp_file, imdates):
    """
    bperp_file (baselines) contains (m: primary (master), s: secondary,
                                     sm: single prime):
          smdate    sdate    bp    dt
        20170302 20170326 130.9  24.0
        20170302 20170314  32.4  12.0

    Old bperp_file contains (m: primary (master), s:secondary,
                             sm: single prime):
        num    mdate    sdate   bp   dt  dt_m_sm dt_s_sm bp_m_sm bp_s_sm
          1 20170218 20170326 96.6 36.0    -12.0    24.0    34.2   130.9
          2 20170302 20170314 32.4 12.0      0.0    12.0     0.0    32.4

    Return: bperp
    """
    bperp = []
    bperp_dict = {}

    ### Determine type of bperp_file; old or not
    with open(bperp_file) as f:
        line = f.readline().split() #list

    if len(line) == 4: ## new format
        bperp_dict[line[0]] = '0.00' ## single prime. unnecessary?
        with open(bperp_file) as f:
            for l in f:
                if len(l.split()) == 4:
                    bperp_dict[l.split()[1]] = l.split()[2]

    else: ## old format
        with open(bperp_file) as f:
            for l in f:
                bperp_dict[l.split()[1]] = l.split()[-2]
                bperp_dict[l.split()[2]] = l.split()[-1]

    for imd in imdates:
        if imd in bperp_dict:
            bperp.append(float(bperp_dict[imd]))
        else: ## If no key exists
            print('ERROR: bperp for {} not found!'.format(imd), file=sys.stderr)
            return False

    return bperp


#%%
def read_geotiff(file, file_ref=None):
    geotiff = gdal.Open(file)

    if file_ref is not None: # Compare size and area
        size = (geotiff.RasterXSize, geotiff.RasterYSize)
        area = geotiff.GetGeoTransform()

        geotiff_ref = gdal.Open(file_ref)
        size_ref = (geotiff_ref.RasterXSize, geotiff_ref.RasterYSize)
        area_ref = geotiff_ref.GetGeoTransform()

        if not (size == size_ref and area == area_ref):
            raise Exception('ERROR: File size or area are not identical between {} and {}'.format(file, file_ref))

    data = geotiff.ReadAsArray()

    return data


#%%
def read_img(file, length, width, dtype=np.float32, endian='little'):
    """
    Read image data into numpy array.
    endian: 'little' or 'big' (not 'little' is regarded as 'big')
    """

    if endian == 'little':
        data = np.fromfile(file, dtype=dtype).reshape((length, width))
    else:
        data = np.fromfile(file, dtype=dtype).byteswap().reshape((length, width))

    return data


#%%
def read_ifg_list(ifg_listfile):
    ifgdates = []
    f = open(ifg_listfile)
    line = f.readline()
    while line:
        ifgd = line.split()[0]
        if ifgd == "#":
            line = f.readline()
            continue # Comment
        else:
            ifgdates.append(ifgd)
            line = f.readline()
    f.close()

    return ifgdates


#%%
def get_param_par(mlipar, field):
    """
    Get parameter from mli.par or dem_par file. Examples of fields are;
     - range_samples
     - azimuth_lines
     - range_looks
     - azimuth_looks
     - range_pixel_spacing (m)
     - azimuth_pixel_spacing (m)
     - radar_frequency  (Hz)
    """
    value = subp.check_output(['grep', field,mlipar]).decode().split()[1].strip()
    return value


#%%
def get_geotiff_info(tif_path):
    """Get geotransform, projection, and shape from a geotiff."""
    try:
        ds = gdal.Open(tif_path)
    except RuntimeError as e:  ## gdal.UseExceptions() is on
        raise ValueError(f"Cannot open {tif_path}") from e
    gt = ds.GetGeoTransform()
    proj = ds.GetProjection()
    shape = (ds.RasterYSize, ds.RasterXSize)
    ds = None
    return gt, proj, shape


#%%
def resample_geotiff(src_path, dst_path, new_gt, new_shape, proj, dtype=None, nodata=None, resample_alg='cubic'):
    """Resample geotiff to new geometry.

    Parameters
    ----------
    resample_alg : str or int
        Resampling algorithm to use. If a string is provided, it is mapped to a
        GDAL resampling constant. The default is 'cubic'.
    """
    try:
        src_ds = gdal.Open(src_path)
    except RuntimeError:  ## gdal.UseExceptions() is on
        return False

    if dtype is None:
        dtype = src_ds.GetRasterBand(1).DataType
    if nodata is None:
        nodata = src_ds.GetRasterBand(1).GetNoDataValue()
        if nodata is None:
            nodata = 0.0
            src_ds.GetRasterBand(1).SetNoDataValue(nodata)

    if isinstance(resample_alg, str):
        alg = resample_alg.strip().lower().replace('-', '').replace('_', '')
        resample_map = {
            'nearest': gdal.GRA_NearestNeighbour,
            'nearestneighbour': gdal.GRA_NearestNeighbour,
            'bilinear': gdal.GRA_Bilinear,
            'cubic': gdal.GRA_Cubic,
            'cubicspline': gdal.GRA_CubicSpline,
            'lanczos': gdal.GRA_Lanczos,
            'average': gdal.GRA_Average,
            'mode': gdal.GRA_Mode,
            'max': gdal.GRA_Max,
            'min': gdal.GRA_Min,
            'med': gdal.GRA_Med,
            'q1': gdal.GRA_Q1,
            'q3': gdal.GRA_Q3,
        }
        resample_alg = resample_map.get(alg, gdal.GRA_NearestNeighbour)

    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(dst_path, new_shape[1], new_shape[0], 1, dtype, options=['COMPRESS=DEFLATE', 'TILED=YES'])
    dst_ds.SetGeoTransform(new_gt)
    dst_ds.SetProjection(proj)
    dst_ds.GetRasterBand(1).SetNoDataValue(nodata)

    gdal.ReprojectImage(src_ds, dst_ds, src_ds.GetProjection(), proj, resample_alg)
    dst_ds.FlushCache()
    dst_ds = None
    src_ds = None
    return True

