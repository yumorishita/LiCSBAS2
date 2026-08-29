#!/usr/bin/env python3
"""
v1.2.0 20260829 Yu Morishita

This script checks if LiCSBAS install is OK or not.

"""

#%% Import
from importlib import import_module
import platform
import shutil
import sys


#%% Main
if __name__ == "__main__":
    flag = True
    ## Keep in sync with LiCSBAS.yml and requirements.txt
    modules = ['astropy',
                'bs4',
                'dateutil',
                'h5py',
                'matplotlib',
                'numpy',
                'pandas',
                'psutil',
                'requests',
                'shapely',
                'statsmodels',
                'threadpoolctl',
                'yaml',
               ]


    py_min = (3, 12)
    print('\nPython version: {}'.format(platform.python_version()))
    pyver = tuple(int(v) for v in platform.python_version_tuple()[:2])
    if pyver < py_min:
        print('  ERROR: must be >= {}.{}'.format(*py_min))
        flag = False
    else:
        print('  OK')


    print('\nCheck required modues and versions')
    for module in modules:
        try:
            imported = import_module(module)
        except Exception as err:
            print('  ERROR: {}'.format(err))
            flag = False
        else:
            ver = getattr(imported, '__version__', 'unknown')
            print('  {}({}) OK'.format(module, ver))


    try:
        imported = import_module('osgeo.gdal')
    except Exception as err:
        print('  ERROR: {}'.format(err))
        flag = False
    else:
        gdal_min = (3, 6)
        gdalver = tuple(int(v) for v in imported.__version__.split('.')[:2])
        ver = imported.__version__
        if gdalver < gdal_min:
            print('  ERROR: gdal ver is {} but must be >= {}.{}'.format(
                ver, *gdal_min))
            flag = False
        else:
            print('  gdal({}) OK'.format(ver))


    print('\nCheck LiCSBAS commands')
    rc = shutil.which('LiCSBAS01_get_geotiff.py')
    if rc is None:
        print('  ERROR: PATH is not set to LiCSBAS commands')
        flag = False
    else:
        print('  OK')


    print('\nCheck LiCSBAS library')
    try:
         imported = import_module('LiCSBAS_io_lib')
    except Exception:
         print('  ERROR: PYTHONPATH is not set to LiCSBAS library')
         flag = False
    else:
         print('  OK')


    if flag:
        print('\nLiCSBAS install is OK\n')
        sys.exit(0)
    else:
        print('\nERROR: LiCSBAS install is NOT OK\n')
        sys.exit(1)

