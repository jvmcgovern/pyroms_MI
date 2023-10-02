#
#
# #####################################################################
# #####################################################################
#
#  Main program
#
#  Build a CROCO initial file using GLORYS12 renanalysis data
#
#
#
# #####################################################################
# #####################################################################
#
#
import glob
import cftime
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from scipy import interpolate
from netCDF4 import date2index as d2i
from netCDF4 import num2date as n2d
from datetime import date, datetime
from calendar import monthrange
import sys
import pandas

# sys.path.insert(0,'')

from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from progressbar import *

if 1 == 1:
    # def main_func():

    #
    # #################### USERS DEFINED VARIABLES ########################
    #

    title = 'Initial file using GLORYS'

    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']

    N = 20
    theta_s = 7.
    theta_b = 0.
    hc = 20.
    vtransform = 2

    Yorig = 1990  # year origin of time : days since Yorig-01-01

    Yini = 2017
    Mini = 1
    Dini = 1

    Ystart = 2013
    Mstart = 3
    Dstart = 1

    Yend = 2013
    Mend = 8
    Dend = 31

    # Get list of files of interest between start and end month for year of interest
    # open each file
    # extract surface temperature, along with the time field
    # initialise the netcdf
    # pass the new array for surface temp and the time field

    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m

    ncg = netcdf(grdname, 'r')
    ncgo = ncg.variables
    lat_rho = np.array(ncgo['lat_rho'][:])
    lon_rho = np.array(ncgo['lon_rho'][:])
    ncg.close()

    crocodir = '/media/dskfour/CROCO_BGC_1p2p1_1314_1017866/'
    # crocodir = '/media/dskfour/CROCO_BGC_1p2p1_1314_1035076/'

    crocofil = 'croco_avg_'
    croco_end = '.nc'
    flist = []
    for mth in range(Mstart, Mend + 1):
        flist = flist + sorted(glob.glob(crocodir + crocofil + 'Y' + str(Ystart) + 'M' + str(mth).zfill(2) + croco_end))

    Tsurfname = crocodir + \
                'croco_Tsurf' + '_Y' + str(Ystart) + '_M' + str(Mstart).zfill(2) + '_M' + str(Mend).zfill(2) + '.nc'

    for fl in range(0, len(flist)):
        ncf = netcdf(flist[fl], 'r')
        ncfo = ncf.variables
        if fl == 0:
            time = np.array(ncfo['time'][:])
            tsurf = np.array(ncfo['temp'][:, ncfo['temp'].shape[1] - 1, :, :])
        else:
            time = np.concatenate((time, np.array(ncfo['time'][:])))
            tsurf = np.concatenate((tsurf, np.array(ncfo['temp'][:, ncfo['temp'].shape[1] - 1, :, :])))
        ncf.close()

    glor.create_surftemp_RCave(Tsurfname, grdname, title, theta_s, theta_b, hc, N, time, vtransform)

    ncmd = netcdf(Tsurfname, 'a')
    ncmdi = ncmd.variables
    ncmdi['temps'][:] = tsurf
    ncmdi['surf_time'][:] = time / 86400
    ncmdi['lat_rho'][:] = lat_rho
    ncmdi['lon_rho'][:] = lon_rho
    ncmd.close()

