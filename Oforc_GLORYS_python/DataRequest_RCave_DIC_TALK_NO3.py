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

    Ystart = 2014
    Mstart = 1
    Dstart = 1

    Yend = 2014
    Mend = 12
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

    crocofil = 'croco_avg_'
    croco_end = '.nc'
    flist = []
    for mth in range(Mstart, Mend + 1):
        flist = flist + sorted(glob.glob(crocodir + crocofil + 'Y' + str(Ystart) + 'M' + str(mth).zfill(2) + croco_end))
        # Tsurfname = crocodir + \
        #             'croco_surf_no3_o2' + '_Y' + str(Ystart) + '_M' + str(Mstart).zfill(2) + \
        #             '_M' + str(Mend).zfill(2) + '.nc'
        # Tsurfname = crocodir + \
        #             'croco_surf_ts' + '_Y' + str(Ystart) + '_M' + str(Mstart).zfill(2) + \
        #             '_M' + str(Mend).zfill(2) + '.nc'
        # Tsurfname = crocodir + \
        #             'croco_bott_no3_o2' + '_Y' + str(Ystart) + '_M' + str(Mstart).zfill(2) + \
        #             '_M' + str(Mend).zfill(2) + '.nc'
        Tsurfname = crocodir + \
                    'croco_bott_ts' + '_Y' + str(Ystart) + '_M' + str(Mstart).zfill(2) + \
                    '_M' + str(Mend).zfill(2) + '.nc'

    for fl in range(0, len(flist)):
        ncf = netcdf(flist[fl], 'r')
        ncfo = ncf.variables
        if fl == 0:
            time = np.array(ncfo['time'][:])
            #
            # temp_surf = np.array(ncfo['temp'][:, ncfo['temp'].shape[1] - 1, :, :])
            # salt_surf = np.array(ncfo['salt'][:, ncfo['salt'].shape[1] - 1, :, :])
            # o2_surf = np.array(ncfo['O2'][:, ncfo['O2'].shape[1] - 1, :, :])
            # no3_surf = np.array(ncfo['NO3'][:, ncfo['NO3'].shape[1] - 1, :, :])
            #
            temp_surf = np.array(ncfo['temp'][:, 0, :, :])
            salt_surf = np.array(ncfo['salt'][:, 0, :, :])
            o2_surf = np.array(ncfo['O2'][:, 0, :, :])
            no3_surf = np.array(ncfo['NO3'][:, 0, :, :])
        else:
            time = np.concatenate((time, np.array(ncfo['time'][:])))
            #
            # temp_surf = np.concatenate((temp_surf, np.array(ncfo['temp'][:, ncfo['temp'].shape[1] - 1, :, :])))
            # salt_surf = np.concatenate((salt_surf, np.array(ncfo['salt'][:, ncfo['salt'].shape[1] - 1, :, :])))
            # o2_surf = np.concatenate((o2_surf, np.array(ncfo['O2'][:, ncfo['O2'].shape[1] - 1, :, :])))
            # no3_surf = np.concatenate((no3_surf, np.array(ncfo['NO3'][:, ncfo['NO3'].shape[1] - 1, :, :])))
            #
            temp_surf = np.concatenate((temp_surf, np.array(ncfo['temp'][:, 0, :, :])))
            salt_surf = np.concatenate((salt_surf, np.array(ncfo['salt'][:, 0, :, :])))
            o2_surf = np.concatenate((o2_surf, np.array(ncfo['O2'][:, 0, :, :])))
            no3_surf = np.concatenate((no3_surf, np.array(ncfo['NO3'][:, 0, :, :])))

        ncf.close()
    temp_surf_abs = temp_surf + 273.15
    # do_c = o2_surf * np.exp(-1 * salt_surf * (0.020573 + ((-12.142 + (2363.1 / temp_surf_abs)) / temp_surf_abs)))
    # do_sat = np.exp(-135.29996 + (1.572288e+05 / temp_surf_abs) - (6.637149e+07 / (temp_surf_abs ** 2)) +
    #                 (1.243678e+10 / (temp_surf_abs ** 3)) -
    #                 (8.621061e+11 / (temp_surf_abs ** 4)) +
    #          (-salt_surf * (0.020573 + (-12.142/temp_surf_abs) + (2363.1/(temp_surf_abs ** 2)))))
    do_c = o2_surf
    do_sat = np.exp(-135.29996 + (1.572288e+05 / temp_surf_abs) - (6.637149e+07 / (temp_surf_abs ** 2)) +
                    (1.243678e+10 / (temp_surf_abs ** 3)) -
                    (8.621061e+11 / (temp_surf_abs ** 4)) +
             (-salt_surf * (0.020573 + (-12.142/temp_surf_abs) + (2363.1/(temp_surf_abs ** 2)))))
    do_sat[temp_surf_abs == 273.15] = np.nan
    do_c[do_c == 0] = np.nan
    dosat_per = (do_c/do_sat)*100

    glor.create_surfbgc_RCave(Tsurfname, grdname, title, theta_s, theta_b, hc, N, time, vtransform)

    ncmd = netcdf(Tsurfname, 'a')
    ncmdi = ncmd.variables
    # ncmdi['NO3'][:] = no3_surf
    # ncmdi['O2'][:] = o2_surf
    # ncmdi['DO_SAT'][:] = o2_surf
    ncmdi['temp'][:] = temp_surf
    ncmdi['salt'][:] = salt_surf
    # ncmdi['surf_time'][:] = time / 86400
    ncmdi['bott_time'][:] = time / 86400
    ncmdi['lat_rho'][:] = lat_rho
    ncmdi['lon_rho'][:] = lon_rho
    ncmd.close()
