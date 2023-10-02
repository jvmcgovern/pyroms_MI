#
######################################################################
######################################################################
#
#  Main program
#
#  Build a CROCO boundary file using GLORYS12 renanalysis data
#
######################################################################
######################################################################
#

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import num2date
from scipy.interpolate import griddata

import glob
from netCDF4 import Dataset as netcdf
from netCDF4 import date2index as d2i
from calendar import monthrange, isleap
from datetime import date, datetime
import PyCO2SYS as pyco2

import sys

# sys.path.insert(0,'/XXX/')

import croco_vgrid as vgrd
import croco_glorys as glor
from interp_Cgrid import *
from progressbar import *
from scipy.spatial import Delaunay

if 1 == 1:
    # def main_func():

    #
    # #################### USERS DEFINED VARIABLES ########################
    #

    title = 'Boundary file using GLORYS'

    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']

    # my_home_dir = '/home/penven/'

    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    bryfiles_dir = '/media/dskthree/UKESM1-0-LL/'
    grdname = crocofiles_dir + 'CELTIC_grd_v5_MSL_h10.nc'  # LAT to MSL included, INFOMAR block averaged hmin = 10m

    N = 20
    theta_s = 7.
    theta_b = 0.
    hc = 50.
    vtransform = 2

    obc = [1, 1, 1, 1]  # open boundaries (1=open , [S E N W])

    time_bry = 0.
    cycle_bry = 0.

    Yorig = 1990  # year origin of time : days since Yorig-01-01

    # Ystart = 2015
    # Mstart = 1
    # Dstart = 1
    #
    # Yend = 2099
    # Mend = 12
    # Dend = 31

    Ystart = 2000
    Mstart = 1
    Dstart = 1

    Yend = 2014
    Mend = 12
    Dend = 31

    ERA5sfold = '/media/dskone/CELTIC/ECMWF/ERA5/ERA5_CELTIC_PISCES_II/'
    ERA5sfile = 'U10M_Y1999M06.nc'

    CMIPfold = '/media/dsksix/UKESM1-0-LL/NATIVE/'

    PHYfold = 'PHY/'
    BGCfold = 'BGC/'
    ATMfold = 'ATM/'

    CMIP6_ATMstr  = ['rsds_rsus', 'rlds', 'tas', 'psl',  'pr',  'huss', 'uas',  'vas']
    CMIP6_ATMlab = ['rsds', 'rlds', 'tas', 'psl', 'pr', 'huss', 'uas', 'vas']
    CMIP6_ATMtres = ['3hr',       '3hr',  'day', 'day',  '3hr', 'day',  'day',  'day']
    ERA5_ATMstr   = ['SSR',       'STRD', 'T2M', 'pmer', 'TP',  'Q',    'U10M', 'V10M']
    atm_ln       = ['surface_net_solar_radiation', 'surface_thermal_radiation_downwards', '2m_temperature',
                    'mean_sea_level_pressure', 'total_precipitation', 'specific_humidity', '10m_u_component_of_wind',
                    '10m_v_component_of_wind']
    atm_units = ['W m-2', 'W m-2', 'K', 'Pa', 'kg m-2 s-1', 'kg kg-1', 'm s-1', 'm s-1']
    ATMpost = '_ERA5g.nc'

    phytracvars = ['thetao', 'so']
    bgctracvars = ['no3', 'po4', 'si', 'o2', 'dfe', 'dissic', 'talk']
    ele = ['zos']
    vec_vars = ['uo', 'vo']

    model = 'UKESM1-0-LL'

    hist = 'historical'
    # proj = 'ssp245'
    proj = 'ssp585'
    exp = 'r1i1p1f2_gn'
    freq = 'Omon'
    #
    # get the sample extent
    #

    nces = netcdf(ERA5sfold + ERA5sfile, 'r')
    nceso = nces.variables
    lon = np.array(nces['lon'][:])
    lat = np.array(nces['lat'][:])
    lonlen = np.shape(lon)[0]
    latlen = np.shape(lat)[0]
    nces.close()
    # Need to list all netcdfs for each of the variables to interpolate.
    # For each file, it will be possible to know what years are covered by it by identifying the year in file name
    # For each year and month, it will be necessary to identify the index of the pre-, present- and post- month in their
    # respective files and pass the files to the interpolation process.
    # Need to firstly sample one of each files to explore that lat/lon data and do the grid weighting

    for v in range(0, len(ERA5_ATMstr)):
        if CMIP6_ATMtres[v] == '3hr':
            tpad = '_????????????-????????????'
            ysta1 = -34
            yend1 = -30
            ysta2 = -21
            yend2 = -17
            hrinit = 1
            hrend = 22
            mint = 30
            tinit = 0.0625
            tendi = 1
            tincr = 0.125

        elif CMIP6_ATMtres[v] == 'day':
            tpad = '_????????-????????'
            ysta1 = -26
            yend1 = -22
            ysta2 = -17
            yend2 = -13
            hrinit = 12
            hrend = 12
            mint = 00
            tinit = 0.5
            tendi = 1
            tincr = 1

        # CMIP6 ATM files
        CMIP6_hist = sorted(glob.glob(CMIPfold + ATMfold + CMIP6_ATMstr[v] + '_' + CMIP6_ATMtres[v] + '_' +
                                      model + '_' + hist +
                                      '_' + exp + tpad + ATMpost))
        CMIP6_proj = sorted(glob.glob(CMIPfold + ATMfold + CMIP6_ATMstr[v] + '_' + CMIP6_ATMtres[v] + '_' +
                                      model + '_' + proj +
                                      '_' + exp + tpad + ATMpost))
        CMIP6files = CMIP6_hist + CMIP6_proj
        files_yrs = np.ones((len(CMIP6files), 2), dtype=int)
        for file in range(0, len(CMIP6files)):
            files_yrs[file, 0] = int(CMIP6files[file][ysta1:yend1])
            files_yrs[file, 1] = int(CMIP6files[file][ysta2:yend2])

        for iyear in range(Ystart, Yend + 1):
            for imonth in range(Mstart, Mend + 1):

                atmname = CMIPfold + ERA5_ATMstr[v] + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'

                print(' ')
                print(' Making atm file: ' + atmname)
                print(' ')

                #
                # Create the CROCO boundary file
                #
                atime = 0

                glor.create_ERA5style_ESMatm(atmname, ERA5_ATMstr[v], atm_ln[v], atm_units[v], atime, latlen, lonlen)

                if imonth is 1:
                    files_3m = [np.argwhere((files_yrs[:, 0] <= iyear-1) & (files_yrs[:, 1] >= iyear-1))[0][0],
                                np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0]]

                elif imonth == 12:
                    files_3m = [np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((files_yrs[:, 0] <= iyear+1) & (files_yrs[:, 1] >= iyear+1))[0][0]]

                else:
                    files_3m = [np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0]]

                files_3m = np.unique(files_3m)

                if imonth == 1:
                    Tinin_pre = datetime(iyear - 1, 12, monthrange(iyear - 1, 12)[1] - 2,
                                         hrinit, mint, 0)
                    Tinin_post = datetime(iyear, imonth + 1, 3, hrend, mint, 0)
                elif imonth == 12:
                    Tinin_pre = datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                         hrinit, mint, 0)
                    Tinin_post = datetime(iyear + 1, 1, 3, hrend, mint, 0)
                else:
                    Tinin_pre = datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                         hrinit, mint, 0)
                    Tinin_post = datetime(iyear, imonth + 1, 3, hrend, mint, 0)

                for fn in range(0, len(files_3m)):
                    # Get time array per file
                    f1 = netcdf(CMIP6files[files_3m[fn]], 'r')
                    f1o = f1.variables
                    try:
                        Tininx_pre = d2i(Tinin_pre, f1o['time'], select='exact', calendar=f1o['time'].calendar)
                    except:
                        Tininx_pre = 0
                    try:
                        Tininx_post = d2i(Tinin_post, f1o['time'], select='exact', calendar=f1o['time'].calendar)
                        # if Tininx_pre == 0:
                        Tininx_post = Tininx_post + 1
                    except:
                        Tininx_post = f1o['time'].shape[0]
                    if fn == 0:
                        var_time = np.array(f1o['time'][Tininx_pre:Tininx_post])
                        var_ext = np.array(f1o[CMIP6_ATMlab[v]][Tininx_pre:Tininx_post, :, :])
                    else:
                        var_time = np.concatenate((var_time, np.array(f1o['time'][Tininx_pre:Tininx_post])))
                        var_ext = np.concatenate((var_ext,
                                                 np.array(f1o[CMIP6_ATMlab[v]][Tininx_pre:Tininx_post, :, :])))
                    f1.close()

                if imonth == 1:
                    Tnew_pre = date.toordinal(datetime(iyear - 1, 12, monthrange(iyear - 1, 12)[1] - 2,
                                                       hrinit, mint, 0)) - \
                               date.toordinal(datetime(Yorig, 1, 1))

                    Tnew_post = date.toordinal(datetime(iyear, imonth + 1, 3, hrend, mint, 0)) - \
                                date.toordinal(datetime(Yorig, 1, 1))
                elif imonth == 12:
                    Tnew_pre = date.toordinal(datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                                       hrinit, mint, 0)) - \
                               date.toordinal(datetime(Yorig, 1, 1))

                    Tnew_post = date.toordinal(datetime(iyear + 1, 1, 3, hrend, mint, 0)) - \
                                date.toordinal(datetime(Yorig, 1, 1))
                else:
                    Tnew_pre = date.toordinal(datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                                       hrinit, mint, 0)) - \
                               date.toordinal(datetime(Yorig, 1, 1))

                    Tnew_post = date.toordinal(datetime(iyear, imonth + 1, 3, hrend, mint, 0)) - \
                                date.toordinal(datetime(Yorig, 1, 1))

                if CMIP6_ATMtres[v] == '3hr':
                    if (imonth == 1) or (imonth == 3) or (imonth == 8):
                        var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :],
                                                  var_ext[16:256, :, :], var_ext[-32:, :, :]))
                    if (imonth == 4) or (imonth == 6):
                        var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :], var_ext[16:280, :, :]))
                    if (imonth == 5) or (imonth == 7) or (imonth == 12) or (imonth == 10):
                        var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :], var_ext[16:288, :, :]))
                    if (imonth == 9) or (imonth == 11):
                        var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :], var_ext[16:280, :, :]))
                    if (imonth == 2):
                        if isleap(iyear) == 0:
                            var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :],
                                                      var_ext[16:240, :, :], var_ext[-24:, :, :]))
                        elif isleap(iyear) == 1:
                            var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :],
                                                      var_ext[16:248, :, :], var_ext[-24:, :, :]))

                if CMIP6_ATMtres[v] == 'day':
                    if (imonth == 1) or (imonth == 3) or (imonth == 8):
                        var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0),
                                                  var_ext[2:34, :, :], var_ext[-3:, :, :]))
                    if (imonth == 4) or (imonth == 6):
                        var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0), var_ext[1:, :, :]))
                    if (imonth == 5) or (imonth == 7) or (imonth == 12):
                        var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0), var_ext[1:, :, :]))
                    if (imonth == 9) or (imonth == 10) or (imonth == 11):
                        var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0), var_ext[1:, :, :]))
                    if (imonth == 2):
                        if isleap(iyear) == 0:
                            var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0),
                                                      var_ext[2:31, :, :], var_ext[-3:, :, :]))
                        elif isleap(iyear) == 1:
                            var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0),
                                                      var_ext[2:32, :, :], var_ext[-3:, :, :]))

                Tnew = np.arange(Tnew_pre + tinit, Tnew_post + tendi, tincr)

                if len(Tnew) != var_ext.shape[0]:
                    print('something')

                #
                # Open the atmospheric file for writing
                #

                ncatm = netcdf(atmname, 'a')
                ncatmo = ncatm.variables
                ncatmo[ERA5_ATMstr[v]][:] = var_ext
                ncatmo['time'][:] = Tnew
                ncatmo['lat'][:] = lat
                ncatmo['lon'][:] = lon

                #
                #  End loop on time
                #
                ncatm.close()
