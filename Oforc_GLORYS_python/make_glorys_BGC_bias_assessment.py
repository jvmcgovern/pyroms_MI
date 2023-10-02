#
#
######################################################################
######################################################################
#
#  Main program
#
#  Build a CROCO initial file using GLORYS12 renanalysis data
#
# #####################################################################
# #####################################################################
#
#

import glob
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter
from netCDF4 import Dataset as netcdf
from netCDF4 import num2date as n2d
from netCDF4 import date2index as d2i
from datetime import date, datetime
from scipy import interpolate
import cmocean.cm as cmo
import numpy.matlib
import PyCO2SYS as pyco2
import time
import pickle
from pyhdf.SD import SD, SDC
import pandas as pd
from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from math import cos, sin, asin, sqrt, radians
from mpl_toolkits.basemap import Basemap

import cartopy.crs as ccrs
import pyproj
import pyepsg
# import cartopy.feature as cfeature
# from scipy.interpolate import griddata
import scipy.io as scio
# from netCDF4 import date2index as d2i
# from netCDF4 import num2date as n2d
# from datetime import date, datetime
# from calendar import monthrange
# import cftime
# import numpy as np
# import numpy.ma as ma
# import croco_vgrid as vgrd
from progressbar import *

# import sys
# sys.path.insert(0,'')

wod_processing = 0
mi_processing = 0
glodap_processing = 0
dailyres = 1
monthlyres = 0

if dailyres == 1:
    npzstr = 'daily'
else:
    npzstr = 'monthly'

# crocodir = '/media/dskfour/CS1KM_19932021_Uncorr/'
# crocodir = '/media/dskfive/CS1KM_19932021_BCd_results/'
crocodir = '/media/dsksix/CS1KM_19932021_BCd_TNuts_results/'

if wod_processing == 1:

    EMODdir = '/media/dskfour/VAL/EMODNET/INSITU/'
    nut = 'NPSO'
    dic = 'DIC'
    talk = 'TALK'

    WODdir = '/media/dskone/VAL/WOD/'
    WODnut = 'WOD_NO3_PO4_Si_O2_OSD_ragged.nc'

    ncwd = netcdf(WODdir + WODnut, 'r')
    ncwdo = ncwd.variables
    z = np.array(ncwdo['z'][:])
    z_row_size = np.array(ncwdo['z_row_size'][:])
    z_row_idx_upr = np.cumsum(z_row_size)
    z_row_idx_lwr = np.insert(z_row_idx_upr[:-1], 0, 0)
    wodtime = np.array(ncwdo['time'][:])  # days since 1770-01-01 00:00:00 UTC
    wodlat = np.array(ncwdo['lat'][:])
    wodlon = np.array(ncwdo['lon'][:])

    wodnit = np.array(ncwdo['Nitrate'][:])
    wodnit_f = np.array(ncwdo['Nitrate_WODflag'][:])
    wodnit_cs = np.array(ncwdo['Nitrate_row_size'][:])
    wodnit_pf = np.array(ncwdo['Nitrate_WODprofileflag'][:])
    wodnit_ibim = np.zeros_like(wodnit)
    wodnit_cs1k = np.zeros_like(wodnit)
    wodnit_niva = np.zeros_like(wodnit)
    wodnit_woa = np.zeros_like(wodnit)
    wodnit_diva = np.zeros_like(wodnit)
    wodnit_diva_L1 = np.zeros_like(wodnit)
    wodnit_diva_L2 = np.zeros_like(wodnit)

    wodpho = np.array(ncwdo['Phosphate'][:])
    wodpho_f = np.array(ncwdo['Phosphate_WODflag'][:])
    wodpho_cs = np.array(ncwdo['Phosphate_row_size'][:])
    wodpho_pf = np.array(ncwdo['Phosphate_WODprofileflag'][:])
    wodpho_ibim = np.zeros_like(wodpho)
    wodpho_cs1k = np.zeros_like(wodpho)
    wodpho_niva = np.zeros_like(wodpho)
    wodpho_woa = np.zeros_like(wodpho)
    wodpho_diva = np.zeros_like(wodpho)
    wodpho_diva_L1 = np.zeros_like(wodpho)
    wodpho_diva_L2 = np.zeros_like(wodpho)

    wodsil = np.array(ncwdo['Silicate'][:])
    wodsil_f = np.array(ncwdo['Silicate_WODflag'][:])
    wodsil_cs = np.array(ncwdo['Silicate_row_size'][:])
    wodsil_pf = np.array(ncwdo['Silicate_WODprofileflag'][:])
    wodsil_ibim = np.zeros_like(wodsil)
    wodsil_cs1k = np.zeros_like(wodsil)
    wodsil_niva = np.zeros_like(wodsil)
    wodsil_woa = np.zeros_like(wodsil)
    wodsil_diva = np.zeros_like(wodsil)
    wodsil_diva_L1 = np.zeros_like(wodsil)
    wodsil_diva_L2 = np.zeros_like(wodsil)

    wodoxy = np.array(ncwdo['Oxygen'][:])
    wodoxy_f = np.array(ncwdo['Oxygen_WODflag'][:])
    wodoxy_cs = np.array(ncwdo['Oxygen_row_size'][:])
    wodoxy_pf = np.array(ncwdo['Oxygen_WODprofileflag'][:])
    wodoxy_ibim = np.zeros_like(wodoxy)
    wodoxy_cs1k = np.zeros_like(wodoxy)
    wodoxy_niva = np.zeros_like(wodoxy)
    wodoxy_woa = np.zeros_like(wodoxy)
    wodoxy_diva = np.zeros_like(wodoxy)
    wodoxy_diva_L1 = np.zeros_like(wodoxy)
    wodoxy_diva_L2 = np.zeros_like(wodoxy)

    woddays, wod_ind, wod_count = np.unique(np.floor(wodtime), return_index=True, return_counts=True)
    latids = np.zeros_like(wodoxy_pf, dtype='int64')
    lonids = np.zeros_like(wodoxy_pf, dtype='int64')
    latidc = np.zeros_like(wodoxy_pf, dtype='int64')
    lonidc = np.zeros_like(wodoxy_pf, dtype='int64')
    latidn = np.zeros_like(wodoxy_pf, dtype='int64')
    lonidn = np.zeros_like(wodoxy_pf, dtype='int64')
    latidw = np.zeros_like(wodoxy_pf, dtype='int64')
    lonidw = np.zeros_like(wodoxy_pf, dtype='int64')
    latidd = np.zeros_like(wodoxy_pf, dtype='int64')
    lonidd = np.zeros_like(wodoxy_pf, dtype='int64')
    land_ibi = np.zeros_like(wodoxy, dtype='bool')
    land_cs1km = np.zeros_like(wodoxy, dtype='bool')
    land_noresm = np.zeros_like(wodoxy, dtype='bool')
    land_woa = np.zeros_like(wodoxy, dtype='bool')
    land_diva = np.zeros_like(wodoxy, dtype='bool')
    land_diva_L1 = np.zeros_like(wodoxy, dtype='bool')
    land_diva_L2 = np.zeros_like(wodoxy, dtype='bool')
    wod_day = np.zeros_like(wodoxy_pf, dtype='int64')
    wod_month = np.zeros_like(wodoxy_pf, dtype='int64')
    wod_year = np.zeros_like(wodoxy_pf, dtype='int64')
    # array for sorting instead of profiles
    wodlat_arr = np.zeros_like(wodoxy)
    wodlon_arr = np.zeros_like(wodoxy)
    wod_day_arr = np.zeros_like(wodoxy, dtype='int64')
    wod_month_arr = np.zeros_like(wodoxy, dtype='int64')
    wod_year_arr = np.zeros_like(wodoxy, dtype='int64')
    #
    time_wdo = n2d(ncwdo['time'][wod_ind], units=ncwdo['time'].units)
    time_wdo_profiles = n2d(ncwdo['time'][:], units=ncwdo['time'].units)
    ncwd.close()

    # Block for IBI processing
    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    PISCES24_prefixm = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_'
    PISCES24_prefixd = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    PISCES24_ending = '_R*_RE01.nc'
    IBI_BIO_file_sample = sorted(glob.glob(PISCES24files_dir +
                                           PISCES24_prefixm +
                                           '*_*' + PISCES24_ending))
    wdibis = netcdf(IBI_BIO_file_sample[0], 'r')
    wdibso = wdibis.variables
    ibilat = np.array(wdibso['latitude'][:])
    ibilon = np.array(wdibso['longitude'][:])
    ibidep = np.array(wdibso['depth'][:])
    wdibis.close()

    for p in range(0, len(wodlat)):
        latids[p] = glor.geo_idx(wodlat[p], ibilat)
        lonids[p] = glor.geo_idx(wodlon[p], ibilon)

    for wd in range(0, len(time_wdo_profiles)):
        wod_day[wd] = time_wdo_profiles[wd].day
        wod_month[wd] = time_wdo_profiles[wd].month
        wod_year[wd] = time_wdo_profiles[wd].year

    for wd in range(0, len(woddays)):
        if monthlyres == 1:

            IBI_BIO_file_m = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixm + str(time_wdo[wd].year) +
                                              str(time_wdo[wd].month).zfill(2) +
                                              '*_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_m[0], 'r')
        if dailyres == 1:

            IBI_BIO_file_d = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixd + str(time_wdo[wd].year) +
                                              str(time_wdo[wd].month).zfill(2) + str(time_wdo[wd].day).zfill(2) +
                                              '_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_d[0], 'r')
        wdibmo = wdibim.variables

        ibinit = np.array(wdibmo['no3'][:])
        ibipho = np.array(wdibmo['po4'][:])
        ibisil = np.array(wdibmo['si'][:])
        ibioxy = np.array(wdibmo['o2'][:])
        for dp in range(0, wod_count[wd]):
            ibinit_smp = ibinit[0, :, latids[wod_ind[wd] + dp], lonids[wod_ind[wd] + dp]]
            ibipho_smp = ibipho[0, :, latids[wod_ind[wd] + dp], lonids[wod_ind[wd] + dp]]
            ibisil_smp = ibisil[0, :, latids[wod_ind[wd] + dp], lonids[wod_ind[wd] + dp]]
            ibioxy_smp = ibioxy[0, :, latids[wod_ind[wd] + dp], lonids[wod_ind[wd] + dp]]
            if (any(ibinit_smp < 0) == True) & (all(ibinit_smp < 0) == False):
                ibinit_smp[ibinit_smp < 0] = ibinit_smp[ibinit_smp >= 0][-1]
                ibipho_smp[ibipho_smp < 0] = ibipho_smp[ibipho_smp >= 0][-1]
                ibisil_smp[ibisil_smp < 0] = ibisil_smp[ibisil_smp >= 0][-1]
                ibioxy_smp[ibioxy_smp < 0] = ibioxy_smp[ibioxy_smp >= 0][-1]
            elif all(ibinit_smp < 0) == True:
                land_ibi[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True
                # land_ibi[wod_ind[wd] + dp] = True
            #
            wodint_n = interpolate.interp1d(ibidep, ibinit_smp, fill_value='extrapolate')
            wodnit_ibim[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_p = interpolate.interp1d(ibidep, ibipho_smp, fill_value='extrapolate')
            wodpho_ibim[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_s = interpolate.interp1d(ibidep, ibisil_smp, fill_value='extrapolate')
            wodsil_ibim[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_o = interpolate.interp1d(ibidep, ibioxy_smp, fill_value='extrapolate')
            wodoxy_ibim[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodlat_arr[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = wodlat[wod_ind[wd]]
            wodlon_arr[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = wodlon[wod_ind[wd]]
            wod_day_arr[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = wod_day[wod_ind[wd]]
            wod_month_arr[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = wod_month[wod_ind[wd]]
            wod_year_arr[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = wod_year[wod_ind[wd]]

    # # Block for CS1KM processing

    croco_finit = 'croco_avg_Y'
    croco_fend = '.nc'
    CS1KM_BIO_file_sample = sorted(glob.glob(crocodir +
                                             croco_finit +
                                             '*' + croco_fend))
    wdcs1kms = netcdf(CS1KM_BIO_file_sample[0], 'r')
    wdcs1kmso = wdcs1kms.variables
    cs1kmlat = np.array(wdcs1kmso['lat_rho'][:])
    cs1kmlon = np.array(wdcs1kmso['lon_rho'][:])
    wdcs1kms.close()

    for p in range(0, len(wodlat)):
        latidc[p] = glor.geo_idx(wodlat[p], cs1kmlat[:, 0])
        lonidc[p] = glor.geo_idx(wodlon[p], cs1kmlon[0, :])

    for wd in range(0, len(woddays)):
        CS1KM_BIO_file = sorted(glob.glob(crocodir +
                                          croco_finit + str(time_wdo[wd].year) + 'M' +
                                          str(time_wdo[wd].month).zfill(2) +
                                          croco_fend))
        wdcs1km = netcdf(CS1KM_BIO_file[0], 'r')
        wdcs1kmo = wdcs1km.variables

        if monthlyres == 1:
            # Monthly read latitude and longitude for indexing
            cs1kmnit = np.mean(np.array(wdcs1kmo['NO3'][:]), 0)
            cs1kmpho = np.mean(np.array(wdcs1kmo['PO4'][:]), 0)
            cs1kmsil = np.mean(np.array(wdcs1kmo['Si'][:]), 0)
            cs1kmoxy = np.mean(np.array(wdcs1kmo['O2'][:]), 0)
            cs1kmz = np.mean(np.array(wdcs1kmo['zeta'][:]), 0)

        if dailyres == 1:
            # Daily read latitude and longitude for indexing
            cs1kmnit = np.array(wdcs1kmo['NO3'][time_wdo[wd].day - 1, :, :, :])
            cs1kmpho = np.array(wdcs1kmo['PO4'][time_wdo[wd].day - 1, :, :, :])
            cs1kmsil = np.array(wdcs1kmo['Si'][time_wdo[wd].day - 1, :, :, :])
            cs1kmoxy = np.array(wdcs1kmo['O2'][time_wdo[wd].day - 1, :, :, :])
            cs1kmz = np.array(wdcs1kmo['zeta'][time_wdo[wd].day - 1, :, :])

        # LAT to MSL included, INFOMAR block averaged hmin = 10m
        grdname = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/' + 'CELTIC_grd_v5_MSL_h10.nc'
        ncg = netcdf(grdname, 'r')
        ncgrd = ncg.variables
        h = np.array(ncgrd['h'][:])
        ncg.close()
        N = 20
        theta_s = 7.
        theta_b = 0.
        hc = 20.
        vtransform = 2
        cs1kmd = vgrd.zlevs(h, cs1kmz, theta_s, theta_b, 50, N, 'r', vtransform) * -1
        for dp in range(0, wod_count[wd]):
            cs1kmnit_smp = cs1kmnit[:, latidc[wod_ind[wd] + dp], lonidc[wod_ind[wd] + dp]]
            cs1kmpho_smp = cs1kmpho[:, latidc[wod_ind[wd] + dp], lonidc[wod_ind[wd] + dp]]
            cs1kmsil_smp = cs1kmsil[:, latidc[wod_ind[wd] + dp], lonidc[wod_ind[wd] + dp]]
            cs1kmoxy_smp = cs1kmoxy[:, latidc[wod_ind[wd] + dp], lonidc[wod_ind[wd] + dp]]
            cs1kmd_smp = cs1kmd[:, latidc[wod_ind[wd] + dp], lonidc[wod_ind[wd] + dp]]
            if (any(cs1kmnit_smp == 0) == True) & (all(cs1kmnit_smp == 0) == False):
                cs1kmnit_smp[cs1kmnit_smp == 0] = cs1kmnit_smp[cs1kmnit_smp > 0][-1]
                cs1kmpho_smp[cs1kmpho_smp == 0] = cs1kmpho_smp[cs1kmpho_smp > 0][-1]
                cs1kmsil_smp[cs1kmsil_smp == 0] = cs1kmsil_smp[cs1kmsil_smp > 0][-1]
                cs1kmoxy_smp[cs1kmoxy_smp == 0] = cs1kmoxy_smp[cs1kmoxy_smp > 0][-1]
            elif all(cs1kmnit_smp < 0) == True:
                # land_cs1km[wod_ind[wd] + dp] = True
                land_cs1km[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True
            #
            wodint_n = interpolate.interp1d(cs1kmd_smp, cs1kmnit_smp, fill_value='extrapolate')
            wodnit_cs1k[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_p = interpolate.interp1d(cs1kmd_smp, cs1kmpho_smp, fill_value='extrapolate')
            wodpho_cs1k[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_s = interpolate.interp1d(cs1kmd_smp, cs1kmsil_smp, fill_value='extrapolate')
            wodsil_cs1k[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_o = interpolate.interp1d(cs1kmd_smp, cs1kmoxy_smp, fill_value='extrapolate')
            wodoxy_cs1k[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])

    # Block for NORESM processing
    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']
    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
    ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
    ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
    ncnoro2 = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[3] + NORESM_ending, 'r')
    ncnornio = ncnorni.variables
    ncnorpoo = ncnorpo.variables
    ncnorsio = ncnorsi.variables
    ncnoro2o = ncnoro2.variables
    noresmlat = np.array(ncnornio['lat'][:])
    noresmlon = np.array(ncnornio['lon'][:])
    depthniva = np.array(ncnornio['depth'][:])
    time_niva = np.array(ncnornio['time'][:])

    for p in range(0, len(wodlat)):
        latidn[p] = glor.geo_idx(wodlat[p], noresmlat)
        lonidn[p] = glor.geo_idx(wodlon[p], noresmlon)

    latidn = latidn - (noresmlat.shape[1] * np.floor(latidn / noresmlat.shape[1]))
    lonidn = lonidn - (noresmlat.shape[1] * np.floor(lonidn / noresmlat.shape[1]))
    latidn = latidn.astype(int)
    lonidn = lonidn.astype(int)

    for wd in range(0, len(woddays)):
        noresmnit = np.array(ncnornio[NIVAvars[0]][:]) * 1000
        noresmpho = np.array(ncnorpoo[NIVAvars[1]][:]) * 1000
        noresmsil = np.array(ncnorsio[NIVAvars[2]][:]) * 1000
        noresmoxy = np.array(ncnoro2o[NIVAvars[3]][:]) * 1000

        for dp in range(0, wod_count[wd]):
            if monthlyres == 1:
                # By month
                Tinin = datetime(time_wdo[wd].year, time_wdo[wd].month, 15)
                Tininx1 = d2i(Tinin, ncnornio['time'], select='exact')
                noresmnit_smp = noresmnit[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                noresmpho_smp = noresmpho[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                noresmsil_smp = noresmsil[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                noresmoxy_smp = noresmoxy[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
            if dailyres == 1:
                # By day
                if time_wdo[wd].day != 15:
                    Tinin = datetime(time_wdo[wd].year, time_wdo[wd].month, time_wdo[wd].day)
                    Tininx1 = d2i(Tinin, ncnornio['time'], select='before')
                    Tininx2 = d2i(Tinin, ncnornio['time'], select='after')
                    if time_wdo[wd].day > 15:
                        tmult = time_wdo[wd].day - 15
                    else:
                        tmult = time_wdo[wd].day + 15

                    # Tidx1
                    noresmnit_1 = noresmnit[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmpho_1 = noresmpho[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmsil_1 = noresmsil[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmoxy_1 = noresmoxy[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]

                    # Tidx2
                    noresmnit_2 = noresmnit[Tininx2, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmpho_2 = noresmpho[Tininx2, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmsil_2 = noresmsil[Tininx2, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmoxy_2 = noresmoxy[Tininx2, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    
                    # Interpolating to the day
                    noresmnit_smp = noresmnit_1 + (((noresmnit_2 - noresmnit_1)/30)*tmult)
                    noresmpho_smp = noresmpho_1 + (((noresmpho_2 - noresmpho_1)/30)*tmult)
                    noresmsil_smp = noresmsil_1 + (((noresmsil_2 - noresmsil_1)/30)*tmult)
                    noresmoxy_smp = noresmoxy_1 + (((noresmoxy_2 - noresmoxy_1)/30)*tmult)
                else:
                    Tinin = datetime(time_wdo[wd].year, time_wdo[wd].month, 15)
                    Tininx1 = d2i(Tinin, ncnornio['time'], select='exact')
                    noresmnit_smp = noresmnit[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmpho_smp = noresmpho[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmsil_smp = noresmsil[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]
                    noresmoxy_smp = noresmoxy[Tininx1, :, latidn[wod_ind[wd] + dp], lonidn[wod_ind[wd] + dp]]

            if (any(noresmnit_smp > 10000) == True) & (all(noresmnit_smp > 10000) == False):
                noresmnit_smp[noresmnit_smp > 10000] = noresmnit_smp[noresmnit_smp < 10000][-1]
                noresmpho_smp[noresmpho_smp > 10000] = noresmpho_smp[noresmpho_smp < 10000][-1]
                noresmsil_smp[noresmsil_smp > 10000] = noresmsil_smp[noresmsil_smp < 10000][-1]
                noresmoxy_smp[noresmoxy_smp > 10000] = noresmoxy_smp[noresmoxy_smp < 10000][-1]
            elif all(noresmnit_smp > 10000) == True:
                # land_noresm[wod_ind[wd] + dp] = True
                land_noresm[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True
            #
            wodint_n = interpolate.interp1d(depthniva, noresmnit_smp, fill_value='extrapolate')
            wodnit_niva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_p = interpolate.interp1d(depthniva, noresmpho_smp, fill_value='extrapolate')
            wodpho_niva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_s = interpolate.interp1d(depthniva, noresmsil_smp, fill_value='extrapolate')
            wodsil_niva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_o = interpolate.interp1d(depthniva, noresmoxy_smp, fill_value='extrapolate')
            wodoxy_niva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])

    ncnorni.close()
    ncnorpo.close()
    ncnorsi.close()
    ncnoro2.close()

    # Block for WOA processing
    valdir = '/media/dskone/VAL/'
    woadir = valdir + 'WOA18/'

    # Common woa elements for monthly climatology
    woa_exn = 'woa18_all_n01_01.nc'
    woa_exo = 'woa18_all_o01_01.nc'
    woan_01 = netcdf(woadir + woa_exn, 'r')
    woao_01 = netcdf(woadir + woa_exo, 'r')
    woan_01o = woan_01.variables
    woao_01o = woao_01.variables
    woa_lat = np.array(woan_01o['lat'][:])
    woa_lon = np.array(woan_01o['lon'][:])
    woa_depth_nps = np.array(woan_01o['depth'][:])
    woa_depth_nps_bnds = np.array(woan_01o['depth_bnds'][:])
    n_01a = np.array(woan_01o['n_an'][:])
    o_01a = np.array(woao_01o['o_an'][:])
    woan_01.close()
    woao_01.close()

    woa_sta = 'woa18_all_'
    woa_end = '_01.nc'

    woa_n = np.tile(n_01a, (12, 1, 1, 1))
    woa_o = np.tile(o_01a, (12, 1, 1, 1))
    woa_n = np.zeros_like(woa_n)
    woa_p = np.zeros_like(woa_n)
    woa_s = np.zeros_like(woa_n)
    woa_o = np.zeros_like(woa_o)

    for mt in range(0, 12):
        # nitrate
        woa_n_mth = netcdf(woadir + woa_sta + 'n' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_n_mth_o = woa_n_mth.variables
        woa_n[mt, :, :, :] = np.array(woa_n_mth_o['n_an'][:])
        woa_depth_nps = np.array(woa_n_mth_o['depth'][:])
        woa_depth_nps_bnds = np.array(woa_n_mth_o['depth_bnds'][:])
        woa_n_mth.close()

        # phosphate
        woa_p_mth = netcdf(woadir + woa_sta + 'p' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_p_mth_o = woa_p_mth.variables
        woa_p[mt, :, :, :] = np.array(woa_p_mth_o['p_an'][:])
        woa_p_mth.close()

        # silicate
        woa_s_mth = netcdf(woadir + woa_sta + 'i' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_s_mth_o = woa_s_mth.variables
        woa_s[mt, :, :, :] = np.array(woa_s_mth_o['i_an'][:])
        woa_s_mth.close()

        # oxygen
        woa_o_mth = netcdf(woadir + woa_sta + 'o' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_o_mth_o = woa_o_mth.variables
        woa_o[mt, :, :, :] = np.array(woa_o_mth_o['o_an'][:])
        woa_depth_ox = np.array(woa_o_mth_o['depth'][:])
        woa_depth_ox_bnds = np.array(woa_o_mth_o['depth_bnds'][:])
        woa_o_mth.close()

    for p in range(0, len(wodlat)):
        latidw[p] = glor.geo_idx(wodlat[p], woa_lat)
        lonidw[p] = glor.geo_idx(wodlon[p], woa_lon)

    for wd in range(0, len(woddays)):
        for dp in range(0, wod_count[wd]):
            if monthlyres == 1:
                # Month
                woanit_smp = woa_n[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                woapho_smp = woa_p[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                woasil_smp = woa_s[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                woaoxy_smp = woa_o[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
            if dailyres == 1:
                # Day
                if time_wdo[wd].day > 15:
                    tmult = time_wdo[wd].day - 15
                    tidx1 = time_wdo[wd].month - 1
                    if time_wdo[wd].month == 12:
                        tidx2 = 0
                    else:
                        tidx2 = time_wdo[wd].month

                    # Tidx1
                    woanit_1 = woa_n[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woapho_1 = woa_p[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woasil_1 = woa_s[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woaoxy_1 = woa_o[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]

                    # Tidx2
                    woanit_2 = woa_n[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woapho_2 = woa_p[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woasil_2 = woa_s[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woaoxy_2 = woa_o[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]

                    # Interpolating to the day
                    woanit_smp = woanit_1 + (((woanit_2 - woanit_1)/30)*tmult)
                    woapho_smp = woapho_1 + (((woapho_2 - woapho_1)/30)*tmult)
                    woasil_smp = woasil_1 + (((woasil_2 - woasil_1)/30)*tmult)
                    woaoxy_smp = woaoxy_1 + (((woaoxy_2 - woaoxy_1)/30)*tmult)

                elif time_wdo[wd].day < 15:
                    tmult = time_wdo[wd].day + 15
                    tidx2 = time_wdo[wd].month - 1
                    if time_wdo[wd].month == 1:
                        tidx1 = 11
                    else:
                        tidx1 = time_wdo[wd].month - 2

                    # Tidx1
                    woanit_1 = woa_n[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woapho_1 = woa_p[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woasil_1 = woa_s[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woaoxy_1 = woa_o[tidx1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]

                    # Tidx2
                    woanit_2 = woa_n[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woapho_2 = woa_p[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woasil_2 = woa_s[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woaoxy_2 = woa_o[tidx2, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]

                    # Interpolating to the day
                    woanit_smp = woanit_1 + (((woanit_2 - woanit_1)/30)*tmult)
                    woapho_smp = woapho_1 + (((woapho_2 - woapho_1)/30)*tmult)
                    woasil_smp = woasil_1 + (((woasil_2 - woasil_1)/30)*tmult)
                    woaoxy_smp = woaoxy_1 + (((woaoxy_2 - woaoxy_1)/30)*tmult)

                else:
                    woanit_smp = woa_n[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woapho_smp = woa_p[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woasil_smp = woa_s[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
                    woaoxy_smp = woa_o[time_wdo[wd].month - 1, :, latidw[wod_ind[wd] + dp], lonidw[wod_ind[wd] + dp]]
            
            if (any(woanit_smp > 10000) == True) & (all(woanit_smp > 10000) == False):
                woanit_smp[woanit_smp > 10000] = woanit_smp[woanit_smp < 10000][-1]
                woapho_smp[woapho_smp > 10000] = woapho_smp[woapho_smp < 10000][-1]
                woasil_smp[woasil_smp > 10000] = woasil_smp[woasil_smp < 10000][-1]
                woaoxy_smp[woaoxy_smp > 10000] = woaoxy_smp[woaoxy_smp < 10000][-1]

            elif all(woanit_smp > 10000) == True:
                # land_woa[wod_ind[wd] + dp] = True
                land_woa[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True
            #
            wodint_n = interpolate.interp1d(woa_depth_nps, woanit_smp, fill_value='extrapolate')
            wodnit_woa[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_p = interpolate.interp1d(woa_depth_nps, woapho_smp, fill_value='extrapolate')
            wodpho_woa[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_s = interpolate.interp1d(woa_depth_nps, woasil_smp, fill_value='extrapolate')
            wodsil_woa[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_o = interpolate.interp1d(woa_depth_ox, woaoxy_smp, fill_value='extrapolate')
            wodoxy_woa[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])

    # Block for EMODNET DIVA processing

    emcdir = '/media/dsktwo/VAL/EMODNET/CLIMATOLOGY/'

    # Common emc elements for monthly climatology
    emc_ex = 'Water_body_dissolved_inorganic_nitrogen_v2021.nc'
    emc_01 = netcdf(emcdir + emc_ex, 'r')
    emc_01o = emc_01.variables
    emc_lat = np.array(emc_01o['lat'][:])
    emc_lon = np.array(emc_01o['lon'][:])
    emc_depth = np.array(emc_01o['depth'][:])
    n_01a = np.array(emc_01o['Water body dissolved inorganic nitrogen (DIN)'][:])
    emc_01.close()

    emc_sta = 'Water_body_'
    emc_end = '_v2021.nc'

    # nitrate
    emc_n_mth = netcdf(emcdir + emc_sta + 'dissolved_inorganic_nitrogen' + emc_end, 'r')
    emc_n_mth_o = emc_n_mth.variables
    emc_n = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)'][:])
    emc_n_L1 = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)_L1'][:])
    emc_n_L2 = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)_L2'][:])
    emc_n_mth.close()

    # phosphate
    emc_p_mth = netcdf(emcdir + emc_sta + 'phosphate' + emc_end, 'r')
    emc_p_mth_o = emc_p_mth.variables
    emc_p = np.array(emc_p_mth_o['Water body phosphate'][:])
    emc_p_L1 = np.array(emc_p_mth_o['Water body phosphate_L1'][:])
    emc_p_L2 = np.array(emc_p_mth_o['Water body phosphate_L2'][:])
    emc_p_mth.close()

    # silicate
    emc_s_mth = netcdf(emcdir + emc_sta + 'silicate' + emc_end, 'r')
    emc_s_mth_o = emc_s_mth.variables
    emc_s = np.array(emc_s_mth_o['Water body silicate'][:])
    emc_s_L1 = np.array(emc_s_mth_o['Water body silicate_L1'][:])
    emc_s_L2 = np.array(emc_s_mth_o['Water body silicate_L2'][:])
    emc_s_mth.close()

    # oxygen
    emc_o_mth = netcdf(emcdir + emc_sta + 'dissolved_oxygen_concentration' + emc_end, 'r')
    emc_o_mth_o = emc_o_mth.variables
    emc_o = np.array(emc_o_mth_o['Water body dissolved oxygen concentration'][:])
    emc_o_L1 = np.array(emc_o_mth_o['Water body dissolved oxygen concentration_L1'][:])
    emc_o_L2 = np.array(emc_o_mth_o['Water body dissolved oxygen concentration_L2'][:])
    emc_o_mth.close()

    for p in range(0, len(wodlat)):
        latidd[p] = glor.geo_idx(wodlat[p], emc_lat)
        lonidd[p] = glor.geo_idx(wodlon[p], emc_lon)

    for wd in range(0, len(woddays)):
        for dp in range(0, wod_count[wd]):
            if monthlyres == 1:
                # Month
                emcnit_smp = emc_n[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcnit_smp_L1 = emc_n_L1[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcnit_smp_L2 = emc_n_L2[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcpho_smp = emc_p[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcpho_smp_L1 = emc_p_L1[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcpho_smp_L2 = emc_p_L2[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcsil_smp = emc_s[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcsil_smp_L1 = emc_s_L1[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcsil_smp_L2 = emc_s_L2[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcoxy_smp = emc_o[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcoxy_smp_L1 = emc_o_L1[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                emcoxy_smp_L2 = emc_o_L2[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
            if dailyres == 1:
                
                # Daily
                if time_wdo[wd].day > 15:
                    tmult = time_wdo[wd].day - 15
                    tidx1 = time_wdo[wd].month - 1
                    if time_wdo[wd].month == 12:
                        tidx2 = 0
                    else:
                        tidx2 = time_wdo[wd].month
    
                    # Tidx1
                    emcnit_1 = emc_n[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_1 = emc_p[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_1 = emc_s[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_1 = emc_o[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    
                    emcnit_1_L1 = emc_n_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_1_L1 = emc_p_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_1_L1 = emc_s_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_1_L1 = emc_o_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]

                    emcnit_1_L2 = emc_n_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_1_L2 = emc_p_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_1_L2 = emc_s_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_1_L2 = emc_o_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
    
                    # Tidx2
                    emcnit_2 = emc_n[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_2 = emc_p[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_2 = emc_s[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_2 = emc_o[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]

                    emcnit_2_L1 = emc_n_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_2_L1 = emc_p_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_2_L1 = emc_s_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_2_L1 = emc_o_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]

                    emcnit_2_L2 = emc_n_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_2_L2 = emc_p_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_2_L2 = emc_s_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_2_L2 = emc_o_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
    
                    # Interpolating to the day
                    emcnit_smp = emcnit_1 + (((emcnit_2 - emcnit_1)/30)*tmult)
                    emcpho_smp = emcpho_1 + (((emcpho_2 - emcpho_1)/30)*tmult)
                    emcsil_smp = emcsil_1 + (((emcsil_2 - emcsil_1)/30)*tmult)
                    emcoxy_smp = emcoxy_1 + (((emcoxy_2 - emcoxy_1)/30)*tmult)

                    emcnit_smp_L1 = emcnit_1_L1 + (((emcnit_2_L1 - emcnit_1_L1)/30)*tmult)
                    emcpho_smp_L1 = emcpho_1_L1 + (((emcpho_2_L1 - emcpho_1_L1)/30)*tmult)
                    emcsil_smp_L1 = emcsil_1_L1 + (((emcsil_2_L1 - emcsil_1_L1)/30)*tmult)
                    emcoxy_smp_L1 = emcoxy_1_L1 + (((emcoxy_2_L1 - emcoxy_1_L1)/30)*tmult)

                    emcnit_smp_L2 = emcnit_1_L2 + (((emcnit_2_L2 - emcnit_1_L2)/30)*tmult)
                    emcpho_smp_L2 = emcpho_1_L2 + (((emcpho_2_L2 - emcpho_1_L2)/30)*tmult)
                    emcsil_smp_L2 = emcsil_1_L2 + (((emcsil_2_L2 - emcsil_1_L2)/30)*tmult)
                    emcoxy_smp_L2 = emcoxy_1_L2 + (((emcoxy_2_L2 - emcoxy_1_L2)/30)*tmult)
    
                elif time_wdo[wd].day < 15:
                    tmult = time_wdo[wd].day + 15
                    tidx2 = time_wdo[wd].month - 1
                    if time_wdo[wd].month == 1:
                        tidx1 = 11
                    else:
                        tidx1 = time_wdo[wd].month - 2
    
                    # Tidx1
                    emcnit_1 = emc_n[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_1 = emc_p[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_1 = emc_s[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_1 = emc_o[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    
                    emcnit_1_L1 = emc_n_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_1_L1 = emc_p_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_1_L1 = emc_s_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_1_L1 = emc_o_L1[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]

                    emcnit_1_L2 = emc_n_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_1_L2 = emc_p_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_1_L2 = emc_s_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_1_L2 = emc_o_L2[tidx1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
    
                    # Tidx2
                    emcnit_2_L1 = emc_n_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_2_L1 = emc_p_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_2_L1 = emc_s_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_2_L1 = emc_o_L1[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]

                    emcnit_2_L2 = emc_n_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_2_L2 = emc_p_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_2_L2 = emc_s_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_2_L2 = emc_o_L2[tidx2, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
    
                    # Interpolating to the day
                    emcnit_smp = emcnit_1 + (((emcnit_2 - emcnit_1)/30)*tmult)
                    emcpho_smp = emcpho_1 + (((emcpho_2 - emcpho_1)/30)*tmult)
                    emcsil_smp = emcsil_1 + (((emcsil_2 - emcsil_1)/30)*tmult)
                    emcoxy_smp = emcoxy_1 + (((emcoxy_2 - emcoxy_1)/30)*tmult)

                    emcnit_smp_L1 = emcnit_1_L1 + (((emcnit_2_L1 - emcnit_1_L1)/30)*tmult)
                    emcpho_smp_L1 = emcpho_1_L1 + (((emcpho_2_L1 - emcpho_1_L1)/30)*tmult)
                    emcsil_smp_L1 = emcsil_1_L1 + (((emcsil_2_L1 - emcsil_1_L1)/30)*tmult)
                    emcoxy_smp_L1 = emcoxy_1_L1 + (((emcoxy_2_L1 - emcoxy_1_L1)/30)*tmult)

                    emcnit_smp_L2 = emcnit_1_L2 + (((emcnit_2_L2 - emcnit_1_L2)/30)*tmult)
                    emcpho_smp_L2 = emcpho_1_L2 + (((emcpho_2_L2 - emcpho_1_L2)/30)*tmult)
                    emcsil_smp_L2 = emcsil_1_L2 + (((emcsil_2_L2 - emcsil_1_L2)/30)*tmult)
                    emcoxy_smp_L2 = emcoxy_1_L2 + (((emcoxy_2_L2 - emcoxy_1_L2)/30)*tmult)
    
                else:
                    emcnit_smp = emc_n[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcpho_smp = emc_p[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcsil_smp = emc_s[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]
                    emcoxy_smp = emc_o[time_wdo[wd].month - 1, :, latidd[wod_ind[wd] + dp], lonidd[wod_ind[wd] + dp]]

            if (any(emcnit_smp > 10000) == True) & (all(emcnit_smp > 10000) == False):
                emcnit_smp[emcnit_smp > 10000] = emcnit_smp[emcnit_smp < 10000][-1]
                emcpho_smp[emcpho_smp > 10000] = emcpho_smp[emcpho_smp < 10000][-1]
                emcsil_smp[emcsil_smp > 10000] = emcsil_smp[emcsil_smp < 10000][-1]
                emcoxy_smp[emcoxy_smp > 10000] = emcoxy_smp[emcoxy_smp < 10000][-1]
            elif all(emcnit_smp > 10000) == True:
                # land_diva[wod_ind[wd] + dp] = True
                land_diva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True

            if (any(emcnit_smp_L1 > 10000) == True) & (all(emcnit_smp_L1 > 10000) == False):
                emcnit_smp_L1[emcnit_smp_L1 > 10000] = emcnit_smp_L1[emcnit_smp_L1 < 10000][-1]
            elif all(emcnit_smp_L1 > 10000) == True:
                # land_diva_L1[wod_ind[wd] + dp] = True
                land_diva_L1[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True

            if (any(emcnit_smp_L2 > 10000) == True) & (all(emcnit_smp_L2 > 10000) == False):
                emcnit_smp_L2[emcnit_smp_L2 > 10000] = emcnit_smp_L2[emcnit_smp_L2 < 10000][-1]
            elif all(emcnit_smp_L2 > 10000) == True:
                # land_diva_L2[wod_ind[wd] + dp] = True
                land_diva_L2[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = True

            if (any(emcpho_smp_L1 > 10000) == True) & (all(emcpho_smp_L1 > 10000) == False):
                emcpho_smp_L1[emcpho_smp_L1 > 10000] = emcpho_smp_L1[emcpho_smp_L1 < 10000][-1]
            if (any(emcpho_smp_L2 > 10000) == True) & (all(emcpho_smp_L2 > 10000) == False):
                emcpho_smp_L2[emcpho_smp_L2 > 10000] = emcpho_smp_L2[emcpho_smp_L2 < 10000][-1]

            if (any(emcsil_smp_L1 > 10000) == True) & (all(emcsil_smp_L1 > 10000) == False):
                emcsil_smp_L1[emcsil_smp_L1 > 10000] = emcsil_smp_L1[emcsil_smp_L1 < 10000][-1]
            if (any(emcsil_smp_L2 > 10000) == True) & (all(emcsil_smp_L2 > 10000) == False):
                emcsil_smp_L2[emcsil_smp_L2 > 10000] = emcsil_smp_L2[emcsil_smp_L2 < 10000][-1]

            if (any(emcoxy_smp_L1 > 10000) == True) & (all(emcoxy_smp_L1 > 10000) == False):
                emcoxy_smp_L1[emcoxy_smp_L1 > 10000] = emcoxy_smp_L1[emcoxy_smp_L1 < 10000][-1]
            if (any(emcoxy_smp_L2 > 10000) == True) & (all(emcoxy_smp_L2 > 10000) == False):
                emcoxy_smp_L2[emcoxy_smp_L2 > 10000] = emcoxy_smp_L2[emcoxy_smp_L2 < 10000][-1]
            #
            wodint_n = interpolate.interp1d(emc_depth, emcnit_smp, fill_value='extrapolate')
            wodnit_diva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_n_L1 = interpolate.interp1d(emc_depth, emcnit_smp_L1, fill_value='extrapolate')
            wodnit_diva_L1[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n_L1(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_n_L2 = interpolate.interp1d(emc_depth, emcnit_smp_L2, fill_value='extrapolate')
            wodnit_diva_L2[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_n_L2(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_p = interpolate.interp1d(emc_depth, emcpho_smp, fill_value='extrapolate')
            wodpho_diva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_p_L1 = interpolate.interp1d(emc_depth, emcpho_smp_L1, fill_value='extrapolate')
            wodpho_diva_L1[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p_L1(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_p_L2 = interpolate.interp1d(emc_depth, emcpho_smp_L2, fill_value='extrapolate')
            wodpho_diva_L2[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_p_L2(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_s = interpolate.interp1d(emc_depth, emcsil_smp, fill_value='extrapolate')
            wodsil_diva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_s_L1 = interpolate.interp1d(emc_depth, emcsil_smp_L1, fill_value='extrapolate')
            wodsil_diva_L1[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s_L1(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_s_L2 = interpolate.interp1d(emc_depth, emcsil_smp_L2, fill_value='extrapolate')
            wodsil_diva_L2[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_s_L2(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            #
            wodint_o = interpolate.interp1d(emc_depth, emcoxy_smp, fill_value='extrapolate')
            wodoxy_diva[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_o_L1 = interpolate.interp1d(emc_depth, emcoxy_smp_L1, fill_value='extrapolate')
            wodoxy_diva_L1[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o_L1(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])
            wodint_o_L2 = interpolate.interp1d(emc_depth, emcoxy_smp_L2, fill_value='extrapolate')
            wodoxy_diva_L2[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
                wodint_o_L2(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])

    wod_bias_file = crocodir + 'wod_bias_' + npzstr + '.npz'
    np.savez(wod_bias_file,
             wodlat=wodlat_arr, wodlon=wodlon_arr,
             z=z,
             wodnit=wodnit, wodnit_f=wodnit_f,
             wodnit_ibim=wodnit_ibim, wodnit_cs1k=wodnit_cs1k, wodnit_niva=wodnit_niva,
             wodnit_woa=wodnit_woa,
             wodnit_diva=wodnit_diva, wodnit_diva_L1=wodnit_diva_L1, wodnit_diva_L2=wodnit_diva_L2,
             wodpho=wodpho, wodpho_f=wodpho_f,
             wodpho_ibim=wodpho_ibim, wodpho_cs1k=wodpho_cs1k, wodpho_niva=wodpho_niva,
             wodpho_woa=wodpho_woa,
             wodpho_diva=wodpho_diva, wodpho_diva_L1=wodpho_diva_L1, wodpho_diva_L2=wodpho_diva_L2,
             wodsil=wodsil, wodsil_f=wodsil_f,
             wodsil_ibim=wodsil_ibim, wodsil_cs1k=wodsil_cs1k, wodsil_niva=wodsil_niva,
             wodsil_woa=wodsil_woa,
             wodsil_diva=wodsil_diva, wodsil_diva_L1=wodsil_diva_L1, wodsil_diva_L2=wodsil_diva_L2,
             wodoxy=wodoxy, wodoxy_f=wodoxy_f,
             wodoxy_ibim=wodoxy_ibim, wodoxy_cs1k=wodoxy_cs1k, wodoxy_niva=wodoxy_niva,
             wodoxy_woa=wodoxy_woa,
             wodoxy_diva=wodoxy_diva, wodoxy_diva_L1=wodoxy_diva_L1, wodoxy_diva_L2=wodoxy_diva_L2,
             latids=latids, lonids=lonids, latidc=latidc, lonidc=lonidc, latidn=latidn, lonidn=lonidn,
             latidw=latidw, lonidw=lonidw, latidd=latidd, lonidd=lonidd,
             land_ibi=land_ibi, land_cs1km=land_cs1km, land_noresm=land_noresm,
             land_diva=land_diva, land_diva_L1=land_diva_L1, land_diva_L2=land_diva_L2,
             wod_day=wod_day_arr, wod_month=wod_month_arr, wod_year=wod_year_arr)

if mi_processing == 1:

    # Winter nutrients from excel
    # insitufile = '/media/dskone/VAL/Insitu_nutrients.xls'
    # insitufile = '/media/dskone/VAL/Insitu_nutrients_to_2021.xls'
    insitufile = '/media/dskone/VAL/Insitu_nutrients_to_2021_tsorted.xls'
    dfis = pd.read_excel(insitufile)
    isyr = np.array(dfis['Year'])
    ismth = np.array(dfis['Month'])
    isday = np.array(dfis['Day'])
    islat = np.array(dfis['Latitude'])
    islon = np.array(dfis['Longitude'])
    isdep = np.array(dfis['Depth'])
    issal = np.array(dfis['salinity'])
    issil = np.array(dfis['silicate'])
    isphs = np.array(dfis['phosphate'])
    isnit = np.array(dfis['nitrate'])
    # itrimidx = np.argwhere((isyr >= Ystart) & (isyr <= Yend) & (islat >= 49) &
    #                        (islat <= 52.95) & (islon >= -10.75) & (islon <= -5.83))
    # itrimidx = np.argwhere((islat >= 49) & (islat <= 52.95) & (islon >= -10.75) & (islon <= -5.83))
    itrimidx = np.argwhere((islat >= 48) & (islat <= 53.95) & (islon >= -11.75) & (islon <= -4.83))
    isyr = isyr[itrimidx][:, 0]
    ismth = ismth[itrimidx][:, 0]
    isday = isday[itrimidx][:, 0]
    islat = islat[itrimidx][:, 0]
    islon = islon[itrimidx][:, 0]
    isdep = isdep[itrimidx][:, 0]
    issal = issal[itrimidx][:, 0]
    issil = issil[itrimidx][:, 0]
    isphs = isphs[itrimidx][:, 0]
    isnit = isnit[itrimidx][:, 0]

    # GLODAPv2
    valdir = '/media/dskone/VAL/'

    EMODdir = '/media/dskfour/VAL/EMODNET/INSITU/'
    nut = 'NPSO'
    dic = 'DIC'
    talk = 'TALK'

    z = isdep
    wodlat = islat
    wodlon = islon

    wodnit = isnit
    wodnit_ibim = np.zeros_like(wodnit)
    wodnit_cs1k = np.zeros_like(wodnit)
    wodnit_niva = np.zeros_like(wodnit)
    wodnit_woa = np.zeros_like(wodnit)
    wodnit_diva = np.zeros_like(wodnit)
    wodnit_diva_L1 = np.zeros_like(wodnit)
    wodnit_diva_L2 = np.zeros_like(wodnit)

    wodpho = isphs
    wodpho_ibim = np.zeros_like(wodpho)
    wodpho_cs1k = np.zeros_like(wodpho)
    wodpho_niva = np.zeros_like(wodpho)
    wodpho_woa = np.zeros_like(wodpho)
    wodpho_diva = np.zeros_like(wodpho)
    wodpho_diva_L1 = np.zeros_like(wodpho)
    wodpho_diva_L2 = np.zeros_like(wodpho)

    wodsil = issil
    wodsil_ibim = np.zeros_like(wodsil)
    wodsil_cs1k = np.zeros_like(wodsil)
    wodsil_niva = np.zeros_like(wodsil)
    wodsil_woa = np.zeros_like(wodsil)
    wodsil_diva = np.zeros_like(wodsil)
    wodsil_diva_L1 = np.zeros_like(wodsil)
    wodsil_diva_L2 = np.zeros_like(wodsil)

    latids = np.zeros_like(wodsil, dtype='int64')
    lonids = np.zeros_like(wodsil, dtype='int64')
    latidc = np.zeros_like(wodsil, dtype='int64')
    lonidc = np.zeros_like(wodsil, dtype='int64')
    latidn = np.zeros_like(wodsil, dtype='int64')
    lonidn = np.zeros_like(wodsil, dtype='int64')
    latidw = np.zeros_like(wodsil, dtype='int64')
    lonidw = np.zeros_like(wodsil, dtype='int64')
    latidd = np.zeros_like(wodsil, dtype='int64')
    lonidd = np.zeros_like(wodsil, dtype='int64')
    land_ibi = np.zeros_like(wodsil, dtype='bool')
    land_cs1km = np.zeros_like(wodsil, dtype='bool')
    land_noresm = np.zeros_like(wodsil, dtype='bool')
    land_woa = np.zeros_like(wodsil, dtype='bool')
    land_diva = np.zeros_like(wodsil, dtype='bool')
    land_diva_L1 = np.zeros_like(wodsil, dtype='bool')
    land_diva_L2 = np.zeros_like(wodsil, dtype='bool')
    wod_day = np.zeros_like(wodsil, dtype='int64')
    wod_month = np.zeros_like(wodsil, dtype='int64')
    wod_year = np.zeros_like(wodsil, dtype='int64')

    # Block for IBI processing
    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    PISCES24_prefixm = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_'
    PISCES24_prefixd = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    PISCES24_ending = '_R*_RE01.nc'
    IBI_BIO_file_sample = sorted(glob.glob(PISCES24files_dir +
                                           PISCES24_prefixm +
                                           '*_*' + PISCES24_ending))
    wdibis = netcdf(IBI_BIO_file_sample[0], 'r')
    wdibso = wdibis.variables
    ibilat = np.array(wdibso['latitude'][:])
    ibilon = np.array(wdibso['longitude'][:])
    ibidep = np.array(wdibso['depth'][:])
    wdibis.close()

    woddays, wod_ind, wod_count = np.unique([isyr, ismth, isday], axis=1, return_index=True, return_counts=True)

    for wd in range(0, woddays.shape[1]):
        wod_day[wd] = woddays[2, wd]
        wod_month[wd] = woddays[1, wd]
        wod_year[wd] = woddays[0, wd]

        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latids = np.zeros_like(prof_ind, dtype='int64')
        lonids = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latids[p] = glor.geo_idx(profs[0, p], ibilat)
            lonids[p] = glor.geo_idx(profs[1, p], ibilon)
            
        if monthlyres == 1:
            IBI_BIO_file_m = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixm + str(wod_year[wd]) +
                                              str(wod_month[wd]).zfill(2) +
                                              '*_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_m[0], 'r')
            
        if dailyres == 1:
            IBI_BIO_file_d = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixd + str(wod_year[wd]) +
                                              str(wod_month[wd]).zfill(2) + str(wod_day[wd]).zfill(2) +
                                              '_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_d[0], 'r')

        wdibmo = wdibim.variables
        # read latitude and longitude for indexing
        ibinit = np.array(wdibmo['no3'][:])
        ibipho = np.array(wdibmo['po4'][:])
        ibisil = np.array(wdibmo['si'][:])
        for dp in range(0, profs.shape[1]):
            ibinit_smp = ibinit[0, :, latids[dp], lonids[dp]]
            ibipho_smp = ibipho[0, :, latids[dp], lonids[dp]]
            ibisil_smp = ibisil[0, :, latids[dp], lonids[dp]]
            if (any(ibinit_smp < 0) == True) & (all(ibinit_smp < 0) == False):
                ibinit_smp[ibinit_smp < 0] = ibinit_smp[ibinit_smp >= 0][-1]
                ibipho_smp[ibipho_smp < 0] = ibipho_smp[ibipho_smp >= 0][-1]
                ibisil_smp[ibisil_smp < 0] = ibisil_smp[ibisil_smp >= 0][-1]
            elif all(ibinit_smp < 0) == True:
                land_ibi[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(ibidep, ibinit_smp, fill_value='extrapolate')
            wodnit_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

            wodint_p = interpolate.interp1d(ibidep, ibipho_smp, fill_value='extrapolate')
            wodpho_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(ibidep, ibisil_smp, fill_value='extrapolate')
            wodsil_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
        wdibim.close()

    # Block for CS1KM processing
    croco_finit = 'croco_avg_Y'
    croco_fend = '.nc'
    CS1KM_BIO_file_sample = sorted(glob.glob(crocodir +
                                             croco_finit +
                                             '*' + croco_fend))
    wdcs1kms = netcdf(CS1KM_BIO_file_sample[0], 'r')
    wdcs1kmso = wdcs1kms.variables
    cs1kmlat = np.array(wdcs1kmso['lat_rho'][:])
    cs1kmlon = np.array(wdcs1kmso['lon_rho'][:])
    wdcs1kms.close()

    for wd in range(0, woddays.shape[1]):

        CS1KM_BIO_file = sorted(glob.glob(crocodir +
                                          croco_finit + str(wod_year[wd]) + 'M' +
                                          str(wod_month[wd]).zfill(2) +
                                          croco_fend))

        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)
        latidc = np.zeros_like(prof_ind, dtype='int64')
        lonidc = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latidc[p] = glor.geo_idx(profs[0, p], cs1kmlat[:, 0])
            lonidc[p] = glor.geo_idx(profs[1, p], cs1kmlon[0, :])

        wdcs1km = netcdf(CS1KM_BIO_file[0], 'r')
        wdcs1kmo = wdcs1km.variables
        
        if monthlyres == 1:
            # Monthly read latitude and longitude for indexing
            cs1kmnit = np.mean(np.array(wdcs1kmo['NO3'][:]), 0)
            cs1kmpho = np.mean(np.array(wdcs1kmo['PO4'][:]), 0)
            cs1kmsil = np.mean(np.array(wdcs1kmo['Si'][:]), 0)
            cs1kmoxy = np.mean(np.array(wdcs1kmo['O2'][:]), 0)
            cs1kmz = np.mean(np.array(wdcs1kmo['zeta'][:]), 0)
            
        if dailyres == 1:
            # Daily read latitude and longitude for indexing
            cs1kmnit = np.array(wdcs1kmo['NO3'][wod_day[wd] - 1, :, :, :])
            cs1kmpho = np.array(wdcs1kmo['PO4'][wod_day[wd] - 1, :, :, :])
            cs1kmsil = np.array(wdcs1kmo['Si'][wod_day[wd] - 1, :, :, :])
            cs1kmoxy = np.array(wdcs1kmo['O2'][wod_day[wd] - 1, :, :, :])
            cs1kmz = np.array(wdcs1kmo['zeta'][wod_day[wd] - 1, :, :])

        # LAT to MSL included, INFOMAR block averaged hmin = 10m
        grdname = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/' + 'CELTIC_grd_v5_MSL_h10.nc'
        ncg = netcdf(grdname, 'r')
        ncgrd = ncg.variables
        h = np.array(ncgrd['h'][:])
        ncg.close()
        N = 20
        theta_s = 7.
        theta_b = 0.
        hc = 20.
        vtransform = 2
        cs1kmd = vgrd.zlevs(h, cs1kmz, theta_s, theta_b, 50, N, 'r', vtransform) * -1

        for dp in range(0, profs.shape[1]):
            cs1kmnit_smp = cs1kmnit[:, latidc[dp], lonidc[dp]]
            cs1kmpho_smp = cs1kmpho[:, latidc[dp], lonidc[dp]]
            cs1kmsil_smp = cs1kmsil[:, latidc[dp], lonidc[dp]]
            # cs1kmoxy_smp = cs1kmoxy[:, latidc[dp], lonidc[dp]]
            cs1kmd_smp = cs1kmd[:, latidc[dp], lonidc[dp]]
            if (any(cs1kmnit_smp == 0) == True) & (all(cs1kmnit_smp == 0) == False):
                cs1kmnit_smp[cs1kmnit_smp == 0] = cs1kmnit_smp[cs1kmnit_smp > 0][-1]
                cs1kmpho_smp[cs1kmpho_smp == 0] = cs1kmpho_smp[cs1kmpho_smp > 0][-1]
                cs1kmsil_smp[cs1kmsil_smp == 0] = cs1kmsil_smp[cs1kmsil_smp > 0][-1]
                # cs1kmoxy_smp[cs1kmoxy_smp == 0] = cs1kmoxy_smp[cs1kmoxy_smp > 0][-1]
            elif all(cs1kmnit_smp < 0) == True:
                land_cs1km[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(cs1kmd_smp, cs1kmnit_smp, fill_value='extrapolate')
            wodnit_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(cs1kmd_smp, cs1kmpho_smp, fill_value='extrapolate')
            wodpho_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(cs1kmd_smp, cs1kmsil_smp, fill_value='extrapolate')
            wodsil_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            # wodint_o = interpolate.interp1d(cs1kmd_smp, cs1kmoxy_smp, fill_value='extrapolate')
            # wodoxy_cs1k[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]] = \
            #     wodint_o(z[z_row_idx_lwr[wod_ind[wd] + dp]:z_row_idx_upr[wod_ind[wd] + dp]])

    # Block for NORESM processing
    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']
    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
    ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
    ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
    ncnornio = ncnorni.variables
    ncnorpoo = ncnorpo.variables
    ncnorsio = ncnorsi.variables
    # ncnoro2o = ncnoro2.variables
    noresmlat = np.array(ncnornio['lat'][:])
    noresmlon = np.array(ncnornio['lon'][:])
    depthniva = np.array(ncnornio['depth'][:])
    time_niva = np.array(ncnornio['time'][:])

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)
        latidn = np.zeros_like(prof_ind, dtype='int64')
        lonidn = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latidn[p] = glor.geo_idx(profs[0, p], noresmlat)
            lonidn[p] = glor.geo_idx(profs[1, p], noresmlon)

        latidn = latidn - (noresmlat.shape[1] * np.floor(latidn / noresmlat.shape[1]))
        lonidn = lonidn - (noresmlat.shape[1] * np.floor(lonidn / noresmlat.shape[1]))
        latidn = latidn.astype(int)
        lonidn = lonidn.astype(int)

        noresmnit = np.array(ncnornio[NIVAvars[0]][:]) * 1000
        noresmpho = np.array(ncnorpoo[NIVAvars[1]][:]) * 1000
        noresmsil = np.array(ncnorsio[NIVAvars[2]][:]) * 1000

        for dp in range(0, profs.shape[1]):
            Tinin = datetime(wod_year[wd], wod_month[wd], 15)
            if wod_year[wd] <= 2020:
                if monthlyres == 1:
                    Tininx1 = d2i(Tinin, ncnornio['time'], select='exact')
                    noresmnit_smp = noresmnit[Tininx1, :, latidn[dp], lonidn[dp]]
                    noresmpho_smp = noresmpho[Tininx1, :, latidn[dp], lonidn[dp]]
                    noresmsil_smp = noresmsil[Tininx1, :, latidn[dp], lonidn[dp]]
                if dailyres == 1:

                    if wod_day[wd] != 15:
                        Tinin = datetime(wod_year[wd], wod_month[wd], wod_day[wd])
                        Tininx1 = d2i(Tinin, ncnornio['time'], select='before')
                        Tininx2 = d2i(Tinin, ncnornio['time'], select='after')
                        if wod_day[wd] > 15:
                            tmult = wod_day[wd] - 15
                        else:
                            tmult = wod_day[wd] + 15
                        # Tidx1
                        noresmnit_1 = noresmnit[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmpho_1 = noresmpho[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmsil_1 = noresmsil[Tininx1, :, latidn[dp], lonidn[dp]]
                        # Tidx2
                        noresmnit_2 = noresmnit[Tininx2, :, latidn[dp], lonidn[dp]]
                        noresmpho_2 = noresmpho[Tininx2, :, latidn[dp], lonidn[dp]]
                        noresmsil_2 = noresmsil[Tininx2, :, latidn[dp], lonidn[dp]]
                        # Interpolating to the day
                        noresmnit_smp = noresmnit_1 + (((noresmnit_2 - noresmnit_1)/30)*tmult)
                        noresmpho_smp = noresmpho_1 + (((noresmpho_2 - noresmpho_1)/30)*tmult)
                        noresmsil_smp = noresmsil_1 + (((noresmsil_2 - noresmsil_1)/30)*tmult)
                    else:
                        Tinin = datetime(wod_year[wd], wod_month[wd], 15)
                        Tininx1 = d2i(Tinin, ncnornio['time'], select='exact')
                        noresmnit_smp = noresmnit[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmpho_smp = noresmpho[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmsil_smp = noresmsil[Tininx1, :, latidn[dp], lonidn[dp]]
                    
                if (any(noresmnit_smp > 10000) == True) & (all(noresmnit_smp > 10000) == False):
                    noresmnit_smp[noresmnit_smp > 10000] = noresmnit_smp[noresmnit_smp < 10000][-1]
                    noresmpho_smp[noresmpho_smp > 10000] = noresmpho_smp[noresmpho_smp < 10000][-1]
                    noresmsil_smp[noresmsil_smp > 10000] = noresmsil_smp[noresmsil_smp < 10000][-1]
                elif all(noresmnit_smp > 10000) == True:
                    land_noresm[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                #
                wodint_n = interpolate.interp1d(depthniva, noresmnit_smp, fill_value='extrapolate')
                wodnit_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                #
                wodint_p = interpolate.interp1d(depthniva, noresmpho_smp, fill_value='extrapolate')
                wodpho_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                #
                wodint_s = interpolate.interp1d(depthniva, noresmsil_smp, fill_value='extrapolate')
                wodsil_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    ncnorni.close()
    ncnorpo.close()
    ncnorsi.close()

    # Block for WOA processing
    valdir = '/media/dskone/VAL/'
    woadir = valdir + 'WOA18/'

    # Common woa elements for monthly climatology
    woa_exn = 'woa18_all_n01_01.nc'
    woa_exo = 'woa18_all_o01_01.nc'
    woan_01 = netcdf(woadir + woa_exn, 'r')
    woao_01 = netcdf(woadir + woa_exo, 'r')
    woan_01o = woan_01.variables
    woao_01o = woao_01.variables
    woa_lat = np.array(woan_01o['lat'][:])
    woa_lon = np.array(woan_01o['lon'][:])
    woa_depth_nps = np.array(woan_01o['depth'][:])
    woa_depth_nps_bnds = np.array(woan_01o['depth_bnds'][:])
    n_01a = np.array(woan_01o['n_an'][:])
    o_01a = np.array(woao_01o['o_an'][:])
    woan_01.close()
    woao_01.close()

    woa_sta = 'woa18_all_'
    woa_end = '_01.nc'

    woa_n = np.tile(n_01a, (12, 1, 1, 1))
    # woa_o = np.tile(o_01a, (12, 1, 1, 1))
    woa_n = np.zeros_like(woa_n)
    woa_p = np.zeros_like(woa_n)
    woa_s = np.zeros_like(woa_n)
    # woa_o = np.zeros_like(woa_o)

    for mt in range(0, 12):
        # nitrate
        woa_n_mth = netcdf(woadir + woa_sta + 'n' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_n_mth_o = woa_n_mth.variables
        woa_n[mt, :, :, :] = np.array(woa_n_mth_o['n_an'][:])
        woa_depth_nps = np.array(woa_n_mth_o['depth'][:])
        woa_depth_nps_bnds = np.array(woa_n_mth_o['depth_bnds'][:])
        woa_n_mth.close()

        # phosphate
        woa_p_mth = netcdf(woadir + woa_sta + 'p' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_p_mth_o = woa_p_mth.variables
        woa_p[mt, :, :, :] = np.array(woa_p_mth_o['p_an'][:])
        woa_p_mth.close()

        # silicate
        woa_s_mth = netcdf(woadir + woa_sta + 'i' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_s_mth_o = woa_s_mth.variables
        woa_s[mt, :, :, :] = np.array(woa_s_mth_o['i_an'][:])
        woa_s_mth.close()

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidw = np.zeros_like(prof_ind, dtype='int64')
        lonidw = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latidw[p] = glor.geo_idx(profs[0, p], woa_lat)
            lonidw[p] = glor.geo_idx(profs[1, p], woa_lon)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                woanit_smp = woa_n[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                woapho_smp = woa_p[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                woasil_smp = woa_s[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
            if dailyres == 1:
                # Day
                if wod_day[wd] > 15:
                    tmult = wod_day[wd] - 15
                    tidx1 = wod_month[wd] - 1
                    if wod_month[wd] == 12:
                        tidx2 = 0
                    else:
                        tidx2 = wod_month[wd]

                    # Tidx1
                    woanit_1 = woa_n[tidx1, :, latidw[dp], lonidw[dp]]
                    woapho_1 = woa_p[tidx1, :, latidw[dp], lonidw[dp]]
                    woasil_1 = woa_s[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    woanit_2 = woa_n[tidx2, :, latidw[dp], lonidw[dp]]
                    woapho_2 = woa_p[tidx2, :, latidw[dp], lonidw[dp]]
                    woasil_2 = woa_s[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    woanit_smp = woanit_1 + (((woanit_2 - woanit_1)/30)*tmult)
                    woapho_smp = woapho_1 + (((woapho_2 - woapho_1)/30)*tmult)
                    woasil_smp = woasil_1 + (((woasil_2 - woasil_1)/30)*tmult)

                elif wod_day[wd] < 15:
                    tmult = wod_day[wd] + 15
                    tidx2 = wod_month[wd] - 1
                    if wod_month[wd] == 1:
                        tidx1 = 11
                    else:
                        tidx1 = wod_month[wd] - 2

                    # Tidx1
                    woanit_1 = woa_n[tidx1, :, latidw[dp], lonidw[dp]]
                    woapho_1 = woa_p[tidx1, :, latidw[dp], lonidw[dp]]
                    woasil_1 = woa_s[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    woanit_2 = woa_n[tidx2, :, latidw[dp], lonidw[dp]]
                    woapho_2 = woa_p[tidx2, :, latidw[dp], lonidw[dp]]
                    woasil_2 = woa_s[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    woanit_smp = woanit_1 + (((woanit_2 - woanit_1)/30)*tmult)
                    woapho_smp = woapho_1 + (((woapho_2 - woapho_1)/30)*tmult)
                    woasil_smp = woasil_1 + (((woasil_2 - woasil_1)/30)*tmult)

                else:
                    woanit_smp = woa_n[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    woapho_smp = woa_p[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    woasil_smp = woa_s[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]

            if (any(woanit_smp > 10000) == True) & (all(woanit_smp > 10000) == False):
                woanit_smp[woanit_smp > 10000] = woanit_smp[woanit_smp < 10000][-1]
                woapho_smp[woapho_smp > 10000] = woapho_smp[woapho_smp < 10000][-1]
                woasil_smp[woasil_smp > 10000] = woasil_smp[woasil_smp < 10000][-1]
            elif all(woanit_smp > 10000) == True:
                land_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(woa_depth_nps, woanit_smp, fill_value='extrapolate')
            wodnit_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(woa_depth_nps, woapho_smp, fill_value='extrapolate')
            wodpho_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(woa_depth_nps, woasil_smp, fill_value='extrapolate')
            wodsil_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for EMODNET DIVA processing

    emcdir = '/media/dsktwo/VAL/EMODNET/CLIMATOLOGY/'

    # Common emc elements for monthly climatology
    emc_ex = 'Water_body_dissolved_inorganic_nitrogen_v2021.nc'
    emc_01 = netcdf(emcdir + emc_ex, 'r')
    emc_01o = emc_01.variables
    emc_lat = np.array(emc_01o['lat'][:])
    emc_lon = np.array(emc_01o['lon'][:])
    emc_depth = np.array(emc_01o['depth'][:])
    n_01a = np.array(emc_01o['Water body dissolved inorganic nitrogen (DIN)'][:])
    emc_01.close()

    emc_sta = 'Water_body_'
    emc_end = '_v2021.nc'

    # nitrate
    emc_n_mth = netcdf(emcdir + emc_sta + 'dissolved_inorganic_nitrogen' + emc_end, 'r')
    emc_n_mth_o = emc_n_mth.variables
    emc_n = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)'][:])
    emc_n_L1 = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)_L1'][:])
    emc_n_L2 = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)_L2'][:])
    emc_n_mth.close()

    # phosphate
    emc_p_mth = netcdf(emcdir + emc_sta + 'phosphate' + emc_end, 'r')
    emc_p_mth_o = emc_p_mth.variables
    emc_p = np.array(emc_p_mth_o['Water body phosphate'][:])
    emc_p_L1 = np.array(emc_p_mth_o['Water body phosphate_L1'][:])
    emc_p_L2 = np.array(emc_p_mth_o['Water body phosphate_L2'][:])
    emc_p_mth.close()

    # silicate
    emc_s_mth = netcdf(emcdir + emc_sta + 'silicate' + emc_end, 'r')
    emc_s_mth_o = emc_s_mth.variables
    emc_s = np.array(emc_s_mth_o['Water body silicate'][:])
    emc_s_L1 = np.array(emc_s_mth_o['Water body silicate_L1'][:])
    emc_s_L2 = np.array(emc_s_mth_o['Water body silicate_L2'][:])
    emc_s_mth.close()

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidd = np.zeros_like(prof_ind, dtype='int64')
        lonidd = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latidd[p] = glor.geo_idx(profs[0, p], emc_lat)
            lonidd[p] = glor.geo_idx(profs[1, p], emc_lon)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                # Monthly
                emcnit_smp = emc_n[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcnit_smp_L1 = emc_n_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcnit_smp_L2 = emc_n_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcpho_smp = emc_p[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcpho_smp_L1 = emc_p_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcpho_smp_L2 = emc_p_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcsil_smp = emc_s[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcsil_smp_L1 = emc_s_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcsil_smp_L2 = emc_s_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                
            if dailyres == 1:
                # Daily
                if wod_day[wd] > 15:
                    tmult = wod_day[wd] - 15
                    tidx1 = wod_month[wd] - 1
                    if wod_month[wd] == 12:
                        tidx2 = 0
                    else:
                        tidx2 = wod_month[wd]
                    # Tidx1
                    emcnit_1 = emc_n[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1 = emc_p[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1 = emc_s[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L1 = emc_n_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L1 = emc_p_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L1 = emc_s_L1[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L2 = emc_n_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L2 = emc_p_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L2 = emc_s_L2[tidx1, :, latidd[dp], lonidd[dp]]

                    # Tidx2
                    emcnit_2 = emc_n[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2 = emc_p[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2 = emc_s[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L1 = emc_n_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L1 = emc_p_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L1 = emc_s_L1[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L2 = emc_n_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L2 = emc_p_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L2 = emc_s_L2[tidx2, :, latidd[dp], lonidd[dp]]

                    # Interpolating to the day
                    emcnit_smp = emcnit_1 + (((emcnit_2 - emcnit_1)/30)*tmult)
                    emcpho_smp = emcpho_1 + (((emcpho_2 - emcpho_1)/30)*tmult)
                    emcsil_smp = emcsil_1 + (((emcsil_2 - emcsil_1)/30)*tmult)

                    emcnit_smp_L1 = emcnit_1_L1 + (((emcnit_2_L1 - emcnit_1_L1)/30)*tmult)
                    emcpho_smp_L1 = emcpho_1_L1 + (((emcpho_2_L1 - emcpho_1_L1)/30)*tmult)
                    emcsil_smp_L1 = emcsil_1_L1 + (((emcsil_2_L1 - emcsil_1_L1)/30)*tmult)

                    emcnit_smp_L2 = emcnit_1_L2 + (((emcnit_2_L2 - emcnit_1_L2)/30)*tmult)
                    emcpho_smp_L2 = emcpho_1_L2 + (((emcpho_2_L2 - emcpho_1_L2)/30)*tmult)
                    emcsil_smp_L2 = emcsil_1_L2 + (((emcsil_2_L2 - emcsil_1_L2)/30)*tmult)

                elif wod_day[wd] < 15:
                    tmult = wod_day[wd] + 15
                    tidx2 = wod_month[wd] - 1
                    if wod_month[wd] == 1:
                        tidx1 = 11
                    else:
                        tidx1 = wod_month[wd] - 2

                    # Tidx1
                    emcnit_1 = emc_n[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1 = emc_p[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1 = emc_s[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L1 = emc_n_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L1 = emc_p_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L1 = emc_s_L1[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L2 = emc_n_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L2 = emc_p_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L2 = emc_s_L2[tidx1, :, latidd[dp], lonidd[dp]]

                    # Tidx2
                    emcnit_2 = emc_n[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2 = emc_p[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2 = emc_s[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L1 = emc_n_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L1 = emc_p_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L1 = emc_s_L1[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L2 = emc_n_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L2 = emc_p_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L2 = emc_s_L2[tidx2, :, latidd[dp], lonidd[dp]]

                    # Interpolating to the day
                    emcnit_smp = emcnit_1 + (((emcnit_2 - emcnit_1)/30)*tmult)
                    emcpho_smp = emcpho_1 + (((emcpho_2 - emcpho_1)/30)*tmult)
                    emcsil_smp = emcsil_1 + (((emcsil_2 - emcsil_1)/30)*tmult)

                    emcnit_smp_L1 = emcnit_1_L1 + (((emcnit_2_L1 - emcnit_1_L1)/30)*tmult)
                    emcpho_smp_L1 = emcpho_1_L1 + (((emcpho_2_L1 - emcpho_1_L1)/30)*tmult)
                    emcsil_smp_L1 = emcsil_1_L1 + (((emcsil_2_L1 - emcsil_1_L1)/30)*tmult)

                    emcnit_smp_L2 = emcnit_1_L2 + (((emcnit_2_L2 - emcnit_1_L2)/30)*tmult)
                    emcpho_smp_L2 = emcpho_1_L2 + (((emcpho_2_L2 - emcpho_1_L2)/30)*tmult)
                    emcsil_smp_L2 = emcsil_1_L2 + (((emcsil_2_L2 - emcsil_1_L2)/30)*tmult)

                else:
                    emcnit_smp = emc_n[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                    emcpho_smp = emc_p[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                    emcsil_smp = emc_s[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]

            if (any(emcnit_smp > 10000) == True) & (all(emcnit_smp > 10000) == False):
                emcnit_smp[emcnit_smp > 10000] = emcnit_smp[emcnit_smp < 10000][-1]
                emcpho_smp[emcpho_smp > 10000] = emcpho_smp[emcpho_smp < 10000][-1]
                emcsil_smp[emcsil_smp > 10000] = emcsil_smp[emcsil_smp < 10000][-1]
            elif all(emcnit_smp > 10000) == True:
                land_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True

            if (any(emcnit_smp_L1 > 10000) == True) & (all(emcnit_smp_L1 > 10000) == False):
                emcnit_smp_L1[emcnit_smp_L1 > 10000] = emcnit_smp_L1[emcnit_smp_L1 < 10000][-1]
            elif all(emcnit_smp_L1 > 10000) == True:
                land_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True

            if (any(emcnit_smp_L2 > 10000) == True) & (all(emcnit_smp_L2 > 10000) == False):
                emcnit_smp_L2[emcnit_smp_L2 > 10000] = emcnit_smp_L2[emcnit_smp_L2 < 10000][-1]
            elif all(emcnit_smp_L2 > 10000) == True:
                land_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True

            if (any(emcpho_smp_L1 > 10000) == True) & (all(emcpho_smp_L1 > 10000) == False):
                emcpho_smp_L1[emcpho_smp_L1 > 10000] = emcpho_smp_L1[emcpho_smp_L1 < 10000][-1]
            if (any(emcpho_smp_L2 > 10000) == True) & (all(emcpho_smp_L2 > 10000) == False):
                emcpho_smp_L2[emcpho_smp_L2 > 10000] = emcpho_smp_L2[emcpho_smp_L2 < 10000][-1]

            if (any(emcsil_smp_L1 > 10000) == True) & (all(emcsil_smp_L1 > 10000) == False):
                emcsil_smp_L1[emcsil_smp_L1 > 10000] = emcsil_smp_L1[emcsil_smp_L1 < 10000][-1]
            if (any(emcsil_smp_L2 > 10000) == True) & (all(emcsil_smp_L2 > 10000) == False):
                emcsil_smp_L2[emcsil_smp_L2 > 10000] = emcsil_smp_L2[emcsil_smp_L2 < 10000][-1]
            #
            wodint_n = interpolate.interp1d(emc_depth, emcnit_smp, fill_value='extrapolate')
            wodnit_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_n_L1 = interpolate.interp1d(emc_depth, emcnit_smp_L1, fill_value='extrapolate')
            wodnit_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_n_L2 = interpolate.interp1d(emc_depth, emcnit_smp_L2, fill_value='extrapolate')
            wodnit_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(emc_depth, emcpho_smp, fill_value='extrapolate')
            wodpho_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_p_L1 = interpolate.interp1d(emc_depth, emcpho_smp_L1, fill_value='extrapolate')
            wodpho_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_p_L2 = interpolate.interp1d(emc_depth, emcpho_smp_L2, fill_value='extrapolate')
            wodpho_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(emc_depth, emcsil_smp, fill_value='extrapolate')
            wodsil_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_s_L1 = interpolate.interp1d(emc_depth, emcsil_smp_L1, fill_value='extrapolate')
            wodsil_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_s_L2 = interpolate.interp1d(emc_depth, emcsil_smp_L2, fill_value='extrapolate')
            wodsil_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    mi_bias_file = crocodir + 'MI_bias_' + npzstr + '.npz'
    np.savez(mi_bias_file,
             wodlat=wodlat, wodlon=wodlon, z=z,
             wodnit=wodnit,
             wodnit_ibim=wodnit_ibim, wodnit_cs1k=wodnit_cs1k, wodnit_niva=wodnit_niva,
             wodnit_woa=wodnit_woa,
             wodnit_diva=wodnit_diva, wodnit_diva_L1=wodnit_diva_L1, wodnit_diva_L2=wodnit_diva_L2,
             wodpho=wodpho,
             wodpho_ibim=wodpho_ibim, wodpho_cs1k=wodpho_cs1k, wodpho_niva=wodpho_niva,
             wodpho_woa=wodpho_woa,
             wodpho_diva=wodpho_diva, wodpho_diva_L1=wodpho_diva_L1, wodpho_diva_L2=wodpho_diva_L2,
             wodsil=wodsil,
             wodsil_ibim=wodsil_ibim, wodsil_cs1k=wodsil_cs1k, wodsil_niva=wodsil_niva,
             wodsil_woa=wodsil_woa,
             wodsil_diva=wodsil_diva, wodsil_diva_L1=wodsil_diva_L1, wodsil_diva_L2=wodsil_diva_L2,
             latids=latids, lonids=lonids, latidc=latidc, lonidc=lonidc, latidn=latidn, lonidn=lonidn,
             latidw=latidw, lonidw=lonidw, latidd=latidd, lonidd=lonidd,
             land_ibi=land_ibi, land_cs1km=land_cs1km, land_noresm=land_noresm,
             land_diva=land_diva, land_diva_L1=land_diva_L1, land_diva_L2=land_diva_L2,
             wod_day=isday, wod_month=ismth, wod_year=isyr)

if glodap_processing == 1:

    # GLODAPv2
    valdir = '/media/dskone/VAL/'
    glodapfil = valdir + 'GLODAPv2.2022_Merged_Master_File.mat'
    mat = scio.loadmat(glodapfil)
    mat2 = pd.DataFrame({'G2year': pd.to_numeric(mat['G2year'][:, 0]),
                         'G2month': pd.to_numeric(mat['G2month'][:, 0]),
                         'G2day': pd.to_numeric(mat['G2day'][:, 0]),
                         'G2latitude': pd.to_numeric(mat['G2latitude'][:, 0]),
                         'G2longitude': pd.to_numeric(mat['G2longitude'][:, 0]),
                         'G2depth': pd.to_numeric(mat['G2depth'][:, 0]),
                         'G2nitrate': pd.to_numeric(mat['G2nitrate'][:, 0]),
                         'G2phosphate': pd.to_numeric(mat['G2phosphate'][:, 0]),
                         'G2oxygen': pd.to_numeric(mat['G2oxygen'][:, 0]),
                         'G2silicate': pd.to_numeric(mat['G2silicate'][:, 0]),
                         'G2fco2': pd.to_numeric(mat['G2fco2'][:, 0]),
                         'G2tco2': pd.to_numeric(mat['G2tco2'][:, 0]),
                         'G2talk': pd.to_numeric(mat['G2talk'][:, 0])})
    mat2.sort_values(by=['G2year', 'G2month', 'G2day', 'G2latitude', 'G2longitude', 'G2depth'],
                     ascending=True, inplace=True)
    gyear = mat2['G2year'].to_numpy()
    gmonth = mat2['G2month'].to_numpy()
    gday = mat2['G2day'].to_numpy()
    glat = mat2['G2latitude'].to_numpy()
    glon = mat2['G2longitude'].to_numpy()
    gdepth = mat2['G2depth'].to_numpy()
    gnitr = mat2['G2nitrate'].to_numpy()
    gphos = mat2['G2phosphate'].to_numpy()
    goxy = mat2['G2oxygen'].to_numpy()
    gsil = mat2['G2silicate'].to_numpy()
    gspco2 = mat2['G2fco2'].to_numpy()
    gdic = mat2['G2tco2'].to_numpy()
    gtalk = mat2['G2talk'].to_numpy()

    trimidx = np.argwhere((glat >= 48) & (glat <= 53.95) & (glon >= -11.75) & (glon <= -4.83))
    gyear = gyear[trimidx][:, 0]
    gmonth = gmonth[trimidx][:, 0]
    gday = gday[trimidx][:, 0]
    glat = glat[trimidx][:, 0]
    glon = glon[trimidx][:, 0]
    gdepth = gdepth[trimidx][:, 0]
    gnitr = gnitr[trimidx][:, 0]
    gphos = gphos[trimidx][:, 0]
    goxy = goxy[trimidx][:, 0]
    gsil = gsil[trimidx][:, 0]
    gspco2 = gspco2[trimidx][:, 0]
    gdic = gdic[trimidx][:, 0]
    gtalk = gtalk[trimidx][:, 0]

    EMODdir = '/media/dskfour/VAL/EMODNET/INSITU/'
    nut = 'NPSO'
    dic = 'DIC'
    talk = 'TALK'

    # WODdir = '/media/dskone/VAL/WOD/'
    # WODnut = 'WOD_NO3_PO4_Si_O2_OSD_ragged.nc'
    #
    # ncwd = netcdf(WODdir + WODnut, 'r')
    # ncwdo = ncwd.variables
    z = gdepth
    wodlat = glat
    wodlon = glon

    wodnit = gnitr
    wodnit_ibim = np.zeros_like(wodnit)
    wodnit_cs1k = np.zeros_like(wodnit)
    wodnit_niva = np.zeros_like(wodnit)
    wodnit_woa = np.zeros_like(wodnit)
    wodnit_diva = np.zeros_like(wodnit)
    wodnit_diva_L1 = np.zeros_like(wodnit)
    wodnit_diva_L2 = np.zeros_like(wodnit)

    wodpho = gphos
    wodpho_ibim = np.zeros_like(wodpho)
    wodpho_cs1k = np.zeros_like(wodpho)
    wodpho_niva = np.zeros_like(wodpho)
    wodpho_woa = np.zeros_like(wodpho)
    wodpho_diva = np.zeros_like(wodpho)
    wodpho_diva_L1 = np.zeros_like(wodpho)
    wodpho_diva_L2 = np.zeros_like(wodpho)

    wodsil = gsil
    wodsil_ibim = np.zeros_like(wodsil)
    wodsil_cs1k = np.zeros_like(wodsil)
    wodsil_niva = np.zeros_like(wodsil)
    wodsil_woa = np.zeros_like(wodsil)
    wodsil_diva = np.zeros_like(wodsil)
    wodsil_diva_L1 = np.zeros_like(wodsil)
    wodsil_diva_L2 = np.zeros_like(wodsil)

    wodoxy = goxy
    wodoxy_ibim = np.zeros_like(wodoxy)
    wodoxy_cs1k = np.zeros_like(wodoxy)
    wodoxy_niva = np.zeros_like(wodoxy)
    wodoxy_woa = np.zeros_like(wodoxy)
    wodoxy_diva = np.zeros_like(wodoxy)
    wodoxy_diva_L1 = np.zeros_like(wodoxy)
    wodoxy_diva_L2 = np.zeros_like(wodoxy)

    woddays, wod_ind, wod_count = np.unique([gyear, gmonth, gday], axis=1, return_index=True, return_counts=True)

    latids = np.zeros_like(wodoxy, dtype='int64')
    lonids = np.zeros_like(wodoxy, dtype='int64')
    latidc = np.zeros_like(wodoxy, dtype='int64')
    lonidc = np.zeros_like(wodoxy, dtype='int64')
    latidn = np.zeros_like(wodoxy, dtype='int64')
    lonidn = np.zeros_like(wodoxy, dtype='int64')
    latidw = np.zeros_like(wodoxy, dtype='int64')
    lonidw = np.zeros_like(wodoxy, dtype='int64')
    latidd = np.zeros_like(wodoxy, dtype='int64')
    lonidd = np.zeros_like(wodoxy, dtype='int64')
    land_ibi = np.zeros_like(wodoxy, dtype='bool')
    land_cs1km = np.zeros_like(wodoxy, dtype='bool')
    land_noresm = np.zeros_like(wodoxy, dtype='bool')
    land_woa = np.zeros_like(wodoxy, dtype='bool')
    land_diva = np.zeros_like(wodoxy, dtype='bool')
    land_diva_L1 = np.zeros_like(wodoxy, dtype='bool')
    land_diva_L2 = np.zeros_like(wodoxy, dtype='bool')
    wod_day = np.zeros_like(wodoxy, dtype='int64')
    wod_month = np.zeros_like(wodoxy, dtype='int64')
    wod_year = np.zeros_like(wodoxy, dtype='int64')

    # Block for IBI processing
    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    PISCES24_prefixm = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_'
    PISCES24_prefixd = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    PISCES24_ending = '_R*_RE01.nc'
    IBI_BIO_file_sample = sorted(glob.glob(PISCES24files_dir +
                                           PISCES24_prefixm +
                                           '*_*' + PISCES24_ending))
    wdibis = netcdf(IBI_BIO_file_sample[0], 'r')
    wdibso = wdibis.variables
    ibilat = np.array(wdibso['latitude'][:])
    ibilon = np.array(wdibso['longitude'][:])
    ibidep = np.array(wdibso['depth'][:])
    wdibis.close()

    for wd in range(0, woddays.shape[1]):
        wod_day[wd] = woddays[2, wd]
        wod_month[wd] = woddays[1, wd]
        wod_year[wd] = woddays[0, wd]

        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latids = np.zeros_like(prof_ind, dtype='int64')
        lonids = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latids[p] = glor.geo_idx(profs[0, p], ibilat)
            lonids[p] = glor.geo_idx(profs[1, p], ibilon)
        if monthlyres == 1:
            IBI_BIO_file_m = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixm + str(wod_year[wd]) +
                                              str(wod_month[wd]).zfill(2) +
                                              '*_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_m[0], 'r')
        if dailyres == 1:
            IBI_BIO_file_d = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixd + str(wod_year[wd]) +
                                              str(wod_month[wd]).zfill(2) + str(wod_day[wd]).zfill(2) +
                                              '_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_d[0], 'r')
        wdibmo = wdibim.variables
        # read latitude and longitude for indexing
        ibinit = np.array(wdibmo['no3'][:])
        ibipho = np.array(wdibmo['po4'][:])
        ibisil = np.array(wdibmo['si'][:])
        ibioxy = np.array(wdibmo['o2'][:])
        for dp in range(0, profs.shape[1]):
            ibinit_smp = ibinit[0, :, latids[dp], lonids[dp]]
            ibipho_smp = ibipho[0, :, latids[dp], lonids[dp]]
            ibisil_smp = ibisil[0, :, latids[dp], lonids[dp]]
            ibioxy_smp = ibioxy[0, :, latids[dp], lonids[dp]]
            if (any(ibinit_smp < 0) == True) & (all(ibinit_smp < 0) == False):
                ibinit_smp[ibinit_smp < 0] = ibinit_smp[ibinit_smp >= 0][-1]
                ibipho_smp[ibipho_smp < 0] = ibipho_smp[ibipho_smp >= 0][-1]
                ibisil_smp[ibisil_smp < 0] = ibisil_smp[ibisil_smp >= 0][-1]
                ibioxy_smp[ibioxy_smp < 0] = ibioxy_smp[ibioxy_smp >= 0][-1]
            elif all(ibinit_smp < 0) == True:
                land_ibi[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(ibidep, ibinit_smp, fill_value='extrapolate')
            wodnit_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(ibidep, ibipho_smp, fill_value='extrapolate')
            wodpho_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(ibidep, ibisil_smp, fill_value='extrapolate')
            wodsil_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_o = interpolate.interp1d(ibidep, ibioxy_smp, fill_value='extrapolate')
            wodoxy_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for CS1KM processing
    croco_finit = 'croco_avg_Y'
    croco_fend = '.nc'
    CS1KM_BIO_file_sample = sorted(glob.glob(crocodir +
                                             croco_finit +
                                             '*' + croco_fend))
    wdcs1kms = netcdf(CS1KM_BIO_file_sample[0], 'r')
    wdcs1kmso = wdcs1kms.variables
    cs1kmlat = np.array(wdcs1kmso['lat_rho'][:])
    cs1kmlon = np.array(wdcs1kmso['lon_rho'][:])
    wdcs1kms.close()

    for wd in range(0, woddays.shape[1]):

        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidc = np.zeros_like(prof_ind, dtype='int64')
        lonidc = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latidc[p] = glor.geo_idx(profs[0, p], cs1kmlat[:, 0])
            lonidc[p] = glor.geo_idx(profs[1, p], cs1kmlon[0, :])

        CS1KM_BIO_file = sorted(glob.glob(crocodir +
                                          croco_finit + str(wod_year[wd]) + 'M' +
                                          str(wod_month[wd]).zfill(2) +
                                          croco_fend))
        wdcs1km = netcdf(CS1KM_BIO_file[0], 'r')
        wdcs1kmo = wdcs1km.variables
        if monthlyres == 1:
            # Monthly read latitude and longitude for indexing
            cs1kmnit = np.mean(np.array(wdcs1kmo['NO3'][:]), 0)
            cs1kmpho = np.mean(np.array(wdcs1kmo['PO4'][:]), 0)
            cs1kmsil = np.mean(np.array(wdcs1kmo['Si'][:]), 0)
            cs1kmoxy = np.mean(np.array(wdcs1kmo['O2'][:]), 0)
            cs1kmz = np.mean(np.array(wdcs1kmo['zeta'][:]), 0)
            
        if dailyres == 1:
            # Daily read latitude and longitude for indexing
            cs1kmnit = np.array(wdcs1kmo['NO3'][wod_day[wd] - 1, :, :, :])
            cs1kmpho = np.array(wdcs1kmo['PO4'][wod_day[wd] - 1, :, :, :])
            cs1kmsil = np.array(wdcs1kmo['Si'][wod_day[wd] - 1, :, :, :])
            cs1kmoxy = np.array(wdcs1kmo['O2'][wod_day[wd] - 1, :, :, :])
            cs1kmz = np.array(wdcs1kmo['zeta'][wod_day[wd] - 1, :, :])

        # LAT to MSL included, INFOMAR block averaged hmin = 10m
        grdname = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/' + 'CELTIC_grd_v5_MSL_h10.nc'
        ncg = netcdf(grdname, 'r')
        ncgrd = ncg.variables
        h = np.array(ncgrd['h'][:])
        ncg.close()
        N = 20
        theta_s = 7.
        theta_b = 0.
        hc = 20.
        vtransform = 2
        cs1kmd = vgrd.zlevs(h, cs1kmz, theta_s, theta_b, 50, N, 'r', vtransform) * -1
        for dp in range(0, profs.shape[1]):
            cs1kmnit_smp = cs1kmnit[:, latidc[dp], lonidc[dp]]
            cs1kmpho_smp = cs1kmpho[:, latidc[dp], lonidc[dp]]
            cs1kmsil_smp = cs1kmsil[:, latidc[dp], lonidc[dp]]
            cs1kmoxy_smp = cs1kmoxy[:, latidc[dp], lonidc[dp]]
            cs1kmd_smp = cs1kmd[:, latidc[dp], lonidc[dp]]
            if (any(cs1kmnit_smp == 0) == True) & (all(cs1kmnit_smp == 0) == False):
                cs1kmnit_smp[cs1kmnit_smp == 0] = cs1kmnit_smp[cs1kmnit_smp > 0][-1]
                cs1kmpho_smp[cs1kmpho_smp == 0] = cs1kmpho_smp[cs1kmpho_smp > 0][-1]
                cs1kmsil_smp[cs1kmsil_smp == 0] = cs1kmsil_smp[cs1kmsil_smp > 0][-1]
                cs1kmoxy_smp[cs1kmoxy_smp == 0] = cs1kmoxy_smp[cs1kmoxy_smp > 0][-1]
            elif all(cs1kmnit_smp < 0) == True:
                land_cs1km[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(cs1kmd_smp, cs1kmnit_smp, fill_value='extrapolate')
            wodnit_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(cs1kmd_smp, cs1kmpho_smp, fill_value='extrapolate')
            wodpho_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(cs1kmd_smp, cs1kmsil_smp, fill_value='extrapolate')
            wodsil_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_o = interpolate.interp1d(cs1kmd_smp, cs1kmoxy_smp, fill_value='extrapolate')
            wodoxy_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for NORESM processing
    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']
    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
    ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
    ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
    ncnoro2 = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[3] + NORESM_ending, 'r')
    ncnornio = ncnorni.variables
    ncnorpoo = ncnorpo.variables
    ncnorsio = ncnorsi.variables
    ncnoro2o = ncnoro2.variables
    noresmlat = np.array(ncnornio['lat'][:])
    noresmlon = np.array(ncnornio['lon'][:])
    depthniva = np.array(ncnornio['depth'][:])
    time_niva = np.array(ncnornio['time'][:])

    for wd in range(0, woddays.shape[1]):

        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidn = np.zeros_like(prof_ind, dtype='int64')
        lonidn = np.zeros_like(prof_ind, dtype='int64')
        for p in range(0, profs.shape[1]):
            latidn[p] = glor.geo_idx(profs[0, p], noresmlat)
            lonidn[p] = glor.geo_idx(profs[1, p], noresmlon)
        latidn = latidn - (noresmlat.shape[1] * np.floor(latidn / noresmlat.shape[1]))
        lonidn = lonidn - (noresmlat.shape[1] * np.floor(lonidn / noresmlat.shape[1]))
        latidn = latidn.astype(int)
        lonidn = lonidn.astype(int)
        noresmnit = np.array(ncnornio[NIVAvars[0]][:]) * 1000
        noresmpho = np.array(ncnorpoo[NIVAvars[1]][:]) * 1000
        noresmsil = np.array(ncnorsio[NIVAvars[2]][:]) * 1000
        noresmoxy = np.array(ncnoro2o[NIVAvars[3]][:]) * 1000

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                Tinin = datetime(gyear[wd], gmonth[wd], 15)
                Tininx1 = d2i(Tinin, ncnornio['time'], select='exact')
                noresmnit_smp = noresmnit[Tininx1, :, latidn[dp], lonidn[dp]]
                noresmpho_smp = noresmpho[Tininx1, :, latidn[dp], lonidn[dp]]
                noresmsil_smp = noresmsil[Tininx1, :, latidn[dp], lonidn[dp]]
                noresmoxy_smp = noresmoxy[Tininx1, :, latidn[dp], lonidn[dp]]
            if dailyres == 1:
                    # By day
                    if gday[wd] != 15:
                        Tinin = datetime(gyear[wd], gmonth[wd], gday[wd])
                        Tininx1 = d2i(Tinin, ncnornio['time'], select='before')
                        Tininx2 = d2i(Tinin, ncnornio['time'], select='after')
                        if gday[wd] > 15:
                            tmult = gday[wd] - 15
                        else:
                            tmult = gday[wd] + 15
                        # Tidx1
                        noresmnit_1 = noresmnit[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmpho_1 = noresmpho[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmsil_1 = noresmsil[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmoxy_1 = noresmoxy[Tininx1, :, latidn[dp], lonidn[dp]]
                        # Tidx2
                        noresmnit_2 = noresmnit[Tininx2, :, latidn[dp], lonidn[dp]]
                        noresmpho_2 = noresmpho[Tininx2, :, latidn[dp], lonidn[dp]]
                        noresmsil_2 = noresmsil[Tininx2, :, latidn[dp], lonidn[dp]]
                        noresmoxy_2 = noresmoxy[Tininx2, :, latidn[dp], lonidn[dp]]
                        # Interpolating to the day
                        noresmnit_smp = noresmnit_1 + (((noresmnit_2 - noresmnit_1)/30)*tmult)
                        noresmpho_smp = noresmpho_1 + (((noresmpho_2 - noresmpho_1)/30)*tmult)
                        noresmsil_smp = noresmsil_1 + (((noresmsil_2 - noresmsil_1)/30)*tmult)
                        noresmoxy_smp = noresmoxy_1 + (((noresmoxy_2 - noresmoxy_1)/30)*tmult)
                    else:
                        Tinin = datetime(gyear[wd], gmonth[wd], 15)
                        Tininx1 = d2i(Tinin, ncnornio['time'], select='exact')
                        noresmnit_smp = noresmnit[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmpho_smp = noresmpho[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmsil_smp = noresmsil[Tininx1, :, latidn[dp], lonidn[dp]]
                        noresmoxy_smp = noresmoxy[Tininx1, :, latidn[dp], lonidn[dp]]

            if (any(noresmnit_smp > 10000) == True) & (all(noresmnit_smp > 10000) == False):
                noresmnit_smp[noresmnit_smp > 10000] = noresmnit_smp[noresmnit_smp < 10000][-1]
                noresmpho_smp[noresmpho_smp > 10000] = noresmpho_smp[noresmpho_smp < 10000][-1]
                noresmsil_smp[noresmsil_smp > 10000] = noresmsil_smp[noresmsil_smp < 10000][-1]
                noresmoxy_smp[noresmoxy_smp > 10000] = noresmoxy_smp[noresmoxy_smp < 10000][-1]
            elif all(noresmnit_smp > 10000) == True:
                land_noresm[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(depthniva, noresmnit_smp, fill_value='extrapolate')
            wodnit_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(depthniva, noresmpho_smp, fill_value='extrapolate')
            wodpho_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(depthniva, noresmsil_smp, fill_value='extrapolate')
            wodsil_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_o = interpolate.interp1d(depthniva, noresmoxy_smp, fill_value='extrapolate')
            wodoxy_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    ncnorni.close()
    ncnorpo.close()
    ncnorsi.close()
    ncnoro2.close()

    # Block for WOA processing
    valdir = '/media/dskone/VAL/'
    woadir = valdir + 'WOA18/'

    # Common woa elements for monthly climatology
    woa_exn = 'woa18_all_n01_01.nc'
    woa_exo = 'woa18_all_o01_01.nc'
    woan_01 = netcdf(woadir + woa_exn, 'r')
    woao_01 = netcdf(woadir + woa_exo, 'r')
    woan_01o = woan_01.variables
    woao_01o = woao_01.variables
    woa_lat = np.array(woan_01o['lat'][:])
    woa_lon = np.array(woan_01o['lon'][:])
    woa_depth_nps = np.array(woan_01o['depth'][:])
    woa_depth_nps_bnds = np.array(woan_01o['depth_bnds'][:])
    n_01a = np.array(woan_01o['n_an'][:])
    o_01a = np.array(woao_01o['o_an'][:])
    woan_01.close()
    woao_01.close()

    woa_sta = 'woa18_all_'
    woa_end = '_01.nc'

    woa_n = np.tile(n_01a, (12, 1, 1, 1))
    woa_o = np.tile(o_01a, (12, 1, 1, 1))
    woa_n = np.zeros_like(woa_n)
    woa_p = np.zeros_like(woa_n)
    woa_s = np.zeros_like(woa_n)
    woa_o = np.zeros_like(woa_o)

    for mt in range(0, 12):
        # nitrate
        woa_n_mth = netcdf(woadir + woa_sta + 'n' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_n_mth_o = woa_n_mth.variables
        woa_n[mt, :, :, :] = np.array(woa_n_mth_o['n_an'][:])
        woa_depth_nps = np.array(woa_n_mth_o['depth'][:])
        woa_depth_nps_bnds = np.array(woa_n_mth_o['depth_bnds'][:])
        woa_n_mth.close()

        # phosphate
        woa_p_mth = netcdf(woadir + woa_sta + 'p' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_p_mth_o = woa_p_mth.variables
        woa_p[mt, :, :, :] = np.array(woa_p_mth_o['p_an'][:])
        woa_p_mth.close()

        # silicate
        woa_s_mth = netcdf(woadir + woa_sta + 'i' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_s_mth_o = woa_s_mth.variables
        woa_s[mt, :, :, :] = np.array(woa_s_mth_o['i_an'][:])
        woa_s_mth.close()

        # oxygen
        woa_o_mth = netcdf(woadir + woa_sta + 'o' + str(mt + 1).zfill(2) + woa_end, 'r')
        woa_o_mth_o = woa_o_mth.variables
        woa_o[mt, :, :, :] = np.array(woa_o_mth_o['o_an'][:])
        woa_depth_ox = np.array(woa_o_mth_o['depth'][:])
        woa_depth_ox_bnds = np.array(woa_o_mth_o['depth_bnds'][:])
        woa_o_mth.close()

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd] + wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd] + wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidw = np.zeros_like(prof_ind, dtype='int64')
        lonidw = np.zeros_like(prof_ind, dtype='int64')

        for p in range(0, profs.shape[1]):
            latidw[p] = glor.geo_idx(profs[0, p], woa_lat)
            lonidw[p] = glor.geo_idx(profs[1, p], woa_lat)
        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                woanit_smp = woa_n[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                woapho_smp = woa_p[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                woasil_smp = woa_s[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                woaoxy_smp = woa_o[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
            if dailyres == 1:
                # Day
                if wod_day[wd] > 15:
                    tmult = wod_day[wd] - 15
                    tidx1 = wod_month[wd] - 1
                    if wod_month[wd] == 12:
                        tidx2 = 0
                    else:
                        tidx2 = wod_month[wd]

                    # Tidx1
                    woanit_1 = woa_n[tidx1, :, latidw[dp], lonidw[dp]]
                    woapho_1 = woa_p[tidx1, :, latidw[dp], lonidw[dp]]
                    woasil_1 = woa_s[tidx1, :, latidw[dp], lonidw[dp]]
                    woaoxy_1 = woa_o[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    woanit_2 = woa_n[tidx2, :, latidw[dp], lonidw[dp]]
                    woapho_2 = woa_p[tidx2, :, latidw[dp], lonidw[dp]]
                    woasil_2 = woa_s[tidx2, :, latidw[dp], lonidw[dp]]
                    woaoxy_2 = woa_o[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    woanit_smp = woanit_1 + (((woanit_2 - woanit_1)/30)*tmult)
                    woapho_smp = woapho_1 + (((woapho_2 - woapho_1)/30)*tmult)
                    woasil_smp = woasil_1 + (((woasil_2 - woasil_1)/30)*tmult)
                    woaoxy_smp = woaoxy_1 + (((woaoxy_2 - woaoxy_1)/30)*tmult)

                elif wod_day[wd] < 15:
                    tmult = wod_day[wd] + 15
                    tidx2 = wod_month[wd] - 1
                    if wod_month[wd] == 1:
                        tidx1 = 11
                    else:
                        tidx1 = wod_month[wd] - 2

                    # Tidx1
                    woanit_1 = woa_n[tidx1, :, latidw[dp], lonidw[dp]]
                    woapho_1 = woa_p[tidx1, :, latidw[dp], lonidw[dp]]
                    woasil_1 = woa_s[tidx1, :, latidw[dp], lonidw[dp]]
                    woaoxy_1 = woa_o[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    woanit_2 = woa_n[tidx2, :, latidw[dp], lonidw[dp]]
                    woapho_2 = woa_p[tidx2, :, latidw[dp], lonidw[dp]]
                    woasil_2 = woa_s[tidx2, :, latidw[dp], lonidw[dp]]
                    woaoxy_2 = woa_o[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    woanit_smp = woanit_1 + (((woanit_2 - woanit_1)/30)*tmult)
                    woapho_smp = woapho_1 + (((woapho_2 - woapho_1)/30)*tmult)
                    woasil_smp = woasil_1 + (((woasil_2 - woasil_1)/30)*tmult)
                    woaoxy_smp = woaoxy_1 + (((woaoxy_2 - woaoxy_1)/30)*tmult)

                else:
                    woanit_smp = woa_n[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    woapho_smp = woa_p[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    woasil_smp = woa_s[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    woaoxy_smp = woa_o[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]

            if (any(woanit_smp > 10000) == True) & (all(woanit_smp > 10000) == False):
                woanit_smp[woanit_smp > 10000] = woanit_smp[woanit_smp < 10000][-1]
                woapho_smp[woapho_smp > 10000] = woapho_smp[woapho_smp < 10000][-1]
                woasil_smp[woasil_smp > 10000] = woasil_smp[woasil_smp < 10000][-1]
                woaoxy_smp[woaoxy_smp > 10000] = woaoxy_smp[woaoxy_smp < 10000][-1]
            elif all(woanit_smp > 10000) == True:
                land_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_n = interpolate.interp1d(woa_depth_nps, woanit_smp, fill_value='extrapolate')
            wodnit_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(woa_depth_nps, woapho_smp, fill_value='extrapolate')
            wodpho_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(woa_depth_nps, woasil_smp, fill_value='extrapolate')
            wodsil_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_o = interpolate.interp1d(woa_depth_ox, woaoxy_smp, fill_value='extrapolate')
            wodoxy_woa[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for EMODNET DIVA processing

    emcdir = '/media/dsktwo/VAL/EMODNET/CLIMATOLOGY/'

    # Common emc elements for monthly climatology
    emc_ex = 'Water_body_dissolved_inorganic_nitrogen_v2021.nc'
    emc_01 = netcdf(emcdir + emc_ex, 'r')
    emc_01o = emc_01.variables
    emc_lat = np.array(emc_01o['lat'][:])
    emc_lon = np.array(emc_01o['lon'][:])
    emc_depth = np.array(emc_01o['depth'][:])
    n_01a = np.array(emc_01o['Water body dissolved inorganic nitrogen (DIN)'][:])
    emc_01.close()

    emc_sta = 'Water_body_'
    emc_end = '_v2021.nc'

    # nitrate
    emc_n_mth = netcdf(emcdir + emc_sta + 'dissolved_inorganic_nitrogen' + emc_end, 'r')
    emc_n_mth_o = emc_n_mth.variables
    emc_n = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)'][:])
    emc_n_L1 = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)_L1'][:])
    emc_n_L2 = np.array(emc_n_mth_o['Water body dissolved inorganic nitrogen (DIN)_L2'][:])
    emc_n_mth.close()

    # phosphate
    emc_p_mth = netcdf(emcdir + emc_sta + 'phosphate' + emc_end, 'r')
    emc_p_mth_o = emc_p_mth.variables
    emc_p = np.array(emc_p_mth_o['Water body phosphate'][:])
    emc_p_L1 = np.array(emc_p_mth_o['Water body phosphate_L1'][:])
    emc_p_L2 = np.array(emc_p_mth_o['Water body phosphate_L2'][:])
    emc_p_mth.close()

    # silicate
    emc_s_mth = netcdf(emcdir + emc_sta + 'silicate' + emc_end, 'r')
    emc_s_mth_o = emc_s_mth.variables
    emc_s = np.array(emc_s_mth_o['Water body silicate'][:])
    emc_s_L1 = np.array(emc_s_mth_o['Water body silicate_L1'][:])
    emc_s_L2 = np.array(emc_s_mth_o['Water body silicate_L2'][:])
    emc_s_mth.close()

    # oxygen
    emc_o_mth = netcdf(emcdir + emc_sta + 'dissolved_oxygen_concentration' + emc_end, 'r')
    emc_o_mth_o = emc_o_mth.variables
    emc_o = np.array(emc_o_mth_o['Water body dissolved oxygen concentration'][:])
    emc_o_L1 = np.array(emc_o_mth_o['Water body dissolved oxygen concentration_L1'][:])
    emc_o_L2 = np.array(emc_o_mth_o['Water body dissolved oxygen concentration_L2'][:])
    emc_o_mth.close()

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd] + wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd] + wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidd = np.zeros_like(prof_ind, dtype='int64')
        lonidd = np.zeros_like(prof_ind, dtype='int64')

        for p in range(0, profs.shape[1]):
            latidd[p] = glor.geo_idx(profs[0, p], emc_lat)
            lonidd[p] = glor.geo_idx(profs[1, p], emc_lon)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                emcnit_smp = emc_n[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcnit_smp_L1 = emc_n_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcnit_smp_L2 = emc_n_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcpho_smp = emc_p[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcpho_smp_L1 = emc_p_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcpho_smp_L2 = emc_p_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcsil_smp = emc_s[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcsil_smp_L1 = emc_s_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcsil_smp_L2 = emc_s_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcoxy_smp = emc_o[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcoxy_smp_L1 = emc_o_L1[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                emcoxy_smp_L2 = emc_o_L2[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]

            if dailyres == 1:
                # Daily
                if wod_day[wd] > 15:
                    tmult = wod_day[wd] - 15
                    tidx1 = wod_month[wd] - 1
                    if wod_month[wd] == 12:
                        tidx2 = 0
                    else:
                        tidx2 = wod_month[wd]

                    # Tidx1
                    emcnit_1 = emc_n[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1 = emc_p[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1 = emc_s[tidx1, :, latidd[dp], lonidd[dp]]
                    emcoxy_1 = emc_o[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L1 = emc_n_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L1 = emc_p_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L1 = emc_s_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcoxy_1_L1 = emc_o_L1[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L2 = emc_n_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L2 = emc_p_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L2 = emc_s_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcoxy_1_L2 = emc_o_L2[tidx1, :, latidd[dp], lonidd[dp]]

                    # Tidx2
                    emcnit_2 = emc_n[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2 = emc_p[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2 = emc_s[tidx2, :, latidd[dp], lonidd[dp]]
                    emcoxy_2 = emc_o[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L1 = emc_n_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L1 = emc_p_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L1 = emc_s_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcoxy_2_L1 = emc_o_L1[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L2 = emc_n_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L2 = emc_p_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L2 = emc_s_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcoxy_2_L2 = emc_o_L2[tidx2, :, latidd[dp], lonidd[dp]]

                    # Interpolating to the day
                    emcnit_smp = emcnit_1 + (((emcnit_2 - emcnit_1)/30)*tmult)
                    emcpho_smp = emcpho_1 + (((emcpho_2 - emcpho_1)/30)*tmult)
                    emcsil_smp = emcsil_1 + (((emcsil_2 - emcsil_1)/30)*tmult)
                    emcoxy_smp = emcoxy_1 + (((emcoxy_2 - emcoxy_1)/30)*tmult)

                    emcnit_smp_L1 = emcnit_1_L1 + (((emcnit_2_L1 - emcnit_1_L1)/30)*tmult)
                    emcpho_smp_L1 = emcpho_1_L1 + (((emcpho_2_L1 - emcpho_1_L1)/30)*tmult)
                    emcsil_smp_L1 = emcsil_1_L1 + (((emcsil_2_L1 - emcsil_1_L1)/30)*tmult)
                    emcoxy_smp_L1 = emcoxy_1_L1 + (((emcoxy_2_L1 - emcoxy_1_L1)/30)*tmult)

                    emcnit_smp_L2 = emcnit_1_L2 + (((emcnit_2_L2 - emcnit_1_L2)/30)*tmult)
                    emcpho_smp_L2 = emcpho_1_L2 + (((emcpho_2_L2 - emcpho_1_L2)/30)*tmult)
                    emcsil_smp_L2 = emcsil_1_L2 + (((emcsil_2_L2 - emcsil_1_L2)/30)*tmult)
                    emcoxy_smp_L2 = emcoxy_1_L2 + (((emcoxy_2_L2 - emcoxy_1_L2)/30)*tmult)

                elif wod_day[wd] < 15:
                    tmult = wod_day[wd] + 15
                    tidx2 = wod_month[wd] - 1
                    if wod_month[wd] == 1:
                        tidx1 = 11
                    else:
                        tidx1 = wod_month[wd] - 2

                    # Tidx1
                    emcnit_1 = emc_n[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1 = emc_p[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1 = emc_s[tidx1, :, latidd[dp], lonidd[dp]]
                    emcoxy_1 = emc_o[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L1 = emc_n_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L1 = emc_p_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L1 = emc_s_L1[tidx1, :, latidd[dp], lonidd[dp]]
                    emcoxy_1_L1 = emc_o_L1[tidx1, :, latidd[dp], lonidd[dp]]

                    emcnit_1_L2 = emc_n_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcpho_1_L2 = emc_p_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcsil_1_L2 = emc_s_L2[tidx1, :, latidd[dp], lonidd[dp]]
                    emcoxy_1_L2 = emc_o_L2[tidx1, :, latidd[dp], lonidd[dp]]

                    # Tidx2
                    emcnit_2 = emc_n[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2 = emc_p[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2 = emc_s[tidx2, :, latidd[dp], lonidd[dp]]
                    emcoxy_2 = emc_o[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L1 = emc_n_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L1 = emc_p_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L1 = emc_s_L1[tidx2, :, latidd[dp], lonidd[dp]]
                    emcoxy_2_L1 = emc_o_L1[tidx2, :, latidd[dp], lonidd[dp]]

                    emcnit_2_L2 = emc_n_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcpho_2_L2 = emc_p_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcsil_2_L2 = emc_s_L2[tidx2, :, latidd[dp], lonidd[dp]]
                    emcoxy_2_L2 = emc_o_L2[tidx2, :, latidd[dp], lonidd[dp]]

                    # Interpolating to the day
                    emcnit_smp = emcnit_1 + (((emcnit_2 - emcnit_1)/30)*tmult)
                    emcpho_smp = emcpho_1 + (((emcpho_2 - emcpho_1)/30)*tmult)
                    emcsil_smp = emcsil_1 + (((emcsil_2 - emcsil_1)/30)*tmult)
                    emcoxy_smp = emcoxy_1 + (((emcoxy_2 - emcoxy_1)/30)*tmult)

                    emcnit_smp_L1 = emcnit_1_L1 + (((emcnit_2_L1 - emcnit_1_L1)/30)*tmult)
                    emcpho_smp_L1 = emcpho_1_L1 + (((emcpho_2_L1 - emcpho_1_L1)/30)*tmult)
                    emcsil_smp_L1 = emcsil_1_L1 + (((emcsil_2_L1 - emcsil_1_L1)/30)*tmult)
                    emcoxy_smp_L1 = emcoxy_1_L1 + (((emcoxy_2_L1 - emcoxy_1_L1)/30)*tmult)

                    emcnit_smp_L2 = emcnit_1_L2 + (((emcnit_2_L2 - emcnit_1_L2)/30)*tmult)
                    emcpho_smp_L2 = emcpho_1_L2 + (((emcpho_2_L2 - emcpho_1_L2)/30)*tmult)
                    emcsil_smp_L2 = emcsil_1_L2 + (((emcsil_2_L2 - emcsil_1_L2)/30)*tmult)
                    emcoxy_smp_L2 = emcoxy_1_L2 + (((emcoxy_2_L2 - emcoxy_1_L2)/30)*tmult)

                else:
                    emcnit_smp = emc_n[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                    emcpho_smp = emc_p[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                    emcsil_smp = emc_s[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]
                    emcoxy_smp = emc_o[wod_month[wd] - 1, :, latidd[dp], lonidd[dp]]

            if (any(emcnit_smp > 10000) == True) & (all(emcnit_smp > 10000) == False):
                emcnit_smp[emcnit_smp > 10000] = emcnit_smp[emcnit_smp < 10000][-1]
                emcpho_smp[emcpho_smp > 10000] = emcpho_smp[emcpho_smp < 10000][-1]
                emcsil_smp[emcsil_smp > 10000] = emcsil_smp[emcsil_smp < 10000][-1]
                emcoxy_smp[emcoxy_smp > 10000] = emcoxy_smp[emcoxy_smp < 10000][-1]
            elif all(emcnit_smp > 10000) == True:
                land_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True

            if (any(emcnit_smp_L1 > 10000) == True) & (all(emcnit_smp_L1 > 10000) == False):
                emcnit_smp_L1[emcnit_smp_L1 > 10000] = emcnit_smp_L1[emcnit_smp_L1 < 10000][-1]
            elif all(emcnit_smp_L1 > 10000) == True:
                land_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True

            if (any(emcnit_smp_L2 > 10000) == True) & (all(emcnit_smp_L2 > 10000) == False):
                emcnit_smp_L2[emcnit_smp_L2 > 10000] = emcnit_smp_L2[emcnit_smp_L2 < 10000][-1]
            elif all(emcnit_smp_L2 > 10000) == True:
                land_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True

            if (any(emcpho_smp_L1 > 10000) == True) & (all(emcpho_smp_L1 > 10000) == False):
                emcpho_smp_L1[emcpho_smp_L1 > 10000] = emcpho_smp_L1[emcpho_smp_L1 < 10000][-1]
            if (any(emcpho_smp_L2 > 10000) == True) & (all(emcpho_smp_L2 > 10000) == False):
                emcpho_smp_L2[emcpho_smp_L2 > 10000] = emcpho_smp_L2[emcpho_smp_L2 < 10000][-1]

            if (any(emcsil_smp_L1 > 10000) == True) & (all(emcsil_smp_L1 > 10000) == False):
                emcsil_smp_L1[emcsil_smp_L1 > 10000] = emcsil_smp_L1[emcsil_smp_L1 < 10000][-1]
            if (any(emcsil_smp_L2 > 10000) == True) & (all(emcsil_smp_L2 > 10000) == False):
                emcsil_smp_L2[emcsil_smp_L2 > 10000] = emcsil_smp_L2[emcsil_smp_L2 < 10000][-1]

            if (any(emcoxy_smp_L1 > 10000) == True) & (all(emcoxy_smp_L1 > 10000) == False):
                emcoxy_smp_L1[emcoxy_smp_L1 > 10000] = emcoxy_smp_L1[emcoxy_smp_L1 < 10000][-1]
            if (any(emcoxy_smp_L2 > 10000) == True) & (all(emcoxy_smp_L2 > 10000) == False):
                emcoxy_smp_L2[emcoxy_smp_L2 > 10000] = emcoxy_smp_L2[emcoxy_smp_L2 < 10000][-1]
            #
            wodint_n = interpolate.interp1d(emc_depth, emcnit_smp, fill_value='extrapolate')
            wodnit_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_n_L1 = interpolate.interp1d(emc_depth, emcnit_smp_L1, fill_value='extrapolate')
            wodnit_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_n_L2 = interpolate.interp1d(emc_depth, emcnit_smp_L2, fill_value='extrapolate')
            wodnit_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_n_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_p = interpolate.interp1d(emc_depth, emcpho_smp, fill_value='extrapolate')
            wodpho_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_p_L1 = interpolate.interp1d(emc_depth, emcpho_smp_L1, fill_value='extrapolate')
            wodpho_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_p_L2 = interpolate.interp1d(emc_depth, emcpho_smp_L2, fill_value='extrapolate')
            wodpho_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(emc_depth, emcsil_smp, fill_value='extrapolate')
            wodsil_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_s_L1 = interpolate.interp1d(emc_depth, emcsil_smp_L1, fill_value='extrapolate')
            wodsil_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_s_L2 = interpolate.interp1d(emc_depth, emcsil_smp_L2, fill_value='extrapolate')
            wodsil_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_o = interpolate.interp1d(emc_depth, emcoxy_smp, fill_value='extrapolate')
            wodoxy_diva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_o_L1 = interpolate.interp1d(emc_depth, emcoxy_smp_L1, fill_value='extrapolate')
            wodoxy_diva_L1[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o_L1(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            wodint_o_L2 = interpolate.interp1d(emc_depth, emcoxy_smp_L2, fill_value='extrapolate')
            wodoxy_diva_L2[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_o_L2(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    glodap_bias_file = crocodir + 'glodap_bias_' + npzstr + '.npz'
    np.savez(glodap_bias_file,
             wodlat=wodlat, wodlon=wodlon, z=z,
             wodnit=wodnit,
             wodnit_ibim=wodnit_ibim, wodnit_cs1k=wodnit_cs1k, wodnit_niva=wodnit_niva,
             wodnit_woa=wodnit_woa,
             wodnit_diva=wodnit_diva, wodnit_diva_L1=wodnit_diva_L1, wodnit_diva_L2=wodnit_diva_L2,
             wodpho=wodpho,
             wodpho_ibim=wodpho_ibim, wodpho_cs1k=wodpho_cs1k, wodpho_niva=wodpho_niva,
             wodpho_woa=wodpho_woa,
             wodpho_diva=wodpho_diva, wodpho_diva_L1=wodpho_diva_L1, wodpho_diva_L2=wodpho_diva_L2,
             wodsil=wodsil,
             wodsil_ibim=wodsil_ibim, wodsil_cs1k=wodsil_cs1k, wodsil_niva=wodsil_niva,
             wodsil_woa=wodsil_woa,
             wodsil_diva=wodsil_diva, wodsil_diva_L1=wodsil_diva_L1, wodsil_diva_L2=wodsil_diva_L2,
             wodoxy=wodoxy,
             wodoxy_ibim=wodoxy_ibim, wodoxy_cs1k=wodoxy_cs1k, wodoxy_niva=wodoxy_niva,
             wodoxy_woa=wodoxy_woa,
             wodoxy_diva=wodoxy_diva, wodoxy_diva_L1=wodoxy_diva_L1, wodoxy_diva_L2=wodoxy_diva_L2,
             latids=latids, lonids=lonids, latidc=latidc, lonidc=lonidc, latidn=latidn, lonidn=lonidn,
             latidw=latidw, lonidw=lonidw, latidd=latidd, lonidd=lonidd,
             land_ibi=land_ibi, land_cs1km=land_cs1km, land_noresm=land_noresm,
             land_diva=land_diva, land_diva_L1=land_diva_L1, land_diva_L2=land_diva_L2,
             wod_day=gday, wod_month=gmonth, wod_year=gyear)

# Load data
# Identify the data of interest and dimension of interest
# e.g. month, time since 1993, depth, spatial.
# timeseries plot (beginning to end, depth binning)
# climatology plot (equivalent depth binning)
# lat/lon plot by depth bins
# Try to include data count for context

# Order to follow:
# Inspect data count to rationalise what is possible and meaningful

wod_plotting = 1
mi_plotting = 1
glodap_plotting = 1

if wod_plotting == 1:
    wod_bias_file = crocodir + 'wod_bias_' + npzstr + '.npz'
    wod_bias = np.load(wod_bias_file)
    wodlat = wod_bias['wodlat']
    wodlon = wod_bias['wodlon']
    z = wod_bias['z']
    wodnit = wod_bias['wodnit']*1.025
    wodnit_f = wod_bias['wodnit_f']
    wodnit_ibim = wod_bias['wodnit_ibim']
    wodnit_cs1k = wod_bias['wodnit_cs1k']
    wodnit_niva = wod_bias['wodnit_niva']
    wodnit_diva = wod_bias['wodnit_diva']
    wodnit_woa = wod_bias['wodnit_woa']
    wodnit_diva_L1 = wod_bias['wodnit_diva_L1']
    wodnit_diva_L2 = wod_bias['wodnit_diva_L2']
    wodpho = wod_bias['wodpho']*1.025
    wodpho_f = wod_bias['wodpho_f']
    wodpho_ibim = wod_bias['wodpho_ibim']
    wodpho_cs1k = wod_bias['wodpho_cs1k']
    wodpho_niva = wod_bias['wodpho_niva']
    wodpho_diva = wod_bias['wodpho_diva']
    wodpho_woa = wod_bias['wodpho_woa']
    wodpho_diva_L1 = wod_bias['wodpho_diva_L1']
    wodpho_diva_L2 = wod_bias['wodpho_diva_L2']
    wodsil = wod_bias['wodsil']*1.025
    wodsil_f = wod_bias['wodsil_f']
    wodsil_ibim = wod_bias['wodsil_ibim']
    wodsil_cs1k = wod_bias['wodsil_cs1k']
    wodsil_niva = wod_bias['wodsil_niva']
    wodsil_diva = wod_bias['wodsil_diva']
    wodsil_woa = wod_bias['wodsil_woa']
    wodsil_diva_L1 = wod_bias['wodsil_diva_L1']
    wodsil_diva_L2 = wod_bias['wodsil_diva_L2']
    wodoxy = wod_bias['wodoxy']*1.025
    wodoxy_f = wod_bias['wodoxy_f']
    wodoxy_ibim = wod_bias['wodoxy_ibim']
    wodoxy_cs1k = wod_bias['wodoxy_cs1k']
    wodoxy_niva = wod_bias['wodoxy_niva']
    wodoxy_diva = wod_bias['wodoxy_diva']
    wodoxy_woa = wod_bias['wodoxy_woa']
    wodoxy_diva_L1 = wod_bias['wodoxy_diva_L1']
    wodoxy_diva_L2 = wod_bias['wodoxy_diva_L2']
    land_ibi = wod_bias['land_ibi']
    land_cs1km = wod_bias['land_cs1km']
    land_noresm = wod_bias['land_noresm']
    land_diva = wod_bias['land_diva']
    land_diva_L1 = wod_bias['land_diva_L1']
    land_diva_L2 = wod_bias['land_diva_L2']
    wod_day = wod_bias['wod_day'].astype(int)
    wod_month = wod_bias['wod_month'].astype(int)
    wod_year = wod_bias['wod_year'].astype(int)

if mi_plotting == 1:
    mi_bias_file = crocodir + 'MI_bias_' + npzstr + '.npz'
    mi_bias = np.load(mi_bias_file)
    wodlat = np.concatenate((wodlat, mi_bias['wodlat']))
    wodlon = np.concatenate((wodlon, mi_bias['wodlon']))
    z = np.concatenate((z, mi_bias['z']))
    wodnit = np.concatenate((wodnit, mi_bias['wodnit']))
    wodnit_ibim = np.concatenate((wodnit_ibim, mi_bias['wodnit_ibim']))
    wodnit_cs1k = np.concatenate((wodnit_cs1k, mi_bias['wodnit_cs1k']))
    wodnit_niva = np.concatenate((wodnit_niva, mi_bias['wodnit_niva']))
    wodnit_diva = np.concatenate((wodnit_diva, mi_bias['wodnit_diva']))
    wodnit_woa = np.concatenate((wodnit_woa, mi_bias['wodnit_woa']))
    wodnit_diva_L1 = np.concatenate((wodnit_diva_L1, mi_bias['wodnit_diva_L1']))
    wodnit_diva_L2 = np.concatenate((wodnit_diva_L2, mi_bias['wodnit_diva_L2']))
    wodpho = np.concatenate((wodpho, mi_bias['wodpho']))
    wodpho_ibim = np.concatenate((wodpho_ibim, mi_bias['wodpho_ibim']))
    wodpho_cs1k = np.concatenate((wodpho_cs1k, mi_bias['wodpho_cs1k']))
    wodpho_niva = np.concatenate((wodpho_niva, mi_bias['wodpho_niva']))
    wodpho_diva = np.concatenate((wodpho_diva, mi_bias['wodpho_diva']))
    wodpho_woa = np.concatenate((wodpho_woa, mi_bias['wodpho_woa']))
    wodpho_diva_L1 = np.concatenate((wodpho_diva_L1, mi_bias['wodpho_diva_L1']))
    wodpho_diva_L2 = np.concatenate((wodpho_diva_L2, mi_bias['wodpho_diva_L2']))
    wodsil = np.concatenate((wodsil, mi_bias['wodsil']))
    wodsil_ibim = np.concatenate((wodsil_ibim, mi_bias['wodsil_ibim']))
    wodsil_cs1k = np.concatenate((wodsil_cs1k, mi_bias['wodsil_cs1k']))
    wodsil_niva = np.concatenate((wodsil_niva, mi_bias['wodsil_niva']))
    wodsil_diva = np.concatenate((wodsil_diva, mi_bias['wodsil_diva']))
    wodsil_woa = np.concatenate((wodsil_woa, mi_bias['wodsil_woa']))
    wodsil_diva_L1 = np.concatenate((wodsil_diva_L1, mi_bias['wodsil_diva_L1']))
    wodsil_diva_L2 = np.concatenate((wodsil_diva_L2, mi_bias['wodsil_diva_L2']))
    land_ibi = np.concatenate((land_ibi, mi_bias['land_ibi']))
    land_cs1km = np.concatenate((land_cs1km, mi_bias['land_cs1km']))
    land_noresm = np.concatenate((land_noresm, mi_bias['land_noresm']))
    land_diva = np.concatenate((land_diva, mi_bias['land_diva']))
    land_diva_L1 = np.concatenate((land_diva_L1, mi_bias['land_diva_L1']))
    land_diva_L2 = np.concatenate((land_diva_L2, mi_bias['land_diva_L2']))
    wod_day = np.concatenate((wod_day, mi_bias['wod_day']))
    wod_month = np.concatenate((wod_month, mi_bias['wod_month']))
    wod_year = np.concatenate((wod_year, mi_bias['wod_year']))

if glodap_plotting == 1:
    glodap_bias_file = crocodir + 'glodap_bias_' + npzstr + '.npz'
    glodap_bias = np.load(glodap_bias_file)
    wodlat = np.concatenate((wodlat, glodap_bias['wodlat']))
    wodlon = np.concatenate((wodlon, glodap_bias['wodlon']))
    z = np.concatenate((z, glodap_bias['z']))
    wodnit = np.concatenate((wodnit, glodap_bias['wodnit']*1.025))
    wodnit_ibim = np.concatenate((wodnit_ibim, glodap_bias['wodnit_ibim']))
    wodnit_cs1k = np.concatenate((wodnit_cs1k, glodap_bias['wodnit_cs1k']))
    wodnit_niva = np.concatenate((wodnit_niva, glodap_bias['wodnit_niva']))
    wodnit_diva = np.concatenate((wodnit_diva, glodap_bias['wodnit_diva']))
    wodnit_woa = np.concatenate((wodnit_woa, glodap_bias['wodnit_woa']))
    wodnit_diva_L1 = np.concatenate((wodnit_diva_L1, glodap_bias['wodnit_diva_L1']))
    wodnit_diva_L2 = np.concatenate((wodnit_diva_L2, glodap_bias['wodnit_diva_L2']))
    wodpho = np.concatenate((wodpho, glodap_bias['wodpho']*1.025))
    wodpho_ibim = np.concatenate((wodpho_ibim, glodap_bias['wodpho_ibim']))
    wodpho_cs1k = np.concatenate((wodpho_cs1k, glodap_bias['wodpho_cs1k']))
    wodpho_niva = np.concatenate((wodpho_niva, glodap_bias['wodpho_niva']))
    wodpho_diva = np.concatenate((wodpho_diva, glodap_bias['wodpho_diva']))
    wodpho_woa = np.concatenate((wodpho_woa, glodap_bias['wodpho_woa']))
    wodpho_diva_L1 = np.concatenate((wodpho_diva_L1, glodap_bias['wodpho_diva_L1']))
    wodpho_diva_L2 = np.concatenate((wodpho_diva_L2, glodap_bias['wodpho_diva_L2']))
    wodsil = np.concatenate((wodsil, glodap_bias['wodsil']*1.025))
    wodsil_ibim = np.concatenate((wodsil_ibim, glodap_bias['wodsil_ibim']))
    wodsil_cs1k = np.concatenate((wodsil_cs1k, glodap_bias['wodsil_cs1k']))
    wodsil_niva = np.concatenate((wodsil_niva, glodap_bias['wodsil_niva']))
    wodsil_diva = np.concatenate((wodsil_diva, glodap_bias['wodsil_diva']))
    wodsil_woa = np.concatenate((wodsil_woa, glodap_bias['wodsil_woa']))
    wodsil_diva_L1 = np.concatenate((wodsil_diva_L1, glodap_bias['wodsil_diva_L1']))
    wodsil_diva_L2 = np.concatenate((wodsil_diva_L2, glodap_bias['wodsil_diva_L2']))
    wodoxy = np.concatenate((wodoxy, glodap_bias['wodoxy']*1.025))
    wodoxy_ibim = np.concatenate((wodoxy_ibim, glodap_bias['wodoxy_ibim']))
    wodoxy_cs1k = np.concatenate((wodoxy_cs1k, glodap_bias['wodoxy_cs1k']))
    wodoxy_niva = np.concatenate((wodoxy_niva, glodap_bias['wodoxy_niva']))
    wodoxy_diva = np.concatenate((wodoxy_diva, glodap_bias['wodoxy_diva']))
    wodoxy_woa = np.concatenate((wodoxy_woa, glodap_bias['wodoxy_woa']))
    wodoxy_diva_L1 = np.concatenate((wodoxy_diva_L1, glodap_bias['wodoxy_diva_L1']))
    wodoxy_diva_L2 = np.concatenate((wodoxy_diva_L2, glodap_bias['wodoxy_diva_L2']))
    land_ibi = np.concatenate((land_ibi, glodap_bias['land_ibi']))
    land_cs1km = np.concatenate((land_cs1km, glodap_bias['land_cs1km']))
    land_noresm = np.concatenate((land_noresm, glodap_bias['land_noresm']))
    land_diva = np.concatenate((land_diva, glodap_bias['land_diva']))
    land_diva_L1 = np.concatenate((land_diva_L1, glodap_bias['land_diva_L1']))
    land_diva_L2 = np.concatenate((land_diva_L2, glodap_bias['land_diva_L2']))
    wod_day = np.concatenate((wod_day, glodap_bias['wod_day']))
    wod_month = np.concatenate((wod_month, glodap_bias['wod_month']))
    wod_year = np.concatenate((wod_year, glodap_bias['wod_year']))
    wodtime = pd.DataFrame({'year': wod_year,
                            'month': wod_month,
                            'day': wod_day})
    tstamp = pd.to_datetime(wodtime)

    def calc_distance(lat1, lon1, lat2, lon2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        c = 2 * asin(sqrt(a))
        km = 6371 * c
        return km

    def nut_bias_plt_stats(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa, dep,
                           wod_yr, wod_mth, wod_dy, otimes,
                           wod_lon, wod_lat, geodist, nutstr):
        ibimfilt = np.argwhere((obs > 0) & (obs < 10000) & (mod_ibim > 0) & (mod_ibim < 10000))[:, 0]
        cs1kfilt = np.argwhere((obs > 0) & (obs < 10000) & (mod_cs1k > 0) & (mod_cs1k < 10000))[:, 0]
        nivafilt = np.argwhere((obs > 0) & (obs < 10000) & (mod_niva > 0) & (mod_niva < 10000))[:, 0]
        divafilt = np.argwhere((obs > 0) & (obs < 10000) & (clim_diva > 0) & (clim_diva < 10000))[:, 0]
        woafilt = np.argwhere((obs > 0) & (obs < 10000) & (clim_woa > 0) & (clim_woa < 10000))[:, 0]
        allfilt = np.argwhere((obs > 0) & (obs < 10000) &
                              (mod_ibim > 0) & (mod_ibim < 10000) &
                              (mod_cs1k > 0) & (mod_cs1k < 10000) &
                              (mod_niva > 0) & (mod_niva < 10000) &
                              (clim_diva > 0) & (clim_diva < 10000) &
                              (clim_woa > 0) & (clim_woa < 10000))[:, 0]
        uppfilt = np.argwhere((dep <= 20))[:, 0]
        midfilt = np.argwhere((dep > 20) & (dep <= 70))[:, 0]
        depfilt = np.argwhere((dep > 70))[:, 0]
        djffilt = np.argwhere((wod_mth == 12) | (wod_mth == 1) | (wod_mth == 2))[:, 0]
        mamfilt = np.argwhere((wod_mth == 3) | (wod_mth == 4) | (wod_mth == 5))[:, 0]
        jjafilt = np.argwhere((wod_mth == 6) | (wod_mth == 7) | (wod_mth == 8))[:, 0]
        sonfilt = np.argwhere((wod_mth == 9) | (wod_mth == 10) | (wod_mth == 11))[:, 0]
        csdfilt = np.argwhere((wod_lat >= 49) & (wod_lat <= 52.95) & (wod_lon >= -10.75) & (wod_lon <= -5.83))[:, 0]

        cs1k_ibi = np.argwhere((obs > 0) & (obs < 10000) &
                               (mod_ibim > 0) & (mod_ibim < 10000) &
                               (mod_cs1k > 0) & (mod_cs1k < 10000))[:, 0]
        cs1kibi = np.intersect1d(cs1k_ibi, csdfilt)
        ibi_niva = np.argwhere((obs > 0) & (obs < 10000) &
                               (mod_ibim > 0) & (mod_ibim < 10000) &
                               (mod_niva > 0) & (mod_niva < 10000))[:, 0]
        ibiniva = np.setdiff1d(ibi_niva, csdfilt)

        woa_diva = np.argwhere((obs > 0) & (obs < 10000) &
                               (clim_woa > 0) & (clim_woa < 10000) &
                               (clim_diva > 0) & (clim_diva < 10000))[:, 0]
        woadiva = np.setdiff1d(woa_diva, csdfilt)

        allfilt_II = np.intersect1d(allfilt, csdfilt)

        # allfilt_t = allfilt
        allfilt_t = allfilt_II

        uppfiltII = np.intersect1d(ibimfilt, uppfilt)
        midfiltII = np.intersect1d(ibimfilt, midfilt)
        depfiltII = np.intersect1d(ibimfilt, depfilt)

        scatr5(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa,
               allfilt_II, nutstr, 'common')
        RMSD_ibim_t, R2_ibim_t, MD_ibim_t, N_ibim_t = get_eans(obs[allfilt_II], mod_ibim[allfilt_II])
        RMSD_cs1k_t, R2_cs1k_t, MD_cs1k_t, N_cs1k_t = get_eans(obs[allfilt_II], mod_cs1k[allfilt_II])
        RMSD_niva_t, R2_niva_t, MD_niva_t, N_niva_t = get_eans(obs[allfilt_II], mod_niva[allfilt_II])
        RMSD_diva_t, R2_diva_t, MD_diva_t, N_diva_t = get_eans(obs[allfilt_II], clim_diva[allfilt_II])
        RMSD_woa_t, R2_woa_t, MD_woa_t, N_woa_t = get_eans(obs[allfilt_II], clim_woa[allfilt_II])

        RMSD_ibim_1dm, R2_ibim_1dm, MD_ibim_1dm, N_ibim_1dm = get_eans(obs[allfilt], mod_ibim[allfilt])
        RMSD_niva_1dm, R2_niva_1dm, MD_niva_1dm, N_niva_1dm = get_eans(obs[allfilt], mod_niva[allfilt])
        RMSD_diva_1dm, R2_diva_1dm, MD_diva_1dm, N_diva_1dm = get_eans(obs[allfilt], clim_diva[allfilt])
        RMSD_woa_1dm, R2_woa_1dm, MD_woa_1dm, N_woa_1dm = get_eans(obs[allfilt], clim_woa[allfilt])

        biscatr(fdir, obs, mod_ibim, mod_cs1k, 'IBI', 'CS1K', cs1kibi, nutstr, 'mmol/m3', 'In-domain', 'Model')
        RMSD_ibi_dom, R2_ibi_dom, MD_ibi_dom, N_ibi_dom = get_eans(obs[cs1kibi], mod_ibim[cs1kibi])
        RMSD_cs1k_dom, R2_cs1k_dom, MD_cs1k_dom, N_cs1k_dom = get_eans(obs[cs1kibi], mod_cs1k[cs1kibi])

        biscatr(fdir, obs, mod_ibim, mod_niva, 'IBI', 'NIVA', ibiniva, nutstr, 'mmol/m3', 'Degree-margin', 'Model')
        RMSD_ibi_bry, R2_ibi_bry, MD_ibi_bry, N_ibi_bry = get_eans(obs[ibiniva], mod_ibim[ibiniva])
        RMSD_niva_bry, R2_niva_bry, MD_niva_bry, N_niva_bry = get_eans(obs[ibiniva], mod_niva[ibiniva])

        biscatr(fdir, obs, clim_woa, clim_diva, 'WOA', 'DIVA', woadiva, nutstr, 'mmol/m3', 'Degree-margin', 'Clim.')
        RMSD_woa_bry, R2_woa_bry, MD_woa_bry, N_woa_bry = get_eans(obs[woadiva], clim_woa[woadiva])
        RMSD_diva_bry, R2_diva_bry, MD_diva_bry, N_diva_bry = get_eans(obs[woadiva], clim_diva[woadiva])

        # Nutrient colormap for spatial plot
        nutcmp = 'tab20c'  # 'bwr_r', 'gist_ncar', 'RdYlBu'
        prod_biases(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa, wodlon, wodlat,
                    nutstr, '(mmol/m3)', 'tab20c', allfilt, cs1kibi)

        obsVsprods(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa, wodlon, wodlat,
                   nutstr, '(mmol/m3)', 'tab20c', allfilt, cs1kibi)

        obs_prod_biases_x9(fdir, obs, mod_ibim, mod_cs1k, mod_niva, wodlon, wodlat,
                           nutstr, '(mmol/m3)', 'tab20c', allfilt, cs1kibi)

        lt_tseries(fdir, otimes, 'bias', obs, mod_ibim, mod_cs1k, mod_niva,
                   allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'bias_models', 'Model')

        lt_tseries(fdir, otimes, 'bias', obs, mod_ibim, clim_woa, clim_diva,
                   allfilt_t, allfilt_t, allfilt_t, 'IBI', 'WOA', 'DIVA', nutstr, '(mmol/m3)', 'bias_iclims', 'IBI+Clims')

        lt_tseries(fdir, otimes, 'bias', obs, mod_niva, clim_woa, clim_diva,
                   allfilt_t, allfilt_t, allfilt_t, 'NIVA', 'WOA', 'DIVA', nutstr, '(mmol/m3)', 'bias_nclims', 'NIVA+Clims')

        lt_tseries(fdir, otimes, 'data', obs, mod_ibim, mod_cs1k, mod_niva,
                   allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'data_models', 'Model')

        lt_tseries(fdir, otimes, 'data', obs, mod_ibim, clim_woa, clim_diva,
                   allfilt_t, allfilt_t, allfilt_t, 'IBI', 'WOA', 'DIVA', nutstr, '(mmol/m3)', 'data_iclims', 'IBI+Clims')

        lt_tseries(fdir, otimes, 'data', obs, mod_niva, clim_woa, clim_diva,
                   allfilt_t, allfilt_t, allfilt_t, 'NIVA', 'WOA', 'DIVA', nutstr, '(mmol/m3)', 'data_nclims', 'NIVA+Clims')

        t_clim(fdir, wod_mth, 'bias', obs, mod_ibim, mod_cs1k, mod_niva,
               allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'bias_models', 'Model')

        t_clim(fdir, wod_mth, 'bias', obs, mod_ibim, clim_woa, clim_diva,
               allfilt_t, allfilt_t, allfilt_t, 'IBI', 'WOA', 'DIVA', nutstr, '(mmol/m3)', 'bias_iclims', 'Model')

        t_clim(fdir, wod_mth, 'data', obs, mod_ibim, mod_cs1k, mod_niva,
               allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'data_models', 'Model')

        t_clim(fdir, wod_mth, 'data', obs, mod_ibim, clim_woa, clim_diva,
               allfilt_t, allfilt_t, allfilt_t, 'IBI', 'WOA', 'DIVA', nutstr, '(mmol/m3)', 'data_iclims', 'Model')

        return (RMSD_ibim_t, R2_ibim_t, MD_ibim_t, N_ibim_t), (RMSD_cs1k_t, R2_cs1k_t, MD_cs1k_t, N_cs1k_t), \
            (RMSD_niva_t, R2_niva_t, MD_niva_t, N_niva_t), (RMSD_diva_t, R2_diva_t, MD_diva_t, N_diva_t), \
            (RMSD_woa_t, R2_woa_t, MD_woa_t, N_woa_t), \
            (RMSD_ibim_1dm, R2_ibim_1dm, MD_ibim_1dm, N_ibim_1dm), (RMSD_niva_1dm, R2_niva_1dm, MD_niva_1dm, N_niva_1dm), \
            (RMSD_diva_1dm, R2_diva_1dm, MD_diva_1dm, N_diva_1dm), (RMSD_woa_1dm, R2_woa_1dm, MD_woa_1dm, N_woa_1dm), \
            (RMSD_ibi_dom, R2_ibi_dom, MD_ibi_dom, N_ibi_dom), (RMSD_cs1k_dom, R2_cs1k_dom, MD_cs1k_dom, N_cs1k_dom), \
            (RMSD_ibi_bry, R2_ibi_bry, MD_ibi_bry, N_ibi_bry), (RMSD_niva_bry, R2_niva_bry, MD_niva_bry, N_niva_bry), \
            (RMSD_woa_bry, R2_woa_bry, MD_woa_bry, N_woa_bry), (RMSD_diva_bry, R2_diva_bry, MD_diva_bry, N_diva_bry)

    def get_eans(obs, mod):
        N = len(obs)
        RMSD = (np.sum((obs - mod) ** 2.) / len(obs)) ** 0.5
        R2 = np.corrcoef(obs, mod)[0, 1]
        MD = np.mean(obs) - np.mean(mod)
        return RMSD, R2, MD, N

    # Scatters
    def scatr5(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa,
               qfilt, nutstr, descr):
        # 5 scatter
        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('In-situ')
        plt.ylabel('model/climatology')
        lims1 = [np.min([obs[qfilt].min(), mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                         clim_diva[qfilt].min(), clim_woa[qfilt].min()]),
                 np.max([obs[qfilt].max(), mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                         clim_diva[qfilt].max(), clim_woa[qfilt].max()])]
        (ax1, ax2, ax3), (ax4, ax5, ax6) = axs
        ax1.scatter(obs[qfilt], mod_ibim[qfilt])
        ax1.grid(visible=True, which='major', color='grey', linestyle='--')
        ax1.set_xlim(lims1)
        ax1.set_ylim(lims1)
        ax1.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        ax2.scatter(obs[qfilt], mod_cs1k[qfilt])
        ax2.grid(visible=True, which='major', color='grey', linestyle='--')
        ax2.set_xlim(lims1)
        ax2.set_ylim(lims1)
        ax2.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        ax3.scatter(obs[qfilt], mod_niva[qfilt])
        ax3.grid(visible=True, which='major', color='grey', linestyle='--')
        ax3.set_xlim(lims1)
        ax3.set_ylim(lims1)
        ax3.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        ax4.scatter(obs[qfilt], clim_diva[qfilt])
        ax4.grid(visible=True, which='major', color='grey', linestyle='--')
        ax4.set_xlim(lims1)
        ax4.set_ylim(lims1)
        ax4.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        ax6.scatter(obs[qfilt], clim_woa[qfilt])
        ax6.grid(visible=True, which='major', color='grey', linestyle='--')
        ax6.set_xlim(lims1)
        ax6.set_ylim(lims1)
        ax6.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        fig.suptitle(nutstr + ' ' + descr + ' In-situ Vs. models/clims.')
        ax1.set_title('IBI', fontsize=10)
        ax2.set_title('CS1KM', fontsize=10)
        ax3.set_title('NIVA', fontsize=10)
        ax4.set_title('DIVA', fontsize=10)
        ax6.set_title('WOA', fontsize=10)
        fig.tight_layout()
        # plt.show()
        plt_tit = fdir + descr + '_' + nutstr + '_scatcomp.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    # Scatters
    def biscatr(fdir, obs, dat1, dat2, lab1, lab2, qfilt, nutstr, nutunits, descr, climormod):
        # 2 scatter
        fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('In-situ, (' + nutunits + ')')
        plt.ylabel(climormod + ', (' + nutunits + ')')
        (ax1, ax2) = axs
        ax1.scatter(obs[qfilt], dat1[qfilt])
        ax1.set_xlim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max()])])
        ax1.set_ylim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max()])])
        lims1 = [np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min()]),
                 np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max()])]
        ax1.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        ax2.scatter(obs[qfilt], dat2[qfilt])
        ax2.set_xlim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max()])])
        ax2.set_ylim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max()])])
        lims2 = [np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min()]),
                 np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max()])]
        ax2.plot(lims2, lims2, 'k--', alpha=0.75, zorder=0)
        fig.suptitle(nutstr + ' ' + descr + ' In-situ Vs. ' + climormod)
        ax1.set_title(lab1, fontsize=10)
        ax2.set_title(lab2, fontsize=10)
        fig.tight_layout()
        # plt.show()
        plt_tit = fdir + descr + '_' + nutstr + lab1 + '_' + lab2 + '_biscatr.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def prod_biases(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa, wodlon, wodlat,
                    nutstr, nutunit, nutcmp, qfilt, cs1kfilt):

        # Spatial map of bias
        # nutstr = 'NO3'
        # nutunit = 'mmol/m3'
        # nutcmp = 'tab20c'  # 'bwr_r', 'gist_ncar', 'RdYlBu'
        csidx = np.intersect1d(qfilt, cs1kfilt)
        # cm = plt.cm.get_cmap(nutcmp)
        cm = plt.cm.get_cmap('seismic')
        nmin = np.min([(obs[qfilt]-mod_ibim[qfilt]).min(),
                       (obs[csidx]-mod_cs1k[csidx]).min(),
                       (obs[qfilt]-mod_niva[qfilt]).min(),
                       (obs[qfilt]-clim_diva[qfilt]).min(),
                       (obs[qfilt]-clim_woa[qfilt]).min()])
        nmax = np.max([(obs[qfilt]-mod_ibim[qfilt]).max(),
                       (obs[csidx]-mod_cs1k[csidx]).max(),
                       (obs[qfilt]-mod_niva[qfilt]).max(),
                       (obs[qfilt]-clim_diva[qfilt]).max(),
                       (obs[qfilt]-clim_woa[qfilt]).max()])
        # if nmin < 0:
        #     bigger = np.max([abs(nmin), abs(nmax)])
        #     if bigger > 10:
        #         bigger2 = np.ceil(bigger/10)*10
        #     elif (bigger < 10) & (bigger > 1):
        #         bigger2 = np.ceil(bigger)
        #     else:
        #         bigger2 = bigger
        #     nmin = bigger2*(-1)
        #     nmax = bigger2
        # else:
        #     nmin = 0
        #     if nmax > 10:
        #         nmax2 = np.ceil(nmax/10)*10
        #     elif (nmax < 10) & (nmax > 1):
        #         nmax2 = np.ceil(nmax)
        #     else:
        #         nmax2 = nmax
        #     nmax = nmax2

        if nmax > 10:
            nmax2 = np.ceil(nmax/10)*10
        elif (nmax < 10) & (nmax > 1):
            nmax2 = np.ceil(nmax)
        else:
            nmax2 = nmax
        nmax = nmax2
        if nmin < 0:
            if abs(nmin) > 10:
                nmin2 = np.ceil(abs(nmin)/10)*10
            elif (abs(nmin) < 10) & (abs(nmin) > 1):
                nmin2 = np.ceil(abs(nmin))
            else:
                nmin2 = nmin
        else:
            if nmin > 10:
                nmin2 = np.floor(nmin/10)*10
            elif (nmin < 10) & (nmin > 1):
                nmin2 = np.floor(nmin)
            else:
                nmin2 = nmin
        nmin = nmin2*(-1)

        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(13, 7.5))
        (ax1, ax2, ax3), (ax4, ax5, ax6) = axs
        ax2.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - mod_ibim[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax3.scatter(wodlon[csidx], wodlat[csidx], c=obs[csidx] - mod_cs1k[csidx],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax4.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - mod_niva[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax5.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - clim_diva[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        im = ax6.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - clim_woa[qfilt],
                         cmap=cm, vmin=nmin, vmax=nmax, s=5)
        fig.suptitle(nutstr + ' ' + 'Bias,' + ' ' + nutunit, fontsize=10)

        ax2.title.set_text('IBI')
        ax3.title.set_text('CS1K')
        ax4.title.set_text('NIVA')
        ax5.title.set_text('DIVA')
        ax6.title.set_text('WOA')
        cbar = plt.colorbar(im, ax=axs)
        # fig.tight_layout()
        plt_tit = fdir + nutstr + '_spat_biases.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def obs_prod_biases_x9(fdir, obs, mod_ibim, mod_cs1k, mod_niva, wodlon, wodlat,
                           nutstr, nutunit, nutcmp, qfilt, cs1kfilt):

        # Spatial map of bias
        # nutstr = 'NO3'
        # nutunit = 'mmol/m3'
        # nutcmp = 'tab20c'  # 'bwr_r', 'gist_ncar', 'RdYlBu'

        csidx = np.intersect1d(qfilt, cs1kfilt)

        # cm = plt.cm.get_cmap(nutcmp)
        # cm = plt.cm.get_cmap('seismic')

        # nmin_bias = np.min([(mod_ibim[csidx]-obs[csidx]).min(),
        #                     (mod_cs1k[csidx]-obs[csidx]).min(),
        #                     (mod_niva[csidx]-obs[csidx]).min()])
        # nmax_bias = np.max([(mod_ibim[csidx]-obs[csidx]).max(),
        #                     (mod_cs1k[csidx]-obs[csidx]).max(),
        #                     (mod_niva[csidx]-obs[csidx]).max()])

        nmin_bias = np.percentile([(mod_ibim[csidx]-obs[csidx]),
                                   (mod_cs1k[csidx]-obs[csidx]),
                                   (mod_niva[csidx]-obs[csidx])], 1)
        nmax_bias = np.percentile([(mod_ibim[csidx]-obs[csidx]),
                                   (mod_cs1k[csidx]-obs[csidx]),
                                   (mod_niva[csidx]-obs[csidx])], 99)

        if nmax_bias > 10:
            nmax2 = np.ceil(nmax_bias/10)*10
        elif (nmax_bias < 10) & (nmax_bias > 1):
            nmax2 = np.ceil(nmax_bias)
        else:
            nmax2 = nmax_bias
        nmax_bias = nmax2
        if nmin_bias < 0:
            if abs(nmin_bias) > 10:
                nmin2 = np.ceil(abs(nmin_bias)/10)*10
            elif (abs(nmin_bias) < 10) & (abs(nmin_bias) > 1):
                nmin2 = np.ceil(abs(nmin_bias))
            else:
                nmin2 = abs(nmin_bias)
            nmin_bias = nmin2*(-1)
        else:
            if nmin_bias > 10:
                nmin2 = np.floor(nmin_bias/10)*10
            elif (nmin_bias < 10) & (nmin_bias > 1):
                nmin2 = np.floor(nmin_bias)
            else:
                nmin2 = nmin_bias
            nmin_bias = nmin2

        if np.sign(nmin_bias) != np.sign(nmax_bias):
            bias_ext = np.max([abs(nmin_bias), abs(nmax_bias)])
            nmin_bias = bias_ext*(-1)
            nmax_bias = bias_ext
            cm_bias = plt.cm.get_cmap('seismic')
        else:
            if (np.sign(nmin_bias) == 1) & (np.sign(nmax_bias) == 1):
                cm_bias = plt.cm.get_cmap('Reds')
            elif (np.sign(nmin_bias) == -1) & (np.sign(nmax_bias) == -1):
                cm_bias = plt.cm.get_cmap('Blues_r')

        # nmin = np.min([(obs[csidx]).min(),
        #                (mod_ibim[csidx]).min(),
        #                (mod_cs1k[csidx]).min(),
        #                (mod_niva[csidx]).min()])
        # nmax = np.max([(obs[csidx]).max(),
        #                (mod_ibim[csidx]).max(),
        #                (mod_cs1k[csidx]).max(),
        #                (mod_niva[csidx]).max()])

        nmin = np.percentile([(obs[csidx]),
                              (mod_ibim[csidx]),
                              (mod_cs1k[csidx]),
                              (mod_niva[csidx])], 1)
        nmax = np.percentile([(obs[csidx]),
                              (mod_ibim[csidx]),
                              (mod_cs1k[csidx]),
                              (mod_niva[csidx])], 99)

        if nmin/nmax < 0.1:
            nmin2 = 0
        else:
            if nmin > 1500:
                nmin2 = np.floor(nmin/50)*50
            elif nmin > 120:
                nmin2 = np.floor(nmin/10)*10
            else:
                nmin2 = np.floor(nmin)
        nmin = nmin2

        if nmax > 10:
            nmax2 = np.ceil(nmax/10)*10
        elif (nmax < 10) & (nmax > 1):
            nmax2 = np.ceil(nmax)
        else:
            nmax2 = nmax
        nmax = nmax2

        # m.drawcoastlines()
        # m.shadedrelief()
        # m.bluemarble()
        # m.drawlsmask(land_color='darkseagreen', ocean_color='gainsboro')

        fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True,
                                figsize=(13, 11.25))
        (ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axs

        ax1.title.set_text('Obs')
        ax1.set_facecolor('gainsboro')
        m1 = Basemap(ax=ax1, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m1.fillcontinents(color='darkseagreen')
        x1, y1 = m1(wodlon[csidx], wodlat[csidx])
        m1.scatter(x1, y1, c=obs[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m1.drawparallels(np.arange(49., 53., 1.), labels=[1, 0, 0, 0], linewidth=0.25)
        m1.drawmeridians(np.arange(-10.75, -5.75, 1.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax2.title.set_text('Obs')
        ax2.set_facecolor('gainsboro')
        m2 = Basemap(ax=ax2, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m2.fillcontinents(color='darkseagreen')
        x2, y2 = m2(wodlon[csidx], wodlat[csidx])
        m2.scatter(x2, y2, c=obs[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m2.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m2.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax3.title.set_text('Obs')
        ax3.set_facecolor('gainsboro')
        m3 = Basemap(ax=ax3, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m3.fillcontinents(color='darkseagreen')
        x3, y3 = m3(wodlon[csidx], wodlat[csidx])
        m3.scatter(x3, y3, c=obs[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m3.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m3.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax4.title.set_text('IBI')
        ax4.set_facecolor('gainsboro')
        m4 = Basemap(ax=ax4, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m4.fillcontinents(color='darkseagreen')
        x4, y4 = m4(wodlon[csidx], wodlat[csidx])
        m4.scatter(x4, y4, c=mod_ibim[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m4.drawparallels(np.arange(49., 53., 1.), labels=[1, 0, 0, 0], linewidth=0.25)
        m4.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax5.title.set_text('NIVA')
        ax5.set_facecolor('gainsboro')
        m5 = Basemap(ax=ax5, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m5.fillcontinents(color='darkseagreen')
        x5, y5 = m5(wodlon[csidx], wodlat[csidx])
        m5.scatter(x5, y5, c=mod_niva[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m5.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m5.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax6.title.set_text('CS1K')
        ax6.set_facecolor('gainsboro')
        m6 = Basemap(ax=ax6, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m6.fillcontinents(color='darkseagreen')
        x6, y6 = m6(wodlon[csidx], wodlat[csidx])
        im1 = m6.scatter(x6, y6, c=mod_niva[csidx],
                         cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m6.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m6.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax7.title.set_text('IBI - Obs')
        ax7.set_facecolor('gainsboro')
        m7 = Basemap(ax=ax7, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m7.fillcontinents(color='darkseagreen')
        x7, y7 = m7(wodlon[csidx], wodlat[csidx])
        m7.scatter(x7, y7, c=mod_ibim[csidx] - obs[csidx],
                   cmap=cm_bias, vmin=nmin_bias, vmax=nmax_bias, s=5)
        m7.drawparallels(np.arange(49., 53., 1.), labels=[1, 0, 0, 0], linewidth=0.25)
        m7.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 1], linewidth=0.25)

        ax8.title.set_text('NIVA - Obs')
        ax8.set_facecolor('gainsboro')
        m8 = Basemap(ax=ax8, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m8.fillcontinents(color='darkseagreen')
        x8, y8 = m8(wodlon[csidx], wodlat[csidx])
        m8.scatter(x8, y8, c=mod_niva[csidx] - obs[csidx],
                   cmap=cm_bias, vmin=nmin_bias, vmax=nmax_bias, s=5)
        m8.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m8.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 1], linewidth=0.25)

        ax9.title.set_text('CS1K - Obs')
        ax9.set_facecolor('gainsboro')
        m9 = Basemap(ax=ax9, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m9.fillcontinents(color='darkseagreen')
        x9, y9 = m9(wodlon[csidx], wodlat[csidx])
        im2 = m9.scatter(x9, y9, c=mod_cs1k[csidx] - obs[csidx],
                         cmap=cm_bias, vmin=nmin_bias, vmax=nmax_bias, s=5)
        m9.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m9.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 1], linewidth=0.25)

        fig.suptitle(nutstr + ' ' + nutunit)

        cbar_upper = fig.colorbar(im1, ax=axs[:2, :], location='right')
        cbar_lower = fig.colorbar(im2, ax=axs[2, :], location='right', aspect=9)

        plt.rcParams['figure.constrained_layout.use'] = True
        plt_tit = fdir + nutstr + '_spat_modobs_bias9.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def obsVsprods(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_diva, clim_woa, wodlon, wodlat,
                   nutstr, nutunit, nutcmp, qfilt, cs1kfilt):
        # cm = plt.cm.get_cmap('bwr_r')
        # cm = plt.cm.get_cmap('gist_ncar')
        # cm = plt.cm.get_cmap('RdYlBu')

        csidx = np.intersect1d(qfilt, cs1kfilt)

        # cm = plt.cm.get_cmap(nutcmp)
        cm = cmo.matter

        nmin = np.min([(obs[qfilt]).min(),
                       (mod_ibim[qfilt]).min(),
                       (mod_cs1k[csidx]).min(),
                       (mod_niva[qfilt]).min(),
                       (clim_diva[qfilt]).min(),
                       (clim_woa[qfilt]).min()])
        nmax = np.max([(obs[qfilt]).max(),
                       (mod_ibim[qfilt]).max(),
                       (mod_cs1k[csidx]).max(),
                       (mod_niva[qfilt]).max(),
                       (clim_diva[qfilt]).max(),
                       (clim_woa[qfilt]).max()])

        if nmax > 10:
            nmax2 = np.ceil(nmax/10)*10
        elif (nmax < 10) & (nmax > 1):
            nmax2 = np.ceil(nmax)
        else:
            nmax2 = nmax
        nmax = nmax2

        # if nmin > 10:
        #     nmin2 = np.floor(nmin/10)*10
        # elif (nmin < 10) & (nmin > 1):
        #     nmin2 = np.floor(nmin)
        # else:
        #     nmin2 = nmin
        # nmin = nmin2

        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(13, 7.5))
        (ax1, ax2, ax3), (ax4, ax5, ax6) = axs
        ax1.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax2.scatter(wodlon[qfilt], wodlat[qfilt], c=mod_ibim[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax3.scatter(wodlon[csidx], wodlat[csidx], c=mod_cs1k[csidx],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax4.scatter(wodlon[qfilt], wodlat[qfilt], c=mod_niva[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax5.scatter(wodlon[qfilt], wodlat[qfilt], c=clim_diva[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        im = ax6.scatter(wodlon[qfilt], wodlat[qfilt], c=clim_woa[qfilt],
                         cmap=cm, vmin=nmin, vmax=nmax, s=5)
        fig.suptitle(nutstr + ', ' + nutunit, fontsize=10)
        ax1.title.set_text('Obs')
        ax2.title.set_text('IBI')
        ax3.title.set_text('CS1K')
        ax4.title.set_text('NIVA')
        ax5.title.set_text('DIVA')
        ax6.title.set_text('WOA')
        cbar = plt.colorbar(im, ax=axs)
        # fig.tight_layout()
        plt_tit = fdir + nutstr + '_spat.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def lt_tseries(fdir, otimes, purp, obs, dat1, dat2, dat3, qflt1, qflt2, qflt3, desc1, desc2, desc3,
                   nutstr, nutunit, aim, climormod):
        # Full timeseries x3 stacked
        fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.ylabel(nutstr + ', ' + nutunit)
        if purp == 'bias':
            axs[0].plot_date(otimes[qflt1], obs[qflt1]-dat1[qflt1])
            axs[0].axhline(y=0, color="black", linestyle="--")
            axs[1].plot_date(otimes[qflt2], obs[qflt2]-dat2[qflt2])
            axs[1].axhline(y=0, color="black", linestyle="--")
            axs[2].plot_date(tstamp[qflt3], obs[qflt3]-dat3[qflt3])
            axs[2].axhline(y=0, color="black", linestyle="--")
            axs[0].title.set_text(desc1)
            axs[1].title.set_text(desc2)
            axs[2].title.set_text(desc3)
            fig.suptitle(nutstr + ', bias timeseries')
        else:  # data
            axs[0].plot_date(otimes[qflt1], dat1[qflt1], label=climormod, markersize=7, marker='s')
            axs[0].plot_date(otimes[qflt1], obs[qflt1], label='Obs', markersize=5, marker='o')
            axs[1].plot_date(otimes[qflt2], dat2[qflt2], label=climormod, markersize=7, marker='s')
            axs[1].plot_date(otimes[qflt2], obs[qflt2], label='Obs', markersize=5, marker='o')
            axs[2].plot_date(tstamp[qflt3], dat3[qflt3], label=climormod, markersize=7, marker='s')
            axs[2].plot_date(tstamp[qflt3], obs[qflt3], label='Obs', markersize=5, marker='o')
            axs[0].title.set_text(desc1)
            axs[1].title.set_text(desc2)
            axs[2].title.set_text(desc3)
            fig.suptitle(nutstr + ', timeseries of comparison by month')
            handles, labels = axs[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', ncol=2)
        axs[0].title.set_text(desc1)
        axs[1].title.set_text(desc2)
        axs[2].title.set_text(desc3)
        plt_tit = fdir + nutstr + '_' + aim + '_tseries.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def t_clim(fdir, tmth, purp, obs, dat1, dat2, dat3, qflt1, qflt2, qflt3, desc1, desc2, desc3,
               nutstr, nutunit, aim, climormod):
        # Climatology x3 stacked
        fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Month')
        plt.ylabel(nutstr + ', ' + nutunit)
        if purp == 'bias':
            axs[0].scatter(tmth[qflt1], obs[qflt1]-dat1[qflt1])
            axs[0].axhline(y=0, color="black", linestyle="--")
            axs[0].title.set_text(desc1)
            axs[1].scatter(tmth[qflt2], obs[qflt2]-dat2[qflt2])
            axs[1].axhline(y=0, color="black", linestyle="--")
            axs[1].title.set_text(desc2)
            axs[2].scatter(tmth[qflt3], obs[qflt3]-dat3[qflt3])
            axs[2].axhline(y=0, color="black", linestyle="--")
            axs[2].title.set_text(desc3)
            fig.suptitle(nutstr + ', bias by month')
        else:  # 'data'
            axs[0].scatter(tmth[qflt1], dat1[qflt1], label=climormod, s=40, marker='s')
            axs[0].scatter(tmth[qflt1], obs[qflt1], label='Obs', s=20, marker='o')
            axs[0].title.set_text(desc1)
            axs[1].scatter(tmth[qflt2], dat2[qflt2], label=climormod, s=40, marker='s')
            axs[1].scatter(tmth[qflt2], obs[qflt2], label='Obs', s=20, marker='o')
            axs[1].title.set_text(desc2)
            axs[2].scatter(tmth[qflt3], dat3[qflt3], label=climormod, s=40, marker='s')
            axs[2].scatter(tmth[qflt3], obs[qflt3], label='Obs', s=20, marker='o')
            axs[2].title.set_text(desc3)
            fig.suptitle(nutstr + ', comparison by month')
            handles, labels = axs[0].get_legend_handles_labels()
            fig.legend(handles, labels, loc='lower center', ncol=2)
        plt_tit = fdir + nutstr + '_' + aim + '_tclim.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    calc_dists = 0
    dist_file = crocodir + 'wod_dist.npz'
    if calc_dists == 1:

        # Read in variable from IBI and do contour to get coastal outline

        # 3 subplots, individual data points aggregated, measure, model and difference
        # separate bar for each plot
        ibigrid = '/media/dskone/CELTIC/CMEMS_IBI/CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_19930101_19930101_R20201201_RE01.nc'
        ibg = netcdf(ibigrid, 'r')
        ibgo = ibg.variables
        ibisal = np.array(ibgo['so'][:])
        ibilat = np.array(ibgo['latitude'][:])
        ibilon = np.array(ibgo['longitude'][:])
        landsal = ibisal[0, 0, :, :]
        landid = np.argwhere((landsal < 0))
        landlat = ibilat[landid[:, 0]]
        landlon = ibilon[landid[:, 1]]

        wod_ll_id1 = np.zeros_like(wodlat, dtype=int)
        wod_dists = np.zeros_like(wodlat)
        lls = np.unique([wodlon, wodlat], axis=1)
        wodkm = np.zeros(lls.shape[1], dtype=lls.dtype)
        dists = np.zeros_like(landlat)
        for p in range(0, lls.shape[1]):
            ll_ind = np.argwhere((wodlon == lls[0, p]) & (wodlat == lls[1, p]))[:, 0]
            for l in range(0, len(landlat)):
                dists[l] = calc_distance(lls[1, p], lls[0, p], landlat[l], landlon[l])
            wodkm[p] = dists.min()
            wod_dists[ll_ind] = wodkm[p]
        np.savez(dist_file,
                 wod_dists=wod_dists)
    else:
        dist_opened = np.load(dist_file)
        wod_dists = dist_opened['wod_dists']

    (RMSD_ibim_t_NO3, R2_ibim_t_NO3, MD_ibim_t_NO3, N_ibim_t_NO3), \
    (RMSD_cs1k_t_NO3, R2_cs1k_t_NO3, MD_cs1k_t_NO3, N_cs1k_t_NO3), \
    (RMSD_niva_t_NO3, R2_niva_t_NO3, MD_niva_t_NO3, N_niva_t_NO3), \
    (RMSD_diva_t_NO3, R2_diva_t_NO3, MD_diva_t_NO3, N_diva_t_NO3), \
    (RMSD_woa_t_NO3, R2_woa_t_NO3, MD_woa_t_NO3, N_woa_t_NO3), \
    (RMSD_ibim_1dm_NO3, R2_ibim_1dm_NO3, MD_ibim_1dm_NO3, N_ibim_1dm_NO3), \
    (RMSD_niva_1dm_NO3, R2_niva_1dm_NO3, MD_niva_1dm_NO3, N_niva_1dm_NO3), \
    (RMSD_diva_1dm_NO3, R2_diva_1dm_NO3, MD_diva_1dm_NO3, N_diva_1dm_NO3), \
    (RMSD_woa_1dm_NO3, R2_woa_1dm_NO3, MD_woa_1dm_NO3, N_woa_1dm_NO3), \
    (RMSD_ibi_dom_NO3, R2_ibi_dom_NO3, MD_ibi_dom_NO3, N_ibi_dom_NO3), \
    (RMSD_cs1k_dom_NO3, R2_cs1k_dom_NO3, MD_cs1k_dom_NO3, N_cs1k_dom_NO3), \
    (RMSD_ibi_bry_NO3, R2_ibi_bry_NO3, MD_ibi_bry_NO3, N_ibi_bry_NO3), \
    (RMSD_niva_bry_NO3, R2_niva_bry_NO3, MD_niva_bry_NO3, N_niva_bry_NO3), \
    (RMSD_woa_bry_NO3, R2_woa_bry_NO3, MD_woa_bry_NO3, N_woa_bry_NO3), \
    (RMSD_diva_bry_NO3, R2_diva_bry_NO3, MD_diva_bry_NO3, N_diva_bry_NO3) = \
        nut_bias_plt_stats(crocodir, wodnit, wodnit_ibim, wodnit_cs1k, wodnit_niva,
                           wodnit_diva*1.025, wodnit_woa*1.025, z,
                           wod_year, wod_month, wod_day, tstamp, wodlon, wodlat, wod_dists, 'NO3')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    no3writer = pd.ExcelWriter(crocodir + 'NO3_stats.xls')

    no3stats_t = pd.DataFrame({'RMSD': [RMSD_ibim_t_NO3, RMSD_cs1k_t_NO3, RMSD_niva_t_NO3,
                                        RMSD_diva_t_NO3, RMSD_woa_t_NO3],
                               'R2': [R2_ibim_t_NO3, R2_cs1k_t_NO3, R2_niva_t_NO3, R2_diva_t_NO3, R2_woa_t_NO3],
                               'MD': [MD_ibim_t_NO3, MD_cs1k_t_NO3, MD_niva_t_NO3, MD_diva_t_NO3, MD_woa_t_NO3],
                               'N': [N_ibim_t_NO3, N_cs1k_t_NO3, N_niva_t_NO3, N_diva_t_NO3, N_woa_t_NO3]})
    no3stats_t.to_excel(no3writer, sheet_name='total_at_cs1k')

    no3stats_1dm = pd.DataFrame({'RMSD': [RMSD_ibim_1dm_NO3, RMSD_niva_1dm_NO3,
                                          RMSD_diva_1dm_NO3, RMSD_woa_1dm_NO3],
                                 'R2': [R2_ibim_1dm_NO3, R2_niva_1dm_NO3, R2_diva_1dm_NO3, R2_woa_1dm_NO3],
                                 'MD': [MD_ibim_1dm_NO3, MD_niva_1dm_NO3, MD_diva_1dm_NO3, MD_woa_1dm_NO3],
                                 'N': [N_ibim_1dm_NO3, N_niva_1dm_NO3, N_diva_1dm_NO3, N_woa_1dm_NO3]})
    no3stats_1dm.to_excel(no3writer, sheet_name='total_1dm')

    no3stats_dom = pd.DataFrame({'RMSD': [RMSD_ibi_dom_NO3, RMSD_cs1k_dom_NO3],
                                'R2': [R2_ibi_dom_NO3, R2_cs1k_dom_NO3],
                                'MD': [MD_ibi_dom_NO3, MD_cs1k_dom_NO3],
                                'N': [N_ibi_dom_NO3, N_cs1k_dom_NO3]})
    no3stats_dom.to_excel(no3writer, sheet_name='domain')

    no3stats_bry = pd.DataFrame({'RMSD': [RMSD_ibi_bry_NO3, RMSD_niva_bry_NO3, RMSD_diva_bry_NO3, RMSD_woa_bry_NO3],
                                'R2': [R2_ibi_bry_NO3, R2_niva_bry_NO3, R2_diva_bry_NO3, R2_woa_bry_NO3],
                                'MD': [MD_ibi_bry_NO3, MD_niva_bry_NO3, MD_diva_bry_NO3, MD_woa_bry_NO3],
                                'N': [N_ibi_bry_NO3, N_niva_bry_NO3, N_diva_bry_NO3, N_woa_bry_NO3]})
    no3stats_bry.to_excel(no3writer, sheet_name='boundary')

    # Close the Pandas Excel writer and output the Excel file.
    no3writer.close()

    (RMSD_ibim_t_PO4, R2_ibim_t_PO4, MD_ibim_t_PO4, N_ibim_t_PO4), \
    (RMSD_cs1k_t_PO4, R2_cs1k_t_PO4, MD_cs1k_t_PO4, N_cs1k_t_PO4), \
    (RMSD_niva_t_PO4, R2_niva_t_PO4, MD_niva_t_PO4, N_niva_t_PO4), \
    (RMSD_diva_t_PO4, R2_diva_t_PO4, MD_diva_t_PO4, N_diva_t_PO4), \
    (RMSD_woa_t_PO4, R2_woa_t_PO4, MD_woa_t_PO4, N_woa_t_PO4), \
    (RMSD_ibim_1dm_PO4, R2_ibim_1dm_PO4, MD_ibim_1dm_PO4, N_ibim_1dm_PO4), \
    (RMSD_niva_1dm_PO4, R2_niva_1dm_PO4, MD_niva_1dm_PO4, N_niva_1dm_PO4), \
    (RMSD_diva_1dm_PO4, R2_diva_1dm_PO4, MD_diva_1dm_PO4, N_diva_1dm_PO4), \
    (RMSD_woa_1dm_PO4, R2_woa_1dm_PO4, MD_woa_1dm_PO4, N_woa_1dm_PO4), \
    (RMSD_ibi_dom_PO4, R2_ibi_dom_PO4, MD_ibi_dom_PO4, N_ibi_dom_PO4), \
    (RMSD_cs1k_dom_PO4, R2_cs1k_dom_PO4, MD_cs1k_dom_PO4, N_cs1k_dom_PO4), \
    (RMSD_ibi_bry_PO4, R2_ibi_bry_PO4, MD_ibi_bry_PO4, N_ibi_bry_PO4), \
    (RMSD_niva_bry_PO4, R2_niva_bry_PO4, MD_niva_bry_PO4, N_niva_bry_PO4), \
    (RMSD_woa_bry_PO4, R2_woa_bry_PO4, MD_woa_bry_PO4, N_woa_bry_PO4), \
    (RMSD_diva_bry_PO4, R2_diva_bry_PO4, MD_diva_bry_PO4, N_diva_bry_PO4) = \
        nut_bias_plt_stats(crocodir, wodpho, wodpho_ibim, wodpho_cs1k, wodpho_niva,
                           wodpho_diva*1.025, wodpho_woa*1.025, z,
                           wod_year, wod_month, wod_day, tstamp, wodlon, wodlat, wod_dists, 'PO4')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    po4writer = pd.ExcelWriter(crocodir + 'PO4_stats.xls')

    po4stats_t = pd.DataFrame({'RMSD': [RMSD_ibim_t_PO4, RMSD_cs1k_t_PO4, RMSD_niva_t_PO4,
                                        RMSD_diva_t_PO4, RMSD_woa_t_PO4],
                              'R2': [R2_ibim_t_PO4, R2_cs1k_t_PO4, R2_niva_t_PO4, R2_diva_t_PO4, R2_woa_t_PO4],
                              'MD': [MD_ibim_t_PO4, MD_cs1k_t_PO4, MD_niva_t_PO4, MD_diva_t_PO4, MD_woa_t_PO4],
                              'N': [N_ibim_t_PO4, N_cs1k_t_PO4, N_niva_t_PO4, N_diva_t_PO4, N_woa_t_PO4]})
    po4stats_t.to_excel(po4writer, sheet_name='total_at_cs1k')

    po4stats_1dm = pd.DataFrame({'RMSD': [RMSD_ibim_1dm_PO4, RMSD_niva_1dm_PO4,
                                          RMSD_diva_1dm_PO4, RMSD_woa_1dm_PO4],
                                 'R2': [R2_ibim_1dm_PO4, R2_niva_1dm_PO4, R2_diva_1dm_PO4, R2_woa_1dm_PO4],
                                 'MD': [MD_ibim_1dm_PO4, MD_niva_1dm_PO4, MD_diva_1dm_PO4, MD_woa_1dm_PO4],
                                 'N': [N_ibim_1dm_PO4, N_niva_1dm_PO4, N_diva_1dm_PO4, N_woa_1dm_PO4]})
    po4stats_1dm.to_excel(po4writer, sheet_name='total_1dm')

    po4stats_dom = pd.DataFrame({'RMSD': [RMSD_ibi_dom_PO4, RMSD_cs1k_dom_PO4],
                                 'R2': [R2_ibi_dom_PO4, R2_cs1k_dom_PO4],
                                 'MD': [MD_ibi_dom_PO4, MD_cs1k_dom_PO4],
                                 'N': [N_ibi_dom_PO4, N_cs1k_dom_PO4]})
    po4stats_dom.to_excel(po4writer, sheet_name='domain')

    po4stats_bry = pd.DataFrame({'RMSD': [RMSD_ibi_bry_PO4, RMSD_niva_bry_PO4, RMSD_diva_bry_PO4, RMSD_woa_bry_PO4],
                                'R2': [R2_ibi_bry_PO4, R2_niva_bry_PO4, R2_diva_bry_PO4, R2_woa_bry_PO4],
                                'MD': [MD_ibi_bry_PO4, MD_niva_bry_PO4, MD_diva_bry_PO4, MD_woa_bry_PO4],
                                'N': [N_ibi_bry_PO4, N_niva_bry_PO4, N_diva_bry_PO4, N_woa_bry_PO4]})
    po4stats_bry.to_excel(po4writer, sheet_name='boundary')

    # Close the Pandas Excel writer and output the Excel file.
    po4writer.close()

    (RMSD_ibim_t_Si, R2_ibim_t_Si, MD_ibim_t_Si, N_ibim_t_Si), \
    (RMSD_cs1k_t_Si, R2_cs1k_t_Si, MD_cs1k_t_Si, N_cs1k_t_Si), \
    (RMSD_niva_t_Si, R2_niva_t_Si, MD_niva_t_Si, N_niva_t_Si), \
    (RMSD_diva_t_Si, R2_diva_t_Si, MD_diva_t_Si, N_diva_t_Si), \
    (RMSD_woa_t_Si, R2_woa_t_Si, MD_woa_t_Si, N_woa_t_Si), \
    (RMSD_ibim_1dm_Si, R2_ibim_1dm_Si, MD_ibim_1dm_Si, N_ibim_1dm_Si), \
    (RMSD_niva_1dm_Si, R2_niva_1dm_Si, MD_niva_1dm_Si, N_niva_1dm_Si), \
    (RMSD_diva_1dm_Si, R2_diva_1dm_Si, MD_diva_1dm_Si, N_diva_1dm_Si), \
    (RMSD_woa_1dm_Si, R2_woa_1dm_Si, MD_woa_1dm_Si, N_woa_1dm_Si), \
    (RMSD_ibi_dom_Si, R2_ibi_dom_Si, MD_ibi_dom_Si, N_ibi_dom_Si), \
    (RMSD_cs1k_dom_Si, R2_cs1k_dom_Si, MD_cs1k_dom_Si, N_cs1k_dom_Si), \
    (RMSD_ibi_bry_Si, R2_ibi_bry_Si, MD_ibi_bry_Si, N_ibi_bry_Si), \
    (RMSD_niva_bry_Si, R2_niva_bry_Si, MD_niva_bry_Si, N_niva_bry_Si), \
    (RMSD_woa_bry_Si, R2_woa_bry_Si, MD_woa_bry_Si, N_woa_bry_Si), \
    (RMSD_diva_bry_Si, R2_diva_bry_Si, MD_diva_bry_Si, N_diva_bry_Si) = \
        nut_bias_plt_stats(crocodir, wodsil, wodsil_ibim, wodsil_cs1k, wodsil_niva,
                           wodsil_diva*1.025, wodsil_woa*1.025, z,
                           wod_year, wod_month, wod_day, tstamp, wodlon, wodlat, wod_dists, 'Si')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    siwriter = pd.ExcelWriter(crocodir + 'Si_stats.xls')

    sistats_t = pd.DataFrame({'RMSD': [RMSD_ibim_t_Si, RMSD_cs1k_t_Si, RMSD_niva_t_Si,
                                       RMSD_diva_t_Si, RMSD_woa_t_Si],
                              'R2': [R2_ibim_t_Si, R2_cs1k_t_Si, R2_niva_t_Si, R2_diva_t_Si, R2_woa_t_Si],
                              'MD': [MD_ibim_t_Si, MD_cs1k_t_Si, MD_niva_t_Si, MD_diva_t_Si, MD_woa_t_Si],
                              'N': [N_ibim_t_Si, N_cs1k_t_Si, N_niva_t_Si, N_diva_t_Si, N_woa_t_Si]})
    sistats_t.to_excel(siwriter, sheet_name='total_at_cs1k')

    sistats_1dm = pd.DataFrame({'RMSD': [RMSD_ibim_1dm_Si, RMSD_niva_1dm_Si,
                                          RMSD_diva_1dm_Si, RMSD_woa_1dm_Si],
                                 'R2': [R2_ibim_1dm_Si, R2_niva_1dm_Si, R2_diva_1dm_Si, R2_woa_1dm_Si],
                                 'MD': [MD_ibim_1dm_Si, MD_niva_1dm_Si, MD_diva_1dm_Si, MD_woa_1dm_Si],
                                 'N': [N_ibim_1dm_Si, N_niva_1dm_Si, N_diva_1dm_Si, N_woa_1dm_Si]})
    sistats_1dm.to_excel(siwriter, sheet_name='total_1dm')

    sistats_dom = pd.DataFrame({'RMSD': [RMSD_ibi_dom_Si, RMSD_cs1k_dom_Si],
                                'R2': [R2_ibi_dom_Si, R2_cs1k_dom_Si],
                                'MD': [MD_ibi_dom_Si, MD_cs1k_dom_Si],
                                'N': [N_ibi_dom_Si, N_cs1k_dom_Si]})
    sistats_dom.to_excel(siwriter, sheet_name='domain')

    sistats_bry = pd.DataFrame({'RMSD': [RMSD_ibi_bry_Si, RMSD_niva_bry_Si, RMSD_diva_bry_Si, RMSD_woa_bry_Si],
                                'R2': [R2_ibi_bry_Si, R2_niva_bry_Si, R2_diva_bry_Si, R2_woa_bry_Si],
                                'MD': [MD_ibi_bry_Si, MD_niva_bry_Si, MD_diva_bry_Si, MD_woa_bry_Si],
                                'N': [N_ibi_bry_Si, N_niva_bry_Si, N_diva_bry_Si, N_woa_bry_Si]})
    sistats_bry.to_excel(siwriter, sheet_name='boundary')

    # Close the Pandas Excel writer and output the Excel file.
    siwriter.close()

    (RMSD_ibim_t_O2, R2_ibim_t_O2, MD_ibim_t_O2, N_ibim_t_O2), \
    (RMSD_cs1k_t_O2, R2_cs1k_t_O2, MD_cs1k_t_O2, N_cs1k_t_O2), \
    (RMSD_niva_t_O2, R2_niva_t_O2, MD_niva_t_O2, N_niva_t_O2), \
    (RMSD_diva_t_O2, R2_diva_t_O2, MD_diva_t_O2, N_diva_t_O2), \
    (RMSD_woa_t_O2, R2_woa_t_O2, MD_woa_t_O2, N_woa_t_O2), \
    (RMSD_ibim_1dm_O2, R2_ibim_1dm_O2, MD_ibim_1dm_O2, N_ibim_1dm_O2), \
    (RMSD_niva_1dm_O2, R2_niva_1dm_O2, MD_niva_1dm_O2, N_niva_1dm_O2), \
    (RMSD_diva_1dm_O2, R2_diva_1dm_O2, MD_diva_1dm_O2, N_diva_1dm_O2), \
    (RMSD_woa_1dm_O2, R2_woa_1dm_O2, MD_woa_1dm_O2, N_woa_1dm_O2), \
    (RMSD_ibi_dom_O2, R2_ibi_dom_O2, MD_ibi_dom_O2, N_ibi_dom_O2), \
    (RMSD_cs1k_dom_O2, R2_cs1k_dom_O2, MD_cs1k_dom_O2, N_cs1k_dom_O2), \
    (RMSD_ibi_bry_O2, R2_ibi_bry_O2, MD_ibi_bry_O2, N_ibi_bry_O2), \
    (RMSD_niva_bry_O2, R2_niva_bry_O2, MD_niva_bry_O2, N_niva_bry_O2), \
    (RMSD_woa_bry_O2, R2_woa_bry_O2, MD_woa_bry_O2, N_woa_bry_O2), \
    (RMSD_diva_bry_O2, R2_diva_bry_O2, MD_diva_bry_O2, N_diva_bry_O2) = \
        nut_bias_plt_stats(crocodir, wodoxy, wodoxy_ibim, wodoxy_cs1k, wodoxy_niva,
                           wodoxy_diva*1.025, wodoxy_woa*1.025, z,
                           wod_year, wod_month, wod_day, tstamp, wodlon, wodlat, wod_dists, 'O2')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    o2writer = pd.ExcelWriter(crocodir + 'O2_stats.xls')

    o2stats_t = pd.DataFrame({'RMSD': [RMSD_ibim_t_O2, RMSD_cs1k_t_O2, RMSD_niva_t_O2,
                                       RMSD_diva_t_O2, RMSD_woa_t_O2],
                              'R2': [R2_ibim_t_O2, R2_cs1k_t_O2, R2_niva_t_O2, R2_diva_t_O2, R2_woa_t_O2],
                              'MD': [MD_ibim_t_O2, MD_cs1k_t_O2, MD_niva_t_O2, MD_diva_t_O2, MD_woa_t_O2],
                              'N': [N_ibim_t_O2, N_cs1k_t_O2, N_niva_t_O2, N_diva_t_O2, N_woa_t_O2]})
    o2stats_t.to_excel(o2writer, sheet_name='total_at_cs1k')

    o2stats_1dm = pd.DataFrame({'RMSD': [RMSD_ibim_1dm_O2, RMSD_niva_1dm_O2,
                                         RMSD_diva_1dm_O2, RMSD_woa_1dm_O2],
                                'R2': [R2_ibim_1dm_O2, R2_niva_1dm_O2, R2_diva_1dm_O2, R2_woa_1dm_O2],
                                'MD': [MD_ibim_1dm_O2, MD_niva_1dm_O2, MD_diva_1dm_O2, MD_woa_1dm_O2],
                                'N': [N_ibim_1dm_O2, N_niva_1dm_O2, N_diva_1dm_O2, N_woa_1dm_O2]})
    o2stats_1dm.to_excel(o2writer, sheet_name='total_1dm')

    o2stats_dom = pd.DataFrame({'RMSD': [RMSD_ibi_dom_O2, RMSD_cs1k_dom_O2],
                                'R2': [R2_ibi_dom_O2, R2_cs1k_dom_O2],
                                'MD': [MD_ibi_dom_O2, MD_cs1k_dom_O2],
                                'N': [N_ibi_dom_O2, N_cs1k_dom_O2]})
    o2stats_dom.to_excel(o2writer, sheet_name='domain')

    o2stats_bry = pd.DataFrame({'RMSD': [RMSD_ibi_bry_O2, RMSD_niva_bry_O2, RMSD_diva_bry_O2, RMSD_woa_bry_O2],
                                'R2': [R2_ibi_bry_O2, R2_niva_bry_O2, R2_diva_bry_O2, R2_woa_bry_O2],
                                'MD': [MD_ibi_bry_O2, MD_niva_bry_O2, MD_diva_bry_O2, MD_woa_bry_O2],
                                'N': [N_ibi_bry_O2, N_niva_bry_O2, N_diva_bry_O2, N_woa_bry_O2]})
    o2stats_bry.to_excel(o2writer, sheet_name='boundary')

    # Close the Pandas Excel writer and output the Excel file.
    o2writer.close()

    print('End')

