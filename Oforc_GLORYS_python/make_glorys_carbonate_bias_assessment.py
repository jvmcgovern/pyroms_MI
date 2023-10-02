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
import numpy.matlib
import PyCO2SYS as pyco2
import time
import pickle
from pyhdf.SD import SD, SDC
import pandas as pd
from interp_Cgrid import *
import cmocean.cm as cmo
import croco_vgrid as vgrd
import croco_glorys as glor
from math import cos, sin, asin, sqrt, radians
from mpl_toolkits.basemap import Basemap

# import cartopy.crs as ccrs
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

mi_processing = 0
glodap_processing = 0
dailyres = 1
monthlyres = 0

if dailyres == 1:
    npzstr = 'daily'
else:
    npzstr = 'monthly'

# hc = 'uncorrected'
hcs = 'bias-corrected'

if hcs == 'uncorrected':
    crocodir = '/media/dskfour/CS1KM_19932021_Uncorr/'
elif hcs == 'bias-corrected':
    # crocodir = '/media/dskfive/CS1KM_19932021_BCd_results/'
    crocodir = '/media/dsksix/CS1KM_19932021_BCd_TNuts_results/'

if mi_processing == 1:

    # Winter nutrients from excel
    # insitufile = '/media/dskone/VAL/Insitu_nutrients.xls'
    insitufile = '/media/dskone/VAL/Insitu_dic_talk_MI_tsorted.xls'
    dfis = pd.read_excel(insitufile)
    isyr = np.array(dfis['Year'])
    ismth = np.array(dfis['Month'])
    isday = np.array(dfis['Day'])
    islat = np.array(dfis['Latitude'])
    islon = np.array(dfis['Longitude'])
    isdep = np.array(dfis['Depth'])
    isdic = np.array(dfis['DIC'])
    istlk = np.array(dfis['TALK'])
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
    isdic = isdic[itrimidx][:, 0]
    istlk = istlk[itrimidx][:, 0]

    # GLODAPv2
    valdir = '/media/dskone/VAL/'

    EMODdir = '/media/dskfour/VAL/EMODNET/INSITU/'
    nut = 'NPSO'
    dic = 'DIC'
    talk = 'TALK'

    z = isdep
    wodlat = islat
    wodlon = islon

    wodtlk = istlk
    wodtlk_ibim = np.zeros_like(wodtlk)
    wodtlk_cs1k = np.zeros_like(wodtlk)
    wodtlk_niva = np.zeros_like(wodtlk)
    wodtlk_brou = np.zeros_like(wodtlk)

    woddic = isdic
    woddic_ibim = np.zeros_like(woddic)
    woddic_cs1k = np.zeros_like(woddic)
    woddic_niva = np.zeros_like(woddic)
    woddic_brou = np.zeros_like(woddic)
    woddic_mobo = np.zeros_like(woddic)

    latids = np.zeros_like(woddic, dtype='int64')
    lonids = np.zeros_like(woddic, dtype='int64')
    latidc = np.zeros_like(woddic, dtype='int64')
    lonidc = np.zeros_like(woddic, dtype='int64')
    latidn = np.zeros_like(woddic, dtype='int64')
    lonidn = np.zeros_like(woddic, dtype='int64')
    latidw = np.zeros_like(woddic, dtype='int64')
    lonidw = np.zeros_like(woddic, dtype='int64')
    latidd = np.zeros_like(woddic, dtype='int64')
    lonidd = np.zeros_like(woddic, dtype='int64')
    land_ibi = np.zeros_like(woddic, dtype='bool')
    land_cs1km = np.zeros_like(woddic, dtype='bool')
    land_noresm = np.zeros_like(woddic, dtype='bool')
    land_brou = np.zeros_like(woddic, dtype='bool')
    land_mobo = np.zeros_like(woddic, dtype='bool')
    wod_day = np.zeros_like(woddic, dtype='int64')
    wod_month = np.zeros_like(woddic, dtype='int64')
    wod_year = np.zeros_like(woddic, dtype='int64')

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
            wdibmo = wdibim.variables

        if dailyres == 1:
            IBI_BIO_file_d = sorted(glob.glob(PISCES24files_dir +
                                              PISCES24_prefixd + str(wod_year[wd]) +
                                              str(wod_month[wd]).zfill(2) + str(wod_day[wd]).zfill(2) +
                                              '_*' + PISCES24_ending))
            wdibim = netcdf(IBI_BIO_file_d[0], 'r')
            wdibmo = wdibim.variables

        # Reanalysis
        glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
        glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
        glorys_mth_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01mav_'
        glorys_ending = '_R*_RE01.nc'

        if monthlyres == 1:
            IBI_PHYmth_files = sorted(glob.glob(glorysfiles_dir +
                                                glorys_mth_prefix + str(wod_year[wd]) +
                                                str(wod_month[wd]).zfill(2) +
                                                '01_*' + glorys_ending))
            wdphym = netcdf(IBI_PHYmth_files[0], 'r')
            wdphmo = wdphym.variables
            ibitmp = np.array(wdphmo['thetao'][:])
            ibisal = np.array(wdphmo['so'][:])

        if dailyres == 1:
            IBI_PHY_files = sorted(glob.glob(glorysfiles_dir +
                                             glorys_prefix + str(wod_year[wd]) +
                                             str(wod_month[wd]).zfill(2) +
                                             str(wod_day[wd]).zfill(2) + '_*' + glorys_ending))
            wdphym = netcdf(IBI_PHY_files[0], 'r')
            wdphmo = wdphym.variables
            ibitmp = np.array(wdphmo['thetao'][:])
            ibisal = np.array(wdphmo['so'][:])
        
        # read latitude and longitude for indexing
        ibiph = np.array(wdibmo['ph'][:])
        ibidic = np.array(wdibmo['dissic'][:])*1000

        for dp in range(0, profs.shape[1]):
            ibiph_smp = ibiph[0, :, latids[dp], lonids[dp]]
            ibidic_smp = ibidic[0, :, latids[dp], lonids[dp]]
            ibitmp_smp = ibitmp[0, :, latids[dp], lonids[dp]]
            ibisal_smp = ibisal[0, :, latids[dp], lonids[dp]]
            if (any(ibidic_smp < 0) == True) & (all(ibidic_smp < 0) == False):
                ibiph_smp[ibiph_smp < 0] = ibiph_smp[ibiph_smp >= 0][-1]
                ibidic_smp[ibidic_smp < 0] = ibidic_smp[ibidic_smp >= 0][-1]
                ibitmp_smp[ibidic_smp < 0] = ibitmp_smp[ibidic_smp >= 0][-1]
                ibisal_smp[ibidic_smp < 0] = ibisal_smp[ibidic_smp >= 0][-1]
                results = pyco2.sys(par1=ibiph_smp, par1_type=3,
                                    par2=ibidic_smp, par2_type=2,
                                    temperature=ibitmp_smp,
                                    salinity=ibisal_smp,
                                    opt_pH_scale=1, opt_k_carbonic=4,
                                    opt_k_bisulfate=1, opt_total_borate=1)
                ibitlk_smp = results["alkalinity"]
            elif all(ibidic_smp < 0) == True:
                land_ibi[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_p = interpolate.interp1d(ibidep, ibitlk_smp, fill_value='extrapolate')
            wodtlk_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(ibidep, ibidic_smp, fill_value='extrapolate')
            woddic_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

        wdibim.close()
        wdphym.close()

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
            # read latitude and longitude for indexing
            cs1kmtlk = np.mean(np.array(wdcs1kmo['TALK'][:]), 0)
            cs1kmdic = np.mean(np.array(wdcs1kmo['DIC'][:]), 0)
            cs1kmz = np.mean(np.array(wdcs1kmo['zeta'][:]), 0)

        if dailyres == 1:
            # read latitude and longitude for indexing
            cs1kmtlk = np.array(wdcs1kmo['TALK'][wod_day[wd] - 1, :, :, :])
            cs1kmdic = np.array(wdcs1kmo['DIC'][wod_day[wd] - 1, :, :, :])
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
            cs1kmtlk_smp = cs1kmtlk[:, latidc[dp], lonidc[dp]]
            cs1kmdic_smp = cs1kmdic[:, latidc[dp], lonidc[dp]]
            cs1kmd_smp = cs1kmd[:, latidc[dp], lonidc[dp]]
            if (any(cs1kmdic_smp == 0) == True) & (all(cs1kmdic_smp == 0) == False):
                cs1kmtlk_smp[cs1kmtlk_smp == 0] = cs1kmtlk_smp[cs1kmtlk_smp > 0][-1]
                cs1kmdic_smp[cs1kmdic_smp == 0] = cs1kmdic_smp[cs1kmdic_smp > 0][-1]
            elif all(cs1kmdic_smp < 0) == True:
                land_cs1km[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_p = interpolate.interp1d(cs1kmd_smp, cs1kmtlk_smp, fill_value='extrapolate')
            wodtlk_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(cs1kmd_smp, cs1kmdic_smp, fill_value='extrapolate')
            woddic_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for NORESM processing
    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']
    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    ncnortlk = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[5] + NORESM_ending, 'r')
    ncnordic = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[4] + NORESM_ending, 'r')
    ncnortlko = ncnortlk.variables
    ncnordico = ncnordic.variables
    noresmlat = np.array(ncnordico['lat'][:])
    noresmlon = np.array(ncnordico['lon'][:])
    depthniva = np.array(ncnordico['depth'][:])
    time_niva = np.array(ncnordico['time'][:])

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

        noresmtlk = np.array(ncnortlko[NIVAvars[5]][:]) * 1000
        noresmdic = np.array(ncnordico[NIVAvars[4]][:]) * 1000

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                Tinin = datetime(wod_year[wd], wod_month[wd], 15)
                Tininx1 = d2i(Tinin, ncnordico['time'], select='exact')
                noresmtlk_smp = noresmtlk[Tininx1, :, latidn[dp], lonidn[dp]]
                noresmdic_smp = noresmdic[Tininx1, :, latidn[dp], lonidn[dp]]
            if dailyres == 1:
                # By day
                if wod_day[wd] != 15:
                    Tinin = datetime(wod_year[wd], wod_month[wd], wod_day[wd])
                    Tininx1 = d2i(Tinin, ncnordico['time'], select='before')
                    Tininx2 = d2i(Tinin, ncnordico['time'], select='after')
                    if wod_day[wd] > 15:
                        tmult = wod_day[wd] - 15
                    else:
                        tmult = wod_day[wd] + 15

                    # Tidx1
                    noresmtlk_1 = noresmtlk[Tininx1, :, latidn[dp], lonidn[dp]]
                    noresmdic_1 = noresmdic[Tininx1, :, latidn[dp], lonidn[dp]]

                    # Tidx2
                    noresmtlk_2 = noresmtlk[Tininx2, :, latidn[dp], lonidn[dp]]
                    noresmdic_2 = noresmdic[Tininx2, :, latidn[dp], lonidn[dp]]
                    
                    # Interpolating to the day
                    noresmtlk_smp = noresmtlk_1 + (((noresmtlk_2 - noresmtlk_1)/30)*tmult)
                    noresmdic_smp = noresmdic_1 + (((noresmdic_2 - noresmdic_1)/30)*tmult)
                else:
                    Tinin = datetime(wod_year[wd], wod_month[wd], 15)
                    Tininx1 = d2i(Tinin, ncnordico['time'], select='exact')
                    noresmtlk_smp = noresmtlk[Tininx1, :, latidn[dp], lonidn[dp]]
                    noresmdic_smp = noresmdic[Tininx1, :, latidn[dp], lonidn[dp]]

            if (any(noresmdic_smp > 10000) == True) & (all(noresmdic_smp > 10000) == False):
                noresmtlk_smp[noresmtlk_smp > 10000] = noresmtlk_smp[noresmtlk_smp < 10000][-1]
                noresmdic_smp[noresmdic_smp > 10000] = noresmdic_smp[noresmdic_smp < 10000][-1]
            elif all(noresmdic_smp > 10000) == True:
                land_noresm[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_p = interpolate.interp1d(depthniva, noresmtlk_smp, fill_value='extrapolate')
            wodtlk_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(depthniva, noresmdic_smp, fill_value='extrapolate')
            woddic_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    ncnortlk.close()
    ncnordic.close()

    # Block for Broullon processing
    valdir = '/media/dskone/VAL/'
    broudir = valdir + 'BroullonClims/'

    # Common woa elements for monthly climatology
    brou_dic = 'TCO2_NNGv2LDEO_climatology_WOA18.nc'
    brou_tlk = 'AT_NNGv2_climatology_WOA18.nc'
    brou_dicv = netcdf(broudir + brou_dic, 'r')
    brou_tlkv = netcdf(broudir + brou_tlk, 'r')
    brou_dicvo = brou_dicv.variables
    brou_tlkvo = brou_tlkv.variables
    brou_lat = np.array(brou_dicvo['latitude'][:])
    brou_lon = np.array(brou_dicvo['longitude'][:])
    brou_depth = np.array(brou_dicvo['depth'][:])
    broudic = np.array(brou_dicvo['TCO2_NNGv2LDEO'][:])
    broutlk = np.array(brou_tlkvo['AT_NNGv2'][:])
    brou_dicv.close()
    brou_tlkv.close()

    woddays, wod_ind, wod_count = np.unique([isyr, ismth, isday], axis=1, return_index=True, return_counts=True)

    for wd in range(0, woddays.shape[1]):
        wod_day[wd] = woddays[2, wd]
        wod_month[wd] = woddays[1, wd]
        wod_year[wd] = woddays[0, wd]

        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidw = np.zeros_like(prof_ind, dtype='int64')
        lonidw = np.zeros_like(prof_ind, dtype='int64')

        for p in range(0, profs.shape[1]):
            latidw[p] = glor.geo_idx(profs[0, p], brou_lat)
            lonidw[p] = glor.geo_idx(profs[1, p], brou_lon)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                broudic_smp = broudic[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                broutlk_smp = broutlk[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
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
                    broudic_1 = broudic[tidx1, :, latidw[dp], lonidw[dp]]
                    broutlk_1 = broutlk[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    broudic_2 = broudic[tidx2, :, latidw[dp], lonidw[dp]]
                    broutlk_2 = broutlk[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    broudic_smp = broudic_1 + (((broudic_2 - broudic_1)/30)*tmult)
                    broutlk_smp = broutlk_1 + (((broutlk_2 - broutlk_1)/30)*tmult)

                elif wod_day[wd] < 15:
                    tmult = wod_day[wd] + 15
                    tidx2 = wod_month[wd] - 1
                    if wod_month[wd] == 1:
                        tidx1 = 11
                    else:
                        tidx1 = wod_month[wd] - 2

                    # Tidx1
                    broudic_1 = broudic[tidx1, :, latidw[dp], lonidw[dp]]
                    broutlk_1 = broutlk[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    broudic_2 = broudic[tidx2, :, latidw[dp], lonidw[dp]]
                    broutlk_2 = broutlk[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    broudic_smp = broudic_1 + (((broudic_2 - broudic_1)/30)*tmult)
                    broutlk_smp = broutlk_1 + (((broutlk_2 - broutlk_1)/30)*tmult)

                else:
                    broudic_smp = broudic[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    broutlk_smp = broutlk[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]

            if (any(np.isnan(broudic_smp)) == True) & (all(np.isnan(broudic_smp)) == False):
                broudic_smp[np.isnan(broudic_smp)] = broudic_smp[np.isnan(broudic_smp) == False][-1]
                broutlk_smp[np.isnan(broutlk_smp)] = broutlk_smp[np.isnan(broutlk_smp) == False][-1]
                #
                wodint_dic = interpolate.interp1d(brou_depth, broudic_smp, fill_value='extrapolate')
                woddic_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_dic(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                #
                wodint_tlk = interpolate.interp1d(brou_depth, broutlk_smp, fill_value='extrapolate')
                wodtlk_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_tlk(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

            elif all(np.isnan(broudic_smp)) == True:
                # land_brou[wod_ind[wd] + dp] = True
                land_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                woddic_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                wodtlk_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
            else:
                #
                wodint_dic = interpolate.interp1d(brou_depth, broudic_smp, fill_value='extrapolate')
                woddic_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_dic(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                #
                wodint_tlk = interpolate.interp1d(brou_depth, broutlk_smp, fill_value='extrapolate')
                wodtlk_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_tlk(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for MOBO_DIC processing
    MOBO_dir = '/media/dskone/VAL/BroullonClims/'
    MOBOfil = 'MOBO_MPI_DIC_2004-2019_Oct2022.nc'

    ncmobodic = netcdf(MOBO_dir + MOBOfil, 'r')
    ncmobodico = ncmobodic.variables
    mobolat = np.array(ncmobodico['lat'][:])
    mobolon = np.array(ncmobodico['lon'][:])
    depthmobo = np.array(ncmobodico['depth'][:])
    time_mobo = np.array(ncmobodico['juld'][:])
    mobodic = np.array(ncmobodico['DIC'][:])
    ncmobodic.close()

    woddays, wod_ind, wod_count = np.unique([isyr, ismth, isday], axis=1, return_index=True, return_counts=True)

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([islat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 islon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)
        latidn = np.zeros_like(prof_ind, dtype='int64')
        lonidn = np.zeros_like(prof_ind, dtype='int64')

        for p in range(0, profs.shape[1]):
            latidn[p] = glor.geo_idx(profs[0, p], mobolat)
            lonidn[p] = glor.geo_idx(profs[1, p], mobolon)

        latidn = latidn.astype(int)
        lonidn = lonidn.astype(int)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                Tyridx = (wod_year[wd] - 2004)*12
                Tyridx = Tyridx - 1
                Tmthidx = wod_month[wd]
                Ttotidx = Tyridx + Tmthidx
                if (Ttotidx >= 0) & (Ttotidx <= mobodic.shape[0] - 1):
                    mobodic_smp = mobodic[Ttotidx, :, latidn[dp], lonidn[dp]]
                    if (any(mobodic_smp < 0) == True) & (all(mobodic_smp < 0) == False):
                        mobodic_smp[mobodic_smp < 0] = mobodic_smp[mobodic_smp > 0][-1]
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                    elif all(mobodic_smp < 0) == True:
                        land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                    #
                    else:
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                else:
                    land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                    woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                    
            if dailyres == 1:
                if wod_day[wd] == 15:
                    Tyridx = (wod_year[wd] - 2004)*12
                    Tyridx = Tyridx - 1
                    Tmthidx = wod_month[wd]
                    Ttotidx = Tyridx + Tmthidx
                    if (Ttotidx >= 0) & (Ttotidx <= mobodic.shape[0] - 1):
                        mobodic_smp = mobodic[Ttotidx, :, latidn[dp], lonidn[dp]]
                if wod_day[wd] < 15:
                    Tyridx = (wod_year[wd] - 2004)*12
                    Tyridx = Tyridx - 1
                    Tmthidx1 = wod_month[wd] - 1
                    Ttotidx1 = Tyridx + Tmthidx1
                    Tmthidx2 = wod_month[wd]
                    Ttotidx2 = Tyridx + Tmthidx2
                    tmult = wod_day[wd] + 15
                    if (Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1):
                        mobodic_1 = mobodic[Ttotidx1, :, latidn[dp], lonidn[dp]]
                    if (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1):
                        mobodic_2 = mobodic[Ttotidx2, :, latidn[dp], lonidn[dp]]
                    if ((Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1) & 
                            (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1)):
                        # Interpolating to the day
                        mobodic_smp = mobodic_1 + (((mobodic_2 - mobodic_1)/30)*tmult)

                if wod_day[wd] > 15:
                    Tyridx = (wod_year[wd] - 2004)*12
                    Tyridx = Tyridx - 1
                    Tmthidx1 = wod_month[wd]
                    Ttotidx1 = Tyridx + Tmthidx1
                    Tmthidx2 = wod_month[wd] + 1
                    Ttotidx2 = Tyridx + Tmthidx2
                    tmult = wod_day[wd] - 15
                    if (Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1):
                        mobodic_1 = mobodic[Ttotidx1, :, latidn[dp], lonidn[dp]]
                    if (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1):
                        mobodic_2 = mobodic[Ttotidx2, :, latidn[dp], lonidn[dp]]
                    if ((Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1) & 
                            (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1)):
                        # Interpolating to the day
                        mobodic_smp = mobodic_1 + (((mobodic_2 - mobodic_1)/30)*tmult)

                    if (any(mobodic_smp < 0) == True) & (all(mobodic_smp < 0) == False):
                        mobodic_smp[mobodic_smp < 0] = mobodic_smp[mobodic_smp > 0][-1]
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                    elif all(mobodic_smp < 0) == True:
                        land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                    #
                    else:
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                else:
                    land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                    woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0

    mi_bias_file = crocodir + 'MIwint_carbon_bias_' + npzstr + '.npz'
    np.savez(mi_bias_file,
             wodlat=wodlat, wodlon=wodlon, z=z,
             wodtlk=wodtlk,
             wodtlk_ibim=wodtlk_ibim, wodtlk_cs1k=wodtlk_cs1k, wodtlk_niva=wodtlk_niva, wodtlk_brou=wodtlk_brou,
             woddic=woddic,
             woddic_ibim=woddic_ibim, woddic_cs1k=woddic_cs1k, woddic_niva=woddic_niva, woddic_brou=woddic_brou,
             woddic_mobo=woddic_mobo,
             latids=latids, lonids=lonids, latidc=latidc, lonidc=lonidc, latidn=latidn, lonidn=lonidn,
             latidw=latidw, lonidw=lonidw, latidd=latidd, lonidd=lonidd,
             land_ibi=land_ibi, land_cs1km=land_cs1km, land_noresm=land_noresm,
             land_brou=land_brou, land_mobo=land_mobo,
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

    z = gdepth
    wodlat = glat
    wodlon = glon

    wodtlk = gtalk
    wodtlk_ibim = np.zeros_like(wodtlk)
    wodtlk_cs1k = np.zeros_like(wodtlk)
    wodtlk_niva = np.zeros_like(wodtlk)
    wodtlk_brou = np.zeros_like(wodtlk)

    woddic = gdic
    woddic_ibim = np.zeros_like(woddic)
    woddic_cs1k = np.zeros_like(woddic)
    woddic_niva = np.zeros_like(woddic)
    woddic_brou = np.zeros_like(woddic)
    woddic_mobo = np.zeros_like(woddic)

    woddays, wod_ind, wod_count = np.unique([gyear, gmonth, gday], axis=1, return_index=True, return_counts=True)

    latids = np.zeros_like(woddic, dtype='int64')
    lonids = np.zeros_like(woddic, dtype='int64')
    latidc = np.zeros_like(woddic, dtype='int64')
    lonidc = np.zeros_like(woddic, dtype='int64')
    latidn = np.zeros_like(woddic, dtype='int64')
    lonidn = np.zeros_like(woddic, dtype='int64')
    latidw = np.zeros_like(woddic, dtype='int64')
    lonidw = np.zeros_like(woddic, dtype='int64')
    latidd = np.zeros_like(woddic, dtype='int64')
    lonidd = np.zeros_like(woddic, dtype='int64')
    land_ibi = np.zeros_like(woddic, dtype='bool')
    land_cs1km = np.zeros_like(woddic, dtype='bool')
    land_noresm = np.zeros_like(woddic, dtype='bool')
    land_brou = np.zeros_like(woddic, dtype='bool')
    land_mobo = np.zeros_like(woddic, dtype='bool')
    wod_day = np.zeros_like(woddic, dtype='int64')
    wod_month = np.zeros_like(woddic, dtype='int64')
    wod_year = np.zeros_like(woddic, dtype='int64')

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

        # Reanalysis
        glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
        glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
        glorys_mth_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01mav_'
        glorys_ending = '_R*_RE01.nc'

        if monthlyres == 1:
            IBI_PHYmth_files = sorted(glob.glob(glorysfiles_dir +
                                                glorys_mth_prefix + str(wod_year[wd]) +
                                                str(wod_month[wd]).zfill(2) +
                                                '01_*' + glorys_ending))
            wdphym = netcdf(IBI_PHYmth_files[0], 'r')
            
        if dailyres == 1:
            IBI_PHY_files = sorted(glob.glob(glorysfiles_dir +
                                             glorys_prefix + str(wod_year[wd]) +
                                             str(wod_month[wd]).zfill(2) + str(wod_day[wd]).zfill(2) +
                                             '_*' + glorys_ending))
            wdphym = netcdf(IBI_PHY_files[0], 'r')

        wdibmo = wdibim.variables
        wdphmo = wdphym.variables
        ibitmp = np.array(wdphmo['thetao'][:])
        ibisal = np.array(wdphmo['so'][:])
        
        # read latitude and longitude for indexing
        ibiph = np.array(wdibmo['ph'][:])
        ibidic = np.array(wdibmo['dissic'][:])*1000

        for dp in range(0, profs.shape[1]):
            ibiph_smp = ibiph[0, :, latids[dp], lonids[dp]]
            ibidic_smp = ibidic[0, :, latids[dp], lonids[dp]]
            ibitmp_smp = ibitmp[0, :, latids[dp], lonids[dp]]
            ibisal_smp = ibisal[0, :, latids[dp], lonids[dp]]
            if (any(ibidic_smp < 0) == True) & (all(ibidic_smp < 0) == False):
                ibiph_smp[ibiph_smp < 0] = ibiph_smp[ibiph_smp >= 0][-1]
                ibidic_smp[ibidic_smp < 0] = ibidic_smp[ibidic_smp >= 0][-1]
                ibitmp_smp[ibidic_smp < 0] = ibitmp_smp[ibidic_smp >= 0][-1]
                ibisal_smp[ibidic_smp < 0] = ibisal_smp[ibidic_smp >= 0][-1]
                results = pyco2.sys(par1=ibiph_smp, par1_type=3,
                                    par2=ibidic_smp, par2_type=2,
                                    temperature=ibitmp_smp,
                                    salinity=ibisal_smp,
                                    opt_pH_scale=1, opt_k_carbonic=4,
                                    opt_k_bisulfate=1, opt_total_borate=1)
                ibitlk_smp = results["alkalinity"]
            elif all(ibidic_smp < 0) == True:
                land_ibi[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_p = interpolate.interp1d(ibidep, ibitlk_smp, fill_value='extrapolate')
            wodtlk_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(ibidep, ibidic_smp, fill_value='extrapolate')
            woddic_ibim[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

        wdibim.close()
        wdphym.close()

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
            # read latitude and longitude for indexing
            cs1kmtlk = np.mean(np.array(wdcs1kmo['TALK'][:]), 0)
            cs1kmdic = np.mean(np.array(wdcs1kmo['DIC'][:]), 0)
            cs1kmz = np.mean(np.array(wdcs1kmo['zeta'][:]), 0)

        if dailyres == 1:
            # read latitude and longitude for indexing
            cs1kmtlk = np.array(wdcs1kmo['TALK'][wod_day[wd] - 1, :, :, :])
            cs1kmdic = np.array(wdcs1kmo['DIC'][wod_day[wd] - 1, :, :, :])
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
            cs1kmtlk_smp = cs1kmtlk[:, latidc[dp], lonidc[dp]]
            cs1kmdic_smp = cs1kmdic[:, latidc[dp], lonidc[dp]]
            cs1kmd_smp = cs1kmd[:, latidc[dp], lonidc[dp]]
            if (any(cs1kmdic_smp == 0) == True) & (all(cs1kmdic_smp == 0) == False):
                cs1kmtlk_smp[cs1kmtlk_smp == 0] = cs1kmtlk_smp[cs1kmtlk_smp > 0][-1]
                cs1kmdic_smp[cs1kmdic_smp == 0] = cs1kmdic_smp[cs1kmdic_smp > 0][-1]
            elif all(cs1kmdic_smp < 0) == True:
                land_cs1km[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_p = interpolate.interp1d(cs1kmd_smp, cs1kmtlk_smp, fill_value='extrapolate')
            wodtlk_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(cs1kmd_smp, cs1kmdic_smp, fill_value='extrapolate')
            woddic_cs1k[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for NORESM processing
    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']
    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    ncnortlk = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[5] + NORESM_ending, 'r')
    ncnordic = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[4] + NORESM_ending, 'r')
    ncnortlko = ncnortlk.variables
    ncnordico = ncnordic.variables
    noresmlat = np.array(ncnordico['lat'][:])
    noresmlon = np.array(ncnordico['lon'][:])
    depthniva = np.array(ncnordico['depth'][:])
    time_niva = np.array(ncnordico['time'][:])

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
        noresmtlk = np.array(ncnortlko[NIVAvars[5]][:]) * 1000
        noresmdic = np.array(ncnordico[NIVAvars[4]][:]) * 1000

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                Tinin = datetime(gyear[wd], gmonth[wd], 15)
                Tininx1 = d2i(Tinin, ncnordico['time'], select='exact')
                noresmtlk_smp = noresmtlk[Tininx1, :, latidn[dp], lonidn[dp]]
                noresmdic_smp = noresmdic[Tininx1, :, latidn[dp], lonidn[dp]]
            if dailyres == 1:
                # By day
                if gday[wd] != 15:
                    Tinin = datetime(gyear[wd], gmonth[wd], gday[wd])
                    Tininx1 = d2i(Tinin, ncnordico['time'], select='before')
                    Tininx2 = d2i(Tinin, ncnordico['time'], select='after')
                    if wod_day[wd] > 15:
                        tmult = gday[wd] - 15
                    else:
                        tmult = gday[wd] + 15

                    # Tidx1
                    noresmtlk_1 = noresmtlk[Tininx1, :, latidn[dp], lonidn[dp]]
                    noresmdic_1 = noresmdic[Tininx1, :, latidn[dp], lonidn[dp]]

                    # Tidx2
                    noresmtlk_2 = noresmtlk[Tininx2, :, latidn[dp], lonidn[dp]]
                    noresmdic_2 = noresmdic[Tininx2, :, latidn[dp], lonidn[dp]]

                    # Interpolating to the day
                    noresmtlk_smp = noresmtlk_1 + (((noresmtlk_2 - noresmtlk_1)/30)*tmult)
                    noresmdic_smp = noresmdic_1 + (((noresmdic_2 - noresmdic_1)/30)*tmult)
                else:
                    Tinin = datetime(gyear[wd], gmonth[wd], 15)
                    Tininx1 = d2i(Tinin, ncnordico['time'], select='exact')
                    noresmtlk_smp = noresmtlk[Tininx1, :, latidn[dp], lonidn[dp]]
                    noresmdic_smp = noresmdic[Tininx1, :, latidn[dp], lonidn[dp]]

            if (any(noresmdic_smp > 10000) == True) & (all(noresmdic_smp > 10000) == False):
                noresmtlk_smp[noresmtlk_smp > 10000] = noresmtlk_smp[noresmtlk_smp < 10000][-1]
                noresmdic_smp[noresmdic_smp > 10000] = noresmdic_smp[noresmdic_smp < 10000][-1]
            elif all(noresmdic_smp > 10000) == True:
                land_noresm[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
            #
            wodint_p = interpolate.interp1d(depthniva, noresmtlk_smp, fill_value='extrapolate')
            wodtlk_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_p(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
            #
            wodint_s = interpolate.interp1d(depthniva, noresmdic_smp, fill_value='extrapolate')
            woddic_niva[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    ncnortlk.close()
    ncnordic.close()

    # Block for Broullon processing
    valdir = '/media/dskone/VAL/'
    broudir = valdir + 'BroullonClims/'

    # Common woa elements for monthly climatology
    brou_dic = 'TCO2_NNGv2LDEO_climatology_WOA18.nc'
    brou_tlk = 'AT_NNGv2_climatology_WOA18.nc'
    brou_dicv = netcdf(broudir + brou_dic, 'r')
    brou_tlkv = netcdf(broudir + brou_tlk, 'r')
    brou_dicvo = brou_dicv.variables
    brou_tlkvo = brou_tlkv.variables
    brou_lat = np.array(brou_dicvo['latitude'][:])
    brou_lon = np.array(brou_dicvo['longitude'][:])
    brou_depth = np.array(brou_dicvo['depth'][:])
    broudic = np.array(brou_dicvo['TCO2_NNGv2LDEO'][:])
    broutlk = np.array(brou_tlkvo['AT_NNGv2'][:])
    brou_dicv.close()
    brou_tlkv.close()

    woddays, wod_ind, wod_count = np.unique([gyear, gmonth, gday], axis=1, return_index=True, return_counts=True)

    for wd in range(0, woddays.shape[1]):
        wod_day[wd] = woddays[2, wd]
        wod_month[wd] = woddays[1, wd]
        wod_year[wd] = woddays[0, wd]

        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)

        latidw = np.zeros_like(prof_ind, dtype='int64')
        lonidw = np.zeros_like(prof_ind, dtype='int64')

        for p in range(0, profs.shape[1]):
            latidw[p] = glor.geo_idx(profs[0, p], brou_lat)
            lonidw[p] = glor.geo_idx(profs[1, p], brou_lon)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                broudic_smp = broudic[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                broutlk_smp = broutlk[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]

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
                    broudic_1 = broudic[tidx1, :, latidw[dp], lonidw[dp]]
                    broutlk_1 = broutlk[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    broudic_2 = broudic[tidx2, :, latidw[dp], lonidw[dp]]
                    broutlk_2 = broutlk[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    broudic_smp = broudic_1 + (((broudic_2 - broudic_1)/30)*tmult)
                    broutlk_smp = broutlk_1 + (((broutlk_2 - broutlk_1)/30)*tmult)

                elif wod_day[wd] < 15:
                    tmult = wod_day[wd] + 15
                    tidx2 = wod_month[wd] - 1
                    if wod_month[wd] == 1:
                        tidx1 = 11
                    else:
                        tidx1 = wod_month[wd] - 2

                    # Tidx1
                    broudic_1 = broudic[tidx1, :, latidw[dp], lonidw[dp]]
                    broutlk_1 = broutlk[tidx1, :, latidw[dp], lonidw[dp]]

                    # Tidx2
                    broudic_2 = broudic[tidx2, :, latidw[dp], lonidw[dp]]
                    broutlk_2 = broutlk[tidx2, :, latidw[dp], lonidw[dp]]

                    # Interpolating to the day
                    broudic_smp = broudic_1 + (((broudic_2 - broudic_1)/30)*tmult)
                    broutlk_smp = broutlk_1 + (((broutlk_2 - broutlk_1)/30)*tmult)

                else:
                    broudic_smp = broudic[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]
                    broutlk_smp = broutlk[wod_month[wd] - 1, :, latidw[dp], lonidw[dp]]

            if (any(np.isnan(broudic_smp)) == True) & (all(np.isnan(broudic_smp)) == False):
                broudic_smp[np.isnan(broudic_smp)] = broudic_smp[np.isnan(broudic_smp) == False][-1]
                broutlk_smp[np.isnan(broutlk_smp)] = broutlk_smp[np.isnan(broutlk_smp) == False][-1]
                #
                wodint_dic = interpolate.interp1d(brou_depth, broudic_smp, fill_value='extrapolate')
                woddic_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_dic(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                #
                wodint_tlk = interpolate.interp1d(brou_depth, broutlk_smp, fill_value='extrapolate')
                wodtlk_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_tlk(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

            elif all(np.isnan(broudic_smp)) == True:
                # land_brou[wod_ind[wd] + dp] = True
                land_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                woddic_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                wodtlk_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
            else:
                #
                wodint_dic = interpolate.interp1d(brou_depth, broudic_smp, fill_value='extrapolate')
                woddic_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_dic(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                #
                wodint_tlk = interpolate.interp1d(brou_depth, broutlk_smp, fill_value='extrapolate')
                wodtlk_brou[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                    wodint_tlk(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    # Block for MOBO_DIC processing
    MOBO_dir = '/media/dskone/VAL/BroullonClims/'
    MOBOfil = 'MOBO_MPI_DIC_2004-2019_Oct2022.nc'

    ncmobodic = netcdf(MOBO_dir + MOBOfil, 'r')
    ncmobodico = ncmobodic.variables
    mobolat = np.array(ncmobodico['lat'][:])
    mobolon = np.array(ncmobodico['lon'][:])
    depthmobo = np.array(ncmobodico['depth'][:])
    time_mobo = np.array(ncmobodico['juld'][:])
    mobodic = np.array(ncmobodico['DIC'][:])
    ncmobodic.close()

    woddays, wod_ind, wod_count = np.unique([gyear, gmonth, gday], axis=1, return_index=True, return_counts=True)

    for wd in range(0, woddays.shape[1]):
        profs, prof_ind, prof_count = np.unique([glat[wod_ind[wd]:wod_ind[wd]+wod_count[wd]],
                                                 glon[wod_ind[wd]:wod_ind[wd]+wod_count[wd]]],
                                                axis=1, return_index=True, return_counts=True)
        latidn = np.zeros_like(prof_ind, dtype='int64')
        lonidn = np.zeros_like(prof_ind, dtype='int64')

        for p in range(0, profs.shape[1]):
            latidn[p] = glor.geo_idx(profs[0, p], mobolat)
            lonidn[p] = glor.geo_idx(profs[1, p], mobolon)

        latidn = latidn.astype(int)
        lonidn = lonidn.astype(int)

        for dp in range(0, profs.shape[1]):
            if monthlyres == 1:
                Tyridx = (wod_year[wd] - 2004)*12
                Tyridx = Tyridx - 1
                Tmthidx = wod_month[wd]
                Ttotidx = Tyridx + Tmthidx
                if (Ttotidx >= 0) & (Ttotidx <= mobodic.shape[0] - 1):
                    mobodic_smp = mobodic[Ttotidx, :, latidn[dp], lonidn[dp]]
                    if (any(mobodic_smp < 0) == True) & (all(mobodic_smp < 0) == False):
                        mobodic_smp[mobodic_smp < 0] = mobodic_smp[mobodic_smp > 0][-1]
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                    elif all(mobodic_smp < 0) == True:
                        land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                    #
                    else:
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                else:
                    land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                    woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0

            if dailyres == 1:
                if wod_day[wd] == 15:
                    Tyridx = (wod_year[wd] - 2004)*12
                    Tyridx = Tyridx - 1
                    Tmthidx = wod_month[wd]
                    Ttotidx = Tyridx + Tmthidx
                    if (Ttotidx >= 0) & (Ttotidx <= mobodic.shape[0] - 1):
                        mobodic_smp = mobodic[Ttotidx, :, latidn[dp], lonidn[dp]]
                if wod_day[wd] < 15:
                    Tyridx = (wod_year[wd] - 2004)*12
                    Tyridx = Tyridx - 1
                    Tmthidx1 = wod_month[wd] - 1
                    Ttotidx1 = Tyridx + Tmthidx1
                    Tmthidx2 = wod_month[wd]
                    Ttotidx2 = Tyridx + Tmthidx2
                    tmult = wod_day[wd] + 15
                    if (Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1):
                        mobodic_1 = mobodic[Ttotidx1, :, latidn[dp], lonidn[dp]]
                    if (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1):
                        mobodic_2 = mobodic[Ttotidx2, :, latidn[dp], lonidn[dp]]
                    if ((Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1) &
                            (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1)):
                        # Interpolating to the day
                        mobodic_smp = mobodic_1 + (((mobodic_2 - mobodic_1)/30)*tmult)

                if wod_day[wd] > 15:
                    Tyridx = (wod_year[wd] - 2004)*12
                    Tyridx = Tyridx - 1
                    Tmthidx1 = wod_month[wd]
                    Ttotidx1 = Tyridx + Tmthidx1
                    Tmthidx2 = wod_month[wd] + 1
                    Ttotidx2 = Tyridx + Tmthidx2
                    tmult = wod_day[wd] - 15
                    if (Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1):
                        mobodic_1 = mobodic[Ttotidx1, :, latidn[dp], lonidn[dp]]
                    if (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1):
                        mobodic_2 = mobodic[Ttotidx2, :, latidn[dp], lonidn[dp]]
                    if ((Ttotidx1 >= 0) & (Ttotidx1 <= mobodic.shape[0] - 1) &
                            (Ttotidx2 >= 0) & (Ttotidx2 <= mobodic.shape[0] - 1)):
                        # Interpolating to the day
                        mobodic_smp = mobodic_1 + (((mobodic_2 - mobodic_1)/30)*tmult)

                    if (any(mobodic_smp < 0) == True) & (all(mobodic_smp < 0) == False):
                        mobodic_smp[mobodic_smp < 0] = mobodic_smp[mobodic_smp > 0][-1]
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])
                    elif all(mobodic_smp < 0) == True:
                        land_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = True
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = 0
                    #
                    else:
                        wodint_s = interpolate.interp1d(depthmobo, mobodic_smp, fill_value='extrapolate')
                        woddic_mobo[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]] = \
                            wodint_s(z[wod_ind[wd] + prof_ind[dp]:wod_ind[wd] + prof_ind[dp] + prof_count[dp]])

    glodap_bias_file = crocodir + 'glodap_carbon_bias_' + npzstr + '.npz'
    np.savez(glodap_bias_file,
             wodlat=wodlat, wodlon=wodlon, z=z,
             wodtlk=wodtlk,
             wodtlk_ibim=wodtlk_ibim, wodtlk_cs1k=wodtlk_cs1k, wodtlk_niva=wodtlk_niva, wodtlk_brou=wodtlk_brou,
             woddic=woddic,
             woddic_ibim=woddic_ibim, woddic_cs1k=woddic_cs1k, woddic_niva=woddic_niva, woddic_brou=woddic_brou,
             woddic_mobo=woddic_mobo,
             latids=latids, lonids=lonids, latidc=latidc, lonidc=lonidc, latidn=latidn, lonidn=lonidn,
             latidw=latidw, lonidw=lonidw, latidd=latidd, lonidd=lonidd,
             land_ibi=land_ibi, land_cs1km=land_cs1km, land_noresm=land_noresm,
             land_brou=land_brou, land_mobo=land_mobo,
             wod_day=gday, wod_month=gmonth, wod_year=gyear)

mi_plotting = 0
glodap_plotting = 0

if mi_plotting == 1:
    mi_bias_file = crocodir + 'MIwint_carbon_bias_' + npzstr + '.npz'
    mi_bias = np.load(mi_bias_file)
    wodlat = mi_bias['wodlat']
    wodlon = mi_bias['wodlon']
    z = mi_bias['z']
    wodtlk = mi_bias['wodtlk']
    wodtlk_ibim = mi_bias['wodtlk_ibim']
    wodtlk_cs1k = mi_bias['wodtlk_cs1k']
    wodtlk_niva = mi_bias['wodtlk_niva']
    wodtlk_brou = mi_bias['wodtlk_brou']
    woddic = mi_bias['woddic']
    woddic_ibim = mi_bias['woddic_ibim']
    woddic_cs1k = mi_bias['woddic_cs1k']
    woddic_niva = mi_bias['woddic_niva']
    woddic_brou = mi_bias['woddic_brou']
    woddic_mobo = mi_bias['woddic_mobo']
    land_ibi = mi_bias['land_ibi']
    land_cs1km = mi_bias['land_cs1km']
    land_noresm = mi_bias['land_noresm']
    land_brou = mi_bias['land_brou']
    land_mobo = mi_bias['land_mobo']
    wod_day = mi_bias['wod_day']
    wod_month = mi_bias['wod_month']
    wod_year = mi_bias['wod_year']

if glodap_plotting == 1:
    glodap_bias_file = crocodir + 'glodap_carbon_bias_' + npzstr + '.npz'
    glodap_bias = np.load(glodap_bias_file)
    wodlat = np.concatenate((wodlat, glodap_bias['wodlat']))
    wodlon = np.concatenate((wodlon, glodap_bias['wodlon']))
    z = np.concatenate((z, glodap_bias['z']))
    wodtlk = np.concatenate((wodtlk, glodap_bias['wodtlk']))
    wodtlk_ibim = np.concatenate((wodtlk_ibim, glodap_bias['wodtlk_ibim']))
    wodtlk_cs1k = np.concatenate((wodtlk_cs1k, glodap_bias['wodtlk_cs1k']))
    wodtlk_niva = np.concatenate((wodtlk_niva, glodap_bias['wodtlk_niva']))
    wodtlk_brou = np.concatenate((wodtlk_brou, glodap_bias['wodtlk_brou']))
    woddic = np.concatenate((woddic, glodap_bias['woddic']))
    woddic_ibim = np.concatenate((woddic_ibim, glodap_bias['woddic_ibim']))
    woddic_cs1k = np.concatenate((woddic_cs1k, glodap_bias['woddic_cs1k']))
    woddic_niva = np.concatenate((woddic_niva, glodap_bias['woddic_niva']))
    woddic_brou = np.concatenate((woddic_brou, glodap_bias['woddic_brou']))
    woddic_mobo = np.concatenate((woddic_mobo, glodap_bias['woddic_mobo']))
    land_ibi = np.concatenate((land_ibi, glodap_bias['land_ibi']))
    land_cs1km = np.concatenate((land_cs1km, glodap_bias['land_cs1km']))
    land_noresm = np.concatenate((land_noresm, glodap_bias['land_noresm']))
    land_brou = np.concatenate((land_brou, glodap_bias['land_brou']))
    land_mobo = np.concatenate((land_mobo, glodap_bias['land_mobo']))
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

    def carb_bias_plt_stats(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, clim_mobo, dep,
                            wod_yr, wod_mth, wod_dy, otimes,
                            wod_lon, wod_lat, geodist, nutstr):
        ibimfilt = np.argwhere((obs > 0) & (obs < 10000) & (mod_ibim > 0) & (mod_ibim < 10000))[:, 0]
        cs1kfilt = np.argwhere((obs > 0) & (obs < 10000) & (mod_cs1k > 0) & (mod_cs1k < 10000))[:, 0]
        nivafilt = np.argwhere((obs > 0) & (obs < 10000) & (mod_niva > 0) & (mod_niva < 10000))[:, 0]
        broufilt = np.argwhere((obs > 0) & (obs < 10000) & (clim_brou > 0) & (clim_brou < 10000))[:, 0]
        mobofilt = np.argwhere((obs > 0) & (obs < 10000) & (clim_mobo > 0) & (clim_mobo < 10000))[:, 0]
        # allfilt = np.argwhere((obs > 0) & (obs < 10000) &
        #                       (mod_ibim > 0) & (mod_ibim < 10000) &
        #                       (mod_cs1k > 0) & (mod_cs1k < 10000) &
        #                       (mod_niva > 0) & (mod_niva < 10000) &
        #                       (clim_brou > 0) & (clim_brou < 10000) &
        #                       (clim_mobo > 0) & (clim_mobo < 10000))[:, 0]
        allfilt = np.argwhere((obs > 0) & (obs < 10000) &
                              (mod_ibim > 0) & (mod_ibim < 10000) &
                              (mod_cs1k > 0) & (mod_cs1k < 10000) &
                              (mod_niva > 0) & (mod_niva < 10000) &
                              (clim_brou > 0) & (clim_brou < 10000))[:, 0]
        uppfilt = np.argwhere((dep <= 20))[:, 0]
        midfilt = np.argwhere((dep > 20) & (dep <= 70))[:, 0]
        depfilt = np.argwhere((dep > 70))[:, 0]
        djffilt = np.argwhere((wod_mth == 12) | (wod_mth == 1) | (wod_mth == 2))[:, 0]
        mamfilt = np.argwhere((wod_mth == 3) | (wod_mth == 4) | (wod_mth == 5))[:, 0]
        jjafilt = np.argwhere((wod_mth == 6) | (wod_mth == 7) | (wod_mth == 8))[:, 0]
        sonfilt = np.argwhere((wod_mth == 9) | (wod_mth == 10) | (wod_mth == 11))[:, 0]
        uppfiltII = np.intersect1d(ibimfilt, uppfilt)
        midfiltII = np.intersect1d(ibimfilt, midfilt)
        depfiltII = np.intersect1d(ibimfilt, depfilt)
        csdfilt = np.argwhere((wod_lat >= 49) & (wod_lat <= 52.95) & (wod_lon >= -10.75) & (wod_lon <= -5.83))[:, 0]

        cs1k_ibi = np.argwhere((obs > 0) & (obs < 10000) &
                               (mod_ibim > 0) & (mod_ibim < 10000) &
                               (mod_cs1k > 0) & (mod_cs1k < 10000))[:, 0]
        cs1kibi = np.intersect1d(cs1k_ibi, csdfilt)

        ibi_niva = np.argwhere((obs > 0) & (obs < 10000) &
                               (mod_ibim > 0) & (mod_ibim < 10000) &
                               (mod_niva > 0) & (mod_niva < 10000))[:, 0]
        ibiniva = np.setdiff1d(ibi_niva, csdfilt)

        brou_mobo_comp = np.argwhere((obs > 0) & (obs < 10000) &
                                     (clim_brou > 0) & (clim_brou < 10000) &
                                     (clim_mobo > 0) & (clim_mobo < 10000))[:, 0]
        ibi_brou_mobo = np.argwhere((obs > 0) & (obs < 10000) &
                                    (mod_ibim > 0) & (mod_ibim < 10000) &
                                    (clim_mobo > 0) & (clim_mobo < 10000) &
                                    (clim_brou > 0) & (clim_brou < 10000))[:, 0]
        brou_mobo = np.argwhere((obs > 0) & (obs < 10000) &
                                (clim_brou > 0) & (clim_brou < 10000))[:, 0]
        broumobo_bry = np.setdiff1d(brou_mobo, csdfilt)

        allfilt_II = np.intersect1d(allfilt, csdfilt)

        # allfilt_t = allfilt
        allfilt_t = allfilt_II

        ibinivabrou = np.intersect1d(ibiniva, broumobo_bry)

        scatr4(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, allfilt_II, nutstr, 'common')
        RMSD_ibim_t, R2_ibim_t, MD_ibim_t, N_ibim_t = get_eans(obs[allfilt_II], mod_ibim[allfilt_II])
        RMSD_cs1k_t, R2_cs1k_t, MD_cs1k_t, N_cs1k_t = get_eans(obs[allfilt_II], mod_cs1k[allfilt_II])
        RMSD_niva_t, R2_niva_t, MD_niva_t, N_niva_t = get_eans(obs[allfilt_II], mod_niva[allfilt_II])
        RMSD_brou_t, R2_brou_t, MD_brou_t, N_brou_t = get_eans(obs[allfilt_II], clim_brou[allfilt_II])

        RMSD_ibim_1dm, R2_ibim_1dm, MD_ibim_1dm, N_ibim_1dm = get_eans(obs[allfilt], mod_ibim[allfilt])
        RMSD_niva_1dm, R2_niva_1dm, MD_niva_1dm, N_niva_1dm = get_eans(obs[allfilt], mod_niva[allfilt])
        RMSD_brou_1dm, R2_brou_1dm, MD_brou_1dm, N_brou_1dm = get_eans(obs[allfilt], clim_brou[allfilt])

        biscatr(fdir, obs, mod_ibim, mod_cs1k, 'IBI', 'CS1K', cs1kibi, nutstr, 'mmol/m3', 'In-domain', 'Model')
        RMSD_ibi_dom, R2_ibi_dom, MD_ibi_dom, N_ibi_dom = get_eans(obs[cs1kibi], mod_ibim[cs1kibi])
        RMSD_cs1k_dom, R2_cs1k_dom, MD_cs1k_dom, N_cs1k_dom = get_eans(obs[cs1kibi], mod_cs1k[cs1kibi])

        biscatr(fdir, obs, clim_brou, clim_mobo, 'BROU', 'MOBO', brou_mobo_comp, nutstr, 'mmol/m3', 'Overlap', 'Clim')
        RMSD_brou_cmn, R2_brou_cmn, MD_brou_cmn, N_brou_cmn = get_eans(obs[brou_mobo_comp], clim_brou[brou_mobo_comp])
        RMSD_mobo_cmn, R2_mobo_cmn, MD_mobo_cmn, N_mobo_cmn = get_eans(obs[brou_mobo_comp], clim_mobo[brou_mobo_comp])

        triscatr(fdir, obs, mod_ibim, mod_niva, clim_brou, 'IBI', 'NIVA', 'BROU', ibinivabrou, nutstr, 'mmol/m3', 'Degree-margin', 'Models+Clim')
        RMSD_ibi_bry, R2_ibi_bry, MD_ibi_bry, N_ibi_bry = get_eans(obs[ibinivabrou], mod_ibim[ibinivabrou])
        RMSD_niva_bry, R2_niva_bry, MD_niva_bry, N_niva_bry = get_eans(obs[ibinivabrou], mod_niva[ibinivabrou])
        RMSD_brou_bry, R2_brou_bry, MD_brou_bry, N_brou_bry = get_eans(obs[ibinivabrou], clim_brou[ibinivabrou])

        # Nutrient colormap for spatial plot
        nutcmp = 'tab20c'  # 'bwr_r', 'gist_ncar', 'RdYlBu'
        prod_biases(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, clim_mobo, wodlon, wodlat,
                    nutstr, '(mmol/m3)', 'tab20c', allfilt, cs1kibi)

        obsVsprods(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, clim_mobo, wodlon, wodlat,
                   nutstr, '(mmol/m3)', 'tab20c', allfilt, cs1kibi)

        obs_prod_biases_x9(fdir, obs, mod_ibim, mod_cs1k, mod_niva, wodlon, wodlat,
                           nutstr, '(mmol/m3)', 'tab20c', allfilt, cs1kibi)

        lt_tseries(fdir, otimes, 'bias', obs, mod_ibim, mod_cs1k, mod_niva,
                   allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'bias_models', 'Model')

        lt_tseries(fdir, otimes, 'bias', obs, mod_ibim, clim_brou, clim_mobo,
                   ibi_brou_mobo, ibi_brou_mobo, ibi_brou_mobo, 'IBI', 'BROU', 'MOBO', nutstr, '(mmol/m3)', 'bias_clims', 'IBI+Clims.')

        lt_tseries(fdir, otimes, 'data', obs, mod_ibim, mod_cs1k, mod_niva,
                   allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'data_models', 'Model')

        lt_tseries(fdir, otimes, 'data', obs, mod_ibim, clim_brou, clim_mobo,
                   ibi_brou_mobo, ibi_brou_mobo, ibi_brou_mobo, 'IBI', 'BROU', 'MOBO', nutstr, '(mmol/m3)', 'data_clims', 'IBI+Clims.')

        t_clim(fdir, wod_mth, 'bias', obs, mod_ibim, mod_cs1k, mod_niva,
               allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'bias_models', 'Model')

        t_clim(fdir, wod_mth, 'bias', obs, mod_ibim, clim_brou, clim_mobo,
               ibi_brou_mobo, ibi_brou_mobo, ibi_brou_mobo, 'IBI', 'BROU', 'MOBO', nutstr, '(mmol/m3)', 'bias_clims', 'IBI+Clims.')

        t_clim(fdir, wod_mth, 'data', obs, mod_ibim, mod_cs1k, mod_niva,
               allfilt_t, allfilt_t, allfilt_t, 'IBI', 'CS1K', 'NIVA', nutstr, '(mmol/m3)', 'data_models', 'Model')

        t_clim(fdir, wod_mth, 'data', obs, mod_ibim, clim_brou, clim_mobo,
               ibi_brou_mobo, ibi_brou_mobo, ibi_brou_mobo, 'IBI', 'BROU', 'MOBO', nutstr, '(mmol/m3)', 'data_clims', 'IBI+Clims.')

        return (RMSD_ibim_t, R2_ibim_t, MD_ibim_t, N_ibim_t), (RMSD_cs1k_t, R2_cs1k_t, MD_cs1k_t, N_cs1k_t), \
            (RMSD_niva_t, R2_niva_t, MD_niva_t, N_niva_t), (RMSD_brou_t, R2_brou_t, MD_brou_t, N_brou_t), \
            (RMSD_ibim_1dm, R2_ibim_1dm, MD_ibim_1dm, N_ibim_1dm), \
            (RMSD_niva_1dm, R2_niva_1dm, MD_niva_1dm, N_niva_1dm), (RMSD_brou_1dm, R2_brou_1dm, MD_brou_1dm, N_brou_1dm), \
            (RMSD_ibi_dom, R2_ibi_dom, MD_ibi_dom, N_ibi_dom), \
            (RMSD_cs1k_dom, R2_cs1k_dom, MD_cs1k_dom, N_cs1k_dom), \
            (RMSD_brou_cmn, R2_brou_cmn, MD_brou_cmn, N_brou_cmn), \
            (RMSD_mobo_cmn, R2_mobo_cmn, MD_mobo_cmn, N_mobo_cmn), \
            (RMSD_ibi_bry, R2_ibi_bry, MD_ibi_bry, N_ibi_bry), (RMSD_niva_bry, R2_niva_bry, MD_niva_bry, N_niva_bry), \
            (RMSD_brou_bry, R2_brou_bry, MD_brou_bry, N_brou_bry)

    def get_eans(obs, mod):
        N = len(obs)
        RMSD = (np.sum((obs - mod) ** 2.) / len(obs)) ** 0.5
        R2 = np.corrcoef(obs, mod)[0, 1]
        MD = np.mean(obs) - np.mean(mod)
        return RMSD, R2, MD, N

    # Scatters
    def scatr4(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, qfilt, nutstr, descr):
        # 5 scatter
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('In-situ')
        plt.ylabel('model/climatology')
        # (ax1, ax2, ax3), (ax4, ax5, ax6) = axs
        (ax1, ax2), (ax3, ax4) = axs
        ax1.scatter(obs[qfilt], mod_ibim[qfilt])
        ax1.grid(visible=True, which='major', color='grey', linestyle='--')
        ax1.set_xlim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        ax1.set_ylim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        lims1 = [np.min([obs[qfilt].min(),
                         mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                         clim_brou[qfilt].min()]),
                 np.max([obs[qfilt].max(),
                         mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                         clim_brou[qfilt].max()])]
        ax1.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)

        ax2.scatter(obs[qfilt], mod_cs1k[qfilt])
        ax2.grid(visible=True, which='major', color='grey', linestyle='--')
        ax2.set_xlim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        ax2.set_ylim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        lims2 = [np.min([obs[qfilt].min(),
                         mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                         clim_brou[qfilt].min()]),
                 np.max([obs[qfilt].max(),
                         mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                         clim_brou[qfilt].max()])]
        ax2.plot(lims2, lims2, 'k--', alpha=0.75, zorder=0)

        ax3.scatter(obs[qfilt], mod_niva[qfilt])
        ax3.grid(visible=True, which='major', color='grey', linestyle='--')
        ax3.set_xlim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        ax3.set_ylim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        lims3 = [np.min([obs[qfilt].min(),
                         mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                         clim_brou[qfilt].min()]),
                 np.max([obs[qfilt].max(),
                         mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                         clim_brou[qfilt].max()])]
        ax3.plot(lims3, lims3, 'k--', alpha=0.75, zorder=0)

        ax4.scatter(obs[qfilt], clim_brou[qfilt])
        ax4.grid(visible=True, which='major', color='grey', linestyle='--')
        ax4.set_xlim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        ax4.set_ylim([np.min([obs[qfilt].min(),
                              mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                              clim_brou[qfilt].min()]),
                      np.max([obs[qfilt].max(),
                              mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                              clim_brou[qfilt].max()])])
        lims4 = [np.min([obs[qfilt].min(),
                         mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
                         clim_brou[qfilt].min()]),
                 np.max([obs[qfilt].max(),
                         mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
                         clim_brou[qfilt].max()])]
        ax4.plot(lims4, lims4, 'k--', alpha=0.75, zorder=0)

        # ax6.scatter(obs[qfilt], clim_mobo[qfilt])
        # ax6.set_xlim([np.min([obs[qfilt].min(),
        #                       mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
        #                       clim_brou[qfilt].min(), clim_mobo[qfilt].min()]),
        #               np.max([obs[qfilt].max(),
        #                       mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
        #                       clim_brou[qfilt].max(), clim_mobo[qfilt].max()])])
        # ax6.set_ylim([np.min([obs[qfilt].min(),
        #                       mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
        #                       clim_brou[qfilt].min(), clim_mobo[qfilt].min()]),
        #               np.max([obs[qfilt].max(),
        #                       mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
        #                       clim_brou[qfilt].max(), clim_mobo[qfilt].max()])])
        # lims6 = [np.min([obs[qfilt].min(),
        #                  mod_ibim[qfilt].min(), mod_cs1k[qfilt].min(), mod_niva[qfilt].min(),
        #                  clim_brou[qfilt].min(), clim_mobo[qfilt].min()]),
        #          np.max([obs[qfilt].max(),
        #                  mod_ibim[qfilt].max(), mod_cs1k[qfilt].max(), mod_niva[qfilt].max(),
        #                  clim_brou[qfilt].max(), clim_mobo[qfilt].max()])]
        # ax6.plot(lims6, lims6, 'k--', alpha=0.75, zorder=0)

        fig.suptitle(nutstr + ' ' + descr + ' In-situ Vs. models/clims.')
        ax1.set_title('IBI', fontsize=10)
        ax2.set_title('CS1KM', fontsize=10)
        ax3.set_title('NIVA', fontsize=10)
        ax4.set_title('BROU', fontsize=10)
        # ax6.set_title('MOBO', fontsize=10)
        fig.tight_layout()
        # plt.show()
        plt_tit = fdir + descr + '_' + nutstr + '_scatcomp.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

        # Scatters
    def biscatr(fdir, obs, dat1, dat2, lab1, lab2, qfilt, nutstr, nutunits, descr, climormod):
        # 2 scatter
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
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


    def triscatr(fdir, obs, dat1, dat2, dat3, lab1, lab2, lab3, qfilt, nutstr, nutunits, descr, climormod):
        # 2 scatter
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('In-situ, (' + nutunits + ')')
        plt.ylabel(climormod + ', (' + nutunits + ')')
        lims1 = [np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                 np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])]
        (ax1, ax2, ax3) = axs
        ax1.scatter(obs[qfilt], dat1[qfilt])
        ax1.set_xlim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])])
        ax1.set_ylim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])])
        ax1.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)

        ax2.scatter(obs[qfilt], dat2[qfilt])
        ax2.set_xlim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])])
        ax2.set_ylim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])])
        ax2.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)

        ax3.scatter(obs[qfilt], dat3[qfilt])
        ax3.set_xlim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])])
        ax3.set_ylim([np.min([obs[qfilt].min(), dat1[qfilt].min(), dat2[qfilt].min(), dat3[qfilt].min()]),
                      np.max([obs[qfilt].max(), dat1[qfilt].max(), dat2[qfilt].max(), dat3[qfilt].max()])])
        ax3.plot(lims1, lims1, 'k--', alpha=0.75, zorder=0)
        fig.suptitle(nutstr + ' ' + descr + ' In-situ Vs. ' + climormod)
        ax1.set_title(lab1, fontsize=10)
        ax2.set_title(lab2, fontsize=10)
        ax3.set_title(lab3, fontsize=10)
        fig.tight_layout()
        # plt.show()
        plt_tit = fdir + descr + '_' + nutstr + lab1 + '_' + lab2 + '_triscatr.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def prod_biases(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, clim_mobo, wodlon, wodlat,
                    nutstr, nutunit, nutcmp, qfilt, cs1kfilt):

        # Spatial map of bias
        # nutstr = 'NO3'
        # nutunit = 'mmol/m3'
        # nutcmp = 'tab20c'  # 'bwr_r', 'gist_ncar', 'RdYlBu'

        # cm = plt.cm.get_cmap(nutcmp)
        cm = plt.cm.get_cmap('seismic')
        nmin = np.min([(obs[qfilt]-mod_ibim[qfilt]).min(),
                       (obs[cs1kfilt]-mod_cs1k[cs1kfilt]).min(),
                       (obs[qfilt]-mod_niva[qfilt]).min(),
                       (obs[qfilt]-clim_brou[qfilt]).min()])
        nmax = np.max([(obs[qfilt]-mod_ibim[qfilt]).max(),
                       (obs[cs1kfilt]-mod_cs1k[cs1kfilt]).max(),
                       (obs[qfilt]-mod_niva[qfilt]).max(),
                       (obs[qfilt]-clim_brou[qfilt]).max()])
        if nmin < 0:
            bigger = np.max([abs(nmin), abs(nmax)])
            if bigger > 10:
                bigger2 = np.ceil(bigger/10)*10
            elif (bigger < 10) & (bigger > 1):
                bigger2 = np.ceil(bigger)
            else:
                bigger2 = bigger
            nmin = bigger2*(-1)
            nmax = bigger2
        else:
            nmin = 0
            if nmax > 10:
                nmax2 = np.ceil(nmax/10)*10
            elif (nmax < 10) & (nmax > 1):
                nmax2 = np.ceil(nmax)
            else:
                nmax2 = nmax
            nmax = nmax2
        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(13, 7.5))
        (ax1, ax2, ax3), (ax4, ax5, ax6) = axs
        # (ax1, ax2), (ax3, ax4) = axs
        ax2.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - mod_ibim[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax3.scatter(wodlon[cs1kfilt], wodlat[cs1kfilt], c=obs[cs1kfilt] - mod_cs1k[cs1kfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax5.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - mod_niva[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        im = ax6.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - clim_brou[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        # im = ax6.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt] - clim_mobo[qfilt],
        #             cmap=cm, vmin=nmin, vmax=nmax, s=5)
        fig.suptitle(nutstr + ' ' + 'Bias,' + ' ' + nutunit, fontsize=10)
        ax2.title.set_text('IBI')
        ax3.title.set_text('CS1K')
        ax5.title.set_text('NIVA')
        ax6.title.set_text('BROU')
        # ax6.title.set_text('MOBO')
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

        ax1.title.set_text('Obs @ IBI')
        ax1.set_facecolor('gainsboro')
        m1 = Basemap(ax=ax1, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m1.fillcontinents(color='darkseagreen')
        x1, y1 = m1(wodlon[csidx], wodlat[csidx])
        m1.scatter(x1, y1, c=obs[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m1.drawparallels(np.arange(49., 53., 1.), labels=[1, 0, 0, 0], linewidth=0.25)
        m1.drawmeridians(np.arange(-10.75, -5.75, 1.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax2.title.set_text('Obs @ NIVA')
        ax2.set_facecolor('gainsboro')
        m2 = Basemap(ax=ax2, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m2.fillcontinents(color='darkseagreen')
        x2, y2 = m2(wodlon[csidx], wodlat[csidx])
        m2.scatter(x2, y2, c=obs[csidx],
                   cmap=cmo.matter, vmin=nmin, vmax=nmax, s=5)
        m2.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m2.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 0], linewidth=0.25)

        ax3.title.set_text('Obs @ CS1K')
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

        ax7.title.set_text('Obs - IBI')
        ax7.set_facecolor('gainsboro')
        m7 = Basemap(ax=ax7, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m7.fillcontinents(color='darkseagreen')
        x7, y7 = m7(wodlon[csidx], wodlat[csidx])
        m7.scatter(x7, y7, c=mod_ibim[csidx] - obs[csidx],
                   cmap=cm_bias, vmin=nmin_bias, vmax=nmax_bias, s=5)
        m7.drawparallels(np.arange(49., 53., 1.), labels=[1, 0, 0, 0], linewidth=0.25)
        m7.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 1], linewidth=0.25)

        ax8.title.set_text('Obs - NIVA')
        ax8.set_facecolor('gainsboro')
        m8 = Basemap(ax=ax8, projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83,
                     resolution='f')
        m8.fillcontinents(color='darkseagreen')
        x8, y8 = m8(wodlon[csidx], wodlat[csidx])
        m8.scatter(x8, y8, c=mod_niva[csidx] - obs[csidx],
                   cmap=cm_bias, vmin=nmin_bias, vmax=nmax_bias, s=5)
        m8.drawparallels(np.arange(49., 53., 1.), labels=[0, 0, 0, 0], linewidth=0.25)
        m8.drawmeridians(np.arange(-10.5, -5.5, 2.), labels=[0, 0, 0, 1], linewidth=0.25)

        ax9.title.set_text('Obs - CS1K')
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

    def obsVsprods(fdir, obs, mod_ibim, mod_cs1k, mod_niva, clim_brou, clim_mobo, wodlon, wodlat,
                   nutstr, nutunit, nutcmp, qfilt, cs1kfilt):
        # cm = plt.cm.get_cmap('bwr_r')
        # cm = plt.cm.get_cmap('gist_ncar')
        # cm = plt.cm.get_cmap('RdYlBu')
        # cm = plt.cm.get_cmap(nutcmp)
        cm = cmo.matter
        nmin = np.min([(obs[qfilt]).min(),
                       (mod_ibim[qfilt]).min(),
                       (mod_cs1k[cs1kfilt]).min(),
                       (mod_niva[qfilt]).min(),
                       (clim_brou[qfilt]).min()])
        nmax = np.max([(obs[qfilt]).max(),
                       (mod_ibim[qfilt]).max(),
                       (mod_cs1k[cs1kfilt]).max(),
                       (mod_niva[qfilt]).max(),
                       (clim_brou[qfilt]).max()])
        nmin = np.floor(nmin/50.)*50.
        nmax = np.ceil(nmax/50.)*50.
        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(13, 7.5))
        (ax1, ax2, ax3), (ax4, ax5, ax6) = axs
        ax1.scatter(wodlon[qfilt], wodlat[qfilt], c=obs[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax2.scatter(wodlon[qfilt], wodlat[qfilt], c=mod_ibim[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax3.scatter(wodlon[cs1kfilt], wodlat[cs1kfilt], c=mod_cs1k[cs1kfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        ax5.scatter(wodlon[qfilt], wodlat[qfilt], c=mod_niva[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        im = ax6.scatter(wodlon[qfilt], wodlat[qfilt], c=clim_brou[qfilt],
                    cmap=cm, vmin=nmin, vmax=nmax, s=5)
        # im = ax6.scatter(wodlon[qfilt], wodlat[qfilt], c=clim_mobo[qfilt],
        #                  cmap=cm, vmin=nmin, vmax=nmax, s=5)
        fig.suptitle(nutstr + ', ' + nutunit, fontsize=10)
        ax1.title.set_text('Obs')
        ax2.title.set_text('IBI')
        ax3.title.set_text('CS1K')
        ax5.title.set_text('NIVA')
        ax6.title.set_text('BROU')
        # ax6.title.set_text('MOBO')
        cbar = plt.colorbar(im, ax=axs)
        # fig.tight_layout()
        plt_tit = fdir + nutstr + '_spat_obsVsprods.jpg'
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def lt_tseries(fdir, otimes, purp, obs, dat1, dat2, dat3, qflt1, qflt2, qflt3, desc1, desc2, desc3,
                   nutstr, nutunit, aim, climormod):
        # Full timeseries x3 stacked
        # fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(8, 8))
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
            fig.suptitle(nutstr + ', bias by month')
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
            fig.suptitle(nutstr + ', comparison by month')
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
        # fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(8, 8))
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

    (RMSD_ibim_t_DIC, R2_ibim_t_DIC, MD_ibim_t_DIC, N_ibim_t_DIC), \
    (RMSD_cs1k_t_DIC, R2_cs1k_t_DIC, MD_cs1k_t_DIC, N_cs1k_t_DIC), \
    (RMSD_niva_t_DIC, R2_niva_t_DIC, MD_niva_t_DIC, N_niva_t_DIC), \
    (RMSD_brou_t_DIC, R2_brou_t_DIC, MD_brou_t_DIC, N_brou_t_DIC), \
    (RMSD_ibim_1dm_DIC, R2_ibim_1dm_DIC, MD_ibim_1dm_DIC, N_ibim_1dm_DIC), \
    (RMSD_niva_1dm_DIC, R2_niva_1dm_DIC, MD_niva_1dm_DIC, N_niva_1dm_DIC), \
    (RMSD_brou_1dm_DIC, R2_brou_1dm_DIC, MD_brou_1dm_DIC, N_brou_1dm_DIC), \
    (RMSD_ibi_dom_DIC, R2_ibi_dom_DIC, MD_ibi_dom_DIC, N_ibi_dom_DIC), \
    (RMSD_cs1k_dom_DIC, R2_cs1k_dom_DIC, MD_cs1k_dom_DIC, N_cs1k_dom_DIC), \
    (RMSD_brou_cmn_DIC, R2_brou_cmn_DIC, MD_brou_cmn_DIC, N_brou_cmn_DIC), \
    (RMSD_mobo_cmn_DIC, R2_mobo_cmn_DIC, MD_mobo_cmn_DIC, N_mobo_cmn_DIC), \
    (RMSD_ibi_bry_DIC, R2_ibi_bry_DIC, MD_ibi_bry_DIC, N_ibi_bry_DIC), \
    (RMSD_niva_bry_DIC, R2_niva_bry_DIC, MD_niva_bry_DIC, N_niva_bry_DIC), \
    (RMSD_brou_bry_DIC, R2_brou_bry_DIC, MD_brou_bry_DIC, N_brou_bry_DIC) = \
    carb_bias_plt_stats(crocodir, woddic*1.025, woddic_ibim, woddic_cs1k, woddic_niva,
                        woddic_brou*1.025, woddic_mobo*1.025, z,
                        wod_year, wod_month, wod_day, tstamp, wodlon, wodlat, wod_dists, 'DIC')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    dicwriter = pd.ExcelWriter(crocodir + 'DIC_' + npzstr + '_' + hcs + 'stats.xls')

    dicstats_t = pd.DataFrame({'RMSD': [RMSD_ibim_t_DIC, RMSD_cs1k_t_DIC,  RMSD_niva_t_DIC,  RMSD_brou_t_DIC],
                               'R2': [R2_ibim_t_DIC, R2_cs1k_t_DIC, R2_niva_t_DIC, R2_brou_t_DIC],
                               'MD': [MD_ibim_t_DIC, MD_cs1k_t_DIC, MD_niva_t_DIC, MD_brou_t_DIC],
                               'N': [N_ibim_t_DIC, N_cs1k_t_DIC, N_niva_t_DIC, N_brou_t_DIC]})
    dicstats_t.to_excel(dicwriter, sheet_name='total_at_cs1k')

    dicstats_1dm = pd.DataFrame({'RMSD': [RMSD_ibim_1dm_DIC, RMSD_niva_1dm_DIC, RMSD_brou_1dm_DIC],
                                 'R2': [R2_ibim_1dm_DIC, R2_niva_1dm_DIC, R2_brou_1dm_DIC],
                                 'MD': [MD_ibim_1dm_DIC, MD_niva_1dm_DIC, MD_brou_1dm_DIC],
                                 'N': [N_ibim_1dm_DIC, N_niva_1dm_DIC, N_brou_1dm_DIC]})
    dicstats_1dm.to_excel(dicwriter, sheet_name='total_1dm')

    dicstats_dom = pd.DataFrame({'RMSD': [RMSD_ibi_dom_DIC, RMSD_cs1k_dom_DIC],
                                 'R2': [R2_ibi_dom_DIC, R2_cs1k_dom_DIC],
                                 'MD': [MD_ibi_dom_DIC, MD_cs1k_dom_DIC],
                                 'N': [N_ibi_dom_DIC, N_cs1k_dom_DIC]})
    dicstats_dom.to_excel(dicwriter, sheet_name='domain')

    dicstats_bry = pd.DataFrame({'RMSD': [RMSD_ibi_bry_DIC, RMSD_niva_bry_DIC, RMSD_brou_bry_DIC],
                                 'R2': [R2_ibi_bry_DIC, R2_niva_bry_DIC, R2_brou_bry_DIC],
                                 'MD': [MD_ibi_bry_DIC, MD_niva_bry_DIC, MD_brou_bry_DIC],
                                 'N': [N_ibi_bry_DIC, N_niva_bry_DIC, N_brou_bry_DIC]})
    dicstats_bry.to_excel(dicwriter, sheet_name='boundary')

    # Close the Pandas Excel writer and output the Excel file.
    dicwriter.close()

    (RMSD_ibim_t_TALK, R2_ibim_t_TALK, MD_ibim_t_TALK, N_ibim_t_TALK), \
    (RMSD_cs1k_t_TALK, R2_cs1k_t_TALK, MD_cs1k_t_TALK, N_cs1k_t_TALK), \
    (RMSD_niva_t_TALK, R2_niva_t_TALK, MD_niva_t_TALK, N_niva_t_TALK), \
    (RMSD_brou_t_TALK, R2_brou_t_TALK, MD_brou_t_TALK, N_brou_t_TALK), \
    (RMSD_ibim_1dm_TALK, R2_ibim_1dm_TALK, MD_ibim_1dm_TALK, N_ibim_1dm_TALK), \
    (RMSD_niva_1dm_TALK, R2_niva_1dm_TALK, MD_niva_1dm_TALK, N_niva_1dm_TALK), \
    (RMSD_brou_1dm_TALK, R2_brou_1dm_TALK, MD_brou_1dm_TALK, N_brou_1dm_TALK), \
    (RMSD_ibi_dom_TALK, R2_ibi_dom_TALK, MD_ibi_dom_TALK, N_ibi_dom_TALK), \
    (RMSD_cs1k_dom_TALK, R2_cs1k_dom_TALK, MD_cs1k_dom_TALK, N_cs1k_dom_TALK), \
    (RMSD_brou_cmn_TALK, R2_brou_cmn_TALK, MD_brou_cmn_TALK, N_brou_cmn_TALK), \
    (RMSD_mobo_cmn_TALK, R2_mobo_cmn_TALK, MD_mobo_cmn_TALK, N_mobo_cmn_TALK), \
    (RMSD_ibi_bry_TALK, R2_ibi_bry_TALK, MD_ibi_bry_TALK, N_ibi_bry_TALK), \
    (RMSD_niva_bry_TALK, R2_niva_bry_TALK, MD_niva_bry_TALK, N_niva_bry_TALK), \
    (RMSD_brou_bry_TALK, R2_brou_bry_TALK, MD_brou_bry_TALK, N_brou_bry_TALK) = \
    carb_bias_plt_stats(crocodir, wodtlk*1.025, wodtlk_ibim, wodtlk_cs1k, wodtlk_niva,
                        wodtlk_brou*1.025, wodtlk_brou*1.025, z,
                        wod_year, wod_month, wod_day, tstamp, wodlon, wodlat, wod_dists, 'TALK')

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    talkwriter = pd.ExcelWriter(crocodir + 'TALK_' + npzstr + '_' + hcs + 'stats.xls')

    talkstats_t = pd.DataFrame({'RMSD': [RMSD_ibim_t_TALK, RMSD_cs1k_t_TALK, RMSD_niva_t_TALK, RMSD_brou_t_TALK],
                                'R2': [R2_ibim_t_TALK, R2_cs1k_t_TALK, R2_niva_t_TALK, R2_brou_t_TALK],
                                'MD': [MD_ibim_t_TALK, MD_cs1k_t_TALK, MD_niva_t_TALK, MD_brou_t_TALK],
                                'N': [N_ibim_t_TALK, N_cs1k_t_TALK, N_niva_t_TALK, N_brou_t_TALK]})
    talkstats_t.to_excel(talkwriter, sheet_name='total_at_cs1k')

    talkstats_1dm = pd.DataFrame({'RMSD': [RMSD_ibim_1dm_TALK, RMSD_niva_1dm_TALK, RMSD_brou_1dm_TALK],
                                 'R2': [R2_ibim_1dm_TALK, R2_niva_1dm_TALK, R2_brou_1dm_TALK],
                                 'MD': [MD_ibim_1dm_TALK, MD_niva_1dm_TALK, MD_brou_1dm_TALK],
                                 'N': [N_ibim_1dm_TALK, N_niva_1dm_TALK, N_brou_1dm_TALK]})
    talkstats_1dm.to_excel(talkwriter, sheet_name='total_1dm')

    talkstats_dom = pd.DataFrame({'RMSD': [RMSD_ibi_dom_TALK, RMSD_cs1k_dom_TALK],
                                  'R2': [R2_ibi_dom_TALK, R2_cs1k_dom_TALK],
                                  'MD': [MD_ibi_dom_TALK, MD_cs1k_dom_TALK],
                                  'N': [N_ibi_dom_TALK, N_cs1k_dom_TALK]})
    talkstats_dom.to_excel(talkwriter, sheet_name='domain')

    talkstats_bry = pd.DataFrame({'RMSD': [RMSD_ibi_bry_TALK, RMSD_niva_bry_TALK, RMSD_brou_bry_TALK],
                                  'R2': [R2_ibi_bry_TALK, R2_niva_bry_TALK, R2_brou_bry_TALK],
                                  'MD': [MD_ibi_bry_TALK, MD_niva_bry_TALK, MD_brou_bry_TALK],
                                  'N': [N_ibi_bry_TALK, N_niva_bry_TALK, N_brou_bry_TALK]})
    talkstats_bry.to_excel(talkwriter, sheet_name='boundary')

    # Close the Pandas Excel writer and output the Excel file.
    talkwriter.close()

    print('end')
