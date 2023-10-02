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
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter
from netCDF4 import Dataset as netcdf
from scipy import interpolate
import numpy.matlib
import PyCO2SYS as pyco2
import time
import pickle
from pyhdf.SD import SD, SDC
import pandas
from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from math import cos, sin, asin, sqrt, radians
import cmocean.cm as cmo

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

data_processing = 1
data_plotting = 0

# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1017866/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1039109/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1042131/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1050210/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1050468/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1051219/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1052836/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1053899/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1054849/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1065017/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1068817/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1082203/'
# crocodir = '/media/dsktwo/CROCO_BGC_1p2p1_1314_1082862/'
# crocodir = '/media/dskfour/CS1KM_19932021_Uncorr/'
crocodir = '/media/dskfive/CS1KM_19932021_BCd_results/'

if data_processing == 1:

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

    Ystart = 2000
    Mstart = 1
    Dstart = 1

    # Yend = 2000
    Yend = 2021
    Mend = 12
    Dend = 31

    # Satellite L3 chlorophyll/PFT processing
    # satchldir = '/media/dsktwo/VAL/SAT_CHL/L3/'
    satchldir = '/media/dskone/VAL/SAT_CHL_GLO/'

    # Sample satellite file name:
    # '20160601_cmems_obs-oc_atl_bgc-plankton_my_l3-multi-1km_P1D.nc'
    # satchlend = '_cmems_obs-oc_atl_bgc-plankton_my_l3-multi-1km_P1D.nc'
    satchlend = '_c3s_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D.nc'
    #   IBI forcing of constituents (Diatom, Nanophytoplankton)
    #   CROCO-PISCES first run

    ncsc = netcdf(satchldir + '20160101' + satchlend, 'r')
    ncsco = ncsc.variables

    sat_chl_lat = np.array(ncsco['latitude'][:])
    sat_chl_lon = np.array(ncsco['longitude'][:])
    # sat_chl_lat = np.array(ncsco['lat'][:])
    # sat_chl_lon = np.array(ncsco['lon'][:])
    # diat = np.array(ncsco['DIATO'][:])
    # nano = np.array(ncsco['NANO'][:])
    # dino = np.array(ncsco['DINO'][:])
    totchl = np.array(ncsco['CHL'][:])
    sat_chl_lat = np.matlib.repmat(sat_chl_lat, totchl.shape[2], 1).transpose()
    sat_chl_lon = np.matlib.repmat(sat_chl_lon, totchl.shape[1], 1)
    # pico = np.array(ncsco['PICO'][:])
    # micr = np.array(ncsco['MICRO'][:])
    # nano = np.array(ncsco['NANO'][:])
    nano_celt = totchl[0, np.argwhere(np.logical_and(sat_chl_lat > 49, sat_chl_lat < 53))[0, 0]:
                          np.argwhere(np.logical_and(sat_chl_lat > 49, sat_chl_lat < 53))[-1, 0],
                np.argwhere(np.logical_and(sat_chl_lon > -10.75, sat_chl_lon < -5.83))[0, 0]:
                np.argwhere(np.logical_and(sat_chl_lon > -10.75, sat_chl_lon < -5.83))[-1, 0]]
    nano_celt_valid = nano_celt[nano_celt != -999]
    ncsc.close()

    crocofil = 'croco_avg_Y2013M01.nc'
    croco_sta = 'croco_avg_'
    croco_end = '.nc'
    croco_nc = netcdf(crocodir + crocofil, 'r')
    croco_nco = croco_nc.variables
    h_crocod = np.array(croco_nco['h'][:])
    lat_crocod = np.array(croco_nco['lat_rho'][:])
    lon_crocod = np.array(croco_nco['lon_rho'][:])
    mask_crocod = np.array(croco_nco['mask_rho'][:])
    s_rho_crocod = np.array(croco_nco['s_rho'][:])
    talk_crocod = np.array(croco_nco['TALK'][:])
    dic_crocod = np.array(croco_nco['DIC'][:])
    nchl_crocod = np.array(croco_nco['NCHL'][:])
    dchl_crocod = np.array(croco_nco['DCHL'][:])
    no3_crocod = np.array(croco_nco['NO3'][:])
    nh4_crocod = np.array(croco_nco['NH4'][:])
    po4_crocod = np.array(croco_nco['PO4'][:])
    micz_crocod = np.array(croco_nco['ZOO'][:])
    mesz_crocod = np.array(croco_nco['MESO'][:])
    o2_crocod = np.array(croco_nco['O2'][:])
    si_crocod = np.array(croco_nco['Si'][:])
    temp_crocod = np.array(croco_nco['temp'][:])
    salt_crocod = np.array(croco_nco['salt'][:])
    croco_nc.close()

    #   IBI forcing
    #   Daily
    ibidird = '/media/dskone/CELTIC/CMEMS_IBI_PISCESvars_012013_202112/daily/'
    ibidex = 'IBI36_cg_1d-m_20170605-20170605_3DT-bgc_hcst_R20170614.nc'
    ibidsta = 'IBI36_cg_1d-m_'
    ibidmid = '_3DT-bgc_hcst_R'
    ibidnc = netcdf(ibidird + ibidex, 'r')
    ibidnco = ibidnc.variables
    alk_id = np.array(ibidnco['alk'][:])
    dchl_id = np.array(ibidnco['dchl'][:])
    nchl_id = np.array(ibidnco['nchl'][:])
    dic_id = np.array(ibidnco['dic'][:])
    no3_id = np.array(ibidnco['no3'][:])
    po4_id = np.array(ibidnco['po4'][:])
    o2_id = np.array(ibidnco['o2'][:])
    lat_id = np.array(ibidnco['nav_lat'][:])
    lon_id = np.array(ibidnco['nav_lon'][:])
    ibidnc.close()

    #   Monthly
    ibidirm = '/media/dskone/CELTIC/CMEMS_IBI_PISCESvars_012013_202112/monthly/'
    ibimex = 'IBI36_cg_1m-m_201301_3DT-bgc_hcst.nc'
    ibimnc = netcdf(ibidirm + ibimex, 'r')
    ibimnco = ibimnc.variables
    alk_im = np.array(ibimnco['alk'][:])
    dchl_im = np.array(ibimnco['dchl'][:])
    nchl_im = np.array(ibimnco['nchl'][:])
    dic_im = np.array(ibimnco['dic'][:])
    no3_im = np.array(ibimnco['no3'][:])
    po4_im = np.array(ibimnco['po4'][:])
    o2_im = np.array(ibimnco['o2'][:])
    ibimnc.close()

    latr = max(abs(lat_crocod[1, 1] - lat_crocod[2, 1]), abs(lat_crocod[1, 1] - lat_crocod[1, 2]))
    lonr = max(abs(lon_crocod[1, 1] - lon_crocod[2, 1]), abs(lon_crocod[1, 1] - lon_crocod[1, 2]))

    # for 1km cmems atl satellite product
    # ltdelta = 0.5 * latr

    # for 4km GLO satellite product
    ltdelta = 0.042
    ltmin = lat_crocod - ltdelta
    ltmax = lat_crocod + ltdelta

    # for 1km cmems atl satellite product
    # lndelta = 0.5 * lonr

    # for 4km GLO satellite product
    lndelta = 0.042
    lnmin = lon_crocod - lndelta
    lnmax = lon_crocod + lndelta

    sat_croco_pkl = crocodir + 'sat_croco_crs_fin.pkl'

    pickled_sat2croco = 1

    if pickled_sat2croco == 0:
        sat_croco_idx = np.empty((lat_crocod.shape[0], lon_crocod.shape[1]), dtype=object)
        for lts in range(0, lat_crocod.shape[0]):
            sta = time.time()
            for lns in range(0, lon_crocod.shape[1]):
                if mask_crocod[lts, lns] == 1.:
                    sat_croco_idx[lts, lns] = np.argwhere(((sat_chl_lat >= ltmin[lts, lns]) &
                                                           (sat_chl_lat <= ltmax[lts, lns])) &
                                                          ((sat_chl_lon >= lnmin[lts, lns]) &
                                                           (sat_chl_lon <= lnmax[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(sat_croco_pkl, 'wb') as f:
            pickle.dump(sat_croco_idx, f)
    else:
        with open(sat_croco_pkl, 'rb') as f:
            sat_croco_idx = pickle.load(f)

    ibi_croco_pkl = crocodir + 'ibi_croco_fin_crs.pkl'

    lati = max(abs(lat_id[1, 1] - lat_id[2, 1]), abs(lat_id[1, 1] - lat_id[1, 2]))
    loni = max(abs(lon_id[1, 1] - lon_id[2, 1]), abs(lon_id[1, 1] - lon_id[1, 2]))

    ltdeltai = 0.5 * lati
    ltmini = lat_id - ltdeltai
    ltmaxi = lat_id + ltdeltai

    lndeltai = 0.5 * loni
    lnmini = lon_id - lndeltai
    lnmaxi = lon_id + lndeltai

    pickled_ibi2croco = 1
    # need to know what ibi cell each croco coordinate falls within,
    if pickled_ibi2croco == 0:
        ibi_croco_idx = np.empty((lat_crocod.shape[0], lon_crocod.shape[1]), dtype=object)
        for lts in range(0, lat_crocod.shape[0]):
            sta = time.time()
            for lns in range(0, lon_crocod.shape[1]):
                if mask_crocod[lts, lns] == 1.:
                    ibi_croco_idx[lts, lns] = np.argwhere(((ltmini <= lat_crocod[lts, lns]) &
                                                           (ltmaxi >= lat_crocod[lts, lns])) &
                                                          ((lnmini <= lon_crocod[lts, lns]) &
                                                           (lnmaxi >= lon_crocod[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(ibi_croco_pkl, 'wb') as f:
            pickle.dump(ibi_croco_idx, f)
    else:
        with open(ibi_croco_pkl, 'rb') as f:
            ibi_croco_idx = pickle.load(f)

    # Organising collocation of satellite temperature and CROCO
    sattemp_croco_pkl = crocodir + 'sattemp_croco_fin_crs.pkl'

    # Satellite L3 temperature processing
    sattempdir = '/media/dsktwo/VAL/SAT_SST/L3/'

    # Sample satellite file name:
    # '20160601_cmems_obs-oc_atl_bgc-plankton_my_l3-multi-1km_P1D.nc'
    # '20131005000000-IFR-L3S_GHRSST-SSTfnd-ODYSSEA-ATL_MY_005_adjusted-v02.0-f01.0.nc'
    sattempend = '000000-IFR-L3S_GHRSST-SSTfnd-ODYSSEA-ATL_MY_005_adjusted-v02.0-f01.0.nc'

    ncst = netcdf(sattempdir + '20131005' + sattempend, 'r')
    ncsto = ncst.variables
    sat_temp_lat = np.array(ncsto['lat'][:])
    sat_temp_lon = np.array(ncsto['lon'][:])
    sattemp = np.array(ncsto['sea_surface_temperature'][:])
    sat_temp_lat = np.matlib.repmat(sat_temp_lat, sattemp.shape[2], 1).transpose()
    sat_temp_lon = np.matlib.repmat(sat_temp_lon, sattemp.shape[1], 1)
    ncst.close()

    latt = max(abs(sat_temp_lat[1, 1] - sat_temp_lat[2, 1]), abs(sat_temp_lat[1, 1] - sat_temp_lat[1, 2]))
    lont = max(abs(sat_temp_lon[1, 1] - sat_temp_lon[2, 1]), abs(sat_temp_lon[1, 1] - sat_temp_lon[1, 2]))

    ltdeltat = 0.5 * latt
    ltmint = sat_temp_lat - ltdeltat
    ltmaxt = sat_temp_lat + ltdeltat

    lndeltat = 0.5 * lont
    lnmint = sat_temp_lon - lndeltat
    lnmaxt = sat_temp_lon + lndeltat

    pickled_sattemp2croco = 1
    # need to know what satellite grid point each croco coordinate falls within,
    if pickled_sattemp2croco == 0:
        sattemp_croco_idx = np.empty((lat_crocod.shape[0], lon_crocod.shape[1]), dtype=object)
        for lts in range(0, lat_crocod.shape[0]):
            sta = time.time()
            for lns in range(0, lon_crocod.shape[1]):
                if mask_crocod[lts, lns] == 1.:
                    sattemp_croco_idx[lts, lns] = np.argwhere(((ltmint <= lat_crocod[lts, lns]) &
                                                               (ltmaxt >= lat_crocod[lts, lns])) &
                                                              ((lnmint <= lon_crocod[lts, lns]) &
                                                               (lnmaxt >= lon_crocod[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(sattemp_croco_pkl, 'wb') as f:
            pickle.dump(sattemp_croco_idx, f)
    else:
        with open(sattemp_croco_pkl, 'rb') as f:
            sattemp_croco_idx = pickle.load(f)

    # list the satellite data that covers the same period as the CROCO model data

    for valyr in range(Ystart, Yend + 1):
        croco_files = sorted(glob.glob(crocodir + croco_sta + 'Y' + str(valyr) + 'M' + '??' + croco_end))
        sat_chl_files = sorted(glob.glob(satchldir + str(valyr) + '????' + satchlend))
        sat_temp_files = sorted(glob.glob(sattempdir + str(valyr) + '????' + sattempend))
        ibi_files = sorted(glob.glob(ibidird + ibidsta + str(valyr) + '????' + '-' + str(valyr) + '????' +
                                     ibidmid + '????????' + '.nc'))

        # if valyr == Ystart:
        #     croco_files = sorted(glob.glob(crocodir + croco_sta + 'Y' + str(valyr) + 'M' + '??' + croco_end))
        #     sat_chl_files = sorted(glob.glob(satchldir + str(valyr) + '????' + satchlend))
        #     sat_temp_files = sorted(glob.glob(sattempdir + str(valyr) + '????' + sattempend))
        #     ibi_files = sorted(glob.glob(ibidird + ibidsta + str(valyr) + '????' + '-' + str(valyr) + '????' +
        #                                  ibidmid + '????????' + '.nc'))
        # else:
        #     croco_files_yr = sorted(glob.glob(crocodir + croco_sta + 'Y' + str(valyr) + 'M' + '??' + croco_end))
        #     sat_chl_files_yr = sorted(glob.glob(satchldir + str(valyr) + '????' + satchlend))
        #     sat_temp_files_yr = sorted(glob.glob(sattempdir + str(valyr) + '????' + sattempend))
        #     ibi_files_yr = sorted(glob.glob(ibidird + ibidsta + str(valyr) + '????' + '-' + str(valyr) + '????' +
        #                                     ibidmid + '????????' + '.nc'))
        #     croco_files = croco_files + croco_files_yr
        #     sat_chl_files = sat_chl_files + sat_chl_files_yr
        #     sat_temp_files = sat_temp_files + sat_temp_files_yr
        #     ibi_files = ibi_files + ibi_files_yr

        # create matrices for croco, satellite and ibi nanophytoplankton and diatoms and temperature
        croco_dchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        croco_nchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))

        sat_mchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_nchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_pchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_mchle = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_nchle = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_pchle = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_mchlb = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_nchlb = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        sat_pchlb = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))

        chl_valmask = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))

        ibi_dchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        ibi_nchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))

        # another as a mask for processing equivalent croco data

        sat_temp = np.zeros((len(sat_temp_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        croco_temp = np.zeros((len(sat_temp_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
        temp_valmask = np.zeros((len(sat_temp_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))

        sat_processing = 0
        # sat_croco_fil = crocodir + 'sat_croco_nd_chl.npz'
        sat_croco_fil = crocodir + 'sat_croco_nd_chl_' + str(valyr) + '.npz'

        if sat_processing == 1:
            # For satellite data in croco domain:
            # inspect satellite data using indexing
            for ifl in np.arange(0, len(sat_chl_files)):
                sta = time.time()
                iflnc = netcdf(sat_chl_files[ifl], 'r')
                iflnco = iflnc.variables
                s_micr = np.array(iflnco['MICRO'][0, :, :])
                s_nano = np.array(iflnco['NANO'][0, :, :])
                s_pico = np.array(iflnco['PICO'][0, :, :])

                s_micre = np.array(iflnco['MICRO_RMSE'][0, :, :])
                s_nanoe = np.array(iflnco['NANO_RMSE'][0, :, :])
                s_picoe = np.array(iflnco['PICO_RMSE'][0, :, :])

                s_micrb = np.array(iflnco['MICRO_BIAS'][0, :, :])
                s_nanob = np.array(iflnco['NANO_BIAS'][0, :, :])
                s_picob = np.array(iflnco['PICO_BIAS'][0, :, :])

                iflnc.close()
                for lts in range(0, lat_crocod.shape[0]):
                    for lns in range(0, lon_crocod.shape[1]):
                        if mask_crocod[lts, lns] == 1.:
                            # Sat data
                            micr_sample = s_micr[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]
                            nano_sample = s_nano[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]
                            pico_sample = s_pico[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]

                            micr_sample = micr_sample[micr_sample > 0]
                            nano_sample = nano_sample[nano_sample > 0]
                            pico_sample = pico_sample[pico_sample > 0]

                            # Sat data RMSE
                            micr_sample_e = s_micre[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]
                            nano_sample_e = s_nanoe[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]
                            pico_sample_e = s_picoe[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]

                            micr_sample_e = micr_sample_e[micr_sample_e > 0]
                            nano_sample_e = nano_sample_e[nano_sample_e > 0]
                            pico_sample_e = pico_sample_e[pico_sample_e > 0]

                            # Sat data BIAS
                            micr_sample_b = s_micrb[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]
                            nano_sample_b = s_nanob[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]
                            pico_sample_b = s_picob[sat_croco_idx[lts, lns][:, 0], sat_croco_idx[lts, lns][:, 1]]

                            micr_sample_b = micr_sample_b[micr_sample_b > 0]
                            nano_sample_b = nano_sample_b[nano_sample_b > 0]
                            pico_sample_b = pico_sample_b[pico_sample_b > 0]

                            if micr_sample.size != 0:
                                sat_mchl[ifl, lts, lns] = np.mean(micr_sample)
                                sat_nchl[ifl, lts, lns] = np.mean(nano_sample)
                                sat_pchl[ifl, lts, lns] = np.mean(pico_sample)
                                sat_mchle[ifl, lts, lns] = np.mean(micr_sample_e)
                                sat_nchle[ifl, lts, lns] = np.mean(nano_sample_e)
                                sat_pchle[ifl, lts, lns] = np.mean(pico_sample_e)
                                sat_mchlb[ifl, lts, lns] = np.mean(micr_sample_b)
                                sat_nchlb[ifl, lts, lns] = np.mean(nano_sample_b)
                                sat_pchlb[ifl, lts, lns] = np.mean(pico_sample_b)
                dne = time.time()
                print(dne - sta)
            np.savez(sat_croco_fil,
                     sat_mchl=sat_mchl, sat_nchl=sat_nchl, sat_pchl=sat_pchl,
                     sat_mchle=sat_mchle, sat_nchle=sat_nchle, sat_pchle=sat_pchle,
                     sat_mchlb=sat_mchlb, sat_nchlb=sat_nchlb, sat_pchlb=sat_pchlb,
                     chl_valmask=chl_valmask, ibi_dchl=ibi_dchl, ibi_nchl=ibi_nchl)
            # Removing this temporarily to allow for croco processing
            # else:
            #     npzfile = np.load(sat_croco_fil)
            #     sat_mchl = npzfile['sat_mchl']
            #     sat_nchl = npzfile['sat_nchl']
            #     sat_pchl = npzfile['sat_pchl']
            #     sat_mchle = npzfile['sat_mchle']
            #     sat_nchle = npzfile['sat_nchle']
            #     sat_pchle = npzfile['sat_pchle']
            #     sat_mchlb = npzfile['sat_mchlb']
            #     sat_nchlb = npzfile['sat_nchlb']
            #     sat_pchlb = npzfile['sat_pchlb']

    # Satellite temperature processing
    sattemp_processing = 0
    sattemp_croco_fil = crocodir + 'sat_croco_temp.npz'

    if sattemp_processing == 1:
        # For satellite data in croco domain:
        # inspect satellite data using indexing
        # for ifl in range(0, len(sat_temp_files)):
        for ifl in progressbar(range(0, len(sat_temp_files))):
            sta = time.time()
            iflnc = netcdf(sat_temp_files[ifl], 'r')
            iflnco = iflnc.variables
            # Read satellite data and convert to celsius
            s_temp = np.array(iflnco['adjusted_sea_surface_temperature'][0, :, :]) - 273.15
            iflnc.close()
            for lts in range(0, lat_crocod.shape[0]):
                for lns in range(0, lon_crocod.shape[1]):
                    if mask_crocod[lts, lns] == 1.:
                        sattemp_sample = s_temp[sattemp_croco_idx[lts, lns][:, 0], sattemp_croco_idx[lts, lns][:, 1]]
                        sattemp_sample = sattemp_sample[sattemp_sample > 0]
                        if sattemp_sample.size != 0:
                            sat_temp[ifl, lts, lns] = np.mean(sattemp_sample)
            dne = time.time()
            print(dne - sta)
        np.savez(sattemp_croco_fil, sat_temp=sat_temp)
    else:
        npzfile = np.load(sattemp_croco_fil)
        sat_temp = npzfile['sat_temp']

    # Croco temperature processing
    crocotemp_processing = 0
    crocotemp_fil = crocodir + 'croco_temp.npz'

    if crocotemp_processing == 1:
        # Process croco files into larger array
        idx_mth = 0
        for cfl in range(0, len(croco_files)):
            cflnc = netcdf(croco_files[cfl], 'r')
            cflnco = cflnc.variables
            s_temp = np.array(cflnco['temp'][:, cflnco['temp'].shape[1] - 1, :, :])
            idx_mth_end = idx_mth + cflnco['temp'].shape[0]
            croco_temp[idx_mth:idx_mth_end, :, :] = s_temp
            idx_mth = idx_mth_end
            cflnc.close()
        # use mask to extract croco data at same location
        temp_valmask[sat_temp > 0] = 1
        croco_temp[temp_valmask == 0] = 0
        np.savez(crocotemp_fil, croco_temp=croco_temp)
    else:
        # use mask to extract croco data at same location
        npzfile = np.load(crocotemp_fil)
        croco_temp = npzfile['croco_temp']

    croco_processing = 0
    croco_fil = crocodir + 'croco_nd_chl.npz'

    if croco_processing == 1:
        # Process croco files into larger array
        idx_mth = 0

        for valyr in range(Ystart, Yend + 1):
            idx_mth = 0
            croco_fil = crocodir + 'croco_nd_chl_' + str(valyr) + '.npz'
            sat_croco_fil = crocodir + 'sat_croco_nd_chl_' + str(valyr) + '.npz'
            #
            croco_files = sorted(glob.glob(crocodir + croco_sta + 'Y' + str(valyr) + 'M' + '??' + croco_end))
            sat_chl_files = sorted(glob.glob(satchldir + str(valyr) + '????' + satchlend))
            sat_temp_files = sorted(glob.glob(sattempdir + str(valyr) + '????' + sattempend))
            ibi_files = sorted(glob.glob(ibidird + ibidsta + str(valyr) + '????' + '-' + str(valyr) + '????' +
                                         ibidmid + '????????' + '.nc'))
            croco_dchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
            croco_nchl = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
            chl_valmask = np.zeros((len(sat_chl_files), dchl_crocod.shape[2], dchl_crocod.shape[3]))
            sta = time.time()
            for cfl in range(0, len(croco_files)):
                cflnc = netcdf(croco_files[cfl], 'r')
                cflnco = cflnc.variables
                # # Extracting surface layer chlorophyll
                # s_diat = np.array(cflnco['DCHL'][:, cflnco['DCHL'].shape[1] - 1, :, :])
                # s_nano = np.array(cflnco['NCHL'][:, cflnco['NCHL'].shape[1] - 1, :, :])
                # Extracting maximum chlorophyll (assumed within the thermocline)
                s_diat = np.max(np.array(cflnco['DCHL'][:, :, :, :]), 1)
                s_nano = np.max(np.array(cflnco['NCHL'][:, :, :, :]), 1)
                idx_mth_end = idx_mth + cflnco['DCHL'].shape[0]
                croco_dchl[idx_mth:idx_mth_end, :, :] = s_diat
                croco_nchl[idx_mth:idx_mth_end, :, :] = s_nano
                idx_mth = idx_mth_end
                cflnc.close()

            dne = time.time()
            # print(dne - sta)
            print('The year', str(valyr), 'took', np.divide(dne - sta, 60), 'mins')

            # use mask to extract croco data at same location
            npzfile = np.load(sat_croco_fil)
            sat_mchl = npzfile['sat_mchl']
            chl_valmask[sat_mchl > 0] = 1
            croco_dchl[chl_valmask == 0] = 0
            croco_nchl[chl_valmask == 0] = 0
            np.savez(croco_fil, croco_dchl=croco_dchl, croco_nchl=croco_nchl)
        else:
            # use mask to extract croco data at same location
            npzfile = np.load(croco_fil)
            croco_dchl = npzfile['croco_dchl']
            croco_nchl = npzfile['croco_nchl']

    # then use indexing to extract data from ibi file for comparison
    ibi_processing = 0
    ibi_fil = crocodir + 'ibi_nd_chl.npz'

    if ibi_processing == 1:
        for ifl in range(0, len(ibi_files)):
            sta = time.time()
            iflnc = netcdf(ibi_files[ifl], 'r')
            iflnco = iflnc.variables
            # s_diat = np.array(iflnco['dchl'][0, 0, :, :])
            # s_nano = np.array(iflnco['nchl'][0, 0, :, :])
            s_diat = np.array(iflnco['dchl'][0, :, :, :])
            s_diat[s_diat > 1.e36] = np.nan
            s_nano = np.array(iflnco['nchl'][0, :, :, :])
            s_nano[s_nano > 1.e36] = np.nan
            s_diat = np.nanmax(s_diat, 0)
            s_nano = np.nanmax(s_nano, 0)
            iflnc.close()
            for lts in range(0, lat_crocod.shape[0]):
                for lns in range(0, lon_crocod.shape[1]):
                    if mask_crocod[lts, lns] == 1.:
                        diat_sample = s_diat[ibi_croco_idx[lts, lns][:, 0], ibi_croco_idx[lts, lns][:, 1]]
                        nano_sample = s_nano[ibi_croco_idx[lts, lns][:, 0], ibi_croco_idx[lts, lns][:, 1]]
                        diat_sample = diat_sample[diat_sample > 0]
                        nano_sample = nano_sample[nano_sample > 0]
                        if diat_sample.size != 0:
                            ibi_dchl[ifl, lts, lns] = np.mean(diat_sample)
                            ibi_nchl[ifl, lts, lns] = np.mean(nano_sample)
            dne = time.time()
            print(dne - sta)

        # use mask to extract croco data at same location
        chl_valmask[sat_dchl > 0] = 1

        for ifl in range(0, len(ibi_files)):
            for lts in range(0, lat_crocod.shape[0]):
                for lns in range(0, lon_crocod.shape[1]):
                    if mask_crocod[lts, lns] == 1. and \
                            chl_valmask[ifl, lts, lns] == 1 and \
                            ibi_dchl[ifl, lts, lns] == 0:
                        if 0 < lts < (lat_crocod.shape[0] - 1) and 0 < lns < (lon_crocod.shape[1] - 1):
                            ibi_dchl[ifl, lts, lns] = (ibi_dchl[ifl, lts - 1, lns - 1]
                                                       + ibi_dchl[ifl, lts - 1, lns]
                                                       + ibi_dchl[ifl, lts - 1, lns + 1]
                                                       + ibi_dchl[ifl, lts, lns - 1]
                                                       + ibi_dchl[ifl, lts, lns + 1]
                                                       + ibi_dchl[ifl, lts + 1, lns - 1]
                                                       + ibi_dchl[ifl, lts + 1, lns]
                                                       + ibi_dchl[ifl, lts + 1, lns + 1]) / 8
                        elif lts == 0 and 0 < lns < (lon_crocod.shape[1] - 1):
                            ibi_dchl[ifl, lts, lns] = (ibi_dchl[ifl, lts, lns - 1]
                                                       + ibi_dchl[ifl, lts, lns + 1]
                                                       + ibi_dchl[ifl, lts + 1, lns - 1]
                                                       + ibi_dchl[ifl, lts + 1, lns]
                                                       + ibi_dchl[ifl, lts + 1, lns + 1]) / 5
                        elif lts == (lat_crocod.shape[0] - 1) and 0 < lns < (lon_crocod.shape[1] - 1):
                            ibi_dchl[ifl, lts, lns] = (ibi_dchl[ifl, lts - 1, lns - 1]
                                                       + ibi_dchl[ifl, lts - 1, lns]
                                                       + ibi_dchl[ifl, lts - 1, lns + 1]
                                                       + ibi_dchl[ifl, lts, lns - 1]
                                                       + ibi_dchl[ifl, lts, lns + 1]) / 5
                        elif 0 < lts < (lat_crocod.shape[0] - 1) and lns == 0:
                            ibi_dchl[ifl, lts, lns] = (ibi_dchl[ifl, lts - 1, lns]
                                                       + ibi_dchl[ifl, lts - 1, lns + 1]
                                                       + ibi_dchl[ifl, lts, lns + 1]
                                                       + ibi_dchl[ifl, lts + 1, lns]
                                                       + ibi_dchl[ifl, lts + 1, lns + 1]) / 5
                        elif 0 < lts < (lat_crocod.shape[0] - 1) and lns == (lon_crocod.shape[1] - 1):
                            ibi_dchl[ifl, lts, lns] = (ibi_dchl[ifl, lts - 1, lns - 1]
                                                       + ibi_dchl[ifl, lts - 1, lns]
                                                       + ibi_dchl[ifl, lts, lns - 1]
                                                       + ibi_dchl[ifl, lts + 1, lns - 1]
                                                       + ibi_dchl[ifl, lts + 1, lns]) / 5

        ibi_dchl[chl_valmask == 0] = 0
        ibi_nchl[chl_valmask == 0] = 0
        np.savez(ibi_fil, ibi_dchl=ibi_dchl, ibi_nchl=ibi_nchl)
    else:
        # use mask to extract croco data at same location
        npzfile = np.load(ibi_fil)
        ibi_dchl = npzfile['ibi_dchl']
        ibi_nchl = npzfile['ibi_nchl']

    # Compare croco data with satellite, IBI
    # Preparing all data for plotting
    for valyr in range(Ystart, Yend + 1):
        # Filtering satellite by year

        # Filtering CROCO by year

        sat_ok = np.logical_and(sat_mchl > 0, sat_mchl < 1000)
        # ibi_ok = np.logical_and(ibi_dchl > 0, ibi_dchl < 1000)
        # ibi_sat_ok = np.logical_and(sat_ok == 1, ibi_ok == 1)
        croco_ok = np.logical_and(croco_dchl > 0, croco_dchl < 1000)
        # isc_ok = np.logical_and(ibi_sat_ok == 1, croco_ok == 1)
        isc_ok = np.logical_and(sat_ok == 1, croco_ok == 1)

        chl_valmask[isc_ok == 1] = 1
        chl_valmask_tot = np.sum(chl_valmask, axis=0)

        # Diatoms (aka Microplankton)
        croco_dchl_tot = np.sum(croco_dchl * isc_ok, axis=0)
        croco_dchl_avg = croco_dchl_tot / chl_valmask_tot
        # ibi_dchl_tot = np.sum(ibi_dchl * isc_ok, axis=0)
        # ibi_dchl_avg = ibi_dchl_tot / chl_valmask_tot
        sat_dchl_tot = np.sum(sat_mchl * isc_ok, axis=0)
        sat_dchl_avg = sat_dchl_tot / chl_valmask_tot

        # Nanophytoplankton
        croco_nchl_tot = np.sum(croco_nchl * isc_ok, axis=0)
        croco_nchl_avg = croco_nchl_tot / chl_valmask_tot
        # ibi_nchl_tot = np.sum(ibi_nchl * isc_ok, axis=0)
        # ibi_nchl_avg = ibi_nchl_tot / chl_valmask_tot
        sat_nchl_tot = np.sum(sat_nchl * isc_ok, axis=0)
        sat_nchl_avg = sat_nchl_tot / chl_valmask_tot



    fig = plt.figure(figsize=(8, 8))
    # ax = plt.contourf(lon_crocod, lat_crocod, croco_nchl_avg, [20,  50])
    ax = plt.contourf(lon_crocod, lat_crocod, croco_nchl_avg, [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
    ax = plt.contourf(lon_crocod, lat_crocod, sat_nchl_avg, [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
    # ax = plt.contourf(lon_crocod, lat_crocod, ibi_nchl_avg, [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])

    valdir = '/media/dskone/VAL/'

    # WOA18 processing
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

    wlat_trm = np.argwhere((np.min(lat_crocod) <= woa_lat) & (woa_lat <= np.max(lat_crocod)))[:, 0]
    wlon_trm = np.argwhere((np.min(lon_crocod) <= woa_lon) & (woa_lon <= np.max(lon_crocod)))[:, 0]

    woa_lat = woa_lat[wlat_trm]
    woa_lon = woa_lon[wlon_trm]
    n_01b = n_01a[:, :, wlat_trm, :]
    n_01b = n_01b[:, :, :, wlon_trm]

    woa_lat = np.matlib.repmat(woa_lat, n_01b.shape[3], 1).transpose()
    woa_lon = np.matlib.repmat(woa_lon, n_01b.shape[2], 1)

    latw = max(abs(woa_lat[1, 1] - woa_lat[2, 1]), abs(woa_lat[1, 1] - woa_lat[1, 2]))
    lonw = max(abs(woa_lon[1, 1] - woa_lon[2, 1]), abs(woa_lon[1, 1] - woa_lon[1, 2]))

    ltdeltaw = 0.5 * latw
    ltminw = woa_lat - ltdeltaw
    ltmaxw = woa_lat + ltdeltaw

    lndeltaw = 0.5 * lonw
    lnminw = woa_lon - lndeltaw
    lnmaxw = woa_lon + lndeltaw

    pickled_croco2woa = 1
    woa_croco_pkl = crocodir + 'croco_woa_fin_crs.pkl'

    # need to know what croco cells are in each woa coordinate
    if pickled_croco2woa == 0:
        croco_woa_idx = np.empty((woa_lat.shape[0], woa_lon.shape[1]), dtype=object)
        for lts in range(0, woa_lat.shape[0]):
            sta = time.time()
            for lns in range(0, woa_lon.shape[1]):
                croco_woa_idx[lts, lns] = np.argwhere(((lat_crocod >= ltminw[lts, lns]) &
                                                       (lat_crocod <= ltmaxw[lts, lns])) &
                                                      ((lon_crocod >= lnminw[lts, lns]) &
                                                       (lon_crocod <= lnmaxw[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(woa_croco_pkl, 'wb') as f:
            pickle.dump(croco_woa_idx, f)
    else:
        with open(woa_croco_pkl, 'rb') as f:
            croco_woa_idx = pickle.load(f)

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
        woa_depth_nps = np.array(woan_01o['depth'][:])
        woa_depth_nps_bnds = np.array(woan_01o['depth_bnds'][:])
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

    # EMODNET climatology processing
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

    clat_trm = np.argwhere((np.min(lat_crocod) <= emc_lat) & (emc_lat <= np.max(lat_crocod)))[:, 0]
    clon_trm = np.argwhere((np.min(lon_crocod) <= emc_lon) & (emc_lon <= np.max(lon_crocod)))[:, 0]

    emc_lat = emc_lat[clat_trm]
    emc_lon = emc_lon[clon_trm]
    n_01b = n_01a[:, :, clat_trm, :]
    n_01b = n_01b[:, :, :, clon_trm]

    emc_lat = np.matlib.repmat(emc_lat, n_01b.shape[3], 1).transpose()
    emc_lon = np.matlib.repmat(emc_lon, n_01b.shape[2], 1)

    latc = max(abs(emc_lat[1, 1] - emc_lat[2, 1]), abs(emc_lat[1, 1] - emc_lat[1, 2]))
    lonc = max(abs(emc_lon[1, 1] - emc_lon[2, 1]), abs(emc_lon[1, 1] - emc_lon[1, 2]))

    ltdeltac = 0.5 * latc
    ltminc = emc_lat - ltdeltac
    ltmaxc = emc_lat + ltdeltac

    lndeltac = 0.5 * lonc
    lnminc = emc_lon - lndeltac
    lnmaxc = emc_lon + lndeltac

    pickled_croco2emc = 1
    emc_croco_pkl = crocodir + 'croco_emc_fin_crs.pkl'

    # need to know what croco cells are in each emc coordinate
    if pickled_croco2emc == 0:
        croco_emc_idx = np.empty((emc_lat.shape[0], emc_lon.shape[1]), dtype=object)
        for lts in range(0, emc_lat.shape[0]):
            sta = time.time()
            for lns in range(0, emc_lon.shape[1]):
                croco_emc_idx[lts, lns] = np.argwhere(((lat_crocod >= ltminc[lts, lns]) &
                                                       (lat_crocod <= ltmaxc[lts, lns])) &
                                                      ((lon_crocod >= lnminc[lts, lns]) &
                                                       (lon_crocod <= lnmaxc[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(emc_croco_pkl, 'wb') as f:
            pickle.dump(croco_emc_idx, f)
    else:
        with open(emc_croco_pkl, 'rb') as f:
            croco_emc_idx = pickle.load(f)

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

    #   NORESM forcing
    noresmdir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    n_nor = 'NORESMOCv1p2_C2CMI_1993_2020_no3no2_reanal.nc'
    n_nor_nc = netcdf(noresmdir + n_nor, 'r')
    n_nor_nco = n_nor_nc.variables
    no3no2_n = np.array(n_nor_nco['no3no2'][:])
    lat_n = np.array(n_nor_nco['lat'][:])
    lon_n = np.array(n_nor_nco['lon'][:])
    n_nor_nc.close()

    p_nor = 'NORESMOCv1p2_C2CMI_1993_2020_po4_reanal.nc'
    p_nor_nc = netcdf(noresmdir + p_nor, 'r')
    p_nor_nco = p_nor_nc.variables
    po4_n = np.array(p_nor_nco['po4'][:])
    p_nor_nc.close()

    o_nor = 'NORESMOCv1p2_C2CMI_1993_2020_o2_reanal.nc'
    o_nor_nc = netcdf(noresmdir + o_nor, 'r')
    o_nor_nco = o_nor_nc.variables
    o2_n = np.array(o_nor_nco['o2'][:])
    o_nor_nc.close()

    s_nor = 'NORESMOCv1p2_C2CMI_1993_2020_si_reanal.nc'
    s_nor_nc = netcdf(noresmdir + s_nor, 'r')
    s_nor_nco = s_nor_nc.variables
    si_n = np.array(s_nor_nco['si'][:])
    s_nor_nc.close()

    pickled_noresm2woa = 1
    woa_noresm_pkl = crocodir + 'noresm_woa_fin_crs.pkl'
    # need to know what noresm cells are in each woa coordinate
    if pickled_noresm2woa == 0:
        noresm_woa_idx = np.empty((woa_lat.shape[0], woa_lon.shape[1]), dtype=object)
        for lts in range(0, woa_lat.shape[0]):
            sta = time.time()
            for lns in range(0, woa_lon.shape[1]):
                noresm_woa_idx[lts, lns] = np.argwhere(((lat_n >= ltminw[lts, lns]) &
                                                        (lat_n <= ltmaxw[lts, lns])) &
                                                       ((lon_n >= lnminw[lts, lns]) &
                                                        (lon_n <= lnmaxw[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(woa_noresm_pkl, 'wb') as f:
            pickle.dump(noresm_woa_idx, f)
    else:
        with open(woa_noresm_pkl, 'rb') as f:
            noresm_woa_idx = pickle.load(f)

    pickled_ibi2woa = 1
    woa_ibi_pkl = crocodir + 'ibi_woa_fin_crs.pkl'
    # need to know what ibi cells are in each woa coordinate
    if pickled_ibi2woa == 0:
        ibi_woa_idx = np.empty((woa_lat.shape[0], woa_lon.shape[1]), dtype=object)
        for lts in range(0, woa_lat.shape[0]):
            sta = time.time()
            for lns in range(0, woa_lon.shape[1]):
                ibi_woa_idx[lts, lns] = np.argwhere(((lat_id >= ltminw[lts, lns]) &
                                                     (lat_id <= ltmaxw[lts, lns])) &
                                                    ((lon_id >= lnminw[lts, lns]) &
                                                     (lon_id <= lnmaxw[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(woa_ibi_pkl, 'wb') as f:
            pickle.dump(ibi_woa_idx, f)
    else:
        with open(woa_ibi_pkl, 'rb') as f:
            ibi_woa_idx = pickle.load(f)

    # Carbon system checks
    socatfil = valdir + 'SOCATv2022_tracks_gridded_monthly.nc'
    soc_nc = netcdf(socatfil, 'r')
    soc_nco = soc_nc.variables
    spco2_socat = np.array(soc_nco['fco2_ave_weighted'][:])
    spco2_socat[spco2_socat < 0] = np.nan
    time_socat = np.array(soc_nco['tmnth'][:])  # days since 1st Jan 1970
    tsta_socat = ((Ystart - 1970) * 12)
    tend_socat = (((Yend - 1970) + 1) * 12)
    xlon_socat = np.array(soc_nco['xlon'][:])
    ylat_socat = np.array(soc_nco['ylat'][:])
    spco2_socat_mth = spco2_socat[tsta_socat:tend_socat, :, :]
    spco2_socat_mth = spco2_socat_mth[:, wlat_trm, :]
    spco2_socat_mth = spco2_socat_mth[:, :, wlon_trm]
    soc_nc.close()

    # WOA trim to Celtic Sea limits
    woa_n = woa_n[:, :, wlat_trm, :]
    woa_n = woa_n[:, :, :, wlon_trm]
    woa_p = woa_p[:, :, wlat_trm, :]
    woa_p = woa_p[:, :, :, wlon_trm]
    woa_s = woa_s[:, :, wlat_trm, :]
    woa_s = woa_s[:, :, :, wlon_trm]
    woa_o = woa_o[:, :, wlat_trm, :]
    woa_o = woa_o[:, :, :, wlon_trm]
    woa_n_mth = np.tile(woa_n, (2, 1, 1, 1))
    woa_p_mth = np.tile(woa_p, (2, 1, 1, 1))
    woa_s_mth = np.tile(woa_s, (2, 1, 1, 1))
    woa_o_mth = np.tile(woa_o, (2, 1, 1, 1))

    # EMODNET Climatology trimmed to Celtic Sea limits
    emc_n = emc_n[:, :, clat_trm, :]
    emc_n = emc_n[:, :, :, clon_trm]
    emc_p = emc_p[:, :, clat_trm, :]
    emc_p = emc_p[:, :, :, clon_trm]
    emc_s = emc_s[:, :, clat_trm, :]
    emc_s = emc_s[:, :, :, clon_trm]
    emc_o = emc_o[:, :, clat_trm, :]
    emc_o = emc_o[:, :, :, clon_trm]

    emc_n_mth = np.tile(emc_n, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_p_mth = np.tile(emc_p, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_s_mth = np.tile(emc_s, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_o_mth = np.tile(emc_o, ((Yend - Ystart) + 1, 1, 1, 1))

    # L1 product
    emc_n_L1 = emc_n_L1[:, :, clat_trm, :]
    emc_n_L1 = emc_n_L1[:, :, :, clon_trm]
    emc_p_L1 = emc_p_L1[:, :, clat_trm, :]
    emc_p_L1 = emc_p_L1[:, :, :, clon_trm]
    emc_s_L1 = emc_s_L1[:, :, clat_trm, :]
    emc_s_L1 = emc_s_L1[:, :, :, clon_trm]
    emc_o_L1 = emc_o_L1[:, :, clat_trm, :]
    emc_o_L1 = emc_o_L1[:, :, :, clon_trm]

    emc_n_mth_L1 = np.tile(emc_n_L1, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_p_mth_L1 = np.tile(emc_p_L1, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_s_mth_L1 = np.tile(emc_s_L1, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_o_mth_L1 = np.tile(emc_o_L1, ((Yend - Ystart) + 1, 1, 1, 1))

    # L2 product
    emc_n_L2 = emc_n_L2[:, :, clat_trm, :]
    emc_n_L2 = emc_n_L2[:, :, :, clon_trm]
    emc_p_L2 = emc_p_L2[:, :, clat_trm, :]
    emc_p_L2 = emc_p_L2[:, :, :, clon_trm]
    emc_s_L2 = emc_s_L2[:, :, clat_trm, :]
    emc_s_L2 = emc_s_L2[:, :, :, clon_trm]
    emc_o_L2 = emc_o_L2[:, :, clat_trm, :]
    emc_o_L2 = emc_o_L2[:, :, :, clon_trm]

    emc_n_mth_L2 = np.tile(emc_n_L2, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_p_mth_L2 = np.tile(emc_p_L2, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_s_mth_L2 = np.tile(emc_s_L2, ((Yend - Ystart) + 1, 1, 1, 1))
    emc_o_mth_L2 = np.tile(emc_o_L2, ((Yend - Ystart) + 1, 1, 1, 1))

    # Eppley-VGPM
    evgpmdir = valdir + 'Eppley-VGPM/'
    evgpmstart = 'eppley.'
    evgpmend = '.hdf'

    # Open file
    eppf = SD(evgpmdir + evgpmstart + '2017001' + evgpmend, SDC.READ)

    # List available SDS datasets
    eppqf = eppf.datasets()

    # Read dataset
    eppley = eppf.select('npp')
    data = eppley[:]
    latsize = np.float64(min(data.shape))
    elat = np.linspace((-90. + (1. / (latsize / 90.))), (90. - (1. / (latsize / 90.))), int(latsize))
    lonsize = np.float64(max(data.shape))
    elon = np.linspace((-180. + (1. / (lonsize / 180.))), (180. - (1. / (lonsize / 180.))), int(lonsize))

    elat_trm = np.argwhere((np.min(lat_crocod) <= elat) & (elat <= np.max(lat_crocod)))[:, 0]
    elon_trm = np.argwhere((np.min(lon_crocod) <= elon) & (elon <= np.max(lon_crocod)))[:, 0]

    elat = elat[elat_trm]
    elon = elon[elon_trm]
    data = data[elat_trm, :]
    data = data[:, elon_trm]

    elat = np.matlib.repmat(elat, data.shape[1], 1).transpose()
    elon = np.matlib.repmat(elon, data.shape[0], 1)

    late = max(abs(elat[1, 1] - elat[2, 1]), abs(elat[1, 1] - elat[1, 2]))
    lone = max(abs(elon[1, 1] - elon[2, 1]), abs(elon[1, 1] - elon[1, 2]))

    ltdeltae = 0.5 * late
    ltmine = elat - ltdeltae
    ltmaxe = elat + ltdeltae

    lndeltae = 0.5 * lone
    lnmine = elon - lndeltae
    lnmaxe = elon + lndeltae

    pickled_croco2eppley = 1
    epp_croco_pkl = crocodir + 'croco_epp_fin_crs.pkl'

    # need to know what croco cells are in each eppley coordinate
    if pickled_croco2eppley == 0:
        croco_epp_idx = np.empty((elat.shape[0], elon.shape[1]), dtype=object)
        for lts in range(0, elat.shape[0]):
            sta = time.time()
            for lns in range(0, elon.shape[1]):
                croco_epp_idx[lts, lns] = np.argwhere(((lat_crocod >= ltmine[lts, lns]) &
                                                       (lat_crocod <= ltmaxe[lts, lns])) &
                                                      ((lon_crocod >= lnmine[lts, lns]) &
                                                       (lon_crocod <= lnmaxe[lts, lns])))
            dne = time.time()
            print(dne - sta)
        with open(epp_croco_pkl, 'wb') as f:
            pickle.dump(croco_epp_idx, f)
    else:
        with open(epp_croco_pkl, 'rb') as f:
            croco_epp_idx = pickle.load(f)

    # Limit insitu data to the spatial and temporal extent of interest
    # Extract file (croco, ibi) programmatically based on date
    # Read respective params as appropriate
    # carry out interpolation of model data vertically to depth of interest

    # Winter nutrients from excel
    insitufile = '/media/dskone/VAL/Insitu_nutrients.xls'
    dfis = pandas.read_excel(insitufile)
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
    issil_cro = np.zeros_like(issil)
    isphs_cro = np.zeros_like(isphs)
    isnit_cro = np.zeros_like(isnit)
    issal_cro = np.zeros_like(issal)
    issil_ibi = np.zeros_like(issil)
    isphs_ibi = np.zeros_like(isphs)
    isnit_ibi = np.zeros_like(isnit)
    issal_ibi = np.zeros_like(issal)

    # for isd in range(0, len(isyr)):
    #     croco_is_files = sorted(glob.glob(crocodir + croco_sta + 'Y' + str(isyr[isd]) + 'M' +
    #                                       str(ismth[isd]).zfill(2) + croco_end))
    #     latidx = glor.geo_idx(islat[isd], lat_crocod[:, 0])
    #     lonidx = glor.geo_idx(islon[isd], lon_crocod[0, :])
    #     croco_is_fil = croco_is_files[0]
    #     croco_is = netcdf(croco_is_fil, 'r')
    #     croco_iso = croco_is.variables
    #     croco_sil = np.array(croco_iso['Si'][isday[isd - 1], :, latidx, lonidx])
    #     croco_no3 = np.array(croco_iso['NO3'][isday[isd - 1], :, latidx, lonidx])
    #     croco_po4 = np.array(croco_iso['PO4'][isday[isd - 1], :, latidx, lonidx])
    #     croco_sal = np.array(croco_iso['salt'][isday[isd - 1], :, latidx, lonidx])
    #     croco_z = np.array(croco_iso['zeta'][isday[isd - 1], latidx, lonidx])
    #     croco_dpths = vgrd.zlevs(h_crocod[latidx, lonidx], croco_z, theta_s, theta_b, 50, N, 'r', vtransform)
    #     no3_is = interpolate.interp1d(croco_dpths * -1, croco_no3, fill_value="extrapolate")
    #     isnit_cro[isd] = no3_is(isdep[isd])
    #     po4_is = interpolate.interp1d(croco_dpths * -1, croco_po4, fill_value="extrapolate")
    #     isphs_cro[isd] = po4_is(isdep[isd])
    #     sil_is = interpolate.interp1d(croco_dpths * -1, croco_sil, fill_value="extrapolate")
    #     issil_cro[isd] = sil_is(isdep[isd])
    #     sal_is = interpolate.interp1d(croco_dpths * -1, croco_sal, fill_value="extrapolate")
    #     issal_cro[isd] = sal_is(isdep[isd])
    #
    #     # extract data at that location, with z param
    #     # interpolate to depth of interest
    #     # save to model and ibi arrays
    #     croco_is.close()

    # GLODAPv2
    glodapfil = valdir + 'GLODAPv2.2022_Merged_Master_File.mat'
    mat = scio.loadmat(glodapfil)
    gyear = mat['G2year'][:, 0]
    gmonth = mat['G2month'][:, 0]
    gday = mat['G2day'][:, 0]
    glat = mat['G2latitude'][:, 0]
    glon = mat['G2longitude'][:, 0]
    gtemp = mat['G2temperature'][:, 0]
    gsalt = mat['G2salinity'][:, 0]
    gdepth = mat['G2depth'][:, 0]
    gnitr = mat['G2nitrate'][:, 0]
    gphos = mat['G2phosphate'][:, 0]
    goxy = mat['G2oxygen'][:, 0]
    gsil = mat['G2silicate'][:, 0]
    gspco2 = mat['G2fco2'][:, 0]
    gdic = mat['G2tco2'][:, 0]
    gtalk = mat['G2talk'][:, 0]
    # gchla = mat['G2chla'][:, 0]
    # trimidx = np.argwhere((gyear >= Ystart) & (gyear <= Yend) & (glat >= 49) &
    #                       (glat <= 52.95) & (glon >= -10.75) & (glon <= -5.83))
    # trimidx = np.argwhere((glat >= 49) & (glat <= 52.95) & (glon >= -10.75) & (glon <= -5.83))
    trimidx = np.argwhere((glat >= 48) & (glat <= 53.95) & (glon >= -11.75) & (glon <= -4.83))
    gyear = gyear[trimidx][:, 0]
    gmonth = gmonth[trimidx][:, 0]
    gday = gday[trimidx][:, 0]
    glat = glat[trimidx][:, 0]
    glon = glon[trimidx][:, 0]
    gtemp = gtemp[trimidx][:, 0]
    gsalt = gsalt[trimidx][:, 0]
    gdepth = gdepth[trimidx][:, 0]
    gnitr = gnitr[trimidx][:, 0]
    gphos = gphos[trimidx][:, 0]
    goxy = goxy[trimidx][:, 0]
    gsil = gsil[trimidx][:, 0]
    gspco2 = gspco2[trimidx][:, 0]
    gdic = gdic[trimidx][:, 0]
    gtalk = gtalk[trimidx][:, 0]
    # gchla = gchla[trimidx][:, 0]

    # # limit data to non-nan, real values - per param basis
    #
    # # Zooplankton weights from Kieran O'Donnell - this data (2016-2022)
    # # won't be valid in current model run (2013-2014)
    # zoofile = '/media/dskone/VAL/WESPAS_Zoo_timeseries_2016_2022.xls'
    # dfz = pandas.read_excel(zoofile)
    # zyr = np.array(dfz['Year'])
    # zmth = np.array(dfz['Month'])
    # zday = np.array(dfz['Day'])
    # zlat = np.array(dfz['Lat'])
    # zlon = np.array(dfz['Long'])
    # zdep = np.array(dfz['Depth'])
    # zsdep = np.array(dfz['SampleDepth'])
    # zdwt = np.array(dfz['DryWeight'])
    # zflt = np.array(dfz['FlowTot'])
    # zlsf = np.array(dfz['LsFiltered'])

    # Landschutzer gridded dataset (NOT climatology)
    lschfil = valdir + 'MPI-SOM_FFN_v2021_NCEI_OCADS.nc'
    lsch_nc = netcdf(lschfil, 'r')
    lsch_nco = lsch_nc.variables
    lsch_spco2 = np.array(lsch_nco['spco2_raw'][:])
    # lsch_time = np.array(lsch_nco['time'][:])  # seconds since 1st January 2000, first date is 15th Jan 1982
    lsch_lat = np.array(lsch_nco['lat'][:])
    lsch_lon = np.array(lsch_nco['lon'][:])
    lsch_date = np.array(lsch_nco['date'][:])
    lsch_spco2[lsch_spco2 > 1000] = np.nan
    lsch_tsta = np.argwhere(lsch_date[:, 0] == Ystart)[0][0]
    lsch_tend = np.argwhere(lsch_date[:, 0] == Yend)[-1][0]
    # lsch_spco2_mth = lsch_spco2[lsch_tsta:lsch_tend, :, :]
    lsch_spco2_mth = lsch_spco2[lsch_tsta:, :, :]
    lsch_spco2_mth[lsch_spco2_mth > 1000] = np.nan
    lsch_spco2_mth = lsch_spco2_mth[0:12 * (Yend + 1 - Ystart), :, :]
    lsch_spco2_mth = lsch_spco2_mth[:, wlat_trm, :]
    lsch_spco2_mth = lsch_spco2_mth[:, :, wlon_trm]
    lsch_nc.close()

    # OceanSODA-ETHZ
    # Time indexing is the same as the Landschutzer
    osefil = valdir + 'OceanSODA-ETHZ/' + 'OceanSODA-ETHZ_GRaCER_v2021a_1982-2020.nc'
    ose_nc = netcdf(osefil, 'r')
    ose_nco = ose_nc.variables
    talk_ose = np.array(ose_nco['talk'][:])
    # talk_ose_mth = talk_ose[lsch_tsta:lsch_tend, :, :]
    talk_ose_mth = talk_ose[lsch_tsta:, :, :]
    talk_ose_mth = talk_ose_mth[0:12 * (Yend + 1 - Ystart), :, :]
    talk_ose_mth = talk_ose_mth[:, wlat_trm, :]
    talk_ose_mth = talk_ose_mth[:, :, wlon_trm]
    dic_ose = np.array(ose_nco['dic'][:])
    # dic_ose_mth = dic_ose[lsch_tsta:lsch_tend, :, :]
    dic_ose_mth = dic_ose[lsch_tsta:, :, :]
    dic_ose_mth = dic_ose_mth[0:12 * (Yend + 1 - Ystart), :, :]
    dic_ose_mth = dic_ose_mth[:, wlat_trm, :]
    dic_ose_mth = dic_ose_mth[:, :, wlon_trm]
    spco2_ose = np.array(ose_nco['spco2'][:])
    # spco2_ose_mth = spco2_ose[lsch_tsta:lsch_tend, :, :]
    spco2_ose_mth = spco2_ose[lsch_tsta:, :, :]
    spco2_ose_mth = spco2_ose_mth[0:12 * (Yend + 1 - Ystart), :, :]
    spco2_ose_mth = spco2_ose_mth[:, wlat_trm, :]
    spco2_ose_mth = spco2_ose_mth[:, :, wlon_trm]
    time_ose = np.array(ose_nco['time'][:])  # days since 1982-01-01
    lat_ose = np.array(ose_nco['lat'][:])
    lon_ose = np.array(ose_nco['lon'][:])
    ose_nc.close()

    # Interested in timeseries and for long term averages
    # Initialise array of croco grid at woa scale, with croco depths
    croco_mth_talk = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_dic = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_no3 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_po4 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_o2 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_si = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_temp = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_salt = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_mth_zeta = np.zeros(((Yend + 1 - Ystart) * 12, lat_crocod.shape[0], lon_crocod.shape[1]))

    croco_clm_talk = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_dic = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_no3 = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_po4 = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_o2 = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_si = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_temp = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_salt = np.zeros((12, no3_crocod.shape[1], lat_crocod.shape[0], lon_crocod.shape[1]))
    croco_clm_zeta = np.zeros((12, lat_crocod.shape[0], lon_crocod.shape[1]))

    # cycle through months
    for mths in range(1, 13):
        croco_flist = sorted(glob.glob(crocodir + croco_sta + 'Y' + '????' + 'M' + str(mths).zfill(2) + croco_end))
        for fil in range(0, len(croco_flist)):
            cr_nc = netcdf(croco_flist[fil], 'r')
            cr_nco = cr_nc.variables
            # get nutrients of interest from netcdf
            talk_cr_mth = np.array(cr_nco['TALK'][:])
            dic_cr_mth = np.array(cr_nco['DIC'][:])
            no3_cr_mth = np.array(cr_nco['NO3'][:])
            po4_cr_mth = np.array(cr_nco['PO4'][:])
            o2_cr_mth = np.array(cr_nco['O2'][:])
            si_cr_mth = np.array(cr_nco['Si'][:])
            temp_cr_mth = np.array(cr_nco['temp'][:])
            salt_cr_mth = np.array(cr_nco['salt'][:])
            zeta_cr_mth = np.array(cr_nco['zeta'][:])
            cr_nc.close()
            # average by month
            croco_mth_talk[(12 * fil) + mths - 1, :, :, :] = np.average(talk_cr_mth, axis=0)
            croco_mth_dic[(12 * fil) + mths - 1, :, :, :] = np.average(dic_cr_mth, axis=0)
            croco_mth_no3[(12 * fil) + mths - 1, :, :, :] = np.average(no3_cr_mth, axis=0)
            croco_mth_po4[(12 * fil) + mths - 1, :, :, :] = np.average(po4_cr_mth, axis=0)
            croco_mth_o2[(12 * fil) + mths - 1, :, :, :] = np.average(o2_cr_mth, axis=0)
            croco_mth_si[(12 * fil) + mths - 1, :, :, :] = np.average(si_cr_mth, axis=0)
            croco_mth_temp[(12 * fil) + mths - 1, :, :, :] = np.average(temp_cr_mth, axis=0)
            croco_mth_salt[(12 * fil) + mths - 1, :, :, :] = np.average(salt_cr_mth, axis=0)
            croco_mth_zeta[(12 * fil) + mths - 1, :, :] = np.average(zeta_cr_mth, axis=0)
            if fil == 0:
                talk_cr = np.average(talk_cr_mth, axis=0)
                dic_cr = np.average(dic_cr_mth, axis=0)
                no3_cr = np.average(no3_cr_mth, axis=0)
                po4_cr = np.average(po4_cr_mth, axis=0)
                o2_cr = np.average(o2_cr_mth, axis=0)
                si_cr = np.average(si_cr_mth, axis=0)
                temp_cr = np.average(temp_cr_mth, axis=0)
                salt_cr = np.average(salt_cr_mth, axis=0)
                zeta_cr = np.average(zeta_cr_mth, axis=0)
            else:
                talk_cr = talk_cr + np.average(talk_cr_mth, axis=0)
                dic_cr = dic_cr + np.average(dic_cr_mth, axis=0)
                no3_cr = no3_cr + np.average(no3_cr_mth, axis=0)
                po4_cr = po4_cr + np.average(po4_cr_mth, axis=0)
                o2_cr = o2_cr + np.average(o2_cr_mth, axis=0)
                si_cr = si_cr + np.average(si_cr_mth, axis=0)
                temp_cr = temp_cr + np.average(temp_cr_mth, axis=0)
                salt_cr = salt_cr + np.average(salt_cr_mth, axis=0)
                zeta_cr = zeta_cr + np.average(zeta_cr_mth, axis=0)
        talk_cr = talk_cr / len(croco_flist)
        dic_cr = dic_cr / len(croco_flist)
        no3_cr = no3_cr / len(croco_flist)
        po4_cr = po4_cr / len(croco_flist)
        o2_cr = o2_cr / len(croco_flist)
        si_cr = si_cr / len(croco_flist)
        temp_cr = temp_cr / len(croco_flist)
        salt_cr = salt_cr / len(croco_flist)
        zeta_cr = zeta_cr / len(croco_flist)
        croco_clm_talk[mths - 1, :, :, :] = talk_cr
        croco_clm_dic[mths - 1, :, :, :] = dic_cr
        croco_clm_no3[mths - 1, :, :, :] = no3_cr
        croco_clm_po4[mths - 1, :, :, :] = po4_cr
        croco_clm_o2[mths - 1, :, :, :] = o2_cr
        croco_clm_si[mths - 1, :, :, :] = si_cr
        croco_clm_temp[mths - 1, :, :, :] = temp_cr
        croco_clm_salt[mths - 1, :, :, :] = salt_cr
        croco_clm_zeta[mths - 1, :, :] = zeta_cr
        talk_cr = []
        dic_cr = []
        no3_cr = []
        po4_cr = []
        o2_cr = []
        si_cr = []
        temp_cr = []
        salt_cr = []
        zeta_cr = []

    # initialise array of croco grid at woa scale, with croco depths
    croco_woa_talk = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_dic = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_no3 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_po4 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_o2 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_si = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_temp = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_salt = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_zeta = np.zeros(((Yend + 1 - Ystart) * 12, woa_lat.shape[0], woa_lon.shape[1]))
    croco_woa_dpth = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))

    # initialise array of croco grid at woa scale, with croco depths
    croco_clm_woa_talk = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_dic = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_no3 = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_po4 = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_o2 = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_si = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_temp = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_salt = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_zeta = np.zeros((12, woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woa_dpth = np.zeros((12, no3_crocod.shape[1], woa_lat.shape[0], woa_lon.shape[1]))

    # block average from croco to woa spatial scale

    # for both monthly gridded data and climatology arrays

    # Establish average bathymetry in woa grid cells
    h_woa = np.zeros((woa_lat.shape[0], woa_lon.shape[1]))
    for wla in range(0, croco_clm_woa_talk.shape[2]):
        for wlo in range(0, croco_clm_woa_talk.shape[3]):
            h_woa[wla, wlo] = np.mean(h_crocod[croco_woa_idx[wla, wlo][:, 0], croco_woa_idx[wla, wlo][:, 1]])

    croco_mth_talk[croco_mth_talk == 0.] = np.nan
    croco_mth_dic[croco_mth_dic == 0.] = np.nan
    croco_mth_no3[croco_mth_no3 == 0.] = np.nan
    croco_mth_po4[croco_mth_po4 == 0.] = np.nan
    croco_mth_o2[croco_mth_o2 == 0.] = np.nan
    croco_mth_si[croco_mth_si == 0.] = np.nan
    croco_mth_temp[croco_mth_temp == 0.] = np.nan
    croco_mth_salt[croco_mth_salt == 0.] = np.nan
    croco_mth_zeta[croco_mth_zeta == 0.] = np.nan

    croco_clm_talk[croco_clm_talk == 0.] = np.nan
    croco_clm_dic[croco_clm_dic == 0.] = np.nan
    croco_clm_no3[croco_clm_no3 == 0.] = np.nan
    croco_clm_po4[croco_clm_po4 == 0.] = np.nan
    croco_clm_o2[croco_clm_o2 == 0.] = np.nan
    croco_clm_si[croco_clm_si == 0.] = np.nan
    croco_clm_temp[croco_clm_temp == 0.] = np.nan
    croco_clm_salt[croco_clm_salt == 0.] = np.nan
    croco_clm_zeta[croco_clm_zeta == 0.] = np.nan

    # loop by month, depth, lat and lon, for each parameter
    for mnt in range(0, croco_woa_talk.shape[0]):
        for dp in range(0, croco_woa_talk.shape[1]):
            for wla in range(0, croco_woa_talk.shape[2]):
                for wlo in range(0, croco_woa_talk.shape[3]):
                    croco_woa_talk[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_talk[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_dic[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_dic[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_no3[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_no3[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_po4[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_po4[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_o2[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_o2[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_si[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_si[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_temp[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_temp[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_salt[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_salt[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_woa_zeta[mnt, wla, wlo] = \
                        np.nanmean(croco_mth_zeta[mnt, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    # at each location get depth in 20 croco layers
                    croco_woa_dpth[mnt, :, wla, wlo] = vgrd.zlevs(h_woa[wla, wlo], croco_woa_zeta[mnt, wla, wlo],
                                                                  theta_s, theta_b, 50, N, 'r', vtransform)

    for mnt in range(0, croco_clm_woa_talk.shape[0]):
        for dp in range(0, croco_clm_woa_talk.shape[1]):
            for wla in range(0, croco_clm_woa_talk.shape[2]):
                for wlo in range(0, croco_clm_woa_talk.shape[3]):
                    croco_clm_woa_talk[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_talk[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_dic[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_dic[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_no3[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_no3[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_po4[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_po4[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_o2[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_o2[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_si[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_si[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_temp[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_temp[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_salt[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_salt[mnt, dp, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    croco_clm_woa_zeta[mnt, wla, wlo] = \
                        np.nanmean(croco_clm_zeta[mnt, croco_woa_idx[wla, wlo][:, 0],
                        croco_woa_idx[wla, wlo][:, 1]])
                    # at each location get depth in 20 croco layers
                    croco_clm_woa_dpth[mnt, :, wla, wlo] = vgrd.zlevs(h_woa[wla, wlo],
                                                                      croco_clm_woa_zeta[mnt, wla, wlo],
                                                                      theta_s, theta_b, 50, N, 'r', vtransform)

    # interpolate onto woa layers (more layers for oxygen)
    # initialise array of croco grid at woa scale, with croco depths
    croco_woad_talk = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_dic = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_no3 = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_po4 = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_o2 = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_ox.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_si = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_temp = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_woad_salt = np.zeros(((Yend + 1 - Ystart) * 12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))

    # initialise array of croco grid at woa scale, with croco depths
    croco_clm_woad_talk = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_dic = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_no3 = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_po4 = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_o2 = np.zeros((12, woa_depth_ox.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_si = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_temp = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))
    croco_clm_woad_salt = np.zeros((12, woa_depth_nps.shape[0], woa_lat.shape[0], woa_lon.shape[1]))

    # for each month and lat/lon, interpolate to woa depths
    # loop by month, depth, lat and lon, for each parameter
    for mnt in range(0, croco_woad_talk.shape[0]):
        for wla in range(0, croco_woad_talk.shape[2]):
            for wlo in range(0, croco_woad_talk.shape[3]):
                talk_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_talk[mnt, :, wla, wlo],
                                              fill_value="extrapolate")
                croco_woad_talk[mnt, :, wla, wlo] = talk_i(woa_depth_nps)
                dic_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_dic[mnt, :, wla, wlo],
                                             fill_value="extrapolate")
                croco_woad_dic[mnt, :, wla, wlo] = dic_i(woa_depth_nps)
                no3_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_no3[mnt, :, wla, wlo],
                                             fill_value="extrapolate")
                croco_woad_no3[mnt, :, wla, wlo] = no3_i(woa_depth_nps)
                po4_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_po4[mnt, :, wla, wlo],
                                             fill_value="extrapolate")
                croco_woad_po4[mnt, :, wla, wlo] = po4_i(woa_depth_nps)
                o2_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_o2[mnt, :, wla, wlo],
                                            fill_value="extrapolate")
                croco_woad_o2[mnt, :, wla, wlo] = o2_i(woa_depth_ox)
                si_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_si[mnt, :, wla, wlo],
                                            fill_value="extrapolate")
                croco_woad_si[mnt, :, wla, wlo] = si_i(woa_depth_nps)
                temp_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_temp[mnt, :, wla, wlo],
                                              fill_value="extrapolate")
                croco_woad_temp[mnt, :, wla, wlo] = temp_i(woa_depth_nps)
                salt_i = interpolate.interp1d(croco_woa_dpth[mnt, :, wla, wlo] * -1, croco_woa_salt[mnt, :, wla, wlo],
                                              fill_value="extrapolate")
                croco_woad_salt[mnt, :, wla, wlo] = salt_i(woa_depth_nps)

    results_croco_woad = pyco2.sys(par1=croco_woad_talk[:, 0, :, :], par1_type=1,
                                   par2=croco_woad_dic[:, 0, :, :], par2_type=2,
                                   temperature=croco_woad_temp[:, 0, :, :],
                                   salinity=croco_woad_salt[:, 0, :, :])
    croco_woad_spco2 = results_croco_woad["pCO2"]

    # for each month and lat/lon, interpolate to woa depths
    # loop by month, depth, lat and lon, for each parameter
    for mnt in range(0, croco_clm_woad_talk.shape[0]):
        for wla in range(0, croco_clm_woad_talk.shape[2]):
            for wlo in range(0, croco_clm_woad_talk.shape[3]):
                talk_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                              croco_clm_woa_talk[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_talk[mnt, :, wla, wlo] = talk_i(woa_depth_nps)
                dic_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                             croco_clm_woa_dic[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_dic[mnt, :, wla, wlo] = dic_i(woa_depth_nps)
                no3_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                             croco_clm_woa_no3[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_no3[mnt, :, wla, wlo] = no3_i(woa_depth_nps)
                po4_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                             croco_clm_woa_po4[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_po4[mnt, :, wla, wlo] = po4_i(woa_depth_nps)
                o2_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                            croco_clm_woa_o2[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_o2[mnt, :, wla, wlo] = o2_i(woa_depth_ox)
                si_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                            croco_clm_woa_si[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_si[mnt, :, wla, wlo] = si_i(woa_depth_nps)
                temp_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                              croco_clm_woa_temp[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_temp[mnt, :, wla, wlo] = temp_i(woa_depth_nps)
                salt_i = interpolate.interp1d(croco_clm_woa_dpth[mnt, :, wla, wlo] * -1,
                                              croco_clm_woa_salt[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_woad_salt[mnt, :, wla, wlo] = salt_i(woa_depth_nps)

    results_croco_clm_woad = pyco2.sys(par1=croco_clm_woad_talk[:, 0, :, :], par1_type=1,
                                       par2=croco_clm_woad_dic[:, 0, :, :], par2_type=2,
                                       temperature=croco_clm_woad_temp[:, 0, :, :],
                                       salinity=croco_clm_woad_salt[:, 0, :, :])
    croco_clm_woad_spco2 = results_croco_clm_woad["pCO2"]

    woa_n[woa_n > 10000.] = np.nan
    woa_n_mth[woa_n_mth > 10000.] = np.nan
    woa_p[woa_p > 10000.] = np.nan
    woa_p_mth[woa_p_mth > 10000.] = np.nan
    woa_s[woa_s > 10000.] = np.nan
    woa_s_mth[woa_s_mth > 10000.] = np.nan
    woa_o[woa_o > 10000.] = np.nan
    woa_o_mth[woa_o_mth > 10000.] = np.nan

    # CROCO nutrient re-gridding to EMODNET climatology grid
    # initialise array of croco grid at emc scale, with croco depths
    croco_emc_talk = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_dic = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_no3 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_po4 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_o2 = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_si = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_temp = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_salt = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_zeta = np.zeros(((Yend + 1 - Ystart) * 12, emc_lat.shape[0], emc_lon.shape[1]))
    croco_emc_dpth = np.zeros(((Yend + 1 - Ystart) * 12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))

    # initialise array of croco grid at emc scale, with croco depths
    croco_clm_emc_talk = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_dic = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_no3 = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_po4 = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_o2 = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_si = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_temp = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_salt = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_zeta = np.zeros((12, emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emc_dpth = np.zeros((12, no3_crocod.shape[1], emc_lat.shape[0], emc_lon.shape[1]))

    # block average from croco to emc spatial scale

    # for both monthly gridded data and climatology arrays

    # Establish average bathymetry in emc grid cells
    h_emc = np.zeros((emc_lat.shape[0], emc_lon.shape[1]))
    for wla in range(0, croco_clm_emc_talk.shape[2]):
        for wlo in range(0, croco_clm_emc_talk.shape[3]):
            h_emc[wla, wlo] = np.mean(h_crocod[croco_emc_idx[wla, wlo][:, 0], croco_emc_idx[wla, wlo][:, 1]])

    croco_mth_talk[croco_mth_talk == 0.] = np.nan
    croco_mth_dic[croco_mth_dic == 0.] = np.nan
    croco_mth_no3[croco_mth_no3 == 0.] = np.nan
    croco_mth_po4[croco_mth_po4 == 0.] = np.nan
    croco_mth_o2[croco_mth_o2 == 0.] = np.nan
    croco_mth_si[croco_mth_si == 0.] = np.nan
    croco_mth_temp[croco_mth_temp == 0.] = np.nan
    croco_mth_salt[croco_mth_salt == 0.] = np.nan
    croco_mth_zeta[croco_mth_zeta == 0.] = np.nan

    croco_clm_talk[croco_clm_talk == 0.] = np.nan
    croco_clm_dic[croco_clm_dic == 0.] = np.nan
    croco_clm_no3[croco_clm_no3 == 0.] = np.nan
    croco_clm_po4[croco_clm_po4 == 0.] = np.nan
    croco_clm_o2[croco_clm_o2 == 0.] = np.nan
    croco_clm_si[croco_clm_si == 0.] = np.nan
    croco_clm_temp[croco_clm_temp == 0.] = np.nan
    croco_clm_salt[croco_clm_salt == 0.] = np.nan
    croco_clm_zeta[croco_clm_zeta == 0.] = np.nan

    # loop by month, depth, lat and lon, for each parameter
    for mnt in range(0, croco_emc_talk.shape[0]):
        for dp in range(0, croco_emc_talk.shape[1]):
            for wla in range(0, croco_emc_talk.shape[2]):
                for wlo in range(0, croco_emc_talk.shape[3]):
                    croco_emc_talk[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_talk[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_dic[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_dic[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_no3[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_no3[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_po4[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_po4[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_o2[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_o2[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_si[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_si[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_temp[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_temp[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_salt[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_mth_salt[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_emc_zeta[mnt, wla, wlo] = \
                        np.nanmean(croco_mth_zeta[mnt, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    # at each location get depth in 20 croco layers
                    croco_emc_dpth[mnt, :, wla, wlo] = vgrd.zlevs(h_emc[wla, wlo], croco_emc_zeta[mnt, wla, wlo],
                                                                  theta_s, theta_b, 50, N, 'r', vtransform)

    for mnt in range(0, croco_clm_emc_talk.shape[0]):
        for dp in range(0, croco_clm_emc_talk.shape[1]):
            for wla in range(0, croco_clm_emc_talk.shape[2]):
                for wlo in range(0, croco_clm_emc_talk.shape[3]):
                    croco_clm_emc_talk[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_talk[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_dic[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_dic[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_no3[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_no3[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_po4[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_po4[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_o2[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_o2[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_si[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_si[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_temp[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_temp[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_salt[mnt, dp, wla, wlo] = \
                        np.nanmean(croco_clm_salt[mnt, dp, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    croco_clm_emc_zeta[mnt, wla, wlo] = \
                        np.nanmean(croco_clm_zeta[mnt, croco_emc_idx[wla, wlo][:, 0],
                        croco_emc_idx[wla, wlo][:, 1]])
                    # at each location get depth in 20 croco layers
                    croco_clm_emc_dpth[mnt, :, wla, wlo] = vgrd.zlevs(h_emc[wla, wlo],
                                                                      croco_clm_emc_zeta[mnt, wla, wlo],
                                                                      theta_s, theta_b, 50, N, 'r', vtransform)

    # interpolate onto emc layers (more layers for oxygen)
    # initialise array of croco grid at emc scale, with croco depths
    croco_emcd_talk = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_dic = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_no3 = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_po4 = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_o2 = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_si = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_temp = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_emcd_salt = np.zeros(((Yend + 1 - Ystart) * 12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))

    # initialise array of croco grid at emc scale, with croco depths
    croco_clm_emcd_talk = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_dic = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_no3 = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_po4 = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_o2 = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_si = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_temp = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))
    croco_clm_emcd_salt = np.zeros((12, emc_depth.shape[0], emc_lat.shape[0], emc_lon.shape[1]))

    # for each month and lat/lon, interpolate to emc depths
    # loop by month, depth, lat and lon, for each parameter
    for mnt in range(0, croco_emcd_talk.shape[0]):
        for wla in range(0, croco_emcd_talk.shape[2]):
            for wlo in range(0, croco_emcd_talk.shape[3]):
                talk_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_talk[mnt, :, wla, wlo],
                                              fill_value="extrapolate")
                croco_emcd_talk[mnt, :, wla, wlo] = talk_i(emc_depth)
                dic_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_dic[mnt, :, wla, wlo],
                                             fill_value="extrapolate")
                croco_emcd_dic[mnt, :, wla, wlo] = dic_i(emc_depth)
                no3_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_no3[mnt, :, wla, wlo],
                                             fill_value="extrapolate")
                croco_emcd_no3[mnt, :, wla, wlo] = no3_i(emc_depth)
                po4_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_po4[mnt, :, wla, wlo],
                                             fill_value="extrapolate")
                croco_emcd_po4[mnt, :, wla, wlo] = po4_i(emc_depth)
                o2_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_o2[mnt, :, wla, wlo],
                                            fill_value="extrapolate")
                croco_emcd_o2[mnt, :, wla, wlo] = o2_i(emc_depth)
                si_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_si[mnt, :, wla, wlo],
                                            fill_value="extrapolate")
                croco_emcd_si[mnt, :, wla, wlo] = si_i(emc_depth)
                temp_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_temp[mnt, :, wla, wlo],
                                              fill_value="extrapolate")
                croco_emcd_temp[mnt, :, wla, wlo] = temp_i(emc_depth)
                salt_i = interpolate.interp1d(croco_emc_dpth[mnt, :, wla, wlo] * -1, croco_emc_salt[mnt, :, wla, wlo],
                                              fill_value="extrapolate")
                croco_emcd_salt[mnt, :, wla, wlo] = salt_i(emc_depth)

    # for each month and lat/lon, interpolate to emc depths
    # loop by month, depth, lat and lon, for each parameter
    for mnt in range(0, croco_clm_emcd_talk.shape[0]):
        for wla in range(0, croco_clm_emcd_talk.shape[2]):
            for wlo in range(0, croco_clm_emcd_talk.shape[3]):
                talk_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                              croco_clm_emc_talk[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_talk[mnt, :, wla, wlo] = talk_i(emc_depth)
                dic_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                             croco_clm_emc_dic[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_dic[mnt, :, wla, wlo] = dic_i(emc_depth)
                no3_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                             croco_clm_emc_no3[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_no3[mnt, :, wla, wlo] = no3_i(emc_depth)
                po4_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                             croco_clm_emc_po4[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_po4[mnt, :, wla, wlo] = po4_i(emc_depth)
                o2_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                            croco_clm_emc_o2[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_o2[mnt, :, wla, wlo] = o2_i(emc_depth)
                si_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                            croco_clm_emc_si[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_si[mnt, :, wla, wlo] = si_i(emc_depth)
                temp_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                              croco_clm_emc_temp[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_temp[mnt, :, wla, wlo] = temp_i(emc_depth)
                salt_i = interpolate.interp1d(croco_clm_emc_dpth[mnt, :, wla, wlo] * -1,
                                              croco_clm_emc_salt[mnt, :, wla, wlo], fill_value="extrapolate")
                croco_clm_emcd_salt[mnt, :, wla, wlo] = salt_i(emc_depth)

    emc_n[emc_n > 10000.] = np.nan
    emc_n_mth[emc_n_mth > 10000.] = np.nan
    emc_p[emc_p > 10000.] = np.nan
    emc_p_mth[emc_p_mth > 10000.] = np.nan
    emc_s[emc_s > 10000.] = np.nan
    emc_s_mth[emc_s_mth > 10000.] = np.nan
    emc_o[emc_o > 10000.] = np.nan
    emc_o_mth[emc_o_mth > 10000.] = np.nan

    ngf = crocodir + 'nutrients_grid.npz'
    nif = crocodir + 'nutrients_insitu.npz'
    cgf = crocodir + 'carbon_grid.npz'

    np.savez(ngf,
             croco_clm_woad_si=croco_clm_woad_si, croco_clm_woad_o2=croco_clm_woad_o2,
             croco_clm_woad_no3=croco_clm_woad_no3, croco_clm_woad_po4=croco_clm_woad_po4,
             croco_woad_si=croco_woad_si, croco_woad_o2=croco_woad_o2,
             croco_woad_no3=croco_woad_no3, croco_woad_po4=croco_woad_po4,
             woa_n=woa_n, woa_p=woa_p, woa_s=woa_s, woa_o=woa_o,
             woa_depth_nps=woa_depth_nps, woa_depth_ox=woa_depth_ox,
             woa_n_mth=woa_n_mth, woa_p_mth=woa_p_mth, woa_s_mth=woa_s_mth, woa_o_mth=woa_o_mth,
             woa_lat=woa_lat, woa_lon=woa_lon, lat_crocod=lat_crocod, lon_crocod=lon_crocod, h_crocod=h_crocod,
             croco_clm_emcd_si=croco_clm_emcd_si, croco_clm_emcd_o2=croco_clm_emcd_o2,
             croco_clm_emcd_no3=croco_clm_emcd_no3, croco_clm_emcd_po4=croco_clm_emcd_po4,
             croco_emcd_si=croco_emcd_si, croco_emcd_o2=croco_emcd_o2,
             croco_emcd_no3=croco_emcd_no3, croco_emcd_po4=croco_emcd_po4,
             emc_n=emc_n, emc_p=emc_p, emc_s=emc_s, emc_o=emc_o,
             emc_n_L1=emc_n_L1, emc_p_L1=emc_p_L1, emc_s_L1=emc_s_L1, emc_o_L1=emc_o_L1,
             emc_n_L2=emc_n_L2, emc_p_L2=emc_p_L2, emc_s_L2=emc_s_L2, emc_o_L2=emc_o_L2,
             emc_depth=emc_depth,
             emc_n_mth=emc_n_mth, emc_p_mth=emc_p_mth, emc_s_mth=emc_s_mth, emc_o_mth=emc_o_mth,
             emc_n_mth_L1=emc_n_mth_L1, emc_p_mth_L1=emc_p_mth_L1, emc_s_mth_L1=emc_s_mth_L1, emc_o_mth_L1=emc_o_mth_L1,
             emc_n_mth_L2=emc_n_mth_L2, emc_p_mth_L2=emc_p_mth_L2, emc_s_mth_L2=emc_s_mth_L2, emc_o_mth_L2=emc_o_mth_L2,
             emc_lat=emc_lat, emc_lon=emc_lon)

    np.savez(cgf,
             croco_clm_woad_spco2=croco_clm_woad_spco2, croco_woad_spco2=croco_woad_spco2,
             spco2_socat_mth=spco2_socat_mth, spco2_ose_mth=spco2_ose_mth, lsch_spco2_mth=lsch_spco2_mth,
             croco_clm_woad_dic=croco_clm_woad_dic, croco_woad_dic=croco_woad_dic,
             croco_clm_woad_talk=croco_clm_woad_talk, croco_woad_talk=croco_woad_talk,
             talk_ose_mth=talk_ose_mth, dic_ose_mth=dic_ose_mth,
             woa_lat=woa_lat, woa_lon=woa_lon, lat_crocod=lat_crocod, lon_crocod=lon_crocod, h_crocod=h_crocod)

    np.savez(nif,
             issil_cro=issil_cro, isphs_cro=isphs_cro, isnit_cro=isnit_cro, issal_cro=issal_cro,
             isyr=isyr, ismth=ismth, isday=isday, islat=islat, islon=islon, isdep=isdep,
             issal=issal, issil=issil, isphs=isphs, isnit=isnit)


if data_plotting == 1:
    # Satellite chlorophyll (diat and nano)
    sat_croco_fil = crocodir + 'sat_croco_nd_chl.npz'
    croco_fil = crocodir + 'croco_nd_chl.npz'
    ibi_fil = crocodir + 'ibi_nd_chl.npz'
    # Satellite temperature (for review of impact of temperature on BGC)
    sattemp_croco_fil = crocodir + 'sat_croco_temp.npz'
    crocotemp_fil = crocodir + 'croco_temp.npz'
    #
    ngf = crocodir + 'nutrients_grid.npz'
    nif = crocodir + 'nutrients_insitu.npz'
    cgf = crocodir + 'carbon_grid.npz'

    # Commenting out to get temperature plot
    # sat_ndfile = np.load(sat_croco_fil)
    # sat_dchl = sat_ndfile['sat_dchl']
    # sat_nchl = sat_ndfile['sat_nchl']
    # sat_dchle = sat_ndfile['sat_dchle']
    # sat_nchle = sat_ndfile['sat_nchle']
    # sat_dchlb = sat_ndfile['sat_dchlb']
    # sat_nchlb = sat_ndfile['sat_nchlb']
    # sat_dchl_corrected = 10 ** (np.log10(sat_dchl) - sat_dchlb)
    # sat_nchl_corrected = 10 ** (np.log10(sat_nchl) - sat_nchlb)
    #
    # croco_ndfile = np.load(croco_fil)
    # croco_dchl = croco_ndfile['croco_dchl']
    # croco_nchl = croco_ndfile['croco_nchl']
    #
    # ibi_ndfile = np.load(ibi_fil)
    # ibi_dchl = ibi_ndfile['ibi_dchl']
    # ibi_nchl = ibi_ndfile['ibi_nchl']

    sat_tempfile = np.load(sattemp_croco_fil)
    sat_temp = sat_tempfile['sat_temp']

    croco_tempfile = np.load(crocotemp_fil)
    croco_temp = croco_tempfile['croco_temp']

    ngfile = np.load(ngf)
    woa_n = ngfile['woa_n']
    woa_p = ngfile['woa_p']
    woa_s = ngfile['woa_s']
    woa_o = ngfile['woa_o']
    woa_n_mth = ngfile['woa_n_mth']
    woa_p_mth = ngfile['woa_p_mth']
    woa_s_mth = ngfile['woa_s_mth']
    woa_o_mth = ngfile['woa_o_mth']
    croco_clm_woad_si = ngfile['croco_clm_woad_si']
    croco_clm_woad_o2 = ngfile['croco_clm_woad_o2']
    croco_clm_woad_no3 = ngfile['croco_clm_woad_no3']
    croco_clm_woad_po4 = ngfile['croco_clm_woad_po4']
    croco_woad_si = ngfile['croco_woad_si']
    croco_woad_o2 = ngfile['croco_woad_o2']
    croco_woad_no3 = ngfile['croco_woad_no3']
    croco_woad_po4 = ngfile['croco_woad_po4']
    croco_clm_woad_si[np.isnan(woa_s)] = np.nan
    croco_clm_woad_o2[np.isnan(woa_o)] = np.nan
    croco_clm_woad_no3[np.isnan(woa_n)] = np.nan
    croco_clm_woad_po4[np.isnan(woa_p)] = np.nan
    croco_woad_si[np.isnan(woa_s_mth)] = np.nan
    croco_woad_o2[np.isnan(woa_o_mth)] = np.nan
    croco_woad_no3[np.isnan(woa_n_mth)] = np.nan
    croco_woad_po4[np.isnan(woa_p_mth)] = np.nan
    woa_lat = ngfile['woa_lat']
    woa_lon = ngfile['woa_lon']
    lat_crocod = ngfile['lat_crocod']
    lon_crocod = ngfile['lon_crocod']
    h_crocod = ngfile['h_crocod']
    woa_depth_nps = ngfile['woa_depth_nps']
    woa_depth_ox = ngfile['woa_depth_ox']
    croco_clm_emcd_si = ngfile['croco_clm_emcd_si']
    croco_clm_emcd_o2 = ngfile['croco_clm_emcd_o2']
    croco_clm_emcd_no3 = ngfile['croco_clm_emcd_no3']
    croco_clm_emcd_po4 = ngfile['croco_clm_emcd_po4']
    croco_emcd_si = ngfile['croco_emcd_si']
    croco_emcd_o2 = ngfile['croco_emcd_o2']
    croco_emcd_no3 = ngfile['croco_emcd_no3']
    croco_emcd_po4 = ngfile['croco_emcd_po4']
    emc_n = ngfile['emc_n']
    emc_p = ngfile['emc_p']
    emc_s = ngfile['emc_s']
    emc_o = ngfile['emc_o']
    emc_depth = ngfile['emc_depth']
    emc_n_mth = ngfile['emc_n_mth']
    emc_p_mth = ngfile['emc_p_mth']
    emc_s_mth = ngfile['emc_s_mth']
    emc_o_mth = ngfile['emc_o_mth']
    croco_clm_emcd_si[np.isnan(emc_s)] = np.nan
    croco_clm_emcd_o2[np.isnan(emc_o)] = np.nan
    croco_clm_emcd_no3[np.isnan(emc_n)] = np.nan
    croco_clm_emcd_po4[np.isnan(emc_p)] = np.nan
    croco_emcd_si[np.isnan(emc_s_mth)] = np.nan
    croco_emcd_o2[np.isnan(emc_o_mth)] = np.nan
    croco_emcd_no3[np.isnan(emc_n_mth)] = np.nan
    croco_emcd_po4[np.isnan(emc_p_mth)] = np.nan
    # Loading L1 and L2
    emc_n_L1 = ngfile['emc_n_L1']
    emc_p_L1 = ngfile['emc_p_L1']
    emc_s_L1 = ngfile['emc_s_L1']
    emc_o_L1 = ngfile['emc_o_L1']
    emc_n_mth_L1 = ngfile['emc_n_mth_L1']
    emc_p_mth_L1 = ngfile['emc_p_mth_L1']
    emc_s_mth_L1 = ngfile['emc_s_mth_L1']
    emc_o_mth_L1 = ngfile['emc_o_mth_L1']
    #
    emc_n_L2 = ngfile['emc_n_L2']
    emc_p_L2 = ngfile['emc_p_L2']
    emc_s_L2 = ngfile['emc_s_L2']
    emc_o_L2 = ngfile['emc_o_L2']
    emc_n_mth_L2 = ngfile['emc_n_mth_L2']
    emc_p_mth_L2 = ngfile['emc_p_mth_L2']
    emc_s_mth_L2 = ngfile['emc_s_mth_L2']
    emc_o_mth_L2 = ngfile['emc_o_mth_L2']
    #
    emc_lat = ngfile['emc_lat']
    emc_lon = ngfile['emc_lon']

    cgfile = np.load(cgf)
    croco_clm_woad_spco2 = cgfile['croco_clm_woad_spco2']
    croco_woad_spco2 = cgfile['croco_woad_spco2']
    spco2_socat_mth = cgfile['spco2_socat_mth']
    spco2_ose_mth = cgfile['spco2_ose_mth']
    lsch_spco2_mth = cgfile['lsch_spco2_mth']
    croco_clm_woad_dic = cgfile['croco_clm_woad_dic']
    croco_woad_dic = cgfile['croco_woad_dic']
    croco_clm_woad_talk = cgfile['croco_clm_woad_talk']
    croco_woad_talk = cgfile['croco_woad_talk']
    talk_ose_mth = cgfile['talk_ose_mth']
    dic_ose_mth = cgfile['dic_ose_mth']
    woa_lat = cgfile['woa_lat']
    woa_lon = cgfile['woa_lon']
    lat_crocod = cgfile['lat_crocod']
    lon_crocod = cgfile['lon_crocod']
    h_crocod = cgfile['h_crocod']

    nifile = np.load(nif)
    issil_cro = nifile['issil_cro']
    isphs_cro = nifile['isphs_cro']
    isnit_cro = nifile['isnit_cro']
    issal_cro = nifile['issal_cro']
    issil = nifile['issil']
    isphs = nifile['isphs']
    isnit = nifile['isnit']
    issal = nifile['issal']
    isyr = nifile['isyr']
    ismth = nifile['ismth']
    isday = nifile['isday']
    islat = nifile['islat']
    islon = nifile['islon']
    isdep = nifile['isdep']

    nps_didx_25m = np.argwhere(woa_depth_nps == 25)[0][0]
    nps_didx_50m = np.argwhere(woa_depth_nps == 50)[0][0]
    nps_didx_100m = np.argwhere(woa_depth_nps == 100)[0][0]

    ox_didx_25m = np.argwhere(woa_depth_ox == 25)[0][0]
    ox_didx_50m = np.argwhere(woa_depth_ox == 50)[0][0]
    ox_didx_100m = np.argwhere(woa_depth_ox == 100)[0][0]

    emc_didx_25m = np.argwhere(emc_depth == 25)[0][0]
    emc_didx_50m = np.argwhere(emc_depth == 50)[0][0]
    emc_didx_100m = np.argwhere(emc_depth == 100)[0][0]

    # Scatter plots from in-situ data
    insit_real = np.argwhere((issil_cro != 0.) & (isnit_cro != 0.) &
                             (isphs_cro != 0.) & (issal_cro != 0.))[:, 0]
    issil_cro = issil_cro[insit_real]
    isphs_cro = isphs_cro[insit_real]
    isnit_cro = isnit_cro[insit_real]
    issal_cro = issal_cro[insit_real]
    issil = issil[insit_real]
    isphs = isphs[insit_real]
    isnit = isnit[insit_real]
    issal = issal[insit_real]
    isyr = isyr[insit_real]
    ismth = ismth[insit_real]
    isday = isday[insit_real]
    islat = islat[insit_real]
    islon = islon[insit_real]
    isdep = isdep[insit_real]


    def is_scatplot(inssil_cro, insphs_cro, insnit_cro, inssal_cro,
                    inssil, insphs, insnit, inssal, plt_tit, descr):
        #
        silmin = np.nanmin((np.nanmin(inssil), np.nanmin(inssil_cro)))
        silmax = np.nanmax((np.nanmax(inssil), np.nanmax(inssil_cro)))
        #
        nitmin = np.nanmin((np.nanmin(insnit), np.nanmin(insnit_cro)))
        nitmax = np.nanmax((np.nanmax(insnit), np.nanmax(insnit_cro)))
        #
        phsmin = np.nanmin((np.nanmin(insphs), np.nanmin(insphs_cro)))
        phsmax = np.nanmax((np.nanmax(insphs), np.nanmax(insphs_cro)))
        #
        salmin = np.nanmin((np.nanmin(inssal), np.nanmin(inssal_cro)))
        salmax = np.nanmax((np.nanmax(inssal), np.nanmax(inssal_cro)))

        sillim = (np.floor(silmin/2.5)*2.5, np.ceil(silmax/2.5)*2.5)
        nitlim = (np.floor(nitmin/2.5)*2.5, np.ceil(nitmax/2.5)*2.5)
        phslim = (np.floor(phsmin/0.1)*0.1, np.ceil(phsmax/0.1)*0.1)
        sallim = (np.floor(salmin/5)*5, np.ceil(salmax/0.5)*0.5)

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('In-situ')
        plt.ylabel('CS1KM')
        (ax1, ax2), (ax3, ax4) = axs
        ax1.scatter(inssil, inssil_cro)
        ax1.set_xlim(sillim)
        ax1.set_ylim(sillim)
        ax1.set_xticks(np.arange(sillim[0], sillim[1], 10))
        ax1.set_yticks(np.arange(sillim[0], sillim[1], 10))
        ax2.scatter(insnit, insnit_cro)
        ax2.set_xlim(nitlim)
        ax2.set_ylim(nitlim)
        ax2.set_xticks(np.arange(nitlim[0], nitlim[1], 5))
        ax2.set_yticks(np.arange(nitlim[0], nitlim[1], 5))
        ax3.scatter(insphs, insphs_cro)
        ax3.set_xlim(phslim)
        ax3.set_ylim(phslim)
        ax3.set_xticks(np.arange(phslim[0], phslim[1], 0.25))
        ax3.set_yticks(np.arange(phslim[0], phslim[1], 0.25))
        ax4.scatter(inssal, inssal_cro)
        ax4.set_xlim(sallim)
        ax4.set_ylim(sallim)
        ax4.set_xticks(np.arange(sallim[0], sallim[1] + 2.5, 2.5))
        ax4.set_yticks(np.arange(sallim[0], sallim[1] + 2.5, 2.5))
        fig.suptitle(descr + ' In-situ Vs. CS1KM')
        ax1.set_title('Si (mmol/m3)', fontsize=10)
        ax2.set_title('NO3 (mmol/m3)', fontsize=10)
        ax3.set_title('PO4 (mmol/m3)', fontsize=10)
        ax4.set_title('Sal (psu)', fontsize=10)
        fig.tight_layout()
        # plt.show()
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    def is_spatplot(ins_mod, ins, inslat, inslon, plt_tit, descr):
        #
        latmin = np.nanmin(inslat)
        latmax = np.nanmax(inslat)
        #
        lonmin = np.nanmin(inslon)
        lonmax = np.nanmax(inslon)
        #
        pmin = np.nanmin((np.nanmin(ins), np.nanmin(ins_mod)))
        pmax = np.nanmax((np.nanmax(ins), np.nanmax(ins_mod)))
        #
        cm = plt.cm.get_cmap('RdYlBu_r')
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Long.')
        plt.ylabel('Lat.')
        (ax1, ax2) = axs
        # In-situ data
        ax1.scatter(inslon, inslat, s=10, c=ins, cmap=cm, vmin=pmin, vmax=pmax, alpha=0.5)
        ax1.set_xlim((lonmin, lonmax))
        ax1.set_ylim((latmin, latmax))
        # Model data
        im = ax2.scatter(inslon, inslat, s=10, c=ins_mod, cmap=cm, vmin=pmin, vmax=pmax, alpha=0.5)
        ax2.set_xlim((lonmin, lonmax))
        ax2.set_ylim((latmin, latmax))
        fig.suptitle(descr + 'In-situ Vs. CS1KM')
        ax1.set_title('In-situ', fontsize=10)
        ax2.set_title('CS1KM', fontsize=10)
        fig.tight_layout()
        plt.colorbar(im, ax=axs)
        # plt.show()
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def tseries_clim_plt(woa_arr, croco_arr, plt_tit, idx_0, idx_25, idx_50, idx_100, pstring):
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Month')
        plt.ylabel('mmol/m3')
        axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_arr[:, idx_0, :, :], 2), 1), 'k',
                       label='Clim')
        axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr[:, idx_0, :, :], 2), 1), 'b-',
                       label='CS1KM')
        axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_arr[:, idx_25, :, :], 2), 1), 'k',
                       label='Clim')
        axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr[:, idx_25, :, :], 2), 1), 'b-',
                       label='CS1KM')
        axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_arr[:, idx_50, :, :], 2), 1), 'k',
                       label='Clim')
        axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr[:, idx_50, :, :], 2), 1), 'b-',
                       label='CS1KM')
        axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_arr[:, idx_100, :, :], 2), 1), 'k',
                       label='Clim')
        axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr[:, idx_100, :, :], 2), 1), 'b-',
                       label='CS1KM')
        fig.suptitle('Monthly Vs.' + pstring)
        axs[0, 0].title.set_text('0m')
        axs[0, 1].title.set_text('25m')
        axs[1, 0].title.set_text('50m')
        axs[1, 1].title.set_text('100m')
        handles, labels = axs[1, 1].get_legend_handles_labels()
        fig.legend(handles, labels, loc='center', ncol=2)
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def tseries_dic_talk_plt(dic_ose, talk_ose, croco_dic, croco_talk, plt_tit):
        croco_dic_ose = np.nan*np.ones_like(croco_dic)
        croco_talk_ose = np.nan*np.ones_like(croco_talk)
        croco_dic_ose[:] = croco_dic
        croco_talk_ose[:] = croco_talk
        croco_dic_ose[np.isnan(dic_ose[:, 0, :, :])] = np.nan
        croco_talk_ose[np.isnan(talk_ose[:, 0, :, :])] = np.nan
        fig, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Month')
        plt.ylabel('mmol/m3')
        axs[0].plot(range(1, 25), np.nanmean(np.nanmean(dic_ose[:, 0, :, :], 2), 1), 'k',
                    label='OSE')
        axs[0].plot(range(1, 25), np.nanmean(np.nanmean(croco_dic_ose, 2), 1), 'b-',
                    label='CS1KM')
        axs[1].plot(range(1, 25), np.nanmean(np.nanmean(talk_ose[:, 0, :, :], 2), 1), 'k',
                    label='OSE')
        axs[1].plot(range(1, 25), np.nanmean(np.nanmean(croco_talk_ose, 2), 1), 'b-',
                    label='CS1KM')
        fig.suptitle('Monthly DIC & TALK Vs. OSE gridded products')
        axs[0].title.set_text('DIC')
        axs[1].title.set_text('TALK')
        handles, labels = axs[1].get_legend_handles_labels()
        axs[0].legend(handles, labels, loc='lower center', ncol=1, fontsize=8)
        axs[1].legend(handles, labels, loc='center right', ncol=2, fontsize=8)
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def tseries_spco2_plt(socat_arr, ose_arr, lsch_arr, croco_arr, plt_tit):
        croco_arr_socat = np.nan*np.ones_like(croco_arr)
        croco_arr_ose = np.nan*np.ones_like(croco_arr)
        croco_arr_lsch = np.nan*np.ones_like(croco_arr)
        croco_arr_socat[:] = croco_arr
        croco_arr_ose[:] = croco_arr
        croco_arr_lsch[:] = croco_arr
        croco_arr_socat[np.isnan(socat_arr)] = np.nan
        croco_arr_ose[np.isnan(ose_arr)] = np.nan
        croco_arr_lsch[np.isnan(ose_arr)] = np.nan
        lsch_arr[np.isnan(ose_arr)] = np.nan
        #
        fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Month')
        plt.ylabel('uatm')
        axs[0].plot(range(1, 25), np.nanmean(np.nanmean(socat_arr, 2), 1), 'k',
                    label='SOCAT')
        axs[0].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr_socat, 2), 1), 'b-',
                    label='CS1KM')
        axs[1].plot(range(1, 25), np.nanmean(np.nanmean(ose_arr, 2), 1), 'k',
                    label='OSE')
        axs[1].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr_ose, 2), 1), 'b-',
                    label='CS1KM')
        axs[2].plot(range(1, 25), np.nanmean(np.nanmean(lsch_arr, 2), 1), 'k',
                    label='LSCH')
        axs[2].plot(range(1, 25), np.nanmean(np.nanmean(croco_arr_lsch, 2), 1), 'b-',
                    label='CS1KM')
        fig.suptitle('Monthly spCO2')
        axs[0].title.set_text('SOCAT')
        axs[1].title.set_text('OSE')
        axs[2].title.set_text('LSCH')
        axs[0].legend(loc='lower center', ncol=2, fontsize=8)
        axs[1].legend(loc='lower center', ncol=2, fontsize=8)
        axs[2].legend(loc='lower center', ncol=2, fontsize=8)
        fig.tight_layout()
        # handles, labels = axs[2].get_legend_handles_labels()
        # fig.legend(handles, labels, loc='center', ncol=2)
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def tseries_CHL_plt(sat_diat, sat_diate, sat_diatb, sat_nano, sat_nanoe, sat_nanob,
                        croc_diat, croc_nano, ibi_diat, ibi_nano, plt_tit):
        fig, axs = plt.subplots(ncols=1, nrows=4, sharex=True, sharey=False, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Day')
        plt.ylabel('mg/m3')
        tlev1 = [0.005, 0.01, 0.02, 0.04, 0.08, 0.15, 0.3, 0.6, 1.2, 2.5, 5, 10]
        tlev1_log = np.log(tlev1)
        tlev1_str = ['0.005', '0.01', '0.02', '0.04', '0.08', '0.15', '0.3', '0.6', '1.2', '2.5', '5', '10']
        mlev = np.arange(15, 730, 30.5)
        mlev_str = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',
                    'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

        # sat_diat_ts = np.nanmean(np.nanmean(sat_diat, 2), 1)
        # sat_diat_ts = np.nanmean(np.nanmean(sat_diat, 2), 1) - np.nanmean(np.nanmean(sat_diatb, 2), 1)
        sat_diat_ts = 10. ** (np.log10(np.nanmean(np.nanmean(sat_diat, 2), 1)) -
                              np.nanmean(np.nanmean(sat_diatb, 2), 1))
        croc_diat_ts = np.nanmean(np.nanmean(croc_diat, 2), 1)
        ibi_diat_ts = np.nanmean(np.nanmean(ibi_diat, 2), 1)
        croc_diat_ts[sat_diat_ts < 0] = np.nan
        ibi_diat_ts[sat_diat_ts < 0] = np.nan
        sat_diat_ts[sat_diat_ts < 0] = np.nan

        # Diatom timeseries
        axs[0].plot(range(1, sat_diat.shape[0] + 1),  np.log(sat_diat_ts), 'k', label='Sat.')
        axs[0].plot(range(1, sat_diat.shape[0] + 1), np.log(croc_diat_ts), 'b-', label='CS1KM')
        axs[0].plot(range(1, sat_diat.shape[0] + 1), np.log(ibi_diat_ts), 'g', label='IBI_MY')
        axs[0].set_yticks(tlev1_log)
        axs[0].set_yticklabels(tlev1_str)
        axs[0].set_xticks(mlev)
        axs[0].set_xticklabels(mlev_str)

        # sat_nano_ts = np.nanmean(np.nanmean(sat_nano, 2), 1)
        # sat_nano_ts = np.nanmean(np.nanmean(sat_nano, 2), 1) - np.nanmean(np.nanmean(sat_nanob, 2), 1)
        sat_nano_ts = 10. ** (np.log10(np.nanmean(np.nanmean(sat_nano, 2), 1)) -
                              np.nanmean(np.nanmean(sat_nanob, 2), 1))
        croc_nano_ts = np.nanmean(np.nanmean(croc_nano, 2), 1)
        ibi_nano_ts = np.nanmean(np.nanmean(ibi_nano, 2), 1)
        croc_nano_ts[sat_nano_ts < 0] = np.nan
        ibi_nano_ts[sat_nano_ts < 0] = np.nan
        sat_nano_ts[sat_nano_ts < 0] = np.nan

        # Nanophytoplankton timeseries
        axs[1].plot(range(1, sat_nano.shape[0] + 1), sat_nano_ts, 'k', label='Sat.')
        axs[1].plot(range(1, sat_nano.shape[0] + 1), croc_nano_ts, 'b-', label='CS1KM')
        axs[1].plot(range(1, sat_nano.shape[0] + 1), ibi_nano_ts, 'g', label='IBI_MY')
        axs[1].set_xticks(mlev)
        axs[1].set_xticklabels(mlev_str)

        # Chl (Nano + Diat) timeseries
        axs[2].plot(range(1, sat_nano.shape[0] + 1), np.log(sat_diat_ts + sat_nano_ts), 'k', label='Sat.')
        axs[2].plot(range(1, sat_nano.shape[0] + 1), np.log(croc_diat_ts + croc_nano_ts), 'b-', label='CS1KM')
        axs[2].plot(range(1, sat_nano.shape[0] + 1), np.log(ibi_diat_ts + ibi_nano_ts), 'g', label='IBI_MY')
        axs[2].set_yticks(tlev1_log)
        axs[2].set_yticklabels(tlev1_str)
        axs[2].set_xticks(mlev)
        axs[2].set_xticklabels(mlev_str)

        # Number of data points
        dcount = np.count_nonzero(np.count_nonzero((croc_diat > 0), 2), 1)
        # dcount[dcount == 0] = np.int(np.nan)
        countmin = np.nanmin(dcount)
        countmax = np.nanmax(dcount)
        # axs[3].plot(range(1, sat_diat.shape[0] + 1), dcount, 'k')
        axs[3].scatter(range(1, sat_diat.shape[0] + 1), dcount)
        fig.suptitle('Satellite Vs. CS1KM')
        axs[0].title.set_text('Diatom')
        axs[1].title.set_text('Nanophyto')
        axs[2].title.set_text('Chl (Diat. + Nano.)')
        axs[3].title.set_text('Data count')
        axs[0].legend(loc='best', ncol=2, fontsize=8)
        axs[1].legend(loc='best', ncol=2, fontsize=8)
        axs[2].legend(loc='best', ncol=2, fontsize=8)
        fig.tight_layout()
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def tseries_temp_plt(sattemp, croctemp, plt_tit):
        croctemp[croctemp == 0.] = np.nan
        sattemp[sattemp == 0.] = np.nan
        rmidx = [1154, 2615, 4076, 5538, 6998, 8459, 9920]
        tdays = np.arange(0, sattemp.shape[0])
        yrdays = np.setdiff1d(tdays, rmidx)
        croct_av = np.nanmean(np.nanmean(croctemp, 2), 1)
        croct_clim = np.nanmean(np.reshape(croct_av[yrdays], (365, 29), order='F'), 1)
        satt_av = np.nanmean(np.nanmean(sattemp, 2), 1)
        satt_clim = np.nanmean(np.reshape(satt_av[yrdays], (365, 29), order='F'), 1)
        fig, axs = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=False, figsize=(8, 8))
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
        plt.xlabel('Day')
        plt.ylabel('$^\circ$ C')
        tlev1 = [9, 10, 11, 12, 13, 14, 15, 16, 17]
        tlev1_str = ['9', '10', '11', '12', '13', '14', '15', '16', '17']
        # mlev = np.arange(15, 730, 30.5)
        # mlev_str = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D',
        #             'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
        mlev = np.arange(15, 365, 30.5)
        mlev_str = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

        # Temperature timeseries
        axs[0].plot(range(1, satt_clim.shape[0] + 1), satt_clim, 'k',
                    label='Sat.')
        axs[0].plot(range(1, croct_clim.shape[0] + 1), croct_clim, 'b-',
                    label='CS1KM')
        axs[0].set_yticks(tlev1)
        axs[0].set_yticklabels(tlev1_str)
        axs[0].set_xticks(mlev)
        axs[0].set_xticklabels(mlev_str)

        # Temperature bias
        dcount = np.count_nonzero(np.count_nonzero((croctemp > 0), 2), 1)
        # countmin = np.nanmin(dcount)
        # countmax = np.nanmax(dcount)
        axs[1].scatter(range(1, croct_clim.shape[0] + 1), croct_clim - satt_clim)
        fig.suptitle('SST: Satellite Vs. CS1KM')
        axs[1].title.set_text('Bias: CS1KM - SST')
        axs[0].legend(loc='best', ncol=2, fontsize=8)
        fig.tight_layout()
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def scat_satCHL_plt(sat_diat, sat_nano, croc_diat, croc_nano, plt_tit):
        tlev = [0.0001, 0.0002, 0.0004, 0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15,
                0.3, 0.6, 1.2, 2.5, 5, 10, 20, 40, 80, 120]
        tlev_log = np.log(tlev)
        tlev_str = ['0.0001', '0.0002', '0.0004', '0.001', '0.002', '0.005', '0.01', '0.02', '0.04', '0.08', '0.15',
                    '0.3', '0.6', '1.2', '2.5', '5', '10', '20', '40', '80', '120']

        # Diatom vars chlorophyll
        tot_sat_diat = sat_diat[sat_diat > 0].flatten()
        tot_croc_diat = croc_diat[croc_diat > 0].flatten()
        spat_sat_diat = np.nanmean(sat_diat, 0).flatten()
        spat_croc_diat = np.nanmean(croc_diat, 0).flatten()
        temp_sat_diat = np.nanmean(np.nanmean(sat_diat, 2), 1).flatten()
        temp_croc_diat = np.nanmean(np.nanmean(croc_diat, 2), 1).flatten()

        # Nanophytoplankton chlorophyll
        tot_sat_nano = sat_nano[sat_nano > 0].flatten()
        tot_croc_nano = croc_nano[croc_nano > 0].flatten()
        spat_sat_nano = np.nanmean(sat_nano, 0).flatten()
        spat_croc_nano = np.nanmean(croc_nano, 0).flatten()
        temp_sat_nano = np.nanmean(np.nanmean(sat_nano, 2), 1).flatten()
        temp_croc_nano = np.nanmean(np.nanmean(croc_nano, 2), 1).flatten()

        # Total chlorophyll
        tot_sat_tchl = tot_sat_diat + tot_sat_nano
        tot_croc_tchl = tot_croc_diat + tot_croc_nano
        spat_sat_tchl = spat_sat_diat + spat_sat_nano
        spat_croc_tchl = spat_croc_diat + spat_croc_nano
        temp_sat_tchl = temp_sat_diat + temp_sat_nano
        temp_croc_tchl = temp_croc_diat + temp_croc_nano
        totmin = np.nanmin((tot_sat_diat, tot_croc_diat,
                            tot_sat_nano, tot_croc_nano,
                            tot_sat_tchl, tot_croc_tchl))
        totmax = np.nanmax((tot_sat_diat, tot_croc_diat,
                            tot_sat_nano, tot_croc_nano,
                            tot_sat_tchl, tot_croc_tchl))
        t1_min_idx = np.argwhere(tlev <= totmin)[:, 0][-1]
        t1_max_idx = np.argwhere(tlev > totmax)[:, 0][0]
        tlev1 = tlev[0:t1_max_idx+1]
        tlev1 = tlev1[t1_min_idx:]
        tlev1_log = tlev_log[0:t1_max_idx+1]
        tlev1_log = tlev1_log[t1_min_idx:]
        tlev1_str = tlev_str[0:t1_max_idx+1]
        tlev1_str = tlev1_str[t1_min_idx:]

        # Spatial Chlorophyll
        spatmin = np.nanmin((spat_sat_diat, spat_croc_diat,
                             spat_sat_nano, spat_croc_nano,
                             spat_sat_tchl, spat_croc_tchl))
        spatmax = np.nanmax((spat_sat_diat, spat_croc_diat,
                             spat_sat_nano, spat_croc_nano,
                             spat_sat_tchl, spat_croc_tchl))
        try:
            t2_min_idx = np.argwhere(tlev <= spatmin)[:, 0][-1]
        except:
            t2_min_idx = 0
        t2_max_idx = np.argwhere(tlev > spatmax)[:, 0][0]
        tlev2 = tlev[0:t2_max_idx+1]
        tlev2 = tlev2[t2_min_idx:]
        tlev2_log = tlev_log[0:t2_max_idx+1]
        tlev2_log = tlev2_log[t2_min_idx:]
        tlev2_str = tlev_str[0:t2_max_idx+1]
        tlev2_str = tlev2_str[t2_min_idx:]

        # Temporal Chlorophyll
        tempmin = np.nanmin((temp_sat_diat, temp_croc_diat,
                             temp_sat_nano, temp_croc_nano,
                             temp_sat_tchl, temp_croc_tchl))
        tempmax = np.nanmax((temp_sat_diat, temp_croc_diat,
                             temp_sat_nano, temp_croc_nano,
                             temp_sat_tchl, temp_croc_tchl))
        try:
            t3_min_idx = np.argwhere(tlev <= tempmin)[:, 0][-1]
        except:
            t3_min_idx = 0
        t3_max_idx = np.argwhere(tlev > tempmax)[:, 0][0]
        tlev3 = tlev[0:t3_max_idx+1]
        tlev3 = tlev3[t3_min_idx:]
        tlev3_log = tlev_log[0:t3_max_idx+1]
        tlev3_log = tlev3_log[t3_min_idx:]
        tlev3_str = tlev_str[0:t3_max_idx+1]
        tlev3_str = tlev3_str[t3_min_idx:]

        # plot for total diat, spat diat, temp diat
        # same for nano and total chlorophyll
        fig1, axs1 = plt.subplots(ncols=3, nrows=1, sharex=True, sharey=True, figsize=(8, 8))
        (ax1, ax2, ax3) = axs1
        ax1.scatter(np.log(tot_sat_diat), np.log(tot_croc_diat))
        ax2.scatter(np.log(tot_sat_nano), np.log(tot_croc_nano))
        ax3.scatter(np.log(tot_sat_tchl), np.log(tot_croc_tchl))
        ax1.set_yticks(tlev1_log)
        ax1.set_yticklabels(tlev1_str)
        ax1.set_xticks(tlev1_log)
        ax1.set_xticklabels(tlev1_str)
        ax2.set_yticks(tlev1_log)
        ax2.set_yticklabels(tlev1_str)
        ax2.set_xticks(tlev1_log)
        ax2.set_xticklabels(tlev1_str)
        ax3.set_yticks(tlev1_log)
        ax3.set_yticklabels(tlev1_str)
        ax3.set_xticks(tlev1_log)
        ax3.set_xticklabels(tlev1_str)
        fig1.suptitle('Daily chlorophyll (all space & time)')
        ax1.title.set_text('Diat. Chl')
        ax2.title.set_text('Nano. Chl')
        ax3.title.set_text('Tot. Chl')
        fig1.tight_layout()
        # plt.show()
        plt.savefig(plt_tit[:-4] + '_total.jpg', dpi='figure', format='jpg')
        fig2, axs2 = plt.subplots(ncols=3, nrows=1, sharex=True, sharey=True, figsize=(8, 8))
        (ax4, ax5, ax6) = axs2
        ax4.scatter(np.log(spat_sat_diat), np.log(spat_croc_diat))
        ax5.scatter(np.log(spat_sat_nano), np.log(spat_croc_nano))
        ax6.scatter(np.log(spat_sat_tchl), np.log(spat_croc_tchl))
        ax4.set_yticks(tlev2_log)
        ax4.set_yticklabels(tlev2_str)
        ax4.set_xticks(tlev2_log)
        ax4.set_xticklabels(tlev2_str)
        ax5.set_yticks(tlev2_log)
        ax5.set_yticklabels(tlev2_str)
        ax5.set_xticks(tlev2_log)
        ax5.set_xticklabels(tlev2_str)
        ax6.set_yticks(tlev2_log)
        ax6.set_yticklabels(tlev2_str)
        ax6.set_xticks(tlev2_log)
        ax6.set_xticklabels(tlev2_str)
        fig2.suptitle('2 yr averaged Chl')
        ax4.title.set_text('Diat. Chl')
        ax5.title.set_text('Nano. Chl')
        ax6.title.set_text('Tot. Chl')
        fig2.tight_layout()
        # plt.show()
        plt.savefig(plt_tit[:-4] + '_2yr_avg.jpg', dpi='figure', format='jpg')
        fig3, axs3 = plt.subplots(ncols=3, nrows=1, sharex=True, sharey=True, figsize=(8, 8))
        (ax7, ax8, ax9) = axs3
        ax7.scatter(np.log(temp_sat_diat), np.log(temp_croc_diat))
        ax8.scatter(np.log(temp_sat_nano), np.log(temp_croc_nano))
        ax9.scatter(np.log(temp_sat_tchl), np.log(temp_croc_tchl))
        ax7.set_yticks(tlev3_log)
        ax7.set_yticklabels(tlev3_str)
        ax7.set_xticks(tlev3_log)
        ax7.set_xticklabels(tlev3_str)
        ax8.set_yticks(tlev3_log)
        ax8.set_yticklabels(tlev3_str)
        ax8.set_xticks(tlev3_log)
        ax8.set_xticklabels(tlev3_str)
        ax9.set_yticks(tlev3_log)
        ax9.set_yticklabels(tlev3_str)
        ax9.set_xticks(tlev3_log)
        ax9.set_xticklabels(tlev3_str)
        fig3.suptitle('Spatially-averaged Chl')
        ax7.title.set_text('Diat. Chl')
        ax8.title.set_text('Nano. Chl')
        ax9.title.set_text('Tot. Chl')
        fig3.tight_layout()
        # plt.show()
        plt.savefig(plt_tit[:-4] + '_spat_avg.jpg', dpi='figure', format='jpg')
        pass


    def spat_CHL_plt(croco_arr, sat_arr, sat_arre, sat_arrb, ibi_arr, clev, clev_str, lon_croco, lat_croco, h_croco,
                     plt_tit, pstring):
        clev_log = np.log(clev)
        croco_arr[croco_arr == 0.] = np.nan
        sat_arr[sat_arr == 0.] = np.nan
        sat_arre[sat_arre == 0.] = np.nan
        sat_arrb[sat_arrb == 0.] = np.nan
        ibi_arr[ibi_arr == 0.] = np.nan
        ibi_arr[ibi_arr >1000.] = np.nan
        cm = plt.cm.get_cmap('RdYlBu_r')
        dat_count = np.count_nonzero((croco_arr > 0), 0)
        countmin = np.min(dat_count[dat_count > 0])
        countmax = np.max(dat_count[dat_count > 0])
        sat_arr = np.nanmean(sat_arr, 0)
        sat_arre = np.nanmean(sat_arre, 0)
        sat_arrb = np.nanmean(sat_arrb, 0)
        # sat_arrII = sat_arr - sat_arrb
        # sat_arrII[sat_arrII < 0] = sat_arr[sat_arrII < 0]
        sat_arrII = 10.**(np.log10(sat_arr) - sat_arrb)
        croco_arr = np.nanmean(croco_arr, 0)
        ibi_arr = np.nanmean(ibi_arr, 0)
        nonnanidx = np.argwhere(croco_arr[:, 0] > 0)[:, 0]
        fullidx = np.arange(0, croco_arr.shape[0])
        nanidx = np.setdiff1d(fullidx, nonnanidx)
        for r in range(0, len(nanidx)):
            sat_arr[nanidx[r], :] = sat_arr[nanidx[r] - 1, :]/2 + sat_arr[nanidx[r] + 1, :]/2
            sat_arre[nanidx[r], :] = sat_arre[nanidx[r] - 1, :]/2 + sat_arre[nanidx[r] + 1, :]/2
            sat_arrb[nanidx[r], :] = sat_arrb[nanidx[r] - 1, :]/2 + sat_arrb[nanidx[r] + 1, :]/2
            sat_arrII[nanidx[r], :] = sat_arrII[nanidx[r] - 1, :]/2 + sat_arrII[nanidx[r] + 1, :]/2
            croco_arr[nanidx[r], :] = croco_arr[nanidx[r] - 1, :]/2 + croco_arr[nanidx[r] + 1, :]/2
            ibi_arr[nanidx[r], :] = ibi_arr[nanidx[r] - 1, :]/2 + ibi_arr[nanidx[r] + 1, :]/2
            dat_count[nanidx[r], :] = dat_count[nanidx[r] - 1, :]/2 + dat_count[nanidx[r] + 1, :]/2
        emin = np.nanmin(sat_arre)
        emax = np.nanmax(sat_arre)
        bmin = np.nanmin(sat_arrb)
        bmax = np.nanmax(sat_arrb)
        nmin = np.nanmin((np.nanpercentile(sat_arr, 0.01), np.nanpercentile(croco_arr, 0.01)))
        nmax = np.nanmax((np.nanpercentile(sat_arr, 99.9), np.nanpercentile(croco_arr, 99.9)))
        c_min_idx = np.argwhere(clev <= nmin)[:, 0][-1]
        c_max_idx = np.argwhere(clev > nmax)[:, 0][0]
        clev = clev[0:c_max_idx+1]
        clev = clev[c_min_idx:]
        clev_log = clev_log[0:c_max_idx+1]
        clev_log = clev_log[c_min_idx:]
        clev_str = clev_str[0:c_max_idx+1]
        clev_str = clev_str[c_min_idx:]

        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8, 8))
        (ax1, ax2, ax3) = axs
        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        ax1.contour(lon_croco, lat_croco, np.log(sat_arrII),
                    levels=clev_log, cmap=cm, vmin=np.log(nmin), vmax=np.log(nmax))
        # ax1.contour(lon_croco, lat_croco, np.log(sat_arr),
        #             levels=clev_log, cmap=cm, vmin=np.log(nmin), vmax=np.log(nmax))
        cxlim = (np.nanmin(lon_croco), np.nanmax(lon_croco))
        cylim = (np.nanmin(lat_croco), np.nanmax(lat_croco))
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)
        ax2.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im = ax2.contour(lon_croco, lat_croco, np.log(croco_arr),
                         levels=clev_log, cmap=cm, vmin=np.log(nmin), vmax=np.log(nmax))
        ax2.set_xlim(cxlim)
        ax2.set_ylim(cylim)
        ax2.axes.get_yaxis().set_visible(False)
        ax3.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im = ax3.contour(lon_croco, lat_croco, np.log(ibi_arr),
                         levels=clev_log, cmap=cm, vmin=np.log(nmin), vmax=np.log(nmax))
        ax3.set_xlim(cxlim)
        ax3.set_ylim(cylim)
        ax3.axes.get_yaxis().set_visible(False)
        fig.suptitle(pstring)
        ax1.title.set_text('Sat.')
        ax2.title.set_text('CS1KM')
        ax3.title.set_text('IBI_MY')
        cbar = plt.colorbar(im, ax=axs, ticks=clev_log)
        cbar.ax.set_yticklabels(clev_str)
        plt.savefig(plt_tit, dpi='figure', format='jpg')

        # Plot satellite error and bias
        figb, axsb = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
        (ax5, ax6) = axsb
        ax5.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im5 = ax5.contour(lon_croco, lat_croco, sat_arre,
                          levels=np.arange(np.ceil(emin*10.)/10., np.ceil(emax*10.)/10., 0.05),
                          cmap=cm, vmin=emin, vmax=emax)
        cxlim = (np.nanmin(lon_croco), np.nanmax(lon_croco))
        cylim = (np.nanmin(lat_croco), np.nanmax(lat_croco))
        ax5.set_xlim(cxlim)
        ax5.set_ylim(cylim)
        ax6.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im6 = ax6.contour(lon_croco, lat_croco, sat_arrb,
                          levels=np.arange(np.floor(bmin*10.)/10., np.ceil(bmax*10.)/10., 0.05),
                          cmap=cm, vmin=bmin, vmax=bmax)
        ax6.set_xlim(cxlim)
        ax6.set_ylim(cylim)
        ax6.axes.get_yaxis().set_visible(False)
        fig.suptitle(pstring)
        ax5.title.set_text('Sat. RMSE')
        ax6.title.set_text('Sat. Bias')
        cbar5 = plt.colorbar(im5, ax=ax5)
        cbar6 = plt.colorbar(im6, ax=ax6)
        # cbar.ax.set_yticklabels(clev_str)
        plt.savefig(plt_tit[:-4] + '_error_bias.jpg', dpi='figure', format='jpg')

        fig1, axs1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im2 = axs1.contour(lon_croco, lat_croco, dat_count, cmap=cm, levels=np.arange(100, 240, 20),
                           vmin=100, vmax=240)
        cbar2 = plt.colorbar(im2, ax=axs1, ticks=np.arange(100, 240, 20))
        plt.savefig(plt_tit[:-4] + '_count.jpg', dpi='figure', format='jpg')

        pass

    def spat_temp_plt(croco_arr, sat_arr, clev, clev_str, lon_croco, lat_croco, h_croco,
                      plt_tit, pstring):
        croco_arr[croco_arr == 0.] = np.nan
        sat_arr[sat_arr == 0.] = np.nan
        # cm = plt.cm.get_cmap('RdYlBu_r')
        cm = cmo.thermal
        dat_count = np.count_nonzero((croco_arr > 0), 0)
        countmin = np.min(dat_count[dat_count > 0])
        countmax = np.max(dat_count[dat_count > 0])
        sat_arr = np.nanmean(sat_arr, 0)
        croco_arr = np.nanmean(croco_arr, 0)
        nonnanidx = np.argwhere(croco_arr[:, 0] > 0)[:, 0]
        fullidx = np.arange(0, croco_arr.shape[0])
        nanidx = np.setdiff1d(fullidx, nonnanidx)
        for r in range(0, len(nanidx)):
            sat_arr[nanidx[r], :] = sat_arr[nanidx[r] - 1, :]/2 + sat_arr[nanidx[r] + 1, :]/2
            croco_arr[nanidx[r], :] = croco_arr[nanidx[r] - 1, :]/2 + croco_arr[nanidx[r] + 1, :]/2
            dat_count[nanidx[r], :] = dat_count[nanidx[r] - 1, :]/2 + dat_count[nanidx[r] + 1, :]/2
        nmin = np.nanmin((np.nanpercentile(sat_arr, 5), np.nanpercentile(croco_arr, 5)))
        nmax = np.nanmax((np.nanpercentile(sat_arr, 95), np.nanpercentile(croco_arr, 95)))
        c_min_idx = np.argwhere(clev <= nmin)[:, 0][-1]
        c_max_idx = np.argwhere(clev > nmax)[:, 0][0]
        clev = clev[0:c_max_idx+1]
        clev = clev[c_min_idx:]
        clev_str = clev_str[0:c_max_idx+1]
        clev_str = clev_str[c_min_idx:]
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
        (ax1, ax2) = axs
        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        ax1.contour(lon_croco, lat_croco, sat_arr,
                    levels=clev, cmap=cm, vmin=nmin, vmax=nmax)
        # ax1.pcolormesh(lon_croco, lat_croco, sat_arr, shading='gouraud', vmin=nmin, vmax=nmax)
        cxlim = (np.nanmin(lon_croco), np.nanmax(lon_croco))
        cylim = (np.nanmin(lat_croco), np.nanmax(lat_croco))
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)
        ax2.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im = ax2.contour(lon_croco, lat_croco, croco_arr,
                         levels=clev, cmap=cm, vmin=nmin, vmax=nmax)
        # im = ax2.pcolormesh(lon_croco, lat_croco, croco_arr, shading='gouraud', vmin=nmin, vmax=nmax)
        ax2.set_xlim(cxlim)
        ax2.set_ylim(cylim)
        ax2.axes.get_yaxis().set_visible(False)
        fig.suptitle(pstring)
        ax1.title.set_text('Sat.')
        ax2.title.set_text('CS1KM')
        cbar = plt.colorbar(im, ax=axs, ticks=clev)
        cbar.ax.set_yticklabels(clev_str)
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        # Separate spatial plot of the number of temperature data plots
        fig1, axs1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im2 = axs1.contour(lon_croco, lat_croco, dat_count, cmap=cmo.solar, levels=np.arange(700, 1400, 100),
                           vmin=700, vmax=1400)
        cbar2 = plt.colorbar(im2, ax=axs1, ticks=np.arange(700, 1400, 100))
        plt.savefig(plt_tit[:-4] + '_count.jpg', dpi='figure', format='jpg')
        pass


    def spat_clim_plt(woa_arr, woa_ln, woa_lt,
                      croco_arr, lon_croco, lat_croco, h_croco,
                      plt_tit, pstring):
        cm = plt.cm.get_cmap('RdYlBu_r')
        nmin = np.nanmin((np.nanmin(woa_arr), np.nanmin(croco_arr)))
        nmax = np.nanmin((np.nanmax(woa_arr), np.nanmax(croco_arr)))
        normalizer = Normalize(nmin, nmax)
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
        (ax1, ax2) = axs
        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        ax1.contour(woa_ln, woa_lt, woa_arr, cmap=cm, vmin=nmin, vmax=nmax)
        cxlim = (np.nanmin(woa_ln), np.nanmax(woa_ln))
        cylim = (np.nanmin(woa_lt), np.nanmax(woa_lt))
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)
        ax2.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im = ax2.contour(woa_ln, woa_lt, croco_arr, cmap=cm, vmin=nmin, vmax=nmax)
        ax2.set_xlim(cxlim)
        ax2.set_ylim(cylim)
        ax2.axes.get_yaxis().set_visible(False)
        fig.suptitle(pstring)
        ax1.title.set_text('Clim')
        ax2.title.set_text('CS1KM')
        fig.tight_layout()
        plt.colorbar(plt.cm.ScalarMappable(norm=normalizer, cmap=cm), ax=axs)
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def spat_spco2_plt(socat_arr, ose_arr, lsch_arr, woa_ln, woa_lt,
                       croco_arr, lon_croco, lat_croco, h_croco,
                       plt_tit):
        croco_arr_socat = np.nan*np.ones_like(croco_arr)
        croco_arr_ose = np.nan*np.ones_like(croco_arr)
        croco_arr_lsch = np.nan*np.ones_like(croco_arr)
        croco_arr_socat[:] = croco_arr
        croco_arr_ose[:] = croco_arr
        croco_arr_lsch[:] = croco_arr
        croco_arr_socat[np.isnan(socat_arr)] = np.nan
        croco_arr_ose[np.isnan(ose_arr)] = np.nan
        croco_arr_lsch[np.isnan(ose_arr)] = np.nan
        lsch_arr[np.isnan(ose_arr)] = np.nan
        croco_arr_socat = np.nanmean(croco_arr_socat, 0)
        croco_arr_ose = np.nanmean(croco_arr_ose, 0)
        croco_arr_lsch = np.nanmean(croco_arr_lsch, 0)
        socat_arr = np.nanmean(socat_arr, 0)
        ose_arr = np.nanmean(ose_arr, 0)
        lsch_arr = np.nanmean(lsch_arr, 0)

        cm = plt.cm.get_cmap('RdYlBu_r')
        nmin = np.nanmin((np.nanmin(socat_arr), np.nanmin(ose_arr), np.nanmin(lsch_arr),
                          np.nanmin(croco_arr_socat), np.nanmin(croco_arr_ose), np.nanmin(croco_arr_lsch)))
        nmin = np.floor(nmin/5)*5
        nmax = np.nanmax((np.nanmax(socat_arr), np.nanmax(ose_arr), np.nanmax(lsch_arr),
                          np.nanmax(croco_arr_socat), np.nanmax(croco_arr_ose), np.nanmax(croco_arr_lsch)))
        nmax = np.ceil(nmax/5)*5
        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 8))
        (ax1, ax2), (ax3, ax4), (ax5, ax6) = axs

        normalizer = Normalize(nmin, nmax)

        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        ax1.contour(woa_ln, woa_lt, socat_arr, cmap=cm, vmin=nmin, vmax=nmax, norm=normalizer)
        cxlim = (np.nanmin(woa_ln), np.nanmax(woa_ln))
        cylim = (np.nanmin(woa_lt), np.nanmax(woa_lt))
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)
        #
        ax2.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im2 = ax2.contour(woa_ln, woa_lt, croco_arr_socat, cmap=cm, vmin=nmin, vmax=nmax, norm=normalizer)
        ax2.set_xlim(cxlim)
        ax2.set_ylim(cylim)
        ax2.axes.get_yaxis().set_visible(False)
        #
        ax3.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im3 = ax3.contour(woa_ln, woa_lt, ose_arr, cmap=cm, vmin=nmin, vmax=nmax, norm=normalizer)
        ax3.set_xlim(cxlim)
        ax3.set_ylim(cylim)
        #
        ax4.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im4 = ax4.contour(woa_ln, woa_lt, croco_arr_ose, cmap=cm, vmin=nmin, vmax=nmax, norm=normalizer)
        ax4.set_xlim(cxlim)
        ax4.set_ylim(cylim)
        ax4.axes.get_yaxis().set_visible(False)
        #
        ax5.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im5 = ax5.contour(woa_ln, woa_lt, lsch_arr, cmap=cm, vmin=nmin, vmax=nmax, norm=normalizer)
        ax5.set_xlim(cxlim)
        ax5.set_ylim(cylim)
        #
        ax6.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im6 = ax6.contour(woa_ln, woa_lt, croco_arr_lsch, cmap=cm, vmin=nmin, vmax=nmax, norm=normalizer)
        ax6.set_xlim(cxlim)
        ax6.set_ylim(cylim)
        ax6.axes.get_yaxis().set_visible(False)
        #
        fig.suptitle('spCO2, (uatm)')
        ax1.title.set_text('SOCAT')
        ax2.title.set_text('CS_SOCAT')
        ax3.title.set_text('OSE')
        ax4.title.set_text('CS_OSE')
        ax5.title.set_text('LSCH')
        ax6.title.set_text('CS_LSCH')
        fig.tight_layout()
        plt.colorbar(plt.cm.ScalarMappable(norm=normalizer, cmap=cm), ax=axs)
        # plt.show()
        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass


    def spat_dic_talk_plt(dic_ose, talk_ose, woa_ln, woa_lt,
                          croco_dic, croco_talk, lon_croco, lat_croco, h_croco,
                          plt_tit):

        croco_dic_ose = np.nan*np.ones((croco_dic.shape[0], croco_dic.shape[2], croco_dic.shape[3]))
        croco_dic_ose[:] = croco_dic[:, 0, :, :]
        croco_talk_ose = np.nan*np.ones((croco_talk.shape[0], croco_talk.shape[2], croco_talk.shape[3]))
        croco_talk_ose[:] = croco_talk[:, 0, :, :]

        croco_dic_ose[np.isnan(dic_ose)] = np.nan
        croco_talk_ose[np.isnan(talk_ose)] = np.nan

        croco_dic_ose = np.nanmean(croco_dic_ose, 0)
        croco_talk_ose = np.nanmean(croco_talk_ose, 0)
        dic_ose = np.nanmean(dic_ose, 0)
        talk_ose = np.nanmean(talk_ose, 0)

        cm = plt.cm.get_cmap('RdYlBu_r')
        nmind = np.nanmin((np.nanmin(croco_dic_ose), np.nanmin(dic_ose)))
        nmind = np.floor(nmind/10)*10
        nmaxd = np.nanmax((np.nanmax(croco_dic_ose), np.nanmax(dic_ose)))
        nmaxd = np.ceil(nmaxd/10)*10
        nmint = np.nanmin((np.nanmin(croco_talk_ose), np.nanmin(talk_ose)))
        nmint = np.floor(nmint/10)*10
        nmaxt = np.nanmax((np.nanmax(croco_talk_ose), np.nanmax(talk_ose)))
        nmaxt = np.ceil(nmaxt/10)*10

        normalizerd = Normalize(nmind, nmaxd)
        normalizert = Normalize(nmint, nmaxt)

        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
        (ax1, ax2), (ax3, ax4) = axs
        ax1.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        ax1.contour(woa_ln, woa_lt, dic_ose, cmap=cm, vmin=nmind, vmax=nmaxd, norm=normalizerd)
        cxlim = (np.nanmin(woa_ln), np.nanmax(woa_ln))
        cylim = (np.nanmin(woa_lt), np.nanmax(woa_lt))
        ax1.set_xlim(cxlim)
        ax1.set_ylim(cylim)
        #
        ax2.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im2 = ax2.contour(woa_ln, woa_lt, croco_dic_ose, cmap=cm, vmin=nmind, vmax=nmaxd, norm=normalizerd)
        ax2.set_xlim(cxlim)
        ax2.set_ylim(cylim)
        ax2.axes.get_yaxis().set_visible(False)
        #
        ax3.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im3 = ax3.contour(woa_ln, woa_lt, talk_ose, cmap=cm, vmin=nmint, vmax=nmaxt, norm=normalizert)
        ax3.set_xlim(cxlim)
        ax3.set_ylim(cylim)
        #
        ax4.contourf(lon_croco, lat_croco, h_croco, [0, 10])
        im4 = ax4.contour(woa_ln, woa_lt, croco_talk_ose, cmap=cm, vmin=nmint, vmax=nmaxt, norm=normalizert)
        ax4.set_xlim(cxlim)
        ax4.set_ylim(cylim)
        ax4.axes.get_yaxis().set_visible(False)
        #
        fig.suptitle('DIC (a, b), TALK (c, d)')
        ax1.title.set_text('(a) OSE')
        ax2.title.set_text('(b) CS1KM')
        ax3.title.set_text('(c) OSE')
        ax4.title.set_text('(d) CS1KM')

        fig.tight_layout()
        plt.colorbar(plt.cm.ScalarMappable(norm=normalizerd, cmap=cm), ax=ax2)
        plt.colorbar(plt.cm.ScalarMappable(norm=normalizert, cmap=cm), ax=ax4)

        plt.savefig(plt_tit, dpi='figure', format='jpg')
        pass

    # Temperature spatial plot
    fig_spat_temp = crocodir + 'temp_spat.jpg'
    clevn = [6, 6.25, 6.5, 6.75,
             7, 7.25, 7.5, 7.75,
             8, 8.25, 8.5, 8.75,
             9, 9.25, 9.5, 9.75,
             10, 10.25, 10.5, 10.75,
             11, 11.25, 11.5, 11.75,
             12, 12.25, 12.5, 12.75,
             13, 13.25, 13.5, 13.75,
             14, 14.25, 14.5, 14.75, 15]

    clevn_str = ['6', '6.25', '6.5', '6.75',
                 '7', '7.25', '7.5', '7.75',
                 '8', '8.25', '8.5', '8.75',
                 '9', '9.25', '9.5', '9.75',
                 '10', '10.25', '10.5', '10.75',
                 '11', '11.25', '11.5', '11.75',
                 '12', '12.25', '12.5', '12.75',
                 '13', '13.25', '13.5', '13.75',
                 '14', '14.25', '14.5', '14.75', '15']
    # clevn_str = ['6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18']
    # spat_temp_plt(croco_temp, sat_temp, clevn, clevn_str, lon_crocod, lat_crocod, h_crocod,
    #               fig_spat_temp, 'Temp $^\circ$ C')

    # Temperature timeseries plot
    fig_temp_temp = crocodir + 'temp_t.jpg'
    tseries_temp_plt(sat_temp, croco_temp, fig_temp_temp)

    # WOA climatology comparisons
    fig_t_n = crocodir + 'woa_n_t.jpg'
    tseries_clim_plt(woa_n_mth, croco_woad_no3, fig_t_n,
                     0, nps_didx_25m, nps_didx_50m, nps_didx_100m, 'WOA NO3')
    fig_t_p = crocodir + 'woa_p_t.jpg'
    tseries_clim_plt(woa_p_mth, croco_woad_po4, fig_t_p,
                     0, nps_didx_25m, nps_didx_50m, nps_didx_100m, 'WOA PO4')
    fig_t_si = crocodir + 'woa_si_t.jpg'
    tseries_clim_plt(woa_s_mth, croco_woad_si, fig_t_si,
                     0, nps_didx_25m, nps_didx_50m, nps_didx_100m, 'WOA Si')
    fig_t_o = crocodir + 'woa_o_t.jpg'
    tseries_clim_plt(woa_o_mth, croco_woad_o2, fig_t_o,
                     0, ox_didx_25m, ox_didx_50m, ox_didx_100m, 'WOA O2')

    # EMODNET climatology comparisons
    fig_t_n = crocodir + 'emc_n_t.jpg'
    tseries_clim_plt(emc_n_mth, croco_emcd_no3, fig_t_n,
                     0, emc_didx_25m, emc_didx_50m, emc_didx_100m, 'EMC NO3')
    fig_t_p = crocodir + 'emc_p_t.jpg'
    tseries_clim_plt(emc_p_mth, croco_emcd_po4, fig_t_p,
                     0, emc_didx_25m, emc_didx_50m, emc_didx_100m, 'EMC PO4')
    fig_t_si = crocodir + 'emc_si_t.jpg'
    tseries_clim_plt(emc_s_mth, croco_emcd_si, fig_t_si,
                     0, emc_didx_25m, emc_didx_50m, emc_didx_100m, 'EMC Si')
    fig_t_o = crocodir + 'emc_o_t.jpg'
    tseries_clim_plt(emc_o_mth, croco_emcd_o2, fig_t_o,
                     0, emc_didx_25m, emc_didx_50m, emc_didx_100m, 'EMC O2')

    # Diatom spatial plot
    fig_spat_diat = crocodir + 'd_chl_spat.jpg'
    clevd = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
             1, 1.5, 2, 3, 4]
    clevd_str = ['0.005', '0.01', '0.02', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9',
                 '1', '1.5', '2', '3', '4']

    for valyr in range(Ystart, Yend + 1):
        fig_spat_diat = crocodir + 'd_chl_spat.jpg'
        croco_fil = crocodir + 'croco_nd_chl_' + str(valyr) + '.npz'
        sat_croco_fil = crocodir + 'sat_croco_nd_chl_' + str(valyr) + '.npz'

        spat_CHL_plt(croco_dchl, sat_dchl, sat_dchle, sat_dchlb, ibi_dchl, clevd, clevd_str,
                     lon_crocod, lat_crocod, h_crocod, fig_spat_diat, 'Diatom')

    # spat_CHL_plt(croco_dchl, sat_dchl_corrected, sat_dchle, sat_dchlb, ibi_dchl, clevd, clevd_str,
    #              lon_crocod, lat_crocod, h_crocod, fig_spat_diat, 'Diatom')

    # Nanophytoplankton spatial plot
    fig_spat_nano = crocodir + 'n_chl_spat.jpg'
    clevn = [0.001, 0.002, 0.003, 0.004, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2,
             0.25, 0.3, 0.35, 0.4, 0.5, 1, 2, 5, 10, 15, 20]
    clevn_str = ['0.001', '0.002', '0.003', '0.004', '0.005', '0.008', '0.01', '0.015', '0.02', '0.03', '0.04',
                 '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.5', '1', '2', '5', '10', '15', '20']
    spat_CHL_plt(croco_nchl, sat_nchl, sat_nchle, sat_nchlb, ibi_nchl, clevn, clevn_str,
                 lon_crocod, lat_crocod, h_crocod, fig_spat_nano, 'Nanophy.')
    # spat_CHL_plt(croco_nchl, sat_nchl_corrected, sat_nchle, sat_nchlb, ibi_nchl, clevn, clevn_str,
    #              lon_crocod, lat_crocod, h_crocod, fig_spat_nano, 'Nanophy.')

    # Combined (diatom + nanophytoplankton) spatial plot
    fig_spat_nd = crocodir + 'nd_chl_spat.jpg'
    clevn = [0.001, 0.002, 0.003, 0.004, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2,
             0.25, 0.3, 0.35, 0.4, 0.5, 1, 2, 5, 10, 15, 20]
    clevn_str = ['0.001', '0.002', '0.003', '0.004', '0.005', '0.008', '0.01', '0.015', '0.02', '0.03', '0.04',
                 '0.05', '0.1', '0.15', '0.2', '0.25', '0.3', '0.35', '0.4', '0.5', '1', '2', '5', '10', '15', '20']
    spat_CHL_plt(croco_nchl + croco_dchl, sat_nchl + sat_dchl, sat_nchle + sat_dchle, sat_nchlb + sat_dchlb,
                 ibi_dchl + ibi_nchl, clevn, clevn_str,
                 lon_crocod, lat_crocod, h_crocod, fig_spat_nd, 'Nano. + Diat. Chl')
    # spat_CHL_plt(croco_nchl + croco_dchl, sat_nchl_corrected + sat_dchl_corrected, sat_nchle + sat_dchle,
    #              sat_nchlb + sat_dchlb, ibi_dchl + ibi_nchl, clevn, clevn_str,
    #              lon_crocod, lat_crocod, h_crocod, fig_spat_nd, 'Nano. + Diat. Chl')

    fig_temp_ndchl = crocodir + 'nano_diat_t.jpg'
    tseries_CHL_plt(sat_dchl, sat_dchle, sat_dchlb, sat_nchl, sat_nchle, sat_nchlb,
                    croco_dchl, croco_nchl, ibi_dchl, ibi_nchl, fig_temp_ndchl)
    # tseries_CHL_plt(sat_dchl_corrected, sat_dchle, sat_dchlb, sat_nchl_corrected, sat_nchle, sat_nchlb,
    #                 croco_dchl, croco_nchl, ibi_dchl, ibi_nchl, fig_temp_ndchl)

    fig_satscat_ndchl = crocodir + 'satCHL_scat.jpg'
    scat_satCHL_plt(sat_dchl, sat_nchl, croco_dchl, croco_nchl, fig_satscat_ndchl)

    fig_spat_spco2 = crocodir + 'clim_spco2_spat.jpg'
    spat_spco2_plt(spco2_socat_mth, spco2_ose_mth, lsch_spco2_mth, woa_lon, woa_lat,
                   croco_woad_spco2, lon_crocod, lat_crocod, h_crocod,
                   fig_spat_spco2)

    fig_t_spco2 = crocodir + 'clim_spco2_t.jpg'
    tseries_spco2_plt(spco2_socat_mth, spco2_ose_mth, lsch_spco2_mth, croco_woad_spco2, fig_t_spco2)

    fig_s_dic_talk = crocodir + 'clim_dic_talk_s.jpg'
    spat_dic_talk_plt(dic_ose_mth, talk_ose_mth, woa_lon, woa_lat,
                      croco_woad_dic, croco_woad_talk, lon_crocod, lat_crocod, h_crocod,
                      fig_s_dic_talk)

    fig_t_dic_talk = crocodir + 'clim_dic_talk_t.jpg'
    tseries_dic_talk_plt(croco_woad_dic, croco_woad_talk, dic_ose_mth, talk_ose_mth, fig_t_dic_talk)

    # NO3 spatial average comparison with WOA climatology
    fig_s_surf_n = crocodir + 'woa_n_s_surf.jpg'
    fig_s_25m_n = crocodir + 'woa_n_s_25m.jpg'
    fig_s_50m_n = crocodir + 'woa_n_s_50m.jpg'
    fig_s_100m_n = crocodir + 'woa_n_s_100m.jpg'
    spat_clim_plt(np.nanmean(woa_n_mth[:, 0, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_no3[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_n, 'WOA average surface NO3')
    spat_clim_plt(np.nanmean(woa_n_mth[:, nps_didx_25m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_no3[:, nps_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_n, 'WOA average 25m NO3')
    spat_clim_plt(np.nanmean(woa_n_mth[:, nps_didx_50m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_no3[:, nps_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_n, 'WOA average 50m NO3')
    spat_clim_plt(np.nanmean(woa_n_mth[:, nps_didx_100m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_no3[:, nps_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_n, 'WOA average 100m NO3')

    # PO4 spatial average comparison with WOA climatology
    fig_s_surf_p = crocodir + 'woa_p_s_surf.jpg'
    fig_s_25m_p = crocodir + 'woa_p_s_25m.jpg'
    fig_s_50m_p = crocodir + 'woa_p_s_50m.jpg'
    fig_s_100m_p = crocodir + 'woa_p_s_100m.jpg'
    spat_clim_plt(np.nanmean(woa_p_mth[:, 0, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_po4[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_p, 'WOA average surface PO4')
    spat_clim_plt(np.nanmean(woa_p_mth[:, nps_didx_25m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_po4[:, nps_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_p, 'WOA average 25m PO4')
    spat_clim_plt(np.nanmean(woa_p_mth[:, nps_didx_50m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_po4[:, nps_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_p, 'WOA average 50m PO4')
    spat_clim_plt(np.nanmean(woa_p_mth[:, nps_didx_100m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_po4[:, nps_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_p, 'WOA average 100m PO4')

    # Si spatial average comparison with WOA climatology
    fig_s_surf_s = crocodir + 'woa_s_s_surf.jpg'
    fig_s_25m_s = crocodir + 'woa_s_s_25m.jpg'
    fig_s_50m_s = crocodir + 'woa_s_s_50m.jpg'
    fig_s_100m_s = crocodir + 'woa_s_s_100m.jpg'
    spat_clim_plt(np.nanmean(woa_s_mth[:, 0, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_si[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_s, 'WOA average surface Si')
    spat_clim_plt(np.nanmean(woa_s_mth[:, nps_didx_25m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_si[:, nps_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_s, 'WOA average 25m Si')
    spat_clim_plt(np.nanmean(woa_s_mth[:, nps_didx_50m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_si[:, nps_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_s, 'WOA average 50m Si')
    spat_clim_plt(np.nanmean(woa_s_mth[:, nps_didx_100m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_si[:, nps_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_s, 'WOA average 100m Si')

    # O2 spatial average comparison with WOA climatology
    fig_s_surf_o = crocodir + 'woa_o_s_surf.jpg'
    fig_s_25m_o = crocodir + 'woa_o_s_25m.jpg'
    fig_s_50m_o = crocodir + 'woa_o_s_50m.jpg'
    fig_s_100m_o = crocodir + 'woa_o_s_100m.jpg'
    spat_clim_plt(np.nanmean(woa_o_mth[:, 0, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_o2[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_o, 'WOA average surface O2')
    spat_clim_plt(np.nanmean(woa_o_mth[:, nps_didx_25m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_o2[:, nps_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_o, 'WOA average 25m O2')
    spat_clim_plt(np.nanmean(woa_o_mth[:, nps_didx_50m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_o2[:, nps_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_o, 'WOA average 50m O2')
    spat_clim_plt(np.nanmean(woa_o_mth[:, nps_didx_100m, :, :], 0), woa_lon, woa_lat,
                  np.nanmean(croco_woad_o2[:, nps_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_o, 'WOA average 100m O2')

    # NO3 spatial average comparison with EMODNET climatology
    fig_s_surf_n = crocodir + 'emc_n_s_surf.jpg'
    fig_s_25m_n = crocodir + 'emc_n_s_25m.jpg'
    fig_s_50m_n = crocodir + 'emc_n_s_50m.jpg'
    fig_s_100m_n = crocodir + 'emc_n_s_100m.jpg'
    nmask = np.isnan(emc_n_mth) | np.isnan(croco_emcd_no3)
    emc_n_mth[nmask] = np.nan
    croco_emcd_no3[nmask] = np.nan
    spat_clim_plt(np.nanmean(emc_n_mth[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_n, 'EMOD average surface NO3')
    spat_clim_plt(np.nanmean(emc_n_mth[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_n, 'EMOD average 25m NO3')
    spat_clim_plt(np.nanmean(emc_n_mth[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_n, 'EMOD average 50m NO3')
    spat_clim_plt(np.nanmean(emc_n_mth[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_n, 'EMOD average 100m NO3')

    # L1 plots
    croco_emcd_no3[emc_n_mth_L1 > 10000] = np.nan
    emc_n_mth_L1[emc_n_mth_L1 > 10000] = np.nan
    nmask_L1 = np.isnan(emc_n_mth_L1) | np.isnan(croco_emcd_no3)
    emc_n_mth_L1[nmask_L1] = np.nan
    croco_emcd_no3[nmask_L1] = np.nan
    spat_clim_plt(np.nanmean(emc_n_mth_L1[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_n[:-10] + 'L1' + fig_s_surf_n[-11:], 'EMOD average surface NO3 (L1)')
    spat_clim_plt(np.nanmean(emc_n_mth_L1[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_n[:-9] + 'L1' + fig_s_25m_n[-10:], 'EMOD average 25m NO3 (L1)')
    spat_clim_plt(np.nanmean(emc_n_mth_L1[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_n[:-9] + 'L1' + fig_s_50m_n[-10:], 'EMOD average 50m NO3 (L1)')
    spat_clim_plt(np.nanmean(emc_n_mth_L1[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_n[:-10] + 'L1' + fig_s_100m_n[-11:], 'EMOD average 100m NO3 (L1)')

    # L2 plots
    croco_emcd_no3[emc_n_mth_L2 > 10000] = np.nan
    emc_n_mth_L2[emc_n_mth_L2 > 10000] = np.nan
    nmask_L2 = np.isnan(emc_n_mth_L2) | np.isnan(croco_emcd_no3)
    emc_n_mth_L2[nmask_L2] = np.nan
    croco_emcd_no3[nmask_L2] = np.nan
    spat_clim_plt(np.nanmean(emc_n_mth_L2[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_n[:-10] + 'L2' + fig_s_surf_n[-11:], 'EMOD average surface NO3 (L2)')
    spat_clim_plt(np.nanmean(emc_n_mth_L2[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_n[:-9] + 'L2' + fig_s_25m_n[-10:], 'EMOD average 25m NO3 (L2)')
    spat_clim_plt(np.nanmean(emc_n_mth_L2[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_n[:-9] + 'L2' + fig_s_50m_n[-10:], 'EMOD average 50m NO3 (L2)')
    spat_clim_plt(np.nanmean(emc_n_mth_L2[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_no3[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_n[:-10] + 'L2' + fig_s_100m_n[-11:], 'EMOD average 100m NO3 (L2)')

    # PO4 spatial average comparison with EMODNET climatology
    fig_s_surf_p = crocodir + 'emc_p_s_surf.jpg'
    fig_s_25m_p = crocodir + 'emc_p_s_25m.jpg'
    fig_s_50m_p = crocodir + 'emc_p_s_50m.jpg'
    fig_s_100m_p = crocodir + 'emc_p_s_100m.jpg'
    pmask = np.isnan(emc_p_mth) | np.isnan(croco_emcd_po4)
    emc_p_mth[pmask] = np.nan
    croco_emcd_po4[pmask] = np.nan
    spat_clim_plt(np.nanmean(emc_p_mth[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_p, 'EMOD average surface PO4')
    spat_clim_plt(np.nanmean(emc_p_mth[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_p, 'EMOD average 25m PO4')
    spat_clim_plt(np.nanmean(emc_p_mth[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_p, 'EMOD average 50m PO4')
    spat_clim_plt(np.nanmean(emc_p_mth[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_p, 'EMOD average 100m PO4')

    # L1 plots
    croco_emcd_po4[emc_p_mth_L1 > 10000] = np.nan
    emc_p_mth_L1[emc_p_mth_L1 > 10000] = np.nan
    pmask_L1 = np.isnan(emc_p_mth_L1) | np.isnan(croco_emcd_po4)
    emc_p_mth_L1[pmask_L1] = np.nan
    croco_emcd_po4[pmask_L1] = np.nan
    spat_clim_plt(np.nanmean(emc_p_mth_L1[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_p[:-10] + 'L1' + fig_s_surf_p[-11:], 'EMOD average surface PO4 (L1)')
    spat_clim_plt(np.nanmean(emc_p_mth_L1[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_p[:-9] + 'L1' + fig_s_25m_p[-10:], 'EMOD average 25m PO4 (L1)')
    spat_clim_plt(np.nanmean(emc_p_mth_L1[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_p[:-9] + 'L1' + fig_s_50m_p[-10:], 'EMOD average 50m PO4 (L1)')
    spat_clim_plt(np.nanmean(emc_p_mth_L1[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_p[:-10] + 'L1' + fig_s_100m_p[-11:], 'EMOD average 100m PO4 (L1)')

    # L2 plots
    croco_emcd_po4[emc_p_mth_L2 > 10000] = np.nan
    emc_p_mth_L2[emc_p_mth_L2 > 10000] = np.nan
    pmask_L2 = np.isnan(emc_p_mth_L2) | np.isnan(croco_emcd_po4)
    emc_p_mth_L2[pmask_L2] = np.nan
    croco_emcd_po4[pmask_L2] = np.nan
    spat_clim_plt(np.nanmean(emc_p_mth_L2[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_p[:-10] + 'L2' + fig_s_surf_p[-11:], 'EMOD average surface PO4 (L2)')
    spat_clim_plt(np.nanmean(emc_p_mth_L2[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_p[:-9] + 'L2' + fig_s_25m_p[-10:], 'EMOD average 25m PO4 (L2)')
    spat_clim_plt(np.nanmean(emc_p_mth_L2[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_p[:-9] + 'L2' + fig_s_50m_p[-10:], 'EMOD average 50m PO4 (L2)')
    spat_clim_plt(np.nanmean(emc_p_mth_L2[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_po4[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_p[:-10] + 'L2' + fig_s_100m_p[-11:], 'EMOD average 100m PO4 (L2)')

    # Si spatial average comparison with EMODNET climatology
    fig_s_surf_s = crocodir + 'emc_s_s_surf.jpg'
    fig_s_25m_s = crocodir + 'emc_s_s_25m.jpg'
    fig_s_50m_s = crocodir + 'emc_s_s_50m.jpg'
    fig_s_100m_s = crocodir + 'emc_s_s_100m.jpg'
    smask = np.isnan(emc_s_mth) | np.isnan(croco_emcd_si)
    emc_s_mth[smask] = np.nan
    croco_emcd_si[smask] = np.nan
    spat_clim_plt(np.nanmean(emc_s_mth[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_s, 'EMOD average surface Si')
    spat_clim_plt(np.nanmean(emc_s_mth[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_s, 'EMOD average 25m Si')
    spat_clim_plt(np.nanmean(emc_s_mth[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_s, 'EMOD average 50m Si')
    spat_clim_plt(np.nanmean(emc_s_mth[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_s, 'EMOD average 100m Si')

    # L1 plots
    croco_emcd_si[emc_s_mth_L1 > 10000] = np.nan
    emc_s_mth_L1[emc_s_mth_L1 > 10000] = np.nan
    smask_L1 = np.isnan(emc_s_mth_L1) | np.isnan(croco_emcd_si)
    emc_s_mth_L1[smask_L1] = np.nan
    croco_emcd_si[smask_L1] = np.nan
    spat_clim_plt(np.nanmean(emc_s_mth_L1[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_s[:-10] + 'L1' + fig_s_surf_s[-11:], 'EMOD average surface Si (L1)')
    spat_clim_plt(np.nanmean(emc_s_mth_L1[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_s[:-9] + 'L1' + fig_s_25m_s[-10:], 'EMOD average 25m Si (L1)')
    spat_clim_plt(np.nanmean(emc_s_mth_L1[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_s[:-9] + 'L1' + fig_s_50m_s[-10:], 'EMOD average 50m Si (L1)')
    spat_clim_plt(np.nanmean(emc_s_mth_L1[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_s[:-10] + 'L1' + fig_s_100m_s[-11:], 'EMOD average 100m Si (L1)')

    # L2 plots
    croco_emcd_si[emc_s_mth_L2 > 10000] = np.nan
    emc_s_mth_L2[emc_s_mth_L2 > 10000] = np.nan
    smask_L2 = np.isnan(emc_s_mth_L2) | np.isnan(croco_emcd_si)
    emc_s_mth_L2[smask_L2] = np.nan
    croco_emcd_si[smask_L2] = np.nan
    spat_clim_plt(np.nanmean(emc_s_mth_L2[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_s[:-10] + 'L2' + fig_s_surf_s[-11:], 'EMOD average surface Si (L2)')
    spat_clim_plt(np.nanmean(emc_s_mth_L2[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_s[:-9] + 'L2' + fig_s_25m_s[-10:], 'EMOD average 25m Si (L2)')
    spat_clim_plt(np.nanmean(emc_s_mth_L2[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_s[:-9] + 'L2' + fig_s_50m_s[-10:], 'EMOD average 50m Si (L2)')
    spat_clim_plt(np.nanmean(emc_s_mth_L2[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_si[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_s[:-10] + 'L2' + fig_s_100m_s[-11:], 'EMOD average 100m Si (L2)')

    # O2 spatial average comparison with EMODNET climatology
    fig_s_surf_o = crocodir + 'emc_o_s_surf.jpg'
    fig_s_25m_o = crocodir + 'emc_o_s_25m.jpg'
    fig_s_50m_o = crocodir + 'emc_o_s_50m.jpg'
    fig_s_100m_o = crocodir + 'emc_o_s_100m.jpg'
    omask = np.isnan(emc_o_mth) | np.isnan(croco_emcd_o2)
    emc_o_mth[omask] = np.nan
    croco_emcd_o2[omask] = np.nan
    spat_clim_plt(np.nanmean(emc_o_mth[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_o, 'EMOD average surface O2')
    spat_clim_plt(np.nanmean(emc_o_mth[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_o, 'EMOD average 25m O2')
    spat_clim_plt(np.nanmean(emc_o_mth[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_o, 'EMOD average 50m O2')
    spat_clim_plt(np.nanmean(emc_o_mth[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_o, 'EMOD average 100m O2')

    # L1 plots
    croco_emcd_o2[emc_o_mth_L1 > 10000] = np.nan
    emc_o_mth_L1[emc_o_mth_L1 > 10000] = np.nan
    omask_L1 = np.isnan(emc_o_mth_L1) | np.isnan(croco_emcd_o2)
    emc_o_mth_L1[omask_L1] = np.nan
    croco_emcd_o2[omask_L1] = np.nan
    spat_clim_plt(np.nanmean(emc_o_mth_L1[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_o[:-10] + 'L1' + fig_s_surf_o[-11:], 'EMOD average surface O2 (L1)')
    spat_clim_plt(np.nanmean(emc_o_mth_L1[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_o[:-9] + 'L1' + fig_s_25m_o[-10:], 'EMOD average 25m O2 (L1)')
    spat_clim_plt(np.nanmean(emc_o_mth_L1[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_o[:-9] + 'L1' + fig_s_50m_o[-10:], 'EMOD average 50m O2 (L1)')
    spat_clim_plt(np.nanmean(emc_o_mth_L1[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_o[:-10] + 'L1' + fig_s_100m_o[-11:], 'EMOD average 100m O2 (L1)')

    # L2 plots
    croco_emcd_o2[emc_o_mth_L2 > 10000] = np.nan
    emc_o_mth_L2[emc_o_mth_L2 > 10000] = np.nan
    omask_L2 = np.isnan(emc_o_mth_L2) | np.isnan(croco_emcd_o2)
    emc_o_mth_L2[omask_L2] = np.nan
    croco_emcd_o2[omask_L2] = np.nan
    spat_clim_plt(np.nanmean(emc_o_mth_L2[:, 0, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, 0, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_surf_o[:-10] + 'L2' + fig_s_surf_o[-11:], 'EMOD average surface O2 (L2)')
    spat_clim_plt(np.nanmean(emc_o_mth_L2[:, emc_didx_25m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_25m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_25m_o[:-9] + 'L2' + fig_s_25m_o[-10:], 'EMOD average 25m O2 (L2)')
    spat_clim_plt(np.nanmean(emc_o_mth_L2[:, emc_didx_50m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_50m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_50m_o[:-9] + 'L2' + fig_s_50m_o[-10:], 'EMOD average 50m O2 (L2)')
    spat_clim_plt(np.nanmean(emc_o_mth_L2[:, emc_didx_100m, :, :], 0), emc_lon, emc_lat,
                  np.nanmean(croco_emcd_o2[:, emc_didx_100m, :, :], 0), lon_crocod, lat_crocod, h_crocod,
                  fig_s_100m_o[:-10] + 'L2' + fig_s_100m_o[-11:], 'EMOD average 100m O2 (L2)')

    # Scatter plots from in-situ data
    insit_real = np.argwhere((issil_cro != 0.) & (isnit_cro != 0.) &
                             (isphs_cro != 0.) & (issal_cro != 0.))[:, 0]
    issil_cro = issil_cro[insit_real]
    isphs_cro = isphs_cro[insit_real]
    isnit_cro = isnit_cro[insit_real]
    issal_cro = issal_cro[insit_real]
    issil = issil[insit_real]
    isphs = isphs[insit_real]
    isnit = isnit[insit_real]
    issal = issal[insit_real]
    isyr = isyr[insit_real]
    ismth = ismth[insit_real]
    isday = isday[insit_real]
    islat = islat[insit_real]
    islon = islon[insit_real]
    isdep = isdep[insit_real]

    fig_scat_alldep = crocodir + 'insitu_alldep.jpg'
    is_scatplot(issil_cro, isphs_cro, isnit_cro, issal_cro, issil, isphs, isnit, issal,
                fig_scat_alldep, '[All depths]')

    # Scatter plots from in-situ data
    insit_3m = np.argwhere(isdep == 3.)[:, 0]
    issil_cro = issil_cro[insit_3m]
    isphs_cro = isphs_cro[insit_3m]
    isnit_cro = isnit_cro[insit_3m]
    issal_cro = issal_cro[insit_3m]
    issil = issil[insit_3m]
    isphs = isphs[insit_3m]
    isnit = isnit[insit_3m]
    issal = issal[insit_3m]
    isyr = isyr[insit_3m]
    ismth = ismth[insit_3m]
    isday = isday[insit_3m]
    islat = islat[insit_3m]
    islon = islon[insit_3m]
    isdep = isdep[insit_3m]

    fig_scat_3m = crocodir + 'insitu_3m.jpg'
    is_scatplot(issil_cro, isphs_cro, isnit_cro, issal_cro, issil, isphs, isnit, issal,
                fig_scat_3m, '[3m]')

    fig_sspat_3m = crocodir + 'insitu_Si_map_3m.jpg'
    is_spatplot(issil_cro, issil, islat, islon, fig_sspat_3m, 'Silicate (mmol/m3)')

    fig_nspat_3m = crocodir + 'insitu_N_map_3m.jpg'
    is_spatplot(isnit_cro, isnit, islat, islon, fig_nspat_3m, 'Nitrate (mmol/m3)')

    fig_pspat_3m = crocodir + 'insitu_P_map_3m.jpg'
    is_spatplot(isphs_cro, isphs, islat, islon, fig_pspat_3m, 'Phosphate (mmol/m3)')

    fig_salspat_3m = crocodir + 'insitu_sal_map_3m.jpg'
    is_spatplot(issal_cro, issal, islat, islon, fig_salspat_3m, 'Salinity (psu)')

    # fig_t_surf_si = crocodir + 'woa_si_t_surf.jpg'
    # fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel('Month')
    # plt.ylabel('mmol/m3')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_s_mth[:, 0, :, :], 2), 1), 'k',
    #                label='WOA Si')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_si[:, 0, :, :], 2), 1), 'b-',
    #                label='CS1KM Si')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_s_mth[:, nps_didx_25m, :, :], 2), 1), 'k',
    #                label='WOA Si')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_si[:, nps_didx_25m, :, :], 2), 1), 'b-',
    #                label='CS1KM Si')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_s_mth[:, nps_didx_50m, :, :], 2), 1), 'k',
    #                label='WOA Si')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_si[:, nps_didx_50m, :, :], 2), 1), 'b-',
    #                label='CS1KM Si')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_s_mth[:, nps_didx_100m, :, :], 2), 1), 'k',
    #                label='WOA Si')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_si[:, nps_didx_100m, :, :], 2), 1), 'b-',
    #                label='CS1KM Si')
    # fig.suptitle('Monthly Si')
    # axs[0, 0].title.set_text('0m')
    # axs[0, 1].title.set_text('25m')
    # axs[1, 0].title.set_text('50m')
    # axs[1, 1].title.set_text('100m')
    # handles, labels = axs[1, 1].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='center', ncol=2)
    # plt.savefig(fig_t_surf_si, dpi='figure', format='jpg')
    #
    # fig_t_surf_n = crocodir + 'woa_n_t_surf.jpg'
    # fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel('Month')
    # plt.ylabel('mmol/m3')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_n_mth[:, 0, :, :], 2), 1), 'k',
    #                label='WOA NO3')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_no3[:, 0, :, :], 2), 1), 'b-',
    #                label='CS1KM NO3')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_n_mth[:, nps_didx_25m, :, :], 2), 1), 'k',
    #                label='WOA NO3')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_no3[:, nps_didx_25m, :, :], 2), 1), 'b-',
    #                label='CS1KM NO3')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_n_mth[:, nps_didx_50m, :, :], 2), 1), 'k',
    #                label='WOA NO3')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_no3[:, nps_didx_50m, :, :], 2), 1), 'b-',
    #                label='CS1KM NO3')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_n_mth[:, nps_didx_100m, :, :], 2), 1), 'k',
    #                label='WOA NO3')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_no3[:, nps_didx_100m, :, :], 2), 1), 'b-',
    #                label='CS1KM NO3')
    # fig.suptitle('Monthly NO3')
    # axs[0, 0].title.set_text('0m')
    # axs[0, 1].title.set_text('25m')
    # axs[1, 0].title.set_text('50m')
    # axs[1, 1].title.set_text('100m')
    # handles, labels = axs[1, 1].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='center', ncol=2)
    # plt.savefig(fig_t_surf_n, dpi='figure', format='jpg')
    #
    # fig_t_surf_p = crocodir + 'woa_p_t_surf.jpg'
    # fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel('Month')
    # plt.ylabel('mmol/m3')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_p_mth[:, 0, :, :], 2), 1), 'k',
    #                label='WOA PO4')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_po4[:, 0, :, :], 2), 1), 'b-',
    #                label='CS1KM PO4')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_p_mth[:, nps_didx_25m, :, :], 2), 1), 'k',
    #                label='WOA PO4')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_po4[:, nps_didx_25m, :, :], 2), 1), 'b-',
    #                label='CS1KM PO4')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_p_mth[:, nps_didx_50m, :, :], 2), 1), 'k',
    #                label='WOA PO4')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_po4[:, nps_didx_50m, :, :], 2), 1), 'b-',
    #                label='CS1KM PO4')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_p_mth[:, nps_didx_100m, :, :], 2), 1), 'k',
    #                label='WOA PO4')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_po4[:, nps_didx_100m, :, :], 2), 1), 'b-',
    #                label='CS1KM PO4')
    # fig.suptitle('Monthly PO4')
    # axs[0, 0].title.set_text('0m')
    # axs[0, 1].title.set_text('25m')
    # axs[1, 0].title.set_text('50m')
    # axs[1, 1].title.set_text('100m')
    # handles, labels = axs[1, 1].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='center', ncol=2)
    # plt.savefig(fig_t_surf_p, dpi='figure', format='jpg')
    #
    # fig_t_surf_o = crocodir + 'woa_o_t_surf.jpg'
    # fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 8))
    # fig.add_subplot(111, frameon=False)
    # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    # plt.xlabel('Month')
    # plt.ylabel('mmol/m3')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_o_mth[:, 0, :, :], 2), 1), 'k',
    #                label='WOA O2')
    # axs[0, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_o2[:, 0, :, :], 2), 1), 'b-',
    #                label='CS1KM O2')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_o_mth[:, nps_didx_25m, :, :], 2), 1), 'k',
    #                label='WOA O2')
    # axs[0, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_o2[:, nps_didx_25m, :, :], 2), 1), 'b-',
    #                label='CS1KM O2')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(woa_o_mth[:, nps_didx_50m, :, :], 2), 1), 'k',
    #                label='WOA O2')
    # axs[1, 0].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_o2[:, nps_didx_50m, :, :], 2), 1), 'b-',
    #                label='CS1KM O2')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(woa_o_mth[:, nps_didx_100m, :, :], 2), 1), 'k',
    #                label='WOA O2')
    # axs[1, 1].plot(range(1, 25), np.nanmean(np.nanmean(croco_woad_o2[:, nps_didx_100m, :, :], 2), 1), 'b-',
    #                label='CS1KM O2')
    # fig.suptitle('Monthly O2')
    # axs[0, 0].title.set_text('0m')
    # axs[0, 1].title.set_text('25m')
    # axs[1, 0].title.set_text('50m')
    # axs[1, 1].title.set_text('100m')
    # handles, labels = axs[1, 1].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='center', ncol=2)
    # plt.savefig(fig_t_surf_o, dpi='figure', format='jpg')

    # cm = plt.cm.get_cmap('RdYlBu_r')
    # nmin = np.nanmin((np.nanmin(np.nanmean(woa_n_mth[:, nps_didx_25m, :, :], 0)),
    #                   np.nanmin(np.nanmean(croco_woad_no3[:, nps_didx_25m, :, :], 0))))
    # nmax = np.nanmin((np.nanmax(np.nanmean(woa_n_mth[:, nps_didx_25m, :, :], 0)),
    #                   np.nanmax(np.nanmean(croco_woad_no3[:, nps_didx_25m, :, :], 0))))
    # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
    # (ax1, ax2) = axs
    # ax1.contourf(lon_crocod, lat_crocod, h_crocod, [0, 10])
    # ax1.contour(woa_lon, woa_lat, np.nanmean(woa_n_mth[:, nps_didx_25m, :, :], 0),
    #             cmap=cm, vmin=nmin, vmax=nmax)
    # cxlim = (np.nanmin(woa_lon), np.nanmax(woa_lon))
    # cylim = (np.nanmin(woa_lat), np.nanmax(woa_lat))
    # ax1.set_xlim(cxlim)
    # ax1.set_ylim(cylim)
    # ax2.contourf(lon_crocod, lat_crocod, h_crocod, [0, 10])
    # im = ax2.contour(woa_lon, woa_lat, np.nanmean(croco_woad_no3[:, nps_didx_25m, :, :], 0),
    #                  cmap=cm, vmin=nmin, vmax=nmax)
    # ax2.set_xlim(cxlim)
    # ax2.set_ylim(cylim)
    # ax2.axes.get_yaxis().set_visible(False)
    # fig.suptitle('24 month average 25m NO3')
    # ax1.title.set_text('WOA')
    # ax2.title.set_text('CS1KM')
    # plt.colorbar(im, ax=axs)
    # plt.savefig(fig_s_25m_n, dpi='figure', format='jpg')
    #
    # fig_s_50m_n = crocodir + 'woa_n_s_50m.jpg'
    # cm = plt.cm.get_cmap('RdYlBu_r')
    # nmin = np.nanmin((np.nanmin(np.nanmean(woa_n_mth[:, nps_didx_50m, :, :], 0)),
    #                   np.nanmin(np.nanmean(croco_woad_no3[:, nps_didx_50m, :, :], 0))))
    # nmax = np.nanmin((np.nanmax(np.nanmean(woa_n_mth[:, nps_didx_50m, :, :], 0)),
    #                   np.nanmax(np.nanmean(croco_woad_no3[:, nps_didx_50m, :, :], 0))))
    # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
    # (ax1, ax2) = axs
    # ax1.contourf(lon_crocod, lat_crocod, h_crocod, [0, 10])
    # ax1.contour(woa_lon, woa_lat, np.nanmean(woa_n_mth[:, nps_didx_50m, :, :], 0),
    #             cmap=cm, vmin=nmin, vmax=nmax)
    # cxlim = (np.nanmin(woa_lon), np.nanmax(woa_lon))
    # cylim = (np.nanmin(woa_lat), np.nanmax(woa_lat))
    # ax1.set_xlim(cxlim)
    # ax1.set_ylim(cylim)
    # ax2.contourf(lon_crocod, lat_crocod, h_crocod, [0, 10])
    # im = ax2.contour(woa_lon, woa_lat, np.nanmean(croco_woad_no3[:, nps_didx_50m, :, :], 0),
    #                  cmap=cm, vmin=nmin, vmax=nmax)
    # ax2.set_xlim(cxlim)
    # ax2.set_ylim(cylim)
    # ax2.axes.get_yaxis().set_visible(False)
    # fig.suptitle('24 month average 50m NO3')
    # ax1.title.set_text('WOA')
    # ax2.title.set_text('CS1KM')
    # plt.colorbar(im, ax=axs)
    # plt.savefig(fig_s_50m_n, dpi='figure', format='jpg')
    #
    # fig_s_100m_n = crocodir + 'woa_n_s_100m.jpg'
    # cm = plt.cm.get_cmap('RdYlBu_r')
    # nmin = np.nanmin((np.nanmin(np.nanmean(woa_n_mth[:, nps_didx_100m, :, :], 0)),
    #                   np.nanmin(np.nanmean(croco_woad_no3[:, nps_didx_100m, :, :], 0))))
    # nmax = np.nanmin((np.nanmax(np.nanmean(woa_n_mth[:, nps_didx_100m, :, :], 0)),
    #                   np.nanmax(np.nanmean(croco_woad_no3[:, nps_didx_100m, :, :], 0))))
    # fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 8))
    # (ax1, ax2) = axs
    # ax1.contourf(lon_crocod, lat_crocod, h_crocod, [0, 10])
    # ax1.contour(woa_lon, woa_lat, np.nanmean(woa_n_mth[:, nps_didx_100m, :, :], 0),
    #             cmap=cm, vmin=nmin, vmax=nmax)
    # cxlim = (np.nanmin(woa_lon), np.nanmax(woa_lon))
    # cylim = (np.nanmin(woa_lat), np.nanmax(woa_lat))
    # ax1.set_xlim(cxlim)
    # ax1.set_ylim(cylim)
    # ax2.contourf(lon_crocod, lat_crocod, h_crocod, [0, 10])
    # im = ax2.contour(woa_lon, woa_lat, np.nanmean(croco_woad_no3[:, nps_didx_100m, :, :], 0),
    #                  cmap=cm, vmin=nmin, vmax=nmax)
    # ax2.set_xlim(cxlim)
    # ax2.set_ylim(cylim)
    # ax2.axes.get_yaxis().set_visible(False)
    # fig.suptitle('24 month average 100m NO3')
    # ax1.title.set_text('WOA')
    # ax2.title.set_text('CS1KM')
    # plt.colorbar(im, ax=axs)
    # plt.savefig(fig_s_100m_n, dpi='figure', format='jpg')


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

