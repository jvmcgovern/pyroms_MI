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
# Task here is to derive temperature and salinity biases with weekly Obs-based product
# Specify year of interest
# Find all files corresponding to that year
# Get date of first obs product, if not 1st Jan then take last product from preceding year
# Get date of last obs product, if not 31st December, take first product from following year


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import num2date
from scipy.interpolate import griddata

import glob
from netCDF4 import Dataset as netcdf
from netCDF4 import date2index as d2i
from calendar import monthrange
from datetime import date, datetime
import numpy.matlib
import time
from datetime import date, datetime
from scipy.interpolate import RegularGridInterpolator

import sys

# sys.path.insert(0,'/XXX/')

import croco_vgrid as vgrd
import croco_glorys as glor
from interp_Cgrid import *
from progressbar import *
from scipy.spatial import Delaunay


# Finding which keys in the netcdf are common to the dictionary of terms
# List is outputted as latitude, longitude, variable, and depth if possible
def keysearch(n, f):
    syek = []  # keys
    seulav = []  # values
    for item in f.keys():
        if item in n.keys():
            syek.append(item)
            seulav.append(f.get(item))
    return syek, seulav


if 1 == 1:
    # def main_func():

    #
    # #################### USERS DEFINED VARIABLES ########################
    #

    title = 'Boundary file using GLORYS'

    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']

    # my_home_dir = '/home/penven/'

    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    # bryname = crocofiles_dir + 'croco_bry_PHYBIO_CELTIC_h6.nc'
    # bryname = crocofiles_dir + 'croco_bry_PHYBIO_CELTIC_h8.nc'
    # bryname = crocofiles_dir + 'croco_bry_MERCATOR_.nc'
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    # grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m
    # grdname = crocofiles_dir + 'CELTIC_grd_v6_MSL_h8_v4.nc'  # LAT to MSL included, hmin = 8m
    grdname = crocofiles_dir + 'croco_grd.nc'

    N = 20
    theta_s = 7.
    theta_b = 0.
    # hc = 20.
    hc = 50.
    vtransform = 2

    obc = [1, 1, 1, 1]  # open boundaries (1=open , [S E N W])

    time_bry = 0.
    cycle_bry = 0.

    Yorig = 1990  # year origin of time : days since Yorig-01-01

    Ystart = 2017
    Mstart = 1
    Dstart = 1

    # Yend = 2019
    # Mend = 12
    # Dend = 31

    Yend = 2017
    Mend = 12
    Dend = 31

    glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    glorys_ending = '_R20201201_RE01.nc'

    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI_PISCESvars_012013_202112/monthly/'
    PISCES24_prefix = 'IBI36_cg_1m-m_'
    PISCES24_ending = '_3DT-bgc_hcst.nc'

    'dataset-armor-3d-rep-weekly_20171018T1200Z_P20201024T1138Z_II.nc'
    obs_dir = '/media/dskone/CELTIC/CMEMS_MYNRT_015_012/'
    obs_prefix = 'dataset-armor-3d-*-weekly_'
    obs_ending = 'T1200Z_P*T*Z.nc'

    glorys_step = 1  # time step between outputs in GLORYS12 [days]

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    ibifile_tot = sorted(glob.glob(glorysfiles_dir + glorys_prefix + '*_*' + glorys_ending))
    nc1 = netcdf(ibifile_tot[0], 'r')
    nc1o = nc1.variables
    lat_ibi = np.array(nc1o['latitude'][:])
    latimin = np.min(lat_ibi)
    latimax = np.max(lat_ibi)
    lon_ibi = np.array(nc1o['longitude'][:])
    lonimin = np.min(lon_ibi)
    lonimax = np.max(lon_ibi)
    ibisamp = np.array(nc1o['thetao'][:])
    ibi_na = nc1o['thetao'].add_offset
    nc1.close()

    obsfile_tot = sorted(glob.glob(obs_dir + obs_prefix + '*' + obs_ending))
    nc0 = netcdf(obsfile_tot[0], 'r')
    nc0o = nc0.variables
    lat_obs = np.array(nc0o['latitude'][:])
    lon_obs = np.array(nc0o['longitude'][:])
    lon_obs = np.concatenate([lon_obs[720:]-360, lon_obs[:720]])

    olonidx = np.argwhere((lon_obs >= lonimin) & (lon_obs <= lonimax))[:, 0]
    olatidx = np.argwhere((lat_obs >= latimin) & (lat_obs <= latimax))[:, 0]

    lon_obs = lon_obs[olonidx]
    lat_obs = lat_obs[olatidx]

    crs_fin_grd = np.empty((ibisamp.shape[1], lat_obs.shape[0], lon_obs.shape[0]), dtype=object)

    latr = abs(lat_obs[0] - lat_obs[1])
    lonr = abs(lon_obs[0] - lon_obs[1])

    ltdelta = 0.5 * latr
    ltmin = lat_obs - ltdelta
    ltmax = lat_obs + ltdelta

    lndelta = 0.5 * lonr
    lnmin = lon_obs - lndelta
    lnmax = lon_obs + lndelta

    lat_ibi = np.matlib.repmat(lat_ibi, lon_ibi.shape[0], 1).transpose()
    lon_ibi = np.matlib.repmat(lon_ibi, lat_ibi.shape[0], 1)

    sta = time.time()
    for lts in range(0, lat_obs.shape[0]):
        for lns in range(0, lon_obs.shape[0]):
            for d in range(0, ibisamp.shape[1]):
                crs_fin_grd[d, lts, lns] = np.argwhere(((lat_ibi >= ltmin[lts]) &
                                                        (lat_ibi <= ltmax[lts])) &
                                                       ((lon_ibi >= lnmin[lns]) &
                                                        (lon_ibi <= lnmax[lns])) &
                                                       (np.squeeze(ibisamp[0, d, :, :]) != ibi_na))
    dne = time.time()
    print('indexing takes', dne - sta, 'seconds')

    for iyear in range(Ystart, Yend + 1):
        obsfiles_pre = sorted(glob.glob(obs_dir + obs_prefix + str(iyear - 1) + '*' + obs_ending))
        obsfiles_main = sorted(glob.glob(obs_dir + obs_prefix + str(iyear) + '*' + obs_ending))
        obsfiles_post = sorted(glob.glob(obs_dir + obs_prefix + str(iyear + 1) + '*' + obs_ending))
        obsfiles = [obsfiles_pre[-1]] + obsfiles_main + [obsfiles_post[0]]

        for wk in range(0, len(obsfiles)-2):
            wkstart = obsfiles[wk][69:77]
            wkend = obsfiles[wk+1][69:77]

            ncostart = netcdf(obsfiles[wk], 'r')
            ncostarto = ncostart.variables
            obs_to_start = np.array(ncostarto['to'][:])
            obs_to_start = np.concatenate([obs_to_start[:, :, :, 720:], obs_to_start[:, :, :, :720]], axis=3)
            obs_to_start = obs_to_start[:, :, olatidx[:, None], olonidx]
            obs_so_start = np.array(ncostarto['so'][:])
            obs_so_start = np.concatenate([obs_so_start[:, :, :, 720:], obs_so_start[:, :, :, :720]], axis=3)
            obs_so_start = obs_so_start[:, :, olatidx[:, None], olonidx]
            obs_to_FillValue = ncostarto['to']._FillValue
            obs_so_FillValue = ncostarto['so']._FillValue
            ncostart.close()

            ncoend = netcdf(obsfiles[wk+1], 'r')
            ncoendo = ncoend.variables
            obs_to_end = np.array(ncoendo['to'][:])
            obs_to_end = np.concatenate([obs_to_end[:, :, :, 720:], obs_to_end[:, :, :, :720]], axis=3)
            obs_to_end = obs_to_end[:, :, olatidx[:, None], olonidx]
            obs_so_end = np.array(ncoendo['so'][:])
            obs_so_end = np.concatenate([obs_so_end[:, :, :, 720:], obs_so_end[:, :, :, :720]], axis=3)
            obs_so_end = obs_so_end[:, :, olatidx[:, None], olonidx]
            ncoend.close()

            ibifilestart = sorted(glob.glob(glorysfiles_dir + glorys_prefix + wkstart + '_' + wkstart + glorys_ending))
            ibifileend = sorted(glob.glob(glorysfiles_dir + glorys_prefix + wkend + '_' + wkend + glorys_ending))

            ncistart = netcdf(ibifilestart[0], 'r')
            ncistarto = ncistart.variables
            ibi_to_start = np.array(ncistarto['thetao'][:])
            ibi_so_start = np.array(ncistarto['so'][:])
            ibi_to_offset = ncistarto['thetao'].add_offset
            ibi_so_offset = ncistarto['so'].add_offset
            ncistart.close()

            nciend = netcdf(ibifileend[0], 'r')
            nciendo = nciend.variables
            ibi_to_end = np.array(nciendo['thetao'][:])
            ibi_so_end = np.array(nciendo['so'][:])
            nciend.close()

            obs_to_start[obs_to_start == obs_to_FillValue] = ibi_to_offset
            obs_to_end[obs_to_end == obs_to_FillValue] = ibi_to_offset
            obs_so_start[obs_so_start == obs_so_FillValue] = ibi_so_offset
            obs_so_end[obs_so_end == obs_so_FillValue] = ibi_so_offset

            ibi_to_start_obsres = np.zeros_like(obs_to_start)
            ibi_to_end_obsres = np.zeros_like(obs_to_end)
            ibi_so_start_obsres = np.zeros_like(obs_so_start)
            ibi_so_end_obsres = np.zeros_like(obs_so_end)
            sta = time.time()
            for lts in range(0, lat_obs.shape[0]):
                for lns in range(0, lon_obs.shape[0]):
                    for d in range(0, ibi_to_start.shape[1]):
                        if crs_fin_grd[d, lts, lns].shape[0] != 0:
                            ibi_to_start_obsres[0, d, lts, lns] = np.mean(ibi_to_start[0, d,
                                                                          crs_fin_grd[d, lts, lns][:, 0],
                                                                          crs_fin_grd[d, lts, lns][:, 1]])
                            ibi_to_end_obsres[0, d, lts, lns] = np.mean(ibi_to_end[0, d,
                                                                                   crs_fin_grd[d, lts, lns][:, 0],
                                                                                   crs_fin_grd[d, lts, lns][:, 1]])
                            ibi_so_start_obsres[0, d, lts, lns] = np.mean(ibi_so_start[0, d,
                                                                          crs_fin_grd[d, lts, lns][:, 0],
                                                                          crs_fin_grd[d, lts, lns][:, 1]])
                            ibi_so_end_obsres[0, d, lts, lns] = np.mean(ibi_so_end[0, d,
                                                                                   crs_fin_grd[d, lts, lns][:, 0],
                                                                                   crs_fin_grd[d, lts, lns][:, 1]])
                        else:
                            ibi_to_start_obsres[0, d, lts, lns] = ibi_to_offset
                            ibi_to_end_obsres[0, d, lts, lns] = ibi_to_offset
                            ibi_so_start_obsres[0, d, lts, lns] = ibi_so_offset
                            ibi_so_end_obsres[0, d, lts, lns] = ibi_so_offset
            dne = time.time()
            print('averaging for 1 week takes', dne-sta, 'seconds')

            ibi_obs_to_mask = np.logical_or(ibi_to_start_obsres == ibi_to_offset, obs_to_start == ibi_to_offset)
            # ibi_to_start_obsres[obs_to_start == ibi_to_offset] = ibi_to_offset
            # ibi_to_end_obsres[obs_to_end == ibi_to_offset] = ibi_to_offset
            ibi_to_start_obsres[ibi_obs_to_mask] = ibi_to_offset
            obs_to_start[ibi_obs_to_mask] = ibi_to_offset
            ibi_to_end_obsres[ibi_obs_to_mask] = ibi_to_offset
            obs_to_end[ibi_obs_to_mask] = ibi_to_offset

            ibi_obs_so_mask = np.logical_or(ibi_so_start_obsres == ibi_so_offset, obs_so_start == ibi_so_offset)
            # ibi_so_start_obsres[obs_so_start == ibi_so_offset] = ibi_so_offset
            # ibi_so_end_obsres[obs_so_end == ibi_so_offset] = ibi_so_offset
            ibi_so_start_obsres[ibi_obs_so_mask] = ibi_so_offset
            obs_so_start[ibi_obs_so_mask] = ibi_so_offset
            ibi_so_end_obsres[ibi_obs_so_mask] = ibi_so_offset
            obs_so_end[ibi_obs_so_mask] = ibi_so_offset

            # to_start_bias = ibi_to_start_obsres - obs_to_start
            # to_end_bias = ibi_to_end_obsres - obs_to_end
            # to_start_bias[to_start_bias == np.nan] = 0
            # to_end_bias[to_end_bias == np.nan] = 0
            to_start_bias = obs_to_start - ibi_to_start_obsres
            to_end_bias = obs_to_end - ibi_to_end_obsres

            # so_start_bias = ibi_so_start_obsres - obs_so_start
            # so_end_bias = ibi_so_end_obsres - obs_so_end
            # so_start_bias[so_start_bias == np.nan] = 0
            # so_end_bias[so_end_bias == np.nan] = 0
            so_start_bias = obs_so_start - ibi_so_start_obsres
            so_end_bias = obs_so_end - ibi_so_end_obsres

            ndays = date.toordinal(date(int(wkend[0:4]), int(wkend[4:6]), int(wkend[6:8]))) - \
                    date.toordinal(date(int(wkstart[0:4]), int(wkstart[4:6]), int(wkstart[6:8])))

            to_bias_week = np.zeros((ndays+1, ibi_to_end.shape[1], ibi_to_end.shape[2], ibi_to_end.shape[3]),
                                    dtype=float)
            so_bias_week = np.zeros((ndays+1, ibi_so_end.shape[1], ibi_so_end.shape[2], ibi_so_end.shape[3]),
                                    dtype=float)
            lat_obs_int = np.append(lat_obs[0]-latr, lat_obs)
            lat_obs_int = np.append(lat_obs_int, lat_obs_int[-1] + latr)
            lon_obs_int = np.append(lon_obs[0]-lonr, lon_obs)
            lon_obs_int = np.append(lon_obs_int, lon_obs_int[-1] + lonr)
            for d in range(0, ibi_to_start.shape[1]):
                # Interpolate temp bias at start of week
                to_start_bias_slice = np.squeeze(to_start_bias[0, d, :, :])
                to_start_bias_slice = np.concatenate((np.expand_dims(to_start_bias_slice[:, 0], axis=1),
                                                      to_start_bias_slice,
                                                      np.expand_dims(to_start_bias_slice[:, -1], axis=1)), axis=1)
                to_start_bias_slice = np.concatenate((np.expand_dims(to_start_bias_slice[0, :], axis=0),
                                                      to_start_bias_slice,
                                                      np.expand_dims(to_start_bias_slice[-1, :], axis=0)), axis=0)

                bath_int_st_t = RegularGridInterpolator((lat_obs_int, lon_obs_int),
                                                        to_start_bias_slice)
                to_st_b_flat = bath_int_st_t((lat_ibi.flat, lon_ibi.flat))
                to_bias_week[0, d, :, :] = np.reshape(to_st_b_flat, lon_ibi.shape)

                # Interpolate temp bias at end of week
                to_end_bias_slice = np.squeeze(to_end_bias[0, d, :, :])
                to_end_bias_slice = np.concatenate((np.expand_dims(to_end_bias_slice[:, 0], axis=1),
                                                    to_end_bias_slice,
                                                    np.expand_dims(to_end_bias_slice[:, -1], axis=1)), axis=1)
                to_end_bias_slice = np.concatenate((np.expand_dims(to_end_bias_slice[0, :], axis=0),
                                                    to_end_bias_slice,
                                                    np.expand_dims(to_end_bias_slice[-1, :], axis=0)), axis=0)

                bath_int_end_t = RegularGridInterpolator((lat_obs_int, lon_obs_int),
                                                         to_end_bias_slice)
                to_end_b_flat = bath_int_end_t((lat_ibi.flat, lon_ibi.flat))
                to_bias_week[7, d, :, :] = np.reshape(to_end_b_flat, lon_ibi.shape)

                # Interpolate salinity bias at start of week
                so_start_bias_slice = np.squeeze(so_start_bias[0, d, :, :])
                so_start_bias_slice = np.concatenate((np.expand_dims(so_start_bias_slice[:, 0], axis=1),
                                                      so_start_bias_slice,
                                                      np.expand_dims(so_start_bias_slice[:, -1], axis=1)), axis=1)
                so_start_bias_slice = np.concatenate((np.expand_dims(so_start_bias_slice[0, :], axis=0),
                                                      so_start_bias_slice,
                                                      np.expand_dims(so_start_bias_slice[-1, :], axis=0)), axis=0)

                bath_int_st_s = RegularGridInterpolator((lat_obs_int, lon_obs_int),
                                                        so_start_bias_slice)
                so_st_b_flat = bath_int_st_s((lat_ibi.flat, lon_ibi.flat))
                so_bias_week[0, d, :, :] = np.reshape(so_st_b_flat, lon_ibi.shape)

                # Interpolate salinity bias at end of week
                so_end_bias_slice = np.squeeze(so_end_bias[0, d, :, :])
                so_end_bias_slice = np.concatenate((np.expand_dims(so_end_bias_slice[:, 0], axis=1),
                                                    so_end_bias_slice,
                                                    np.expand_dims(so_end_bias_slice[:, -1], axis=1)), axis=1)
                so_end_bias_slice = np.concatenate((np.expand_dims(so_end_bias_slice[0, :], axis=0),
                                                    so_end_bias_slice,
                                                    np.expand_dims(so_end_bias_slice[-1, :], axis=0)), axis=0)

                bath_int_end_s = RegularGridInterpolator((lat_obs_int, lon_obs_int),
                                                         so_end_bias_slice)
                so_end_b_flat = bath_int_end_s((lat_ibi.flat, lon_ibi.flat))
                so_bias_week[7, d, :, :] = np.reshape(so_end_b_flat, lon_ibi.shape)
            # Interpolate temporally to different days of intervening week

            to_bias_delta = (to_bias_week[7, :, :, :] - to_bias_week[0, :, :, :])/7
            so_bias_delta = (so_bias_week[7, :, :, :] - so_bias_week[0, :, :, :])/7

            for dst in range(1, 7):
                to_bias_week[dst, :, :, :] = to_bias_week[0, :, :, :] + (dst*to_bias_delta)
                so_bias_week[dst, :, :, :] = so_bias_week[0, :, :, :] + (dst*so_bias_delta)

            # List of 8 files
            setstart = ibifile_tot.index(ibifilestart[0])
            setend = ibifile_tot.index(ibifileend[0])
            ibifile_set = ibifile_tot[setstart:setend+1]

            varlist = ['thetao_b', 'thetao_bc', 'so_b', 'so_bc']

            # If wk is 1, cycle through all files
            # Otherwise, only do latter 7 files
            for bc in range(0, len(ibifile_set)-1):
                bcfile = ibifile_set[bc]
                nctest = netcdf(bcfile, 'r+')
                nctesto = nctest.variables
                tname, _ = keysearch(nctest.variables,
                                     {'thetao_b': 'thetao', 'thetao_bc': 'thetao', 'so_b': 'so', 'so_bc': 'so'})
                #
                if 'thetao_b' not in tname:
                    nc_thetao_b = nctest.createVariable('thetao_b', np.float64,
                                                        ('time', 'depth', 'latitude', 'longitude',),
                                                        fill_value=1.e+37, complevel=1)
                    nc_thetao_b.long_name = 'Temperature bias'
                    nc_thetao_b.units = 'degrees_C'
                    nc_thetao_b.standard_name = 'sea_water_potential_temperature_bias'
                    nc_thetao_b.add_offset = 10.
                    nc_thetao_b.scale_factor = 0.001
                    nc_thetao_b.unit_long = 'degrees_C'
                    nc_thetao_b.valid_min = -12000
                    nc_thetao_b.valid_max = 22000
                    nc_thetao_b[0, :, :, :] = np.squeeze(to_bias_week[bc, :, :, :])
                else:
                    nctest['thetao_b'][0, :, :, :] = np.squeeze(to_bias_week[bc, :, :, :])

                #
                if 'thetao_bc' not in tname:
                    thetao = np.array(nctest['thetao'][:])
                    nc_thetao_bc = nctest.createVariable('thetao_bc', np.float64,
                                                         ('time', 'depth', 'latitude', 'longitude',),
                                                         fill_value=1.e+37, complevel=1)
                    nc_thetao_bc.long_name = 'Temperature bias corrected'
                    nc_thetao_bc.units = 'degrees_C'
                    nc_thetao_bc.standard_name = 'sea_water_potential_temperature_bias_corrected'
                    nc_thetao_bc.add_offset = 10.
                    nc_thetao_bc.scale_factor = 0.001
                    nc_thetao_bc.unit_long = 'degrees_C'
                    nc_thetao_bc.valid_min = -12000
                    nc_thetao_bc.valid_max = 22000
                    nc_thetao_bc[0, :, :, :] = np.squeeze(thetao[0, :, :, :] + to_bias_week[bc, :, :, :])
                else:
                    thetao = np.array(nctest['thetao'][:])
                    nctest['thetao_bc'][0, :, :, :] = np.squeeze(thetao[0, :, :, :] + to_bias_week[bc, :, :, :])

                if 'so_b' not in tname:
                    nc_so_b = nctest.createVariable('so_b', np.float64,
                                                    ('time', 'depth', 'latitude', 'longitude',),
                                                    fill_value=1.e+37, complevel=1)
                    nc_so_b.long_name = 'Salinity bias'
                    nc_so_b.units = '1e-3'
                    nc_so_b.standard_name = 'sea_water_salinity_bias'
                    nc_so_b.add_offset = 20.
                    nc_so_b.scale_factor = 0.001
                    nc_so_b.unit_long = 'Practical Salinity Unit'
                    nc_so_b.valid_min = -20000
                    nc_so_b.valid_max = 20000
                    nc_so_b[0, :, :, :] = np.squeeze(so_bias_week[bc, :, :, :])
                else:
                    nctest['so_b'][0, :, :, :] = np.squeeze(so_bias_week[bc, :, :, :])

                if 'so_bc' not in tname:
                    so = np.array(nctest['so'][:])
                    nc_so_bc = nctest.createVariable('so_bc', np.float64,
                                                     ('time', 'depth', 'latitude', 'longitude',),
                                                     fill_value=1.e+37, complevel=1)
                    nc_so_bc.long_name = 'Salinity bias corrected'
                    nc_so_bc.units = '1e-3'
                    nc_so_bc.standard_name = 'sea_water_salinity_bias_corrected'
                    nc_so_bc.add_offset = 20.
                    nc_so_bc.scale_factor = 0.001
                    nc_so_bc.unit_long = 'Practical Salinity Unit'
                    nc_so_bc.valid_min = -20000
                    nc_so_bc.valid_max = 20000
                    nc_so_bc[0, :, :, :] = np.squeeze(so[0, :, :, :] + so_bias_week[bc, :, :, :])
                else:
                    so = np.array(nctest['so'][:])
                    nctest['so_bc'][0, :, :, :] = np.squeeze(so[0, :, :, :] + so_bias_week[bc, :, :, :])
                nctest.close()
                print(bcfile, '==> bias corrected')

        # Need to spatially interpolate coarsely gridded bias onto CMEMS IBI grid, layer by layer
        # Then need to temporally interpolate for intervening 7 days
        # Resultant data can be added to the existing array or saved as a new, separate variable
        # How about:
        # thetao, thetao_b, thetao_bc
        # so, so_b, so_bc

        # # for imonth in range(Mstart, Mend + 1):
        #
        #     bryname = crocofiles_dir + 'croco_bry_MERCATOR_hc50m' + '_Y' + str(iyear) + 'M' \
        #               + str(imonth).zfill(2) + '.nc'
        #
        #     dpm = monthrange(iyear, imonth)
        #
        #     print(' ')
        #     print(' Making boundary file: ' + bryname)
        #     print(' ')
        #     print(' Title: ' + title)
        #
        #     #
        #     # Create the CROCO boundary file
        #     #
        #     glor.create_bryfile_PISCES_NORESM(bryname, grdname, title, obc,
        #                                       theta_s, theta_b, hc, N,
        #                                       time_bry, cycle_bry, vtransform)
        #
        #     #
        #     # get the CROCO grid
        #     #
        #     ncg = netcdf(grdname, 'r')
        #     ncgrd = ncg.variables
        #     lon_rho = np.array(ncgrd['lon_rho'][:])
        #     lat_rho = np.array(ncgrd['lat_rho'][:])
        #     lon_u = np.array(ncgrd['lon_u'][:])
        #     lat_u = np.array(ncgrd['lat_u'][:])
        #     lon_v = np.array(ncgrd['lon_v'][:])
        #     lat_v = np.array(ncgrd['lat_v'][:])
        #     h = np.array(ncgrd['h'][:])
        #     mask = np.array(ncgrd['mask_rho'][:])
        #     angle = np.array(ncgrd['angle'][:])
        #     [M, L] = np.shape(lon_rho)
        #     ncg.close()
        #
        #     #
        #     # Open the CROCO boundary file for writing
        #     #
        #     ncbr = netcdf(bryname, 'a')
        #     ncbry = ncbr.variables
        #
        #     #
        #     # Get the Delaunay triangulations for each boundary
        #     #
        #
        #     print(' ')
        #     print(' Get the Delaunay triangulations for each boundary')
        #     print(' ')
        #
        #     #
        #     # Get the first GLORYS file name (for the positions of variables)
        #     #
        #
        #     if imonth is 1:
        #         mercator_PHY_files_pre = sorted(glob.glob(glorysfiles_dir +
        #                                                   glorys_prefix + str(iyear-1) +
        #                                                   str(12).zfill(2) + '31_' + str(iyear-1) +
        #                                                   str(12).zfill(2) + '31' + glorys_ending))
        #         # Get the time in days since Yorig, 1, 1
        #         Tstart = date.toordinal(date(iyear-1, 12, 31)) - date.toordinal(date(Yorig, 1, 1))
        #         Tstart = Tstart + 0.5  # 12H
        #     else:
        #         mercator_PHY_files_pre = sorted(glob.glob(glorysfiles_dir +
        #                                                   glorys_prefix + str(iyear) +
        #                                                   str(imonth-1).zfill(2) +
        #                                                   str(monthrange(iyear, imonth-1)[1]) +
        #                                                   '_' + str(iyear) + str(imonth-1).zfill(2) +
        #                                                   str(monthrange(iyear, imonth-1)[1]) + glorys_ending))
        #         # Get the time in days since Yorig, 1, 1
        #         Tstart = date.toordinal(date(iyear, imonth-1,
        #                                      monthrange(iyear, imonth-1)[1])) - date.toordinal(date(Yorig, 1, 1))
        #         Tstart = Tstart + 0.5  # 12H
        #
        #     mercator_PHY_files_main = sorted(glob.glob(glorysfiles_dir +
        #                                                glorys_prefix + str(iyear) +
        #                                                str(imonth).zfill(2) + '??_' + str(iyear) +
        #                                                str(imonth).zfill(2)
        #                                                + '??_R*_RE01.nc'))
        #     if imonth is 12:
        #         if len(mercator_PHY_files_main) is dpm[1]:
        #             mercator_PHY_files_post1 = sorted(glob.glob(glorysfiles_dir +
        #                                               glorys_prefix + str(iyear+1) +
        #                                               str(1).zfill(2) + str(1).zfill(2) + '_'
        #                                               + str(iyear+1) +
        #                                               str(1).zfill(2) + str(1).zfill(2) + glorys_ending))
        #             mercator_PHY_files_post2 = sorted(glob.glob(glorysfiles_dir +
        #                                               glorys_prefix + str(iyear+1) +
        #                                               str(1).zfill(2) + str(2).zfill(2) + '_'
        #                                               + str(iyear+1) +
        #                                               str(1).zfill(2) + str(2).zfill(2) + glorys_ending))
        #             mercator_PHY_files = \
        #                 mercator_PHY_files_pre + mercator_PHY_files_main + mercator_PHY_files_post1\
        #                 + mercator_PHY_files_post2
        #
        #             Tend = date.toordinal(date(iyear+1, 1, 2)) - date.toordinal(date(Yorig, 1, 1))
        #             Tend = Tend + 0.5  # 12H
        #         else:
        #             mercator_PHY_files = mercator_PHY_files_pre + mercator_PHY_files_main
        #             Tend = date.toordinal(date(iyear, imonth,
        #                                        len(mercator_PHY_files)-1)) - date.toordinal(date(Yorig, 1, 1))
        #             Tend = Tend + 0.5  # 12H
        #     else:
        #         mercator_PHY_files_post1 = sorted(glob.glob(glorysfiles_dir +
        #                                                     glorys_prefix + str(iyear) +
        #                                                     str(imonth+1).zfill(2) + str(1).zfill(2) +
        #                                                     '_' + str(iyear) + str(imonth+1).zfill(2) +
        #                                                     str(1).zfill(2) + glorys_ending))
        #         mercator_PHY_files_post2 = sorted(glob.glob(glorysfiles_dir +
        #                                                     glorys_prefix + str(iyear) +
        #                                                     str(imonth+1).zfill(2) + str(2).zfill(2) +
        #                                                     '_' + str(iyear) + str(imonth+1).zfill(2) +
        #                                                     str(2).zfill(2) + glorys_ending))
        #         mercator_PHY_files = \
        #             mercator_PHY_files_pre + mercator_PHY_files_main + mercator_PHY_files_post1 \
        #             + mercator_PHY_files_post2
        #         Tend = date.toordinal(date(iyear, imonth+1, 2)) - date.toordinal(date(Yorig, 1, 1))
        #         Tend = Tend + 0.5  # 12H
        #
        #     if imonth is 1:
        #         mercator_BIO_files_pre = sorted(glob.glob(PISCES24files_dir +
        #                                                   PISCES24_prefix + str(iyear-1) +
        #                                                   str(12).zfill(2) + PISCES24_ending))
        #     else:
        #         mercator_BIO_files_pre = sorted(glob.glob(PISCES24files_dir +
        #                                                   PISCES24_prefix + str(iyear) +
        #                                                   str(imonth-1).zfill(2) + PISCES24_ending))
        #
        #     mercator_BIO_files_main = sorted(glob.glob(PISCES24files_dir +
        #                                                PISCES24_prefix + str(iyear) +
        #                                                str(imonth).zfill(2) + PISCES24_ending))
        #     if imonth is 12:
        #         mercator_BIO_files_post = sorted(glob.glob(PISCES24files_dir +
        #                                                    PISCES24_prefix + str(iyear+1) + str(1).zfill(2) +
        #                                                    PISCES24_ending))
        #         mercator_BIO_files = \
        #             mercator_BIO_files_pre + mercator_BIO_files_main + mercator_BIO_files_post
        #     else:
        #         mercator_BIO_files_post = sorted(glob.glob(PISCES24files_dir +
        #                                                    PISCES24_prefix + str(iyear) +
        #                                                    str(imonth+1).zfill(2) + PISCES24_ending))
        #         mercator_BIO_files = \
        #             mercator_BIO_files_pre + mercator_BIO_files_main + mercator_BIO_files_post
        #
        #     Tbry_str = "%06d" % Tstart
        #
        #     glorysname = mercator_PHY_files[0]
        #
        #     print(' OPEN : ' + glorysname)
        #
        #     NORESMnamebgc_ex = NORESMfiles_dir + NORESM_prefix + 'no3no2' + NORESM_ending
        #
        #     print(NORESMnamebgc_ex)
        #
        #     PISCESnamebgc = PISCES24files_dir + PISCES24_prefix + str(Ystart) + str(Mstart).zfill(2) + PISCES24_ending
        #
        #     print(PISCESnamebgc)
        #
        #     #
        #     # open the first GLORYS file
        #     #
        #     ncgl = netcdf(glorysname, 'r')
        #     ncglo = ncgl.variables
        #     depthg = np.array(ncglo['depth'][:])
        #
        #     ncpis = netcdf(PISCESnamebgc, 'r')
        #     ncpiso = ncpis.variables
        #     depthp = np.array(ncpiso['deptht'][:])
        #
        #     # Nitrate from NORESM
        #     ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
        #     ncnornio = ncnorni.variables
        #     depthn = np.array(ncnornio['depth'][:])
        #
        #     [Nz] = np.shape(depthg)
        #
        #     dl = 1
        #     #
        #
        #     if obc[0] == 1:
        #         #
        #         #  Southern boundary
        #         #
        #
        #         print(' ')
        #         print(' Southern Boundary')
        #         print(' ')
        #
        #         lon_south = lon_rho[0:2, :]
        #         lat_south = lat_rho[0:2, :]
        #         h_south = h[0:2, :]
        #         angle_south = angle[0:2, :]
        #
        #         (LonT_south, LatT_south, iminT_south, imaxT_south, jminT_south, jmaxT_south, elemT_south, coefT_south,
        #          LonU_south, LatU_south, iminU_south, imaxU_south, jminU_south, jmaxU_south, elemU_south, coefU_south,
        #          LonV_south, LatV_south, iminV_south, imaxV_south, jminV_south, jmaxV_south, elemV_south, coefV_south,
        #          LonN_south, LatN_south, iminN_south, imaxN_south, jminN_south, jmaxN_south, elemN_south, coefN_south,
        #          LonP_south, LatP_south, iminP_south, imaxP_south, jminP_south, jmaxP_south, elemP_south, coefP_south) \
        #             = glor.get_delaunay_bry_PISCES_NORESM(lon_south, lat_south, dl, ncglo, ncpiso, ncnornio)
        #
        #     if obc[1] == 1:
        #         #
        #         #  Eastern boundary
        #         #
        #
        #         print(' ')
        #         print(' Eastern Boundary')
        #         print(' ')
        #
        #         lon_east = lon_rho[:, -2:]
        #         lat_east = lat_rho[:, -2:]
        #         h_east = h[:, -2:]
        #         angle_east = angle[:, -2:]
        #
        #         (LonT_east, LatT_east, iminT_east, imaxT_east, jminT_east, jmaxT_east, elemT_east, coefT_east,
        #          LonU_east, LatU_east, iminU_east, imaxU_east, jminU_east, jmaxU_east, elemU_east, coefU_east,
        #          LonV_east, LatV_east, iminV_east, imaxV_east, jminV_east, jmaxV_east, elemV_east, coefV_east,
        #          LonN_east, LatN_east, iminN_east, imaxN_east, jminN_east, jmaxN_east, elemN_east, coefN_east,
        #          LonP_east, LatP_east, iminP_east, imaxP_east, jminP_east, jmaxP_east, elemP_east, coefP_east) \
        #             = glor.get_delaunay_bry_PISCES_NORESM(lon_east, lat_east, dl, ncglo, ncpiso, ncnornio)
        #
        #     if obc[2] == 1:
        #         #
        #         #  Northern boundary
        #         #
        #
        #         print(' ')
        #         print(' Northern Boundary')
        #         print(' ')
        #
        #         lon_north = lon_rho[-2:, :]
        #         lat_north = lat_rho[-2:, :]
        #         h_north = h[-2:, :]
        #         angle_north = angle[-2:, :]
        #
        #         (LonT_north, LatT_north, iminT_north, imaxT_north, jminT_north, jmaxT_north, elemT_north, coefT_north,
        #          LonU_north, LatU_north, iminU_north, imaxU_north, jminU_north, jmaxU_north, elemU_north, coefU_north,
        #          LonV_north, LatV_north, iminV_north, imaxV_north, jminV_north, jmaxV_north, elemV_north, coefV_north,
        #          LonN_north, LatN_north, iminN_north, imaxN_north, jminN_north, jmaxN_north, elemN_north, coefN_north,
        #          LonP_north, LatP_north, iminP_north, imaxP_north, jminP_north, jmaxP_north, elemP_north, coefP_north) \
        #             = glor.get_delaunay_bry_PISCES_NORESM(lon_north, lat_north, dl, ncglo, ncpiso, ncnornio)
        #
        #     if obc[3] == 1:
        #         #
        #         #  Western boundary
        #         #
        #
        #         print(' ')
        #         print(' Western Boundary')
        #         print(' ')
        #
        #         lon_west = lon_rho[:, 0:2]
        #         lat_west = lat_rho[:, 0:2]
        #         h_west = h[:, 0:2]
        #         angle_west = angle[:, 0:2]
        #
        #         (LonT_west, LatT_west, iminT_west, imaxT_west, jminT_west, jmaxT_west, elemT_west, coefT_west,
        #          LonU_west, LatU_west, iminU_west, imaxU_west, jminU_west, jmaxU_west, elemU_west, coefU_west,
        #          LonV_west, LatV_west, iminV_west, imaxV_west, jminV_west, jmaxV_west, elemV_west, coefV_west,
        #          LonN_west, LatN_west, iminN_west, imaxN_west, jminN_west, jmaxN_west, elemN_west, coefN_west,
        #          LonP_west, LatP_west, iminP_west, imaxP_west, jminP_west, jmaxP_west, elemP_west, coefP_west) \
        #             = glor.get_delaunay_bry_PISCES_NORESM(lon_west, lat_west, dl, ncglo, ncpiso, ncnornio)
        #
        #     #
        #     #  Close the GLORYS netcdf file
        #     #
        #
        #     ncgl.close()
        #     ncnorni.close()
        #     ncpis.close()
        #
        #     #
        #     # Get the GLORYS file name from the date (only one time step per file)
        #     #
        #     # open the first GLORYS file
        #     #
        #     # Do the interpolations for each boundary
        #     #
        #
        #     print(' ')
        #     print(' Do the interpolations for each boundary')
        #     print(' ')
        #
        #     tndx_bry = -1
        #
        #     #
        #     # Loop on files
        #     #
        #     Tbries = np.arange(Tstart, Tend + 1, glorys_step)
        #
        #     PISCES_NORESM_interpd = 1
        #
        #     for mpf in np.arange(0, Tend - Tstart + 1, glorys_step):
        #         Tbry = Tbries[int(mpf)]
        #         print(Tbry)
        #
        #         tndx_bry = tndx_bry + 1
        #         # PHYSICS GETS A DAILY TIMESTEP
        #         ncbry['bry_time'][tndx_bry] = Tbry
        #         ncbry['tclm_time'][tndx_bry] = Tbry
        #         ncbry['temp_time'][tndx_bry] = Tbry
        #         ncbry['sclm_time'][tndx_bry] = Tbry
        #         ncbry['salt_time'][tndx_bry] = Tbry
        #         ncbry['uclm_time'][tndx_bry] = Tbry
        #         ncbry['vclm_time'][tndx_bry] = Tbry
        #         ncbry['v2d_time'][tndx_bry] = Tbry
        #         ncbry['v3d_time'][tndx_bry] = Tbry
        #         ncbry['ssh_time'][tndx_bry] = Tbry
        #         ncbry['zeta_time'][tndx_bry] = Tbry
        #
        #         if PISCES_NORESM_interpd is 1:
        #             # BGC VARIABLES GETS MONTHLY TIMESTEP (3 timepoints: prior, present and next month)
        #             ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
        #             ncnornio = ncnorni.variables
        #             Tinin = datetime(iyear, imonth, 15)
        #             Tininxcur = d2i(Tinin, ncnornio['time'], select='exact')
        #             Tininxpri = Tininxcur - 1
        #             Tininxpos = Tininxcur + 1
        #             Tininxn = [Tininxpri, Tininxcur, Tininxpos]
        #
        #             if imonth == 1:
        #                 Tinin1 = date.toordinal(date(iyear - 1, 12, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #                 Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #                 Tinin3 = date.toordinal(date(iyear, imonth + 1, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #             elif imonth == 12:
        #                 Tinin1 = date.toordinal(date(iyear, imonth - 1, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #                 Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #                 Tinin3 = date.toordinal(date(iyear + 1, 1, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #             else:
        #                 Tinin1 = date.toordinal(date(iyear, imonth - 1, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #                 Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #                 Tinin3 = date.toordinal(date(iyear, imonth + 1, 15)) - \
        #                          date.toordinal(date(Yorig, 1, 1))
        #
        #             ncnorni.close()
        #
        #             ncpispri = netcdf(mercator_BIO_files[0], 'r')
        #             ncpisprio = ncpispri.variables
        #             Tinip1 = ncpisprio['time_centered'][:][0]/86400 + date.toordinal(date(1950, 1, 1)) - \
        #                      date.toordinal(date(Yorig, 1, 1))
        #             ncpiscur = netcdf(mercator_BIO_files[1], 'r')
        #             ncpiscuro = ncpiscur.variables
        #             Tinip2 = ncpiscuro['time_centered'][:][0]/86400 + date.toordinal(date(1950, 1, 1)) - \
        #                      date.toordinal(date(Yorig, 1, 1))
        #             ncpispos = netcdf(mercator_BIO_files[2], 'r')
        #             ncpisposo = ncpispos.variables
        #             Tinip3 = ncpisposo['time_centered'][:][0]/86400 + date.toordinal(date(1950, 1, 1)) - \
        #                      date.toordinal(date(Yorig, 1, 1))
        #
        #             # NORESM variables get their timestamps
        #             ncbry['no3_time'][0] = Tinin1
        #             ncbry['po4_time'][0] = Tinin1
        #             ncbry['si_time'][0] = Tinin1
        #             ncbry['o2_time'][0] = Tinin1
        #             ncbry['dic_time'][0] = Tinin1
        #             ncbry['alkalini_time'][0] = Tinin1
        #             ncbry['no3_time'][1] = Tinin2
        #             ncbry['po4_time'][1] = Tinin2
        #             ncbry['si_time'][1] = Tinin2
        #             ncbry['o2_time'][1] = Tinin2
        #             ncbry['dic_time'][1] = Tinin2
        #             ncbry['alkalini_time'][1] = Tinin2
        #             ncbry['no3_time'][2] = Tinin3
        #             ncbry['po4_time'][2] = Tinin3
        #             ncbry['si_time'][2] = Tinin3
        #             ncbry['o2_time'][2] = Tinin3
        #             ncbry['dic_time'][2] = Tinin3
        #             ncbry['alkalini_time'][2] = Tinin3
        #
        #             # PISCES variables get their timestamps
        #             ncbry['fer_time'][0] = Tinip1
        #             ncbry['nh4_time'][0] = Tinip1
        #             ncbry['caco3_time'][0] = Tinip1
        #             ncbry['poc_time'][0] = Tinip1
        #             ncbry['phy_time'][0] = Tinip1
        #             ncbry['zoo_time'][0] = Tinip1
        #             ncbry['doc_time'][0] = Tinip1
        #             ncbry['phy2_time'][0] = Tinip1
        #             ncbry['zoo2_time'][0] = Tinip1
        #             ncbry['bsi_time'][0] = Tinip1
        #             ncbry['bfe_time'][0] = Tinip1
        #             ncbry['goc_time'][0] = Tinip1
        #             ncbry['sfe_time'][0] = Tinip1
        #             ncbry['dfe_time'][0] = Tinip1
        #             ncbry['dsi_time'][0] = Tinip1
        #             ncbry['nfe_time'][0] = Tinip1
        #             ncbry['nchl_time'][0] = Tinip1
        #             ncbry['dchl_time'][0] = Tinip1
        #             ncbry['fer_time'][1] = Tinip2
        #             ncbry['nh4_time'][1] = Tinip2
        #             ncbry['caco3_time'][1] = Tinip2
        #             ncbry['poc_time'][1] = Tinip2
        #             ncbry['phy_time'][1] = Tinip2
        #             ncbry['zoo_time'][1] = Tinip2
        #             ncbry['doc_time'][1] = Tinip2
        #             ncbry['phy2_time'][1] = Tinip2
        #             ncbry['zoo2_time'][1] = Tinip2
        #             ncbry['bsi_time'][1] = Tinip2
        #             ncbry['bfe_time'][1] = Tinip2
        #             ncbry['goc_time'][1] = Tinip2
        #             ncbry['sfe_time'][1] = Tinip2
        #             ncbry['dfe_time'][1] = Tinip2
        #             ncbry['dsi_time'][1] = Tinip2
        #             ncbry['nfe_time'][1] = Tinip2
        #             ncbry['nchl_time'][1] = Tinip2
        #             ncbry['dchl_time'][1] = Tinip2
        #             ncbry['fer_time'][2] = Tinip3
        #             ncbry['nh4_time'][2] = Tinip3
        #             ncbry['caco3_time'][2] = Tinip3
        #             ncbry['poc_time'][2] = Tinip3
        #             ncbry['phy_time'][2] = Tinip3
        #             ncbry['zoo_time'][2] = Tinip3
        #             ncbry['doc_time'][2] = Tinip3
        #             ncbry['phy2_time'][2] = Tinip3
        #             ncbry['zoo2_time'][2] = Tinip3
        #             ncbry['bsi_time'][2] = Tinip3
        #             ncbry['bfe_time'][2] = Tinip3
        #             ncbry['goc_time'][2] = Tinip3
        #             ncbry['sfe_time'][2] = Tinip3
        #             ncbry['dfe_time'][2] = Tinip3
        #             ncbry['dsi_time'][2] = Tinip3
        #             ncbry['nfe_time'][2] = Tinip3
        #             ncbry['nchl_time'][2] = Tinip3
        #             ncbry['dchl_time'][2] = Tinip3
        #
        #         #
        #         # Get the first GLORYS file name (for the positions of variables)
        #         #
        #         Tbry_str = "%06d" % Tbry
        #
        #         glorysname = mercator_PHY_files[int(mpf)]
        #         print(' OPEN : ' + glorysname)
        #         ncgl = netcdf(glorysname, 'r')
        #         ncglo = ncgl.variables
        #
        #         if obc[0] == 1:
        #             #
        #             #  Southern boundary
        #             #
        #             print(' ')
        #             print(' Soutern Boundary')
        #             print(' ')
        #
        #             ncbry = glor.interp_bry_PISCES_NORESM('s', PISCES_NORESM_interpd, Tininxn,
        #                                                   ncglo, ncpisprio, ncpiscuro, ncpisposo,
        #                                                   NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
        #                                                   tndx_glo, ncbry, tndx_bry, h_south, theta_s, theta_b,
        #                                                   hc, N, vtransform, Nzgoodmin,
        #                                                   depthg, depthp, depthn, angle_south,
        #                                                   LonT_south, LatT_south, iminT_south, imaxT_south,
        #                                                   jminT_south, jmaxT_south, elemT_south, coefT_south,
        #                                                   LonU_south, LatU_south, iminU_south, imaxU_south,
        #                                                   jminU_south, jmaxU_south, elemU_south, coefU_south,
        #                                                   LonV_south, LatV_south, iminV_south, imaxV_south,
        #                                                   jminV_south, jmaxV_south, elemV_south, coefV_south,
        #                                                   LonN_south, LatN_south, iminN_south, imaxN_south,
        #                                                   jminN_south, jmaxN_south, elemN_south, coefN_south,
        #                                                   LonP_south, LatP_south, iminP_south, imaxP_south,
        #                                                   jminP_south, jmaxP_south, elemP_south, coefP_south)
        #
        #         if obc[1] == 1:
        #             #
        #             #  Eastern boundary
        #             #
        #             print(' ')
        #             print(' Eastern Boundary')
        #             print(' ')
        #
        #             ncbry = glor.interp_bry_PISCES_NORESM('e', PISCES_NORESM_interpd, Tininxn,
        #                                                   ncglo, ncpisprio, ncpiscuro, ncpisposo,
        #                                                   NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
        #                                                   tndx_glo, ncbry, tndx_bry, h_east, theta_s, theta_b,
        #                                                   hc, N, vtransform, Nzgoodmin,
        #                                                   depthg, depthp, depthn, angle_east,
        #                                                   LonT_east, LatT_east, iminT_east, imaxT_east,
        #                                                   jminT_east, jmaxT_east, elemT_east, coefT_east,
        #                                                   LonU_east, LatU_east, iminU_east, imaxU_east,
        #                                                   jminU_east, jmaxU_east, elemU_east, coefU_east,
        #                                                   LonV_east, LatV_east, iminV_east, imaxV_east,
        #                                                   jminV_east, jmaxV_east, elemV_east, coefV_east,
        #                                                   LonN_east, LatN_east, iminN_east, imaxN_east,
        #                                                   jminN_east, jmaxN_east, elemN_east, coefN_east,
        #                                                   LonP_east, LatP_east, iminP_east, imaxP_east,
        #                                                   jminP_east, jmaxP_east, elemP_east, coefP_east)
        #
        #         if obc[2] == 1:
        #             #
        #             #  Northern boundary
        #             #
        #             print(' ')
        #             print(' Northern Boundary')
        #             print(' ')
        #
        #             ncbry = glor.interp_bry_PISCES_NORESM('n', PISCES_NORESM_interpd, Tininxn,
        #                                                   ncglo, ncpisprio, ncpiscuro, ncpisposo,
        #                                                   NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
        #                                                   tndx_glo, ncbry, tndx_bry, h_north, theta_s, theta_b,
        #                                                   hc, N, vtransform, Nzgoodmin,
        #                                                   depthg, depthp, depthn, angle_north,
        #                                                   LonT_north, LatT_north, iminT_north, imaxT_north,
        #                                                   jminT_north, jmaxT_north, elemT_north, coefT_north,
        #                                                   LonU_north, LatU_north, iminU_north, imaxU_north,
        #                                                   jminU_north, jmaxU_north, elemU_north, coefU_north,
        #                                                   LonV_north, LatV_north, iminV_north, imaxV_north,
        #                                                   jminV_north, jmaxV_north, elemV_north, coefV_north,
        #                                                   LonN_north, LatN_north, iminN_north, imaxN_north,
        #                                                   jminN_north, jmaxN_north, elemN_north, coefN_north,
        #                                                   LonP_north, LatP_north, iminP_north, imaxP_north,
        #                                                   jminP_north, jmaxP_north, elemP_north, coefP_north)
        #
        #         if obc[3] == 1:
        #             #
        #             #  Western boundary
        #             #
        #             print(' ')
        #             print(' Western Boundary')
        #             print(' ')
        #
        #             ncbry = glor.interp_bry_PISCES_NORESM('w', PISCES_NORESM_interpd, Tininxn,
        #                                                   ncglo, ncpisprio, ncpiscuro, ncpisposo,
        #                                                   NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
        #                                                   tndx_glo, ncbry, tndx_bry, h_west, theta_s, theta_b,
        #                                                   hc, N, vtransform, Nzgoodmin,
        #                                                   depthg, depthp, depthn, angle_west,
        #                                                   LonT_west, LatT_west, iminT_west, imaxT_west,
        #                                                   jminT_west, jmaxT_west, elemT_west, coefT_west,
        #                                                   LonU_west, LatU_west, iminU_west, imaxU_west,
        #                                                   jminU_west, jmaxU_west, elemU_west, coefU_west,
        #                                                   LonV_west, LatV_west, iminV_west, imaxV_west,
        #                                                   jminV_west, jmaxV_west, elemV_west, coefV_west,
        #                                                   LonN_west, LatN_west, iminN_west, imaxN_west,
        #                                                   jminN_west, jmaxN_west, elemN_west, coefN_west,
        #                                                   LonP_west, LatP_west, iminP_west, imaxP_west,
        #                                                   jminP_west, jmaxP_west, elemP_west, coefP_west)
        #
        #         if PISCES_NORESM_interpd is 1:
        #             PISCES_NORESM_interpd = 0
        #             ncpispri.close()
        #             ncpiscur.close()
        #             ncpispos.close()
        #
        #         ncgl.close()
        #     #
        #     #  End loop on time
        #     #
        #     ncbr.close()
