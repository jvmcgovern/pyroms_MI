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
from calendar import monthrange
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

    crocofiles_bcd_dir = '/media/dskfour/CS1KM_19932021_BCd/'
    crocofiles_ubcd_dir = '/media/dskfour/CS1KM_19932021_Uncorr/'
    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    # bryname = crocofiles_dir + 'croco_bry_PHYBIO_CELTIC_h6.nc'
    # bryname = crocofiles_dir + 'croco_bry_PHYBIO_CELTIC_h8.nc'
    # bryname = crocofiles_dir + 'croco_bry_MERCATOR_.nc'
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    # grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m
    # grdname = crocofiles_dir + 'CELTIC_grd_v6_MSL_h8_v4.nc'  # LAT to MSL included, hmin = 8m
    grdname = crocofiles_dir + 'CELTIC_grd_v5_MSL_h10.nc'  # LAT to MSL included, INFOMAR block averaged hmin = 10m
    # grdname = crocofiles_dir + 'croco_grd.nc'

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

    Ystart = 2021
    Mstart = 12
    Dstart = 1

    Yend = 2021
    Mend = 12
    Dend = 31

    # CMEMS_IBI MY (Reanalysis data)
    glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # # Reanalysis
    # glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    # glorys_ending = '_R*_RE01.nc'
    # Near realtime
    glorys_prefix = 'CMEMS_v6r1_IBI_PHY_NRT_NL_01dav_'
    glorys_ending = '_R*_AN*.nc'

    # CMEMS_IBI MY (monthly average for TALK derivation from temp/sal/DIC/pH)
    # 'CELTIC/CMEMS_IBI/CMEMS_v5r1_IBI_PHY_MY_PdE_01mav_20051201_20051231_R20201201_RE01.nc'
    # # Reanalysis
    # glorys_mth_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01mav_'
    # Near realtime
    glorys_mth_prefix = 'CMEMS_v6r1_IBI_PHY_NRT_NL_01mav_'

    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    # PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI_PISCESvars_012013_202112/monthly/'
    # PISCES24_prefix = 'IBI36_cg_1m-m_'
    # PISCES24_ending = '_3DT-bgc_hcst.nc'

    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # # Reanalysis
    # PISCES24_prefix = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_'
    # PISCES24_ending = '_R*_RE01.nc'
    # Near realtime
    PISCES24_prefix = 'CMEMS_v7r1_IBI_BIO_NRT_NL_01mav_'
    PISCES24_ending = '_R*_AN*.nc'

    # 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_20190201_20190228_R20201201_RE01.nc'

    glorys_step = 1  # time step between outputs in GLORYS12 [days]

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    #
    # END USERS DEFINED VARIABLES
    #

    for iyear in range(Ystart, Yend + 1):
        for imonth in range(Mstart, Mend + 1):

            # bryname = crocofiles_bcd_dir + 'croco_bry_MERCATOR_hc50m_MSL_ng' + '_Y' + str(iyear) + 'M' \
            #           + str(imonth).zfill(2) + '.nc'

            bryname = crocofiles_bcd_dir + 'croco_bry_nrt_bc_MERCATOR_hc50m_MSL_ng' + '_Y' + str(iyear) + 'M' \
                      + str(imonth).zfill(2) + '.nc'

            # bryname = crocofiles_bcd_dir + 'croco_bry_nrt_MERCATOR_hc50m_MSL_ng' + '_Y' + str(iyear) + 'M' \
            #           + str(imonth).zfill(2) + '.nc'

            # bry_bc_name = crocofiles_bc_dir + 'croco_bry_bias2remove_MERCATOR_hc50m_MSL_ng_Y1993M' \
            #               + str(imonth).zfill(2) + '.nc'

            bry_ubc_name = crocofiles_ubcd_dir + 'croco_bry_nrt_MERCATOR_hc50m_MSL_ng' + '_Y' + str(iyear) + 'M' \
                           + str(imonth).zfill(2) + '.nc'

            print(' ')
            print(' Correcting boundary file: ' + bryname)
            print(' ')

            #
            # Open the CROCO boundary file for processing
            #
            # Uncorrected
            ncbr = netcdf(bryname, 'r+')
            ncbry = ncbr.variables

            # Bias-correction data to remove from uncorrected data
            ncbrbc = netcdf(bry_ubc_name, 'r+')
            ncbrbcy = ncbrbc.variables

            # NO3_south_BCd = ncbry['NO3_south'][:] - ncbrbcy['NO3_south'][:]
            # NO3_south_BCd[NO3_south_BCd < 0] = 0
            #
            # NO3_north_BCd = ncbry['NO3_north'][:] - ncbrbcy['NO3_north'][:]
            # NO3_north_BCd[NO3_north_BCd < 0] = 0
            #
            # NO3_east_BCd = ncbry['NO3_east'][:] - ncbrbcy['NO3_east'][:]
            # NO3_east_BCd[NO3_east_BCd < 0] = 0
            #
            # NO3_west_BCd = ncbry['NO3_west'][:] - ncbrbcy['NO3_west'][:]
            # NO3_west_BCd[NO3_west_BCd < 0] = 0
            #
            # ncbry['NO3_south'][:] = NO3_south_BCd
            # ncbry['NO3_north'][:] = NO3_north_BCd
            # ncbry['NO3_east'][:] = NO3_east_BCd
            # ncbry['NO3_west'][:] = NO3_west_BCd
            #
            # PO4_south_BCd = ncbry['PO4_south'][:] - ncbrbcy['PO4_south'][:]
            # PO4_south_BCd[PO4_south_BCd < 0] = 0
            # PO4_north_BCd = ncbry['PO4_north'][:] - ncbrbcy['PO4_north'][:]
            # PO4_north_BCd[PO4_north_BCd < 0] = 0
            # PO4_east_BCd = ncbry['PO4_east'][:] - ncbrbcy['PO4_east'][:]
            # PO4_east_BCd[PO4_east_BCd < 0] = 0
            # PO4_west_BCd = ncbry['PO4_west'][:] - ncbrbcy['PO4_west'][:]
            # PO4_west_BCd[PO4_west_BCd < 0] = 0
            #
            # ncbry['PO4_south'][:] = PO4_south_BCd
            # ncbry['PO4_north'][:] = PO4_north_BCd
            # ncbry['PO4_east'][:] = PO4_east_BCd
            # ncbry['PO4_west'][:] = PO4_west_BCd
            #
            # Si_south_BCd = ncbry['Si_south'][:] - ncbrbcy['Si_south'][:]
            # Si_south_BCd[Si_south_BCd < 0] = 0
            # Si_north_BCd = ncbry['Si_north'][:] - ncbrbcy['Si_north'][:]
            # Si_north_BCd[Si_north_BCd < 0] = 0
            # Si_east_BCd = ncbry['Si_east'][:] - ncbrbcy['Si_east'][:]
            # Si_east_BCd[Si_east_BCd < 0] = 0
            # Si_west_BCd = ncbry['Si_west'][:] - ncbrbcy['Si_west'][:]
            # Si_west_BCd[Si_west_BCd < 0] = 0
            #
            # ncbry['Si_south'][:] = Si_south_BCd
            # ncbry['Si_north'][:] = Si_north_BCd
            # ncbry['Si_east'][:] = Si_east_BCd
            # ncbry['Si_west'][:] = Si_west_BCd

            # O2_south_BCd = ncbry['O2_south'][:] - ncbrbcy['O2_south'][:]
            # O2_south_BCd[O2_south_BCd < 0] = 0
            # O2_north_BCd = ncbry['O2_north'][:] - ncbrbcy['O2_north'][:]
            # O2_north_BCd[O2_north_BCd < 0] = 0
            # O2_east_BCd = ncbry['O2_east'][:] - ncbrbcy['O2_east'][:]
            # O2_east_BCd[O2_east_BCd < 0] = 0
            # O2_west_BCd = ncbry['O2_west'][:] - ncbrbcy['O2_west'][:]
            # O2_west_BCd[O2_west_BCd < 0] = 0
            #
            # ncbry['O2_south'][:] = O2_south_BCd
            # ncbry['O2_north'][:] = O2_north_BCd
            # ncbry['O2_east'][:] = O2_east_BCd
            # ncbry['O2_west'][:] = O2_west_BCd

            ncbry['O2_south'][:] = ncbrbcy['O2_south'][:]
            ncbry['O2_north'][:] = ncbrbcy['O2_north'][:]
            ncbry['O2_east'][:] = ncbrbcy['O2_east'][:]
            ncbry['O2_west'][:] = ncbrbcy['O2_west'][:]

            # DIC_south_BCd = ncbry['DIC_south'][:] - ncbrbcy['DIC_south'][:]
            # DIC_south_BCd[DIC_south_BCd < 0] = 0
            # DIC_north_BCd = ncbry['DIC_north'][:] - ncbrbcy['DIC_north'][:]
            # DIC_north_BCd[DIC_north_BCd < 0] = 0
            # DIC_east_BCd = ncbry['DIC_east'][:] - ncbrbcy['DIC_east'][:]
            # DIC_east_BCd[DIC_east_BCd < 0] = 0
            # DIC_west_BCd = ncbry['DIC_west'][:] - ncbrbcy['DIC_west'][:]
            # DIC_west_BCd[DIC_west_BCd < 0] = 0
            #
            # ncbry['DIC_south'][:] = DIC_south_BCd
            # ncbry['DIC_north'][:] = DIC_north_BCd
            # ncbry['DIC_east'][:] = DIC_east_BCd
            # ncbry['DIC_west'][:] = DIC_west_BCd
            #
            # TALK_south_BCd = ncbry['TALK_south'][:] - ncbrbcy['TALK_south'][:]
            # TALK_south_BCd[TALK_south_BCd < 0] = 0
            # TALK_north_BCd = ncbry['TALK_north'][:] - ncbrbcy['TALK_north'][:]
            # TALK_north_BCd[TALK_north_BCd < 0] = 0
            # TALK_east_BCd = ncbry['TALK_east'][:] - ncbrbcy['TALK_east'][:]
            # TALK_east_BCd[TALK_east_BCd < 0] = 0
            # TALK_west_BCd = ncbry['TALK_west'][:] - ncbrbcy['TALK_west'][:]
            # TALK_west_BCd[TALK_west_BCd < 0] = 0
            #
            # ncbry['TALK_south'][:] = TALK_south_BCd
            # ncbry['TALK_north'][:] = TALK_north_BCd
            # ncbry['TALK_east'][:] = TALK_east_BCd
            # ncbry['TALK_west'][:] = TALK_west_BCd

            ncbr.close()
            ncbrbc.close()
