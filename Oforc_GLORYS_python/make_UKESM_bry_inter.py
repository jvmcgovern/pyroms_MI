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

    # # Multiyear/reanalysis
    # Ystart = 1993
    # Mstart = 1
    # Dstart = 1
    #
    # Yend = 2021
    # Mend = 12
    # Dend = 31

    # # Near realtime option
    # Ystart = 2022
    # Mstart = 1
    # Dstart = 1
    #
    # Yend = 2022
    # Mend = 1
    # Dend = 31
    #
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

    vec_post = '_0fillval_unrotated_fixedgrid_IE_1degreg.nc'
    phytrac_post = '_IE_1degreg.nc'
    bgctrac_post = '_II_0p5.nc'

    CMIPfold = '/media/dskthree/UKESM1-0-LL/NATIVE/'

    PHYfold = 'PHY/'
    BGCfold = 'BGC/'

    phytracvars = ['thetao', 'so']
    bgctracvars = ['no3', 'po4', 'si', 'o2', 'dfe', 'dissic', 'talk']
    ele = ['zos']
    vec_vars = ['uo', 'vo']

    hist = 'historical'
    proj = 'ssp245'
    # proj = 'ssp585'
    exp = 'r1i1p1f2_gn'
    freq = 'Omon'
    model = 'UKESM1-0-LL'

    # Need to list all netcdfs for each of the variables to interpolate.
    # For each file, it will be possible to know what years are covered by it by identifying the year in file name
    # For each year and month, it will be necessary to identify the index of the pre-, present- and post- month in their
    # respective files and pass the files to the interpolation process.
    # Need to firstly sample one of each files to explore that lat/lon data and do the grid weighting

    # U velocities
    ufiles_hist = sorted(glob.glob(CMIPfold + PHYfold + vec_vars[0] + '_' + freq + '_' + model + '_' + hist +
                                   '_' + exp + '_??????-??????' + vec_post))
    ufiles_proj = sorted(glob.glob(CMIPfold + PHYfold + vec_vars[0] + '_' + freq + '_' + model + '_' + proj +
                                   '_' + exp + '_??????-??????' + vec_post))
    ufiles = ufiles_hist + ufiles_proj
    ufiles_yrs = np.ones((len(ufiles), 2), dtype=int)
    for file in range(0, len(ufiles)):
        ufiles_yrs[file, 0] = int(ufiles[file][-56:-52])
        ufiles_yrs[file, 1] = int(ufiles[file][-49:-45])

    # V velocities
    vfiles_hist = sorted(glob.glob(CMIPfold + PHYfold + vec_vars[1] + '_' + freq + '_' + model + '_' + hist +
                                   '_' + exp + '_??????-??????' + vec_post))
    vfiles_proj = sorted(glob.glob(CMIPfold + PHYfold + vec_vars[1] + '_' + freq + '_' + model + '_' + proj +
                                   '_' + exp + '_??????-??????' + vec_post))
    vfiles = vfiles_hist + vfiles_proj
    vfiles_yrs = np.ones((len(vfiles), 2), dtype=int)
    for file in range(0, len(vfiles)):
        vfiles_yrs[file, 0] = int(vfiles[file][-56:-52])
        vfiles_yrs[file, 1] = int(vfiles[file][-49:-45])

    # SSH
    zfiles_hist = sorted(glob.glob(CMIPfold + PHYfold + ele[0] + '_' + freq + '_' + model + '_' + hist +
                                   '_' + exp + '_??????-??????' + phytrac_post))
    zfiles_proj = sorted(glob.glob(CMIPfold + PHYfold + ele[0] + '_' + freq + '_' + model + '_' + proj +
                                   '_' + exp + '_??????-??????' + phytrac_post))
    zfiles = zfiles_hist + zfiles_proj
    zfiles_yrs = np.ones((len(zfiles), 2), dtype=int)
    for file in range(0, len(zfiles)):
        zfiles_yrs[file, 0] = int(zfiles[file][-27:-23])
        zfiles_yrs[file, 1] = int(zfiles[file][-20:-16])

    # Temperature
    tempfiles_hist = sorted(glob.glob(CMIPfold + PHYfold + phytracvars[0] + '_' + freq + '_' + model + '_' + hist +
                                      '_' + exp + '_??????-??????' + phytrac_post))
    tempfiles_proj = sorted(glob.glob(CMIPfold + PHYfold + phytracvars[0] + '_' + freq + '_' + model + '_' + proj +
                                      '_' + exp + '_??????-??????' + phytrac_post))
    tempfiles = tempfiles_hist + tempfiles_proj
    tempfiles_yrs = np.ones((len(tempfiles), 2), dtype=int)
    for file in range(0, len(tempfiles)):
        tempfiles_yrs[file, 0] = int(tempfiles[file][-27:-23])
        tempfiles_yrs[file, 1] = int(tempfiles[file][-20:-16])

    # Salinity
    saltfiles_hist = sorted(glob.glob(CMIPfold + PHYfold + phytracvars[1] + '_' + freq + '_' + model + '_' + hist +
                                      '_' + exp + '_??????-??????' + phytrac_post))
    saltfiles_proj = sorted(glob.glob(CMIPfold + PHYfold + phytracvars[1] + '_' + freq + '_' + model + '_' + proj +
                                      '_' + exp + '_??????-??????' + phytrac_post))
    saltfiles = saltfiles_hist + saltfiles_proj
    saltfiles_yrs = np.ones((len(saltfiles), 2), dtype=int)
    for file in range(0, len(saltfiles)):
        saltfiles_yrs[file, 0] = int(saltfiles[file][-27:-23])
        saltfiles_yrs[file, 1] = int(saltfiles[file][-20:-16])

    # Nitrate
    nitfiles_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[0] + '_' + freq + '_' + model + '_' + hist +
                                     '_' + exp + '_??????-??????' + bgctrac_post))
    nitfiles_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[0] + '_' + freq + '_' + model + '_' + proj +
                                     '_' + exp + '_??????-??????' + bgctrac_post))
    nitfiles = nitfiles_hist + nitfiles_proj
    nitfiles_yrs = np.ones((len(nitfiles), 2), dtype=int)
    for file in range(0, len(nitfiles)):
        nitfiles_yrs[file, 0] = int(nitfiles[file][-23:-19])
        nitfiles_yrs[file, 1] = int(nitfiles[file][-16:-12])

    # Phosphate
    po4files_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[1] + '_' + freq + '_' + model + '_' + hist +
                                     '_' + exp + '_??????-??????' + bgctrac_post))
    po4files_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[1] + '_' + freq + '_' + model + '_' + proj +
                                     '_' + exp + '_??????-??????' + bgctrac_post))
    po4files = po4files_hist + po4files_proj
    po4files_yrs = np.ones((len(po4files), 2), dtype=int)
    for file in range(0, len(po4files)):
        po4files_yrs[file, 0] = int(po4files[file][-23:-19])
        po4files_yrs[file, 1] = int(po4files[file][-16:-12])

    # Silicate
    sifiles_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[2] + '_' + freq + '_' + model + '_' + hist +
                                    '_' + exp + '_??????-??????' + bgctrac_post))
    sifiles_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[2] + '_' + freq + '_' + model + '_' + proj +
                                    '_' + exp + '_??????-??????' + bgctrac_post))
    sifiles = sifiles_hist + sifiles_proj
    sifiles_yrs = np.ones((len(sifiles), 2), dtype=int)
    for file in range(0, len(sifiles)):
        sifiles_yrs[file, 0] = int(sifiles[file][-23:-19])
        sifiles_yrs[file, 1] = int(sifiles[file][-16:-12])

    # Oxygen
    o2files_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[3] + '_' + freq + '_' + model + '_' + hist +
                                    '_' + exp + '_??????-??????' + bgctrac_post))
    o2files_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[3] + '_' + freq + '_' + model + '_' + proj +
                                    '_' + exp + '_??????-??????' + bgctrac_post))
    o2files = o2files_hist + o2files_proj
    o2files_yrs = np.ones((len(o2files), 2), dtype=int)
    for file in range(0, len(o2files)):
        o2files_yrs[file, 0] = int(o2files[file][-23:-19])
        o2files_yrs[file, 1] = int(o2files[file][-16:-12])

    # Iron
    fefiles_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[4] + '_' + freq + '_' + model + '_' + hist +
                                    '_' + exp + '_??????-??????' + bgctrac_post))
    fefiles_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[4] + '_' + freq + '_' + model + '_' + proj +
                                    '_' + exp + '_??????-??????' + bgctrac_post))
    fefiles = fefiles_hist + fefiles_proj
    fefiles_yrs = np.ones((len(fefiles), 2), dtype=int)
    for file in range(0, len(fefiles)):
        fefiles_yrs[file, 0] = int(fefiles[file][-23:-19])
        fefiles_yrs[file, 1] = int(fefiles[file][-16:-12])

    # DIC
    dicfiles_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[5] + '_' + freq + '_' + model + '_' + hist +
                                     '_' + exp + '_??????-??????' + bgctrac_post))
    dicfiles_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[5] + '_' + freq + '_' + model + '_' + proj +
                                     '_' + exp + '_??????-??????' + bgctrac_post))
    dicfiles = dicfiles_hist + dicfiles_proj
    dicfiles_yrs = np.ones((len(dicfiles), 2), dtype=int)
    for file in range(0, len(dicfiles)):
        dicfiles_yrs[file, 0] = int(dicfiles[file][-23:-19])
        dicfiles_yrs[file, 1] = int(dicfiles[file][-16:-12])

    # TALK
    talkfiles_hist = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[6] + '_' + freq + '_' + model + '_' + hist +
                                      '_' + exp + '_??????-??????' + bgctrac_post))
    talkfiles_proj = sorted(glob.glob(CMIPfold + BGCfold + bgctracvars[6] + '_' + freq + '_' + model + '_' + proj +
                                      '_' + exp + '_??????-??????' + bgctrac_post))
    talkfiles = talkfiles_hist + talkfiles_proj
    talkfiles_yrs = np.ones((len(talkfiles), 2), dtype=int)
    for file in range(0, len(talkfiles)):
        talkfiles_yrs[file, 0] = int(talkfiles[file][-23:-19])
        talkfiles_yrs[file, 1] = int(talkfiles[file][-16:-12])

    #
    # get the CROCO grid
    #
    ncg = netcdf(grdname, 'r')
    ncgrd = ncg.variables
    lon_rho = np.array(ncgrd['lon_rho'][:])
    lat_rho = np.array(ncgrd['lat_rho'][:])
    lon_u = np.array(ncgrd['lon_u'][:])
    lat_u = np.array(ncgrd['lat_u'][:])
    lon_v = np.array(ncgrd['lon_v'][:])
    lat_v = np.array(ncgrd['lat_v'][:])
    h = np.array(ncgrd['h'][:])
    mask = np.array(ncgrd['mask_rho'][:])
    angle = np.array(ncgrd['angle'][:])
    [M, L] = np.shape(lon_rho)
    ncg.close()

    ncll = netcdf(ufiles[0], 'r')
    ncllo = ncll.variables
    lon_esm = np.array(ncllo['lon'][:])
    lat_esm = np.array(ncllo['lat'][:])
    ncll.close()

    udfiles = sorted(glob.glob(CMIPfold + PHYfold + vec_vars[0] + '_' + freq + '_' + model + '_' + proj +
                               '_' + exp + '_??????-??????.nc'))

    ncd = netcdf(udfiles[0], 'r')
    ncdo = ncd.variables
    depth_esm = np.array(ncdo['lev'][:])
    ncd.close()

    #
    # Get the Delaunay triangulations for each boundary
    #

    print(' ')
    print(' Get the Delaunay triangulations for each boundary')
    print(' ')

    dl = 2
    #

    if obc[0] == 1:
        #
        #  Southern boundary
        #

        print(' ')
        print(' Southern Boundary')
        print(' ')

        lon_south = lon_rho[0:2, :]
        lat_south = lat_rho[0:2, :]
        h_south = h[0:2, :]
        angle_south = angle[0:2, :]

        (LonT_south, LatT_south, iminT_south, imaxT_south, jminT_south, jmaxT_south, elemT_south, coefT_south) \
            = glor.get_delaunay_bry_ESM(lon_south, lat_south, dl, lon_esm, lat_esm)

    if obc[1] == 1:
        #
        #  Eastern boundary
        #

        print(' ')
        print(' Eastern Boundary')
        print(' ')

        lon_east = lon_rho[:, -2:]
        lat_east = lat_rho[:, -2:]
        h_east = h[:, -2:]
        angle_east = angle[:, -2:]

        (LonT_east, LatT_east, iminT_east, imaxT_east, jminT_east, jmaxT_east, elemT_east, coefT_east) \
            = glor.get_delaunay_bry_ESM(lon_east, lat_east, dl, lon_esm, lat_esm)

    if obc[2] == 1:
        #
        #  Northern boundary
        #

        print(' ')
        print(' Northern Boundary')
        print(' ')

        lon_north = lon_rho[-2:, :]
        lat_north = lat_rho[-2:, :]
        h_north = h[-2:, :]
        angle_north = angle[-2:, :]

        (LonT_north, LatT_north, iminT_north, imaxT_north, jminT_north, jmaxT_north, elemT_north, coefT_north) \
            = glor.get_delaunay_bry_ESM(lon_north, lat_north, dl, lon_esm, lat_esm)

    if obc[3] == 1:
        #
        #  Western boundary
        #

        print(' ')
        print(' Western Boundary')
        print(' ')

        lon_west = lon_rho[:, 0:2]
        lat_west = lat_rho[:, 0:2]
        h_west = h[:, 0:2]
        angle_west = angle[:, 0:2]

        (LonT_west, LatT_west, iminT_west, imaxT_west, jminT_west, jmaxT_west, elemT_west, coefT_west) \
            = glor.get_delaunay_bry_ESM(lon_west, lat_west, dl, lon_esm, lat_esm)

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

            bryname = bryfiles_dir + 'croco_bry_' + model + '_' + proj + '_hc50m_MSL_ng' + \
                      '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'

            #
            # Get the first GLORYS file name (for the positions of variables)
            #

            if imonth is 1:
                zfiles_3m = [np.argwhere((zfiles_yrs[:, 0] <= iyear-1) & (zfiles_yrs[:, 1] >= iyear-1))[0][0],
                             np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0]]
                tempfiles_3m = [np.argwhere((tempfiles_yrs[:, 0] <= iyear-1) & (tempfiles_yrs[:, 1] >= iyear-1))[0][0],
                                np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0]]
                saltfiles_3m = [np.argwhere((saltfiles_yrs[:, 0] <= iyear-1) & (saltfiles_yrs[:, 1] >= iyear-1))[0][0],
                                np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0]]
                ufiles_3m = [np.argwhere((ufiles_yrs[:, 0] <= iyear-1) & (ufiles_yrs[:, 1] >= iyear-1))[0][0],
                             np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0]]
                vfiles_3m = [np.argwhere((vfiles_yrs[:, 0] <= iyear-1) & (vfiles_yrs[:, 1] >= iyear-1))[0][0],
                             np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0]]
                nitfiles_3m = [np.argwhere((nitfiles_yrs[:, 0] <= iyear-1) & (nitfiles_yrs[:, 1] >= iyear-1))[0][0],
                               np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0]]
                po4files_3m = [np.argwhere((po4files_yrs[:, 0] <= iyear-1) & (po4files_yrs[:, 1] >= iyear-1))[0][0],
                               np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0]]
                sifiles_3m = [np.argwhere((sifiles_yrs[:, 0] <= iyear-1) & (sifiles_yrs[:, 1] >= iyear-1))[0][0],
                              np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0]]
                o2files_3m = [np.argwhere((o2files_yrs[:, 0] <= iyear-1) & (o2files_yrs[:, 1] >= iyear-1))[0][0],
                              np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0]]
                fefiles_3m = [np.argwhere((fefiles_yrs[:, 0] <= iyear-1) & (fefiles_yrs[:, 1] >= iyear-1))[0][0],
                              np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0]]
                dicfiles_3m = [np.argwhere((dicfiles_yrs[:, 0] <= iyear-1) & (dicfiles_yrs[:, 1] >= iyear-1))[0][0],
                               np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0]]
                talkfiles_3m = [np.argwhere((talkfiles_yrs[:, 0] <= iyear-1) & (talkfiles_yrs[:, 1] >= iyear-1))[0][0],
                                np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0]]
            elif imonth == 12:
                zfiles_3m = [np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((zfiles_yrs[:, 0] <= iyear+1) & (zfiles_yrs[:, 1] >= iyear+1))[0][0]]
                tempfiles_3m = [np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((tempfiles_yrs[:, 0] <= iyear+1) & (tempfiles_yrs[:, 1] >= iyear+1))[0][0]]
                saltfiles_3m = [np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((saltfiles_yrs[:, 0] <= iyear+1) & (saltfiles_yrs[:, 1] >= iyear+1))[0][0]]
                ufiles_3m = [np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((ufiles_yrs[:, 0] <= iyear+1) & (ufiles_yrs[:, 1] >= iyear+1))[0][0]]
                vfiles_3m = [np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((vfiles_yrs[:, 0] <= iyear+1) & (vfiles_yrs[:, 1] >= iyear+1))[0][0]]
                nitfiles_3m = [np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((nitfiles_yrs[:, 0] <= iyear+1) & (nitfiles_yrs[:, 1] >= iyear+1))[0][0]]
                po4files_3m = [np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((po4files_yrs[:, 0] <= iyear+1) & (po4files_yrs[:, 1] >= iyear+1))[0][0]]
                sifiles_3m = [np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((sifiles_yrs[:, 0] <= iyear+1) & (sifiles_yrs[:, 1] >= iyear+1))[0][0]]
                o2files_3m = [np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((o2files_yrs[:, 0] <= iyear+1) & (o2files_yrs[:, 1] >= iyear+1))[0][0]]
                fefiles_3m = [np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((fefiles_yrs[:, 0] <= iyear+1) & (fefiles_yrs[:, 1] >= iyear+1))[0][0]]
                dicfiles_3m = [np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((dicfiles_yrs[:, 0] <= iyear+1) & (dicfiles_yrs[:, 1] >= iyear+1))[0][0]]
                talkfiles_3m = [np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((talkfiles_yrs[:, 0] <= iyear+1) & (talkfiles_yrs[:, 1] >= iyear+1))[0][0]]
            else:
                zfiles_3m = [np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((zfiles_yrs[:, 0] <= iyear) & (zfiles_yrs[:, 1] >= iyear))[0][0]]
                tempfiles_3m = [np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((tempfiles_yrs[:, 0] <= iyear) & (tempfiles_yrs[:, 1] >= iyear))[0][0]]
                saltfiles_3m = [np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((saltfiles_yrs[:, 0] <= iyear) & (saltfiles_yrs[:, 1] >= iyear))[0][0]]
                ufiles_3m = [np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((ufiles_yrs[:, 0] <= iyear) & (ufiles_yrs[:, 1] >= iyear))[0][0]]
                vfiles_3m = [np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0],
                             np.argwhere((vfiles_yrs[:, 0] <= iyear) & (vfiles_yrs[:, 1] >= iyear))[0][0]]
                nitfiles_3m = [np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((nitfiles_yrs[:, 0] <= iyear) & (nitfiles_yrs[:, 1] >= iyear))[0][0]]
                po4files_3m = [np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((po4files_yrs[:, 0] <= iyear) & (po4files_yrs[:, 1] >= iyear))[0][0]]
                sifiles_3m = [np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((sifiles_yrs[:, 0] <= iyear) & (sifiles_yrs[:, 1] >= iyear))[0][0]]
                o2files_3m = [np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((o2files_yrs[:, 0] <= iyear) & (o2files_yrs[:, 1] >= iyear))[0][0]]
                fefiles_3m = [np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0],
                              np.argwhere((fefiles_yrs[:, 0] <= iyear) & (fefiles_yrs[:, 1] >= iyear))[0][0]]
                dicfiles_3m = [np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0],
                               np.argwhere((dicfiles_yrs[:, 0] <= iyear) & (dicfiles_yrs[:, 1] >= iyear))[0][0]]
                talkfiles_3m = [np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0],
                                np.argwhere((talkfiles_yrs[:, 0] <= iyear) & (talkfiles_yrs[:, 1] >= iyear))[0][0]]

            print(' ')
            print(' Making boundary file: ' + bryname)
            print(' ')
            print(' Title: ' + title)

            #
            # Create the CROCO boundary file
            #
            glor.create_bryfile_CMIP6(bryname, grdname, title, obc,
                                      theta_s, theta_b, hc, N,
                                      time_bry, 0, 0, cycle_bry, vtransform)

            #
            # Open the CROCO boundary file for writing
            #

            ncbr = netcdf(bryname, 'a')
            ncbry = ncbr.variables

            #
            # Get the GLORYS file name from the date (only one time step per file)
            #
            # open the first GLORYS file
            #
            # Do the interpolations for each boundary
            #

            print(' ')
            print(' Do the interpolations for each boundary')
            print(' ')

            if imonth == 1:
                Tinin1 = date.toordinal(date(iyear - 1, 12, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
                Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
                Tinin3 = date.toordinal(date(iyear, imonth + 1, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
            elif imonth == 12:
                Tinin1 = date.toordinal(date(iyear, imonth - 1, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
                Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
                Tinin3 = date.toordinal(date(iyear + 1, 1, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
            else:
                Tinin1 = date.toordinal(date(iyear, imonth - 1, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
                Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))
                Tinin3 = date.toordinal(date(iyear, imonth + 1, 15)) - \
                         date.toordinal(date(Yorig, 1, 1))

            # PHYSICS
            ncbry['bry_time'][0] = np.float64(Tinin1)
            ncbry['temp_time'][0] = np.float64(Tinin1)
            ncbry['salt_time'][0] = np.float64(Tinin1)
            ncbry['uclm_time'][0] = np.float64(Tinin1)
            ncbry['vclm_time'][0] = np.float64(Tinin1)
            ncbry['v2d_time'][0] = np.float64(Tinin1)
            ncbry['v3d_time'][0] = np.float64(Tinin1)
            ncbry['ssh_time'][0] = np.float64(Tinin1)
            ncbry['zeta_time'][0] = np.float64(Tinin1)

            ncbry['bry_time'][1] = np.float64(Tinin2)
            ncbry['temp_time'][1] = np.float64(Tinin2)
            ncbry['salt_time'][1] = np.float64(Tinin2)
            ncbry['uclm_time'][1] = np.float64(Tinin2)
            ncbry['vclm_time'][1] = np.float64(Tinin2)
            ncbry['v2d_time'][1] = np.float64(Tinin2)
            ncbry['v3d_time'][1] = np.float64(Tinin2)
            ncbry['ssh_time'][1] = np.float64(Tinin2)
            ncbry['zeta_time'][1] = np.float64(Tinin2)

            ncbry['bry_time'][2] = np.float64(Tinin3)
            ncbry['temp_time'][2] = np.float64(Tinin3)
            ncbry['salt_time'][2] = np.float64(Tinin3)
            ncbry['uclm_time'][2] = np.float64(Tinin3)
            ncbry['vclm_time'][2] = np.float64(Tinin3)
            ncbry['v2d_time'][2] = np.float64(Tinin3)
            ncbry['v3d_time'][2] = np.float64(Tinin3)
            ncbry['ssh_time'][2] = np.float64(Tinin3)
            ncbry['zeta_time'][2] = np.float64(Tinin3)

            # BGC variables get their timestamps
            ncbry['no3_time'][0] = np.float64(Tinin1)
            ncbry['po4_time'][0] = np.float64(Tinin1)
            ncbry['si_time'][0] = np.float64(Tinin1)
            ncbry['o2_time'][0] = np.float64(Tinin1)
            ncbry['fer_time'][0] = np.float64(Tinin1)
            ncbry['dic_time'][0] = np.float64(Tinin1)
            ncbry['talk_time'][0] = np.float64(Tinin1)
            ncbry['no3_time'][1] = np.float64(Tinin2)
            ncbry['po4_time'][1] = np.float64(Tinin2)
            ncbry['si_time'][1] = np.float64(Tinin2)
            ncbry['o2_time'][1] = np.float64(Tinin2)
            ncbry['fer_time'][1] = np.float64(Tinin2)
            ncbry['dic_time'][1] = np.float64(Tinin2)
            ncbry['talk_time'][1] = np.float64(Tinin2)
            ncbry['no3_time'][2] = np.float64(Tinin3)
            ncbry['po4_time'][2] = np.float64(Tinin3)
            ncbry['si_time'][2] = np.float64(Tinin3)
            ncbry['o2_time'][2] = np.float64(Tinin3)
            ncbry['fer_time'][2] = np.float64(Tinin3)
            ncbry['dic_time'][2] = np.float64(Tinin3)
            ncbry['talk_time'][2] = np.float64(Tinin3)

            if obc[0] == 1:
                #
                #  Southern boundary
                #
                print(' ')
                print(' Southern Boundary')
                print(' ')

                ncbry = glor.interp_bry_ESM('s',
                                            zfiles_3m, tempfiles_3m, saltfiles_3m, ufiles_3m, vfiles_3m,
                                            nitfiles_3m, po4files_3m, sifiles_3m, o2files_3m, fefiles_3m,
                                            dicfiles_3m, talkfiles_3m,
                                            zfiles, tempfiles, saltfiles, ufiles, vfiles,
                                            nitfiles, po4files, sifiles, o2files, fefiles,
                                            dicfiles, talkfiles,
                                            iyear, imonth,
                                            ncbry, h_south, theta_s, theta_b,
                                            hc, N, vtransform, Nzgoodmin,
                                            depth_esm, angle_south,
                                            LonT_south, LatT_south, iminT_south, imaxT_south,
                                            jminT_south, jmaxT_south, elemT_south, coefT_south)

            if obc[1] == 1:
                #
                #  Eastern boundary
                #
                print(' ')
                print(' Eastern Boundary')
                print(' ')

                ncbry = glor.interp_bry_ESM('e',
                                            zfiles_3m, tempfiles_3m, saltfiles_3m, ufiles_3m, vfiles_3m,
                                            nitfiles_3m, po4files_3m, sifiles_3m, o2files_3m, fefiles_3m,
                                            dicfiles_3m, talkfiles_3m,
                                            zfiles, tempfiles, saltfiles, ufiles, vfiles,
                                            nitfiles, po4files, sifiles, o2files, fefiles,
                                            dicfiles, talkfiles,
                                            iyear, imonth,
                                            ncbry, h_east, theta_s, theta_b,
                                            hc, N, vtransform, Nzgoodmin,
                                            depth_esm, angle_east,
                                            LonT_east, LatT_east, iminT_east, imaxT_east,
                                            jminT_east, jmaxT_east, elemT_east, coefT_east)

            if obc[2] == 1:
                #
                #  Northern boundary
                #
                print(' ')
                print(' Northern Boundary')
                print(' ')

                ncbry = glor.interp_bry_ESM('n',
                                            zfiles_3m, tempfiles_3m, saltfiles_3m, ufiles_3m, vfiles_3m,
                                            nitfiles_3m, po4files_3m, sifiles_3m, o2files_3m, fefiles_3m,
                                            dicfiles_3m, talkfiles_3m,
                                            zfiles, tempfiles, saltfiles, ufiles, vfiles,
                                            nitfiles, po4files, sifiles, o2files, fefiles,
                                            dicfiles, talkfiles,
                                            iyear, imonth,
                                            ncbry, h_north, theta_s, theta_b,
                                            hc, N, vtransform, Nzgoodmin,
                                            depth_esm, angle_north,
                                            LonT_north, LatT_north, iminT_north, imaxT_north,
                                            jminT_north, jmaxT_north, elemT_north, coefT_north)

            if obc[3] == 1:
                #
                #  Western boundary
                #
                print(' ')
                print(' Western Boundary')
                print(' ')

                ncbry = glor.interp_bry_ESM('w',
                                            zfiles_3m, tempfiles_3m, saltfiles_3m, ufiles_3m, vfiles_3m,
                                            nitfiles_3m, po4files_3m, sifiles_3m, o2files_3m, fefiles_3m,
                                            dicfiles_3m, talkfiles_3m,
                                            zfiles, tempfiles, saltfiles, ufiles, vfiles,
                                            nitfiles, po4files, sifiles, o2files, fefiles,
                                            dicfiles, talkfiles,
                                            iyear, imonth,
                                            ncbry, h_west, theta_s, theta_b,
                                            hc, N, vtransform, Nzgoodmin,
                                            depth_esm, angle_west,
                                            LonT_west, LatT_west, iminT_west, imaxT_west,
                                            jminT_west, jmaxT_west, elemT_west, coefT_west)
            #
            #  End loop on time
            #
            ncbr.close()
