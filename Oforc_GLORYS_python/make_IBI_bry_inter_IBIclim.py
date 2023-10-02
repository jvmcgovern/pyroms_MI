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

    # Ystart = 1993
    # Mstart = 1
    # Dstart = 1
    #
    # Yend = 2021
    # Mend = 12
    # Dend = 31

    Ystart = 1993
    Mstart = 1
    Dstart = 1

    Yend = 1993
    Mend = 12
    Dend = 31

    # CMEMS_IBI MY (Reanalysis data)
    glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # Reanalysis
    glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    glorys_ending = '_R*_RE01.nc'
    # # Near realtime
    # glorys_prefix = 'CMEMS_v6r1_IBI_PHY_NRT_NL_01dav_'
    # glorys_ending = '_R*_AN*.nc'

    # CMEMS_IBI MY (monthly average for TALK derivation from temp/sal/DIC/pH)
    # 'CELTIC/CMEMS_IBI/CMEMS_v5r1_IBI_PHY_MY_PdE_01mav_20051201_20051231_R20201201_RE01.nc'
    # Reanalysis
    glorys_mth_prefix = 'PHY_CMEMS_v5r1_IBI_PHY_MY_PdE_01mav_'
    # # Near realtime
    # glorys_mth_prefix = 'CMEMS_v6r1_IBI_PHY_NRT_NL_01mav_'

    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    # PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI_PISCESvars_012013_202112/monthly/'
    # PISCES24_prefix = 'IBI36_cg_1m-m_'
    # PISCES24_ending = '_3DT-bgc_hcst.nc'
    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # Reanalysis
    PISCES24_prefix = 'NPSO_TALK_CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_'
    PISCES24_ending = '_R*_RE01_WOAres.nc'
    # # Near realtime
    # PISCES24_prefix = 'CMEMS_v7r1_IBI_BIO_NRT_NL_01mav_'
    # PISCES24_ending = '_R*_AN*.nc'

    # 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_20190201_20190228_R20201201_RE01.nc'

    glorys_step = 1  # time step between outputs in GLORYS12 [days]

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    #
    # END USERS DEFINED VARIABLES
    #

    for iyear in range(Ystart, Yend + 1):
        for imonth in range(Mstart, Mend + 1):

            # bryname = crocofiles_dir + 'croco_bry_MERCATOR_hc50m_MSL_ng' + '_Y' + str(iyear) + 'M' \
            #           + str(imonth).zfill(2) + '.nc'

            bryname = crocofiles_dir + 'croco_bry_IBICLIM_MERCATOR_hc50m_MSL_ng' + '_Y' + str(iyear) + 'M' \
                      + str(imonth).zfill(2) + '.nc'

            dpm = monthrange(iyear, imonth)

            print(' ')
            print(' Making boundary file: ' + bryname)
            print(' ')
            print(' Title: ' + title)

            #
            # Create the CROCO boundary file
            #
            glor.create_bryfile_BGC_PHY(bryname, grdname, title, obc,
                                        theta_s, theta_b, hc, N,
                                        time_bry, 0, 0, cycle_bry, vtransform)

            #
            # Open the CROCO boundary file for writing
            #

            ncbr = netcdf(bryname, 'a')
            ncbry = ncbr.variables

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

            #
            # Get the Delaunay triangulations for each boundary
            #

            print(' ')
            print(' Get the Delaunay triangulations for each boundary')
            print(' ')

            #
            # Get the first GLORYS file name (for the positions of variables)
            #

            if imonth is 1:
                mercator_PHY_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                          glorys_prefix + str(iyear - 1) +
                                                          str(12).zfill(2) + '31_' + str(iyear - 1) +
                                                          str(12).zfill(2) + '31' + glorys_ending))
                mercator_PHYmth_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                             glorys_mth_prefix +
                                                             str(12).zfill(2) + '_' +
                                                             str(12).zfill(2) + '_R20201201_RE01_WOAres.nc'))
                # Get the time in days since Yorig, 1, 1
                Tstart = date.toordinal(date(iyear - 1, 12, 31)) - date.toordinal(date(Yorig, 1, 1))
                Tstart = Tstart + 0.5  # 12H
            else:
                mercator_PHY_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                          glorys_prefix + str(iyear) +
                                                          str(imonth - 1).zfill(2) +
                                                          str(monthrange(iyear, imonth - 1)[1]) +
                                                          '_' + str(iyear) + str(imonth - 1).zfill(2) +
                                                          str(monthrange(iyear, imonth - 1)[1]) + glorys_ending))
                mercator_PHYmth_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                             glorys_mth_prefix +
                                                             str(imonth - 1).zfill(2) +
                                                             '_' + str(imonth - 1).zfill(2) +
                                                             '_R20201201_RE01_WOAres.nc'))
                # Get the time in days since Yorig, 1, 1
                Tstart = date.toordinal(date(iyear, imonth - 1,
                                             monthrange(iyear, imonth - 1)[1])) - date.toordinal(date(Yorig, 1, 1))
                Tstart = Tstart + 0.5  # 12H

            mercator_PHY_files_main = sorted(glob.glob(glorysfiles_dir +
                                                       glorys_prefix + str(iyear) +
                                                       str(imonth).zfill(2) + '??_' + str(iyear) +
                                                       str(imonth).zfill(2)
                                                       + '??' + glorys_ending))

            mercator_PHYmth_files_main = sorted(glob.glob(glorysfiles_dir +
                                                          glorys_mth_prefix +
                                                          str(imonth).zfill(2) + '_' +
                                                          str(imonth).zfill(2) + '_R20201201_RE01_WOAres.nc'))

            if imonth is 12:
                if len(mercator_PHY_files_main) is dpm[1]:
                    mercator_PHY_files_post1 = sorted(glob.glob(glorysfiles_dir +
                                                                glorys_prefix + str(iyear + 1) +
                                                                str(1).zfill(2) + str(1).zfill(2) + '_'
                                                                + str(iyear + 1) +
                                                                str(1).zfill(2) + str(1).zfill(2) + glorys_ending))
                    mercator_PHY_files_post2 = sorted(glob.glob(glorysfiles_dir +
                                                                glorys_prefix + str(iyear + 1) +
                                                                str(1).zfill(2) + str(2).zfill(2) + '_'
                                                                + str(iyear + 1) +
                                                                str(1).zfill(2) + str(2).zfill(2) + glorys_ending))
                    mercator_PHY_files = \
                        mercator_PHY_files_pre + mercator_PHY_files_main + mercator_PHY_files_post1 \
                        + mercator_PHY_files_post2

                    mercator_PHYmth_files_post = sorted(glob.glob(glorysfiles_dir +
                                                                  glorys_mth_prefix +
                                                                  str(1).zfill(2) + '_' + str(1).zfill(2) +
                                                                  '_R20201201_RE01_WOAres.nc'))

                    Tend = date.toordinal(date(iyear + 1, 1, 2)) - date.toordinal(date(Yorig, 1, 1))
                    Tend = Tend + 0.5  # 12H
                else:
                    mercator_PHY_files = mercator_PHY_files_pre + mercator_PHY_files_main
                    Tend = date.toordinal(date(iyear, imonth,
                                               len(mercator_PHY_files) - 1)) - date.toordinal(date(Yorig, 1, 1))
                    Tend = Tend + 0.5  # 12H
            else:
                mercator_PHY_files_post1 = sorted(glob.glob(glorysfiles_dir +
                                                            glorys_prefix + str(iyear) +
                                                            str(imonth + 1).zfill(2) + str(1).zfill(2) +
                                                            '_' + str(iyear) + str(imonth + 1).zfill(2) +
                                                            str(1).zfill(2) + glorys_ending))
                mercator_PHY_files_post2 = sorted(glob.glob(glorysfiles_dir +
                                                            glorys_prefix + str(iyear) +
                                                            str(imonth + 1).zfill(2) + str(2).zfill(2) +
                                                            '_' + str(iyear) + str(imonth + 1).zfill(2) +
                                                            str(2).zfill(2) + glorys_ending))
                if len(mercator_PHY_files_pre) == 0:
                    mercator_PHY_files = [mercator_PHY_files_main[0]] + mercator_PHY_files_main + \
                                         mercator_PHY_files_post1 + mercator_PHY_files_post2
                else:
                    mercator_PHY_files = \
                        mercator_PHY_files_pre + mercator_PHY_files_main + mercator_PHY_files_post1 \
                        + mercator_PHY_files_post2

                mercator_PHYmth_files_post = sorted(glob.glob(glorysfiles_dir +
                                                              glorys_prefix + str(iyear) +
                                                              str(imonth + 1).zfill(2) + '_' +
                                                              str(imonth + 1).zfill(2) +
                                                              '_R20201201_RE01_WOAres.nc'))

                Tend = date.toordinal(date(iyear, imonth + 1, 2)) - date.toordinal(date(Yorig, 1, 1))
                Tend = Tend + 0.5  # 12H

            if imonth is 1:
                mercator_BIO_files_pre = sorted(glob.glob(PISCES24files_dir +
                                                          PISCES24_prefix +
                                                          str(12).zfill(2) + '_' + str(12).zfill(2) + PISCES24_ending))
            else:
                mercator_BIO_files_pre = sorted(glob.glob(PISCES24files_dir +
                                                          PISCES24_prefix +
                                                          str(imonth - 1).zfill(2) + '_' +
                                                          str(imonth - 1).zfill(2) + PISCES24_ending))

            mercator_BIO_files_main = sorted(glob.glob(PISCES24files_dir +
                                                       PISCES24_prefix +
                                                       str(imonth).zfill(2) + '_' +
                                                       str(imonth).zfill(2) + PISCES24_ending))
            if imonth is 12:
                mercator_BIO_files_post = sorted(glob.glob(PISCES24files_dir +
                                                           PISCES24_prefix + str(1).zfill(2) + '_' + str(1).zfill(2) +
                                                           PISCES24_ending))
                mercator_BIO_files = \
                    mercator_BIO_files_pre + mercator_BIO_files_main + mercator_BIO_files_post
            else:
                mercator_BIO_files_post = sorted(glob.glob(PISCES24files_dir +
                                                           PISCES24_prefix +
                                                           str(imonth + 1).zfill(2) + '_' +
                                                           str(imonth + 1).zfill(2) + PISCES24_ending))
            if len(mercator_BIO_files_pre) == 0:
                mercator_BIO_files = mercator_BIO_files_main + mercator_BIO_files_main
            else:
                mercator_BIO_files = mercator_BIO_files_pre + mercator_BIO_files_main

            if len(mercator_BIO_files_post) == 0:
                mercator_BIO_files = mercator_BIO_files + mercator_BIO_files_main
            else:
                mercator_BIO_files = mercator_BIO_files + mercator_BIO_files_post

            if len(mercator_PHYmth_files_pre) == 0:
                mercator_PHYmth_files = mercator_PHYmth_files_main + mercator_PHYmth_files_main
            else:
                mercator_PHYmth_files = mercator_PHYmth_files_pre + mercator_PHYmth_files_main

            if len(mercator_PHYmth_files_post) == 0:
                mercator_PHYmth_files = mercator_PHYmth_files + mercator_PHYmth_files_main
            else:
                mercator_PHYmth_files = mercator_PHYmth_files + mercator_PHYmth_files_post

            Tbry_str = "%06d" % Tstart

            glorysname = mercator_PHY_files[0]

            print(' OPEN : ' + glorysname)

            NORESMnamebgc_ex = NORESMfiles_dir + NORESM_prefix + 'no3no2' + NORESM_ending

            print(NORESMnamebgc_ex)

            PISCESnamebgc = mercator_BIO_files[0]

            print(PISCESnamebgc)

            #
            # open the first GLORYS file
            #
            ncgl = netcdf(glorysname, 'r')
            ncglo = ncgl.variables
            depthg = np.array(ncglo['depth'][:])

            ncpis = netcdf(PISCESnamebgc, 'r')
            ncpiso = ncpis.variables
            depthp = np.array(ncpiso['depth'][:])

            # Nitrate from NORESM
            ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
            ncnornio = ncnorni.variables
            depthn = np.array(ncnornio['depth'][:])

            [Nz] = np.shape(depthg)

            dlg = 1
            dlp = 2
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

                (LonT_south, LatT_south, iminT_south, imaxT_south, jminT_south, jmaxT_south, elemT_south, coefT_south,
                 LonU_south, LatU_south, iminU_south, imaxU_south, jminU_south, jmaxU_south, elemU_south, coefU_south,
                 LonV_south, LatV_south, iminV_south, imaxV_south, jminV_south, jmaxV_south, elemV_south, coefV_south,
                 LonP_south, LatP_south, iminP_south, imaxP_south, jminP_south, jmaxP_south, elemP_south, coefP_south) \
                    = glor.get_delaunay_bry_IBIclim(lon_south, lat_south, dlg, dlp, ncglo, ncpiso)

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

                (LonT_east, LatT_east, iminT_east, imaxT_east, jminT_east, jmaxT_east, elemT_east, coefT_east,
                 LonU_east, LatU_east, iminU_east, imaxU_east, jminU_east, jmaxU_east, elemU_east, coefU_east,
                 LonV_east, LatV_east, iminV_east, imaxV_east, jminV_east, jmaxV_east, elemV_east, coefV_east,
                 LonP_east, LatP_east, iminP_east, imaxP_east, jminP_east, jmaxP_east, elemP_east, coefP_east) \
                    = glor.get_delaunay_bry_IBIclim(lon_east, lat_east, dlg, dlp, ncglo, ncpiso)

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

                (LonT_north, LatT_north, iminT_north, imaxT_north, jminT_north, jmaxT_north, elemT_north, coefT_north,
                 LonU_north, LatU_north, iminU_north, imaxU_north, jminU_north, jmaxU_north, elemU_north, coefU_north,
                 LonV_north, LatV_north, iminV_north, imaxV_north, jminV_north, jmaxV_north, elemV_north, coefV_north,
                 LonP_north, LatP_north, iminP_north, imaxP_north, jminP_north, jmaxP_north, elemP_north, coefP_north) \
                    = glor.get_delaunay_bry_IBIclim(lon_north, lat_north, dlg, dlp, ncglo, ncpiso)

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

                (LonT_west, LatT_west, iminT_west, imaxT_west, jminT_west, jmaxT_west, elemT_west, coefT_west,
                 LonU_west, LatU_west, iminU_west, imaxU_west, jminU_west, jmaxU_west, elemU_west, coefU_west,
                 LonV_west, LatV_west, iminV_west, imaxV_west, jminV_west, jmaxV_west, elemV_west, coefV_west,
                 LonP_west, LatP_west, iminP_west, imaxP_west, jminP_west, jmaxP_west, elemP_west, coefP_west) \
                    = glor.get_delaunay_bry_IBIclim(lon_west, lat_west, dlg, dlp, ncglo, ncpiso)

            #
            #  Close the GLORYS netcdf file
            #

            ncgl.close()
            ncnorni.close()
            ncpis.close()

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

            tndx_bry = -1

            #
            # Loop on files
            #
            Tbries = np.arange(Tstart, Tend + 1, glorys_step)

            PISCES_NORESM_interpd = 1

            for mpf in np.arange(0, Tend - Tstart + 1, glorys_step):
                Tbry = Tbries[int(mpf)]
                print(Tbry)

                tndx_bry = tndx_bry + 1

                if PISCES_NORESM_interpd is 0:
                    # PHYSICS GETS A DAILY TIMESTEP
                    ncbry['bry_time'][tndx_bry] = Tbry
                    ncbry['temp_time'][tndx_bry] = Tbry
                    ncbry['salt_time'][tndx_bry] = Tbry
                    ncbry['uclm_time'][tndx_bry] = Tbry
                    ncbry['vclm_time'][tndx_bry] = Tbry
                    ncbry['v2d_time'][tndx_bry] = Tbry
                    ncbry['v3d_time'][tndx_bry] = Tbry
                    ncbry['ssh_time'][tndx_bry] = Tbry
                    ncbry['zeta_time'][tndx_bry] = Tbry

                if PISCES_NORESM_interpd is 1:
                    # BGC VARIABLES GETS MONTHLY TIMESTEP (3 timepoints: prior, present and next month)
                    ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
                    ncnornio = ncnorni.variables
                    Tinin = datetime(iyear, imonth, 15)
                    # Tininxcur = d2i(Tinin, ncnornio['time'], select='exact')
                    # Tininxpri = Tininxcur - 1
                    # Tininxpos = Tininxcur + 1
                    # Tininxn = [Tininxpri, Tininxcur, Tininxpos]
                    Tininxn = [0, 1, 2]

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
                    ncnorni.close()
                    n_bry = [np.float64(Tinin1), np.float64(Tinin2), np.float64(Tinin3)]

                    ncpispri = netcdf(mercator_BIO_files[0], 'r')
                    ncpisprio = ncpispri.variables
                    # Tinip1 = ncpisprio['time_centered'][:][0]/86400 + date.toordinal(date(1950, 1, 1)) - \
                    #          date.toordinal(date(Yorig, 1, 1))
                    Tinip1 = ncpisprio['time'][:][0] / 24 + date.toordinal(date(1950, 1, 1)) - \
                             date.toordinal(date(Yorig, 1, 1))
                    ncpiscur = netcdf(mercator_BIO_files[1], 'r')
                    ncpiscuro = ncpiscur.variables
                    # Tinip2 = ncpiscuro['time_centered'][:][0]/86400 + date.toordinal(date(1950, 1, 1)) - \
                    #          date.toordinal(date(Yorig, 1, 1))
                    Tinip2 = ncpiscuro['time'][:][0] / 24 + date.toordinal(date(1950, 1, 1)) - \
                             date.toordinal(date(Yorig, 1, 1))
                    if len(mercator_BIO_files_pre) == 0:
                        Tinip1 = np.float64(Tinip2) - 32
                    ncpispos = netcdf(mercator_BIO_files[2], 'r')
                    ncpisposo = ncpispos.variables
                    # Tinip3 = ncpisposo['time_centered'][:][0]/86400 + date.toordinal(date(1950, 1, 1)) - \
                    #          date.toordinal(date(Yorig, 1, 1))
                    Tinip3 = ncpisposo['time'][:][0] / 24 + date.toordinal(date(1950, 1, 1)) - \
                             date.toordinal(date(Yorig, 1, 1))

                    ncphympri = netcdf(mercator_PHYmth_files[0], 'r')
                    ncphymprio = ncphympri.variables

                    ncphymcur = netcdf(mercator_PHYmth_files[1], 'r')
                    ncphymcuro = ncphymcur.variables

                    ncphympos = netcdf(mercator_PHYmth_files[2], 'r')
                    ncphymposo = ncphympos.variables

                    p_bry = [Tinin1, Tinin2, Tinin3]

                    # PHYSICS GETS A DAILY TIMESTEP
                    ncbry['bry_time'][tndx_bry] = Tbry
                    ncbry['temp_time'][tndx_bry] = Tbry
                    ncbry['salt_time'][tndx_bry] = Tbry
                    ncbry['uclm_time'][tndx_bry] = Tbry
                    ncbry['vclm_time'][tndx_bry] = Tbry
                    ncbry['v2d_time'][tndx_bry] = Tbry
                    ncbry['v3d_time'][tndx_bry] = Tbry
                    ncbry['ssh_time'][tndx_bry] = Tbry
                    ncbry['zeta_time'][tndx_bry] = Tbry

                    # BGC variables get their timestamps
                    ncbry['no3_time'][0] = np.float64(Tinin1)
                    ncbry['po4_time'][0] = np.float64(Tinin1)
                    ncbry['si_time'][0] = np.float64(Tinin1)
                    ncbry['o2_time'][0] = np.float64(Tinin1)
                    ncbry['dic_time'][0] = np.float64(Tinin1)
                    ncbry['talk_time'][0] = np.float64(Tinin1)
                    ncbry['no3_time'][1] = np.float64(Tinin2)
                    ncbry['po4_time'][1] = np.float64(Tinin2)
                    ncbry['si_time'][1] = np.float64(Tinin2)
                    ncbry['o2_time'][1] = np.float64(Tinin2)
                    ncbry['dic_time'][1] = np.float64(Tinin2)
                    ncbry['talk_time'][1] = np.float64(Tinin2)
                    ncbry['no3_time'][2] = np.float64(Tinin3)
                    ncbry['po4_time'][2] = np.float64(Tinin3)
                    ncbry['si_time'][2] = np.float64(Tinin3)
                    ncbry['o2_time'][2] = np.float64(Tinin3)
                    ncbry['dic_time'][2] = np.float64(Tinin3)
                    ncbry['talk_time'][2] = np.float64(Tinin3)
                    ncbry['fer_time'][0] = np.float64(Tinin1)
                    ncbry['nh4_time'][0] = np.float64(Tinin1)
                    ncbry['fer_time'][1] = np.float64(Tinin2)
                    ncbry['nh4_time'][1] = np.float64(Tinin2)
                    ncbry['fer_time'][2] = np.float64(Tinin3)
                    ncbry['nh4_time'][2] = np.float64(Tinin3)
                    ncbry['ph_time'][0] = np.float64(Tinin1)
                    ncbry['ph_time'][1] = np.float64(Tinin2)
                    ncbry['ph_time'][2] = np.float64(Tinin3)
                    # Monthly average salinity and temperature for CO2SYS
                    ncbry['salm_time'][0] = np.float64(Tinin1)
                    ncbry['tmpm_time'][0] = np.float64(Tinin1)
                    ncbry['salm_time'][1] = np.float64(Tinin2)
                    ncbry['tmpm_time'][1] = np.float64(Tinin2)
                    ncbry['salm_time'][2] = np.float64(Tinin3)
                    ncbry['tmpm_time'][2] = np.float64(Tinin3)

                #
                # Get the first GLORYS file name (for the positions of variables)
                #
                Tbry_str = "%06d" % Tbry

                glorysname = mercator_PHY_files[int(mpf)]
                print(' OPEN : ' + glorysname)
                ncgl = netcdf(glorysname, 'r')
                ncglo = ncgl.variables

                if obc[0] == 1:
                    #
                    #  Southern boundary
                    #
                    print(' ')
                    print(' Southern Boundary')
                    print(' ')

                    ncbry = glor.interp_bry_IBIclim('s', PISCES_NORESM_interpd, Tininxn,
                                                    ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                    ncphymprio, ncphymcuro, ncphymposo,
                                                    tndx_glo, ncbry, tndx_bry, h_south, theta_s, theta_b,
                                                    hc, N, vtransform, Nzgoodmin,
                                                    depthg, depthp, angle_south,
                                                    LonT_south, LatT_south, iminT_south, imaxT_south,
                                                    jminT_south, jmaxT_south, elemT_south, coefT_south,
                                                    LonU_south, LatU_south, iminU_south, imaxU_south,
                                                    jminU_south, jmaxU_south, elemU_south, coefU_south,
                                                    LonV_south, LatV_south, iminV_south, imaxV_south,
                                                    jminV_south, jmaxV_south, elemV_south, coefV_south,
                                                    LonP_south, LatP_south, iminP_south, imaxP_south,
                                                    jminP_south, jmaxP_south, elemP_south, coefP_south)

                if obc[1] == 1:
                    #
                    #  Eastern boundary
                    #
                    print(' ')
                    print(' Eastern Boundary')
                    print(' ')

                    ncbry = glor.interp_bry_IBIclim('e', PISCES_NORESM_interpd, Tininxn,
                                                    ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                    ncphymprio, ncphymcuro, ncphymposo,
                                                    tndx_glo, ncbry, tndx_bry, h_east, theta_s, theta_b,
                                                    hc, N, vtransform, Nzgoodmin,
                                                    depthg, depthp, angle_east,
                                                    LonT_east, LatT_east, iminT_east, imaxT_east,
                                                    jminT_east, jmaxT_east, elemT_east, coefT_east,
                                                    LonU_east, LatU_east, iminU_east, imaxU_east,
                                                    jminU_east, jmaxU_east, elemU_east, coefU_east,
                                                    LonV_east, LatV_east, iminV_east, imaxV_east,
                                                    jminV_east, jmaxV_east, elemV_east, coefV_east,
                                                    LonP_east, LatP_east, iminP_east, imaxP_east,
                                                    jminP_east, jmaxP_east, elemP_east, coefP_east)

                if obc[2] == 1:
                    #
                    #  Northern boundary
                    #
                    print(' ')
                    print(' Northern Boundary')
                    print(' ')

                    ncbry = glor.interp_bry_IBIclim('n', PISCES_NORESM_interpd, Tininxn,
                                                    ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                    ncphymprio, ncphymcuro, ncphymposo,
                                                    tndx_glo, ncbry, tndx_bry, h_north, theta_s, theta_b,
                                                    hc, N, vtransform, Nzgoodmin,
                                                    depthg, depthp, angle_north,
                                                    LonT_north, LatT_north, iminT_north, imaxT_north,
                                                    jminT_north, jmaxT_north, elemT_north, coefT_north,
                                                    LonU_north, LatU_north, iminU_north, imaxU_north,
                                                    jminU_north, jmaxU_north, elemU_north, coefU_north,
                                                    LonV_north, LatV_north, iminV_north, imaxV_north,
                                                    jminV_north, jmaxV_north, elemV_north, coefV_north,
                                                    LonP_north, LatP_north, iminP_north, imaxP_north,
                                                    jminP_north, jmaxP_north, elemP_north, coefP_north)

                if obc[3] == 1:
                    #
                    #  Western boundary
                    #
                    print(' ')
                    print(' Western Boundary')
                    print(' ')

                    ncbry = glor.interp_bry_IBIclim('w', PISCES_NORESM_interpd, Tininxn,
                                                    ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                    ncphymprio, ncphymcuro, ncphymposo,
                                                    tndx_glo, ncbry, tndx_bry, h_west, theta_s, theta_b,
                                                    hc, N, vtransform, Nzgoodmin,
                                                    depthg, depthp, angle_west,
                                                    LonT_west, LatT_west, iminT_west, imaxT_west,
                                                    jminT_west, jmaxT_west, elemT_west, coefT_west,
                                                    LonU_west, LatU_west, iminU_west, imaxU_west,
                                                    jminU_west, jmaxU_west, elemU_west, coefU_west,
                                                    LonV_west, LatV_west, iminV_west, imaxV_west,
                                                    jminV_west, jmaxV_west, elemV_west, coefV_west,
                                                    LonP_west, LatP_west, iminP_west, imaxP_west,
                                                    jminP_west, jmaxP_west, elemP_west, coefP_west)

                if PISCES_NORESM_interpd is 1:
                    PISCES_NORESM_interpd = 0

                ncgl.close()
            # Code here will take temp, salinity, DIC and pH and derive TALK from it

            # South
            results_south = pyco2.sys(par1=ncbry['PH_south'][:], par1_type=3,
                                      par2=ncbry['DIC_south'][:], par2_type=2,
                                      temperature=ncbry['tmpm_south'][:],
                                      salinity=ncbry['salm_south'][:],
                                      opt_pH_scale=1, opt_k_carbonic=4,
                                      opt_k_bisulfate=1, opt_total_borate=1)
            ncbry['TALK_south'][:] = results_south["alkalinity"]

            # North
            results_north = pyco2.sys(par1=ncbry['PH_north'][:], par1_type=3,
                                      par2=ncbry['DIC_north'][:], par2_type=2,
                                      temperature=ncbry['tmpm_north'][:],
                                      salinity=ncbry['salm_north'][:],
                                      opt_pH_scale=1, opt_k_carbonic=4,
                                      opt_k_bisulfate=1, opt_total_borate=1)
            ncbry['TALK_north'][:] = results_north["alkalinity"]

            # East
            results_east = pyco2.sys(par1=ncbry['PH_east'][:], par1_type=3,
                                     par2=ncbry['DIC_east'][:], par2_type=2,
                                     temperature=ncbry['tmpm_east'][:],
                                     salinity=ncbry['salm_east'][:],
                                     opt_pH_scale=1, opt_k_carbonic=4,
                                     opt_k_bisulfate=1, opt_total_borate=1)
            ncbry['TALK_east'][:] = results_east["alkalinity"]

            # West
            results_west = pyco2.sys(par1=ncbry['PH_west'][:], par1_type=3,
                                     par2=ncbry['DIC_west'][:], par2_type=2,
                                     temperature=ncbry['tmpm_west'][:],
                                     salinity=ncbry['salm_west'][:],
                                     opt_pH_scale=1, opt_k_carbonic=4,
                                     opt_k_bisulfate=1, opt_total_borate=1)
            ncbry['TALK_west'][:] = results_west["alkalinity"]

            #
            #  End loop on time
            #
            ncbr.close()
