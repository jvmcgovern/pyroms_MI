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
import glob
from netCDF4 import Dataset as netcdf
from netCDF4 import num2date
from scipy.interpolate import griddata

from datetime import date
from calendar import monthrange

import sys

# sys.path.insert(0,'/XXX/')

from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from progressbar import *
from scipy.spatial import Delaunay

if 1 == 1:
    # def main_func():

    #
    # #################### USERS DEFINED VARIABLES ########################
    #

    title = 'Boundary file using GLORYS'

    # my_home_dir = '/home/penven/'

    # crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/NEATL/'
    # bryname = crocofiles_dir + 'croco_bry_PHYBIO_CELTIC_h6.nc'
    # bryname = crocofiles_dir + 'croco_bry_PHYBIO_CELTIC_h8.nc'
    # bryname = crocofiles_dir + 'croco_bry_MERCATOR_.nc'
    # grdname = crocofiles_dir + 'croco_grd.nc'
    grdname = crocofiles_dir + 'ne_atlantic_r4_v2.nc'
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    # grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m

    N = 40
    theta_s = 4.
    theta_b = 0.85
    hc = 10.
    vtransform = 1

    obc = [1, 1, 1, 1]  # open boundaries (1=open , [S E N W])

    time_bry = 0.
    cycle_bry = 0.

    Yorig = 1968  # year origin of time
    Morig = 5  # month origin of time
    Dorig = 23  # day origin of time

    Ystart = 2014
    Mstart = 12
    Dstart = 1

    Yend = 2014
    Mend = 12
    Dend = 31

    # glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    # glorys_ending = '_R20201201_RE01.nc'

    glorysfiles_dir = '/media/dskthree/CMEMS_GLO/'
    # ex 'mercatorglorys12v1_gl12_mean_20150208_R20150211_crop.nc'
    glorys_prefix = 'mercatorglorys12v1_gl12_mean_'
    glorys_ending = '_hazem.nc'
    # glorys_bgc_prefix = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    glorys_bgc_prefix = 'mercatorfreebiorys2v4_global_mean_'
    glorys_bgc_ending = '.nc'

    glorys_step = 1  # time step between outputs in GLORYS12 [days]

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    #
    # END USERS DEFINED VARIABLES
    #

    for iyear in range(Ystart, Yend + 1):
        for imonth in range(Mstart, Mend + 1):

            bryname = crocofiles_dir + 'croco_bry_MERCATOR' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(
                2) + '.nc'

            dpm = monthrange(iyear, imonth)

            print(' ')
            print(' Making boundary file: ' + bryname)
            print(' ')
            print(' Title: ' + title)

            #
            # Create the CROCO boundary file
            #

            glor.create_bryfile(bryname, grdname, title, obc,
                                theta_s, theta_b, hc, N,
                                time_bry, cycle_bry, vtransform)

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
            # Open the CROCO boundary file for writing
            #

            ncbr = netcdf(bryname, 'a')
            ncbry = ncbr.variables

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
                                                          glorys_prefix + str(iyear-1) +
                                                          str(12).zfill(2) + '31' + '_R*.nc'))
                # Get the time in days since Yorig, Morig, Dorig
                Tstart = date.toordinal(date(iyear-1, 12, 31)) - date.toordinal(date(Yorig, Morig, Dorig))
                Tstart = Tstart + 0.5  # 12H
            else:
                mercator_PHY_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                          glorys_prefix + str(iyear) +
                                                          str(imonth-1).zfill(2) + str(monthrange(iyear, imonth-1)[1]) +
                                                          '_R*.nc'))
                # Get the time in days since Yorig, Morig, Dorig
                Tstart = date.toordinal(date(iyear, imonth-1, monthrange(iyear, imonth-1)[1])) - \
                         date.toordinal(date(Yorig, Morig, Dorig))
                Tstart = Tstart + 0.5  # 12H

            mercator_PHY_files_main = sorted(glob.glob(glorysfiles_dir +
                                                       glorys_prefix + str(iyear) +
                                                       str(imonth).zfill(2) + '??' + '_R*.nc'))
            if imonth is 12:
                if len(mercator_PHY_files_main) is dpm[1]:
                    mercator_PHY_files_post1 = sorted(glob.glob(glorysfiles_dir +
                                                                glorys_prefix + str(iyear+1) +
                                                                str(1).zfill(2) + str(1).zfill(2) + '_R*.nc'))
                    mercator_PHY_files_post2 = sorted(glob.glob(glorysfiles_dir +
                                                                glorys_prefix + str(iyear+1) +
                                                                str(1).zfill(2) + str(2).zfill(2) + '_R*.nc'))
                    mercator_PHY_files = \
                        mercator_PHY_files_pre + mercator_PHY_files_main + mercator_PHY_files_post1\
                        + mercator_PHY_files_post2

                    Tend = date.toordinal(date(iyear+1, 1, 2)) - date.toordinal(date(Yorig, Morig, Dorig))
                    Tend = Tend + 0.5  # 12H
                else:
                    mercator_PHY_files = mercator_PHY_files_pre + mercator_PHY_files_main
                    Tend = date.toordinal(date(iyear, imonth,
                                               len(mercator_PHY_files)-1)) - date.toordinal(date(Yorig, Morig, Dorig))
                    Tend = Tend + 0.5  # 12H
            else:
                mercator_PHY_files_post1 = sorted(glob.glob(glorysfiles_dir +
                                                            glorys_prefix + str(iyear) +
                                                            str(imonth+1).zfill(2) + str(1).zfill(2) +
                                                            '_R*.nc'))
                mercator_PHY_files_post2 = sorted(glob.glob(glorysfiles_dir +
                                                            glorys_prefix + str(iyear) +
                                                            str(imonth+1).zfill(2) + str(2).zfill(2) +
                                                            '_R*.nc'))
                mercator_PHY_files = \
                    mercator_PHY_files_pre + mercator_PHY_files_main + mercator_PHY_files_post1 \
                    + mercator_PHY_files_post2
                Tend = date.toordinal(date(iyear, imonth+1, 2)) - date.toordinal(date(Yorig, Morig, Dorig))
                Tend = Tend + 0.5  # 12H

            if imonth is 1:
                mercator_BIO_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                          glorys_bgc_prefix + str(iyear-1) +
                                                          str(12).zfill(2) + '.nc'))
            else:
                mercator_BIO_files_pre = sorted(glob.glob(glorysfiles_dir +
                                                          glorys_bgc_prefix + str(iyear) +
                                                          str(imonth-1).zfill(2) + '.nc'))

            mercator_BIO_files_main = sorted(glob.glob(glorysfiles_dir +
                                                       glorys_bgc_prefix + str(iyear) + str(imonth).zfill(2) +
                                                       '.nc'))
            if imonth is 12:
                mercator_BIO_files_post = sorted(glob.glob(glorysfiles_dir +
                                                 glorys_bgc_prefix + str(iyear+1) +
                                                 str(1).zfill(2) + '.nc'))

                mercator_BIO_files = \
                    mercator_BIO_files_pre + mercator_BIO_files_main + mercator_BIO_files_post

            else:
                mercator_BIO_files_post = sorted(glob.glob(glorysfiles_dir +
                                                           glorys_bgc_prefix + str(iyear) +
                                                           str(imonth+1).zfill(2) + '.nc'))

                mercator_BIO_files = \
                    mercator_BIO_files_pre + mercator_BIO_files_main + mercator_BIO_files_post

            Tbry_str = "%06d" % Tstart

            glorysname = mercator_PHY_files[0]

            print(' OPEN : ' + glorysname)

            glorysnamebgc = mercator_BIO_files[0]

            print(' OPEN : ' + glorysname)

            #
            # open the first GLORYS file
            #

            ncgl = netcdf(glorysname, 'r')
            ncglo = ncgl.variables
            depthp = np.array(ncglo['depth'][:])

            ncbgc = netcdf(glorysnamebgc, 'r')
            ncbgco = ncbgc.variables
            depthb = np.array(ncbgco['depth'][:])

            [Nzp] = np.shape(depthp)
            [Nzb] = np.shape(depthb)

            dl = 1
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
                    = glor.get_delaunay_bry_GLO_BGC_PHY(lon_south, lat_south, dl, ncglo, ncbgco)

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
                    = glor.get_delaunay_bry_GLO_BGC_PHY(lon_east, lat_east, dl, ncglo, ncbgco)

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
                    = glor.get_delaunay_bry_GLO_BGC_PHY(lon_north, lat_north, dl, ncglo, ncbgco)

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
                    = glor.get_delaunay_bry_GLO_BGC_PHY(lon_west, lat_west, dl, ncglo, ncbgco)

            #
            #  Close the GLORYS netcdf file
            #

            ncgl.close()
            ncbgc.close()

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

            PISCES_interpd = 1

            for mpf in np.arange(0, Tend - Tstart + 1, glorys_step):
                Tbry = Tbries[int(mpf)]
                print(Tbry)

                tndx_bry = tndx_bry + 1

                if PISCES_interpd is 0:
                    # PHYSICS GETS A DAILY TIMESTEP
                    ncbry['bry_time'][tndx_bry] = Tbry
                    ncbry['tclm_time'][tndx_bry] = Tbry
                    ncbry['temp_time'][tndx_bry] = Tbry
                    ncbry['sclm_time'][tndx_bry] = Tbry
                    ncbry['salt_time'][tndx_bry] = Tbry
                    ncbry['uclm_time'][tndx_bry] = Tbry
                    ncbry['vclm_time'][tndx_bry] = Tbry
                    ncbry['v2d_time'][tndx_bry] = Tbry
                    ncbry['v3d_time'][tndx_bry] = Tbry
                    ncbry['ssh_time'][tndx_bry] = Tbry
                    ncbry['zeta_time'][tndx_bry] = Tbry

                if PISCES_interpd is 1:
                    # BGC VARIABLES GETS MONTHLY TIMESTEP (3 timepoints: prior, present and next month)
                    if imonth == 1:
                        Tinin1 = date.toordinal(date(iyear - 1, 12, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                        Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                        Tinin3 = date.toordinal(date(iyear, imonth + 1, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                    elif imonth == 12:
                        Tinin1 = date.toordinal(date(iyear, imonth - 1, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                        Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                        Tinin3 = date.toordinal(date(iyear + 1, 1, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                    else:
                        Tinin1 = date.toordinal(date(iyear, imonth - 1, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                        Tinin2 = date.toordinal(date(iyear, imonth, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                        Tinin3 = date.toordinal(date(iyear, imonth + 1, 15)) - \
                                 date.toordinal(date(Yorig, Morig, Dorig))
                    n_bry = [np.float64(Tinin1), np.float64(Tinin2), np.float64(Tinin3)]

                    ncpispri = netcdf(mercator_BIO_files[0], 'r')
                    ncpisprio = ncpispri.variables

                    Tinip1 = ncpisprio['time'][:][0] / 24 + date.toordinal(date(1950, 1, 1)) - \
                             date.toordinal(date(Yorig, Morig, Dorig))
                    ncpiscur = netcdf(mercator_BIO_files[1], 'r')
                    ncpiscuro = ncpiscur.variables

                    Tinip2 = ncpiscuro['time'][:][0] / 24 + date.toordinal(date(1950, 1, 1)) - \
                             date.toordinal(date(Yorig, Morig, Dorig))
                    if len(mercator_BIO_files_pre) == 0:
                        Tinip1 = np.float64(Tinip2) - 32
                    ncpispos = netcdf(mercator_BIO_files[2], 'r')
                    ncpisposo = ncpispos.variables

                    Tinip3 = ncpisposo['time'][:][0] / 24 + date.toordinal(date(1950, 1, 1)) - \
                             date.toordinal(date(Yorig, Morig, Dorig))

                    # PHYSICS GETS A DAILY TIMESTEP
                    ncbry['bry_time'][tndx_bry] = Tbry
                    ncbry['tclm_time'][tndx_bry] = Tbry
                    ncbry['temp_time'][tndx_bry] = Tbry
                    ncbry['sclm_time'][tndx_bry] = Tbry
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
                    ncbry['fe_time'][0] = np.float64(Tinin1)
                    ncbry['no3_time'][1] = np.float64(Tinin2)
                    ncbry['po4_time'][1] = np.float64(Tinin2)
                    ncbry['si_time'][1] = np.float64(Tinin2)
                    ncbry['o2_time'][1] = np.float64(Tinin2)
                    ncbry['dic_time'][1] = np.float64(Tinin2)
                    ncbry['talk_time'][1] = np.float64(Tinin2)
                    ncbry['fe_time'][1] = np.float64(Tinin2)
                    ncbry['no3_time'][2] = np.float64(Tinin3)
                    ncbry['po4_time'][2] = np.float64(Tinin3)
                    ncbry['si_time'][2] = np.float64(Tinin3)
                    ncbry['o2_time'][2] = np.float64(Tinin3)
                    ncbry['dic_time'][2] = np.float64(Tinin3)
                    ncbry['talk_time'][2] = np.float64(Tinin3)
                    ncbry['fe_time'][2] = np.float64(Tinin3)

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
                    print(' Soutern Boundary')
                    print(' ')

                    ncbry = glor.interp_bry_GLO_BGC_PHY('s', PISCES_interpd,
                                                        ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                        tndx_glo, ncbry, tndx_bry, h_south, theta_s, theta_b,
                                                        hc, N, vtransform, Nzgoodmin,
                                                        depthp, depthb, angle_south,
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

                    ncbry = glor.interp_bry_GLO_BGC_PHY('e', PISCES_interpd,
                                                        ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                        tndx_glo, ncbry, tndx_bry, h_east, theta_s, theta_b,
                                                        hc, N, vtransform, Nzgoodmin,
                                                        depthp, depthb, angle_east,
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

                    ncbry = glor.interp_bry_GLO_BGC_PHY('n', PISCES_interpd,
                                                        ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                        tndx_glo, ncbry, tndx_bry, h_north, theta_s, theta_b,
                                                        hc, N, vtransform, Nzgoodmin,
                                                        depthp, depthb, angle_north,
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

                    ncbry = glor.interp_bry_GLO_BGC_PHY('w', PISCES_interpd,
                                                        ncglo, ncpisprio, ncpiscuro, ncpisposo,
                                                        tndx_glo, ncbry, tndx_bry, h_west, theta_s, theta_b,
                                                        hc, N, vtransform, Nzgoodmin,
                                                        depthp, depthb, angle_west,
                                                        LonT_west, LatT_west, iminT_west, imaxT_west,
                                                        jminT_west, jmaxT_west, elemT_west, coefT_west,
                                                        LonU_west, LatU_west, iminU_west, imaxU_west,
                                                        jminU_west, jmaxU_west, elemU_west, coefU_west,
                                                        LonV_west, LatV_west, iminV_west, imaxV_west,
                                                        jminV_west, jmaxV_west, elemV_west, coefV_west,
                                                        LonP_west, LatP_west, iminP_west, imaxP_west,
                                                        jminP_west, jmaxP_west, elemP_west, coefP_west)

                if PISCES_interpd is 1:
                    PISCES_interpd = 0
                    ncpispri.close()
                    ncpiscur.close()
                    ncpispos.close()

                ncgl.close()

            #
            #  End loop on time
            #

            ncbr.close()
