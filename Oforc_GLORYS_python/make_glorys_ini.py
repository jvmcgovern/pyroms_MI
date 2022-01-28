#
#
######################################################################
######################################################################
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
#

import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset as netcdf
from scipy.interpolate import griddata

from datetime import date
import sys

# sys.path.insert(0,'')

from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from progressbar import *

if 1 == 1:
    # def main_func():

    #
    #
    #
    # #################### USERS DEFINED VARIABLES ########################
    #
    #
    #

    title = 'Initial file using GLORYS'

    # my_home_dir = '/home/penven/'

    # crocofiles_dir = my_home_dir + 'SWAG/Run_TEST/CROCO_FILES/'
    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    # ininame = crocofiles_dir + 'croco_ini_PHY_CELTIC.nc'
    ininame = crocofiles_dir + 'croco_ini_PHYBIO_CELTIC.nc'
    grdname = crocofiles_dir + 'croco_grd.nc'

    N = 20
    theta_s = 7.
    theta_b = 0.
    hc = 20.
    vtransform = 2

    Yorig = 1990  # year origin of time : days since Yorig-01-01

    Yini = 2005
    Mini = 1
    Dini = 1

    # Date in form YYYYMMDD
    ini_date = '20050101'

    # glorysfiles_dir = my_home_dir + 'SWAG/GLORYS_FILES/'
    glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    glorys_bgc_prefix = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    glorys_ending = '_R20201201_RE01.nc'
    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation
    comp_delaunay = 1  # 1: compute delaunay triangulations - 0: use saved matrices (for debugging)

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)
    tndx_ini = 0  # time index in the CROCO initial file (should be only 0 here)

    #
    #
    # ################## END USERS DEFINED VARIABLES ######################
    #
    #

    #
    # Get the time in days since Yorig,1,1
    #

    date_str = (Yini, Mini, Dini)

    Tini = date.toordinal(date(Yini, Mini, Dini)) - date.toordinal(date(Yorig, 1, 1))
    Tini = Tini + 0.5  # 12H
    Tini_str = "%06d" % Tini

    #
    # Get the GLRYS file name from the date (only one time step per file)
    #

    glorysname = glorysfiles_dir + glorys_prefix + ini_date + '_' + ini_date + glorys_ending

    print(glorysname)

    glorysnamebgc = glorysfiles_dir + glorys_bgc_prefix + ini_date + '_' + ini_date + glorys_ending

    print(glorysnamebgc)

    #
    # Title
    #

    print(' ')
    print(' Making initial file: ' + ininame)
    print(' ')
    print(' Title: ' + title)

    #
    # Initial file
    #

    glor.create_inifile(ininame, grdname, title, theta_s, theta_b, hc, N, Tini, vtransform)

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

    lonmin = np.min(lon_rho)
    lonmax = np.max(lon_rho)
    latmin = np.min(lat_rho)
    latmax = np.max(lat_rho)

    #
    # Open the CROCO initial file for writing
    #

    ncin = netcdf(ininame, 'a')
    ncini = ncin.variables

    #
    # get a GLORYS subgrid
    #

    ncphy = netcdf(glorysname, 'r')
    ncphyo = ncphy.variables
    depth = np.array(ncphyo['depth'][:])

    ncbgc = netcdf(glorysnamebgc, 'r')
    ncbgco = ncbgc.variables

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncphyo['latitude'][:])
    lonT = np.array(ncphyo['longitude'][:])

    iminT = glor.geo_idx(lonmin - 1, lonT)
    imaxT = glor.geo_idx(lonmax + 1, lonT)
    jminT = glor.geo_idx(latmin - 1, latT)
    jmaxT = glor.geo_idx(latmax + 1, latT)

    lonT = lonT[iminT:imaxT]
    latT = latT[jminT:jmaxT]
    (LonT, LatT) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncphyo['latitude'][:])
    lonU = np.array(ncphyo['longitude'][:])

    iminU = glor.geo_idx(lonmin - 1, lonU)
    imaxU = glor.geo_idx(lonmax + 1, lonU)
    jminU = glor.geo_idx(latmin - 1, latU)
    jmaxU = glor.geo_idx(latmax + 1, latU)

    lonU = lonU[iminU:imaxU]
    latU = latU[jminU:jmaxU]
    (LonU, LatU) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncphyo['latitude'][:])
    lonV = np.array(ncphyo['longitude'][:])

    iminV = glor.geo_idx(lonmin - 1, lonV)
    imaxV = glor.geo_idx(lonmax + 1, lonV)
    jminV = glor.geo_idx(latmin - 1, latV)
    jmaxV = glor.geo_idx(latmax + 1, latV)

    lonV = lonV[iminV:imaxV]
    latV = latV[jminV:jmaxV]
    (LonV, LatV) = np.meshgrid(lonV, latV)

    [Nz] = np.shape(depth)

    #
    # Horizontal and vertical interp/extrapolations
    #
    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        #
        # Compute the Delaunay triangulation matrices (long but only done once)
        #
        # (u and v are interpolated on croco rho_points because we may need to rotate them)
        #

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT, coefT] = glor.get_tri_coef(LonT, LatT, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefT, axis=2)
        coefT = coefT / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU, coefU] = glor.get_tri_coef(LonU, LatU, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefU, axis=2)
        coefU = coefU / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV, coefV] = glor.get_tri_coef(LonV, LatV, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefV, axis=2)
        coefV = coefV / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs.npz', coefT=coefT, elemT=elemT,
                 coefU=coefU, elemU=elemU, coefV=coefV, elemV=elemV)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs.npz')
        coefT = data['coefT']
        elemT = data['elemT']
        coefU = data['coefU']
        elemU = data['elemU']
        coefV = data['coefV']
        elemV = data['elemV']

    print('Delaunay triangulation done')

    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta, NzGood) = glor.interp_tracers(ncphyo, 'zos', tndx_glo, -1, iminT, imaxT, jminT, jmaxT, LonT, LatT, coefT,
                                         elemT)

    ncini['zeta'][tndx_ini, :, :] = zeta

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h, zeta, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h, zeta, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp = glor.interp3d(ncphyo, 'thetao', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['temp'][tndx_ini, :, :, :] = temp

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt = glor.interp3d(ncphyo, 'so', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['salt'][tndx_ini, :, :, :] = salt


    #
    #
    # 3: Nitrate
    #
    #

    print('Interpolate Nitrate...')

    no3 = glor.interp3d(ncbgco, 'no3', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['NO3'][tndx_ini, :, :, :] = no3


    #
    #
    # 3: Ammonium
    #
    #

    print('Interpolate Ammonium...')

    nh4 = glor.interp3d(ncbgco, 'nh4', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['NH4'][tndx_ini, :, :, :] = nh4


    #
    #
    # 3: Orthophosphate
    #
    #

    print('Interpolate Orthophosphate...')

    po4 = glor.interp3d(ncbgco, 'po4', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['PO4'][tndx_ini, :, :, :] = po4


    #
    #
    # 3: Silicate
    #
    #

    print('Interpolate Silicate...')

    si = glor.interp3d(ncbgco, 'si', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['Si'][tndx_ini, :, :, :] = si


    #
    #
    # 3: Iron
    #
    #

    print('Interpolate Iron...')

    fe = glor.interp3d(ncbgco, 'fe', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['FER'][tndx_ini, :, :, :] = fe


    #
    #
    # 3: Dissolved Oxygen
    #
    #

    print('Interpolate Dissolved Oxygen...')

    o2 = glor.interp3d(ncbgco, 'o2', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['O2'][tndx_ini, :, :, :] = o2

    #
    #
    # 3: Dissolved Oxygen
    #
    #

    print('Interpolate Dissolved Oxygen...')

    dissic = glor.interp3d(ncbgco, 'dissic', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['DIC'][tndx_ini, :, :, :] = dissic


    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle)
    sina = np.sin(angle)

    [u, v] = glor.interp3d_uv(ncphyo, tndx_glo, Nzgoodmin, depth, z_rho, cosa, sina,
                              iminU, imaxU, jminU, jmaxU, LonU, LatU, coefU, elemU,
                              iminV, imaxV, jminV, jmaxV, LonV, LatV, coefV, elemV)

    ncini['u'][tndx_ini, :, :, :] = u
    ncini['v'][tndx_ini, :, :, :] = v

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar, h0) = vgrd.vintegr(u, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h)
    (vbar, h0) = vgrd.vintegr(v, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h)

    ncini['ubar'][tndx_ini, :, :] = ubar
    ncini['vbar'][tndx_ini, :, :] = vbar

    ncphy.close()
    ncin.close()

    print('Figure')
    plt.contourf(lon_rho, lat_rho, temp[N - 1, :, :])
    plt.show()

#  return

#
# End
#
######################################################################

# main_func()
