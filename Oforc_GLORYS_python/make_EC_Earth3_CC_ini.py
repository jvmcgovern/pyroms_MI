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
# #####################################################################
# #####################################################################
#
# PHY variables
#
# 'vo_Omon_EC-Earth3-CC_ssp245_r1i1p1f1_gn_201501-201512.nc'
# variable_frq_esm_case_exp_gn_period.nc

import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset as netcdf
from scipy.interpolate import griddata

from datetime import date
import sys

# sys.path.insert(0,'')

from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_EC_Earth3_CC as ece3cc
from progressbar import *

# ESM name
esm = 'EC-Earth3-CC'

# Data frequency
frq = 'Omon'  # monthly
# frq = 'Oyr'   # annual

# Scenario case
# case = 'esm-hist'  # or 'historical'
case = 'ssp245'
# case = 'ssp585'

# Experiment ID
exp = 'r1i1p1f1'

# Variable name in file
variables_phy = ['zos', 'thetao', 'so', 'uo', 'vo']

# #
# # BGC variables
# #
# frq = 'Omon'
# frq = 'Oyr'
#
variables_bgc = ['no3', 'zmicro', 'zmeso', 'talk', 'spco2', 'fgco2', 'si', 'po4', 'phydiat', 'ph', 'o2',  'nh4',
                 'dissoc', 'dissic', 'dfe', 'chl', 'calc']

vari_bgc_ssp_Omon = ['no3', 'nh4', 'si', 'dfe', 'dissic', 'dissoc', 'spco2', 'fgco2']

vari_bgc_ssp_Oyr = ['po4', 'o2', 'zmicro', 'zmeso', 'chl', 'calc', 'talk', 'phydiat', 'ph']

# #
# # ATM variables
# #
# frq = 'Amon'
# frq = 'day'
#
# variables_atm = {'fco2nat', 'hfls', 'hfss', 'hur', 'hurs', 'hursmax', 'hursmin', 'pr', 'psl', 'rlds', 'rldscs',
# 'rlus', 'rlut', 'rlutcs', 'rsds', 'rsdscs', 'rsdt', 'rsus', 'rsuscs', 'rsut', 'rsutcs', 'sfcWind', 'sfcWindmax',
# 'ta', 'tas', 'tasmax', 'tasmin', 'tauu', 'tauv', 'tos', 'ts', 'ua', 'uas', 'va', 'vas'}

if 1 == 1:
    # def main_func():

    #
    # #################### USERS DEFINED VARIABLES ########################
    #

    title = 'Initial CROCO-PISCES file from EC-Earth3-CC'

    # Position of different forcing files for scenario modelling
    # crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    crocofiles_dir = '/media/dskthree/EC_Earth3_CC/INTERP/HIST/'
    # crocofiles_dir = '/media/dskthree/EC_Earth3_CC/INTERP/SSP245/'
    # crocofiles_dir = '/media/dskthree/EC_Earth3_CC/INTERP/SSP585/'

    # ininame = crocofiles_dir + 'croco_ini_HIST_1984_CELTIC.nc'  # was hmin = 20m
    # ininame = crocofiles_dir + 'croco_ini_HIST_1984_CELTIC_h6.nc'  # was hmin = 6m
    # ininame = crocofiles_dir + 'croco_ini_HIST_1984_CELTIC_h8.nc'  # was hmin = 8m
    ininame = crocofiles_dir + 'croco_ini_HIST_1984_CELTIC_h8.nc'  # was hmin = 8m
    # grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m

    N = 20
    theta_s = 7.
    theta_b = 0.
    hc = 20.
    vtransform = 2

    Yorig = 1984  # year origin of time : days since Yorig-01-01

    # Yini = 2005
    Yini = 2015
    Mini = 1
    Dini = 1

    # period = str(Yini).zfill(4) + str(Mini).zfill(2)
    period = str(Yini).zfill(4) + '01' + '-' + str(Yini).zfill(4) + '12'

    # Date in form YYYYMMDD
    # ini_date = '20050101'
    ini_date = '19840101'

    # glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    # glorys_bgc_prefix = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    # glorys_ending = '_R20201201_RE01.nc'

    esmphy_dir = '/media/dskthree/EC_Earth3_CC/NATIVE/PHY/'
    esmbgc_dir = '/media/dskthree/EC_Earth3_CC/NATIVE/BGC/'

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation
    comp_delaunay = 1  # 1: compute delaunay triangulations - 0: use saved matrices (for debugging)

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)
    tndx_ini = 0  # time index in the CROCO initial file (should be only 0 here)

    #
    # ################## END USERS DEFINED VARIABLES ######################
    #

    #
    # Get the time in days since Yorig, 1, 1
    #

    date_str = (Yini, Mini, Dini)

    Tini = date.toordinal(date(Yini, Mini, Dini)) - date.toordinal(date(Yorig, 1, 1))
    Tini = Tini + 0.5  # 12H
    Tini_str = "%06d" % Tini

    #
    # Get the GLORYS file name from the date (only one time step per file)
    #
    #
    # glorysname = glorysfiles_dir + glorys_prefix + ini_date + '_' + ini_date + glorys_ending

    esmname = esmphy_dir + variables_phy[1] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
              + period + '.nc'
    #
    # print(esmname)
    #
    # glorysnamebgc = glorysfiles_dir + glorys_bgc_prefix + ini_date + '_' + ini_date + glorys_ending

    esmbgcname = esmbgc_dir + variables_bgc[0] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
                 + period + '.nc'
    #
    # print(esmbgcname)

    #
    # Title
    #

    print(' ')
    print(' Making initial file: ' + ininame)
    print(' ')
    print(' Title: ' + title)

    #
    # Check what BGC variables are available first at desired resolution
    #   Can then tailor init file to not include variables that have not been disseminated via ESGF
    #

    #
    # Initial file
    #

    ece3cc.create_inifile(ininame, grdname, title, theta_s, theta_b, hc, N, Tini, vtransform)

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

    ncphy = netcdf(esmname, 'r')
    ncphyo = ncphy.variables

    depth = np.array(ncphyo['lev'][:])

    ncbgc = netcdf(esmbgcname, 'r')
    ncbgco = ncbgc.variables

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncphyo['latitude'][:])
    lonT = np.array(ncphyo['longitude'][:])

    iminT = ece3cc.geo_idx(lonmin - 1, lonT)
    imaxT = ece3cc.geo_idx(lonmax + 1, lonT)
    jminT = ece3cc.geo_idx(latmin - 1, latT)
    jmaxT = ece3cc.geo_idx(latmax + 1, latT)

    lonT = lonT[iminT:imaxT]
    latT = latT[jminT:jmaxT]
    (LonT, LatT) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncphyo['latitude'][:])
    lonU = np.array(ncphyo['longitude'][:])

    iminU = ece3cc.geo_idx(lonmin - 1, lonU)
    imaxU = ece3cc.geo_idx(lonmax + 1, lonU)
    jminU = ece3cc.geo_idx(latmin - 1, latU)
    jmaxU = ece3cc.geo_idx(latmax + 1, latU)

    lonU = lonU[iminU:imaxU]
    latU = latU[jminU:jmaxU]
    (LonU, LatU) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncphyo['latitude'][:])
    lonV = np.array(ncphyo['longitude'][:])

    iminV = ece3cc.geo_idx(lonmin - 1, lonV)
    imaxV = ece3cc.geo_idx(lonmax + 1, lonV)
    jminV = ece3cc.geo_idx(latmin - 1, latV)
    jmaxV = ece3cc.geo_idx(latmax + 1, latV)

    lonV = lonV[iminV:imaxV]
    latV = latV[jminV:jmaxV]
    (LonV, LatV) = np.meshgrid(lonV, latV)

    ncphy.close()

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
        [elemT, coefT] = ece3cc.get_tri_coef(LonT, LatT, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefT, axis=2)
        coefT = coefT / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU, coefU] = ece3cc.get_tri_coef(LonU, LatU, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefU, axis=2)
        coefU = coefU / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV, coefV] = ece3cc.get_tri_coef(LonV, LatV, lon_rho, lat_rho, 1)
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

    esmnamez = esmphy_dir + '_' + variables_phy[0] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
              + period + '.nc'
    print(esmnamez)
    ncphyz = netcdf(esmnamez, 'r')
    ncphyzo = ncphyz.variables

    print('Interpolate SSH...')

    (zeta, NzGood) = ece3cc.interp_tracers(ncphyzo, 'zos', tndx_glo, -1, iminT, imaxT, jminT, jmaxT, LonT, LatT, coefT,
                                           elemT)

    ncini['zeta'][tndx_ini, :, :] = zeta

    ncphyz.close()

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

    esmnamet = esmphy_dir + '_' + variables_phy[1] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
              + period + '.nc'
    print(esmnamet)

    ncphyt = netcdf(esmnamet, 'r')
    ncphyto = ncphyt.variables

    print('Interpolate Temperature...')

    temp = ece3cc.interp3d(ncphyto, 'thetao', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                           coefT, elemT)

    ncini['temp'][tndx_ini, :, :, :] = temp

    ncphyt.close()

    #
    #
    # 3: Salinity
    #
    #

    esmnames = esmphy_dir + '_' + variables_phy[2] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnames)

    ncphys = netcdf(esmnames, 'r')
    ncphyso = ncphys.variables

    print('Interpolate Salinity...')

    salt = ece3cc.interp3d(ncphyso, 'so', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                           coefT, elemT)

    ncini['salt'][tndx_ini, :, :, :] = salt

    ncphys.close()

    #
    #
    # 3: Nitrate
    #
    #

    esmnamen1 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[0] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnamen1)

    ncbgcn1 = netcdf(esmnamen1, 'r')
    ncbgcn1o = ncbgcn1.variables

    print('Interpolate Nitrate...')

    no3 = ece3cc.interp3d(ncbgcn1o, 'no3', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                          coefT, elemT)

    ncini['NO3'][tndx_ini, :, :, :] = no3

    ncbgcn1.close()

    #
    #
    # 3: Ammonium
    #
    #

    esmnamen2 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[1] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnamen2)

    ncbgcn2 = netcdf(esmnamen2, 'r')
    ncbgcn2o = ncbgcn2.variables

    print('Interpolate Ammonium...')

    nh4 = ece3cc.interp3d(ncbgcn2o, 'nh4', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                          coefT, elemT)

    ncini['NH4'][tndx_ini, :, :, :] = nh4

    ncbgcn2.close()

    # #
    # #
    # # 3: Orthophosphate
    # #
    # #
    #
    # esmnamep1 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[2] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
    #            + period + '.nc'
    # print(esmnamep1)
    # ncbgcp1 = netcdf(esmnamep1, 'r')
    # ncbgcp1o = ncbgcp1.variables
    #
    # print('Interpolate Orthophosphate...')
    #
    # po4 = ece3cc.interp3d(ncbgcp1o, 'po4', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
    #                       coefT, elemT)
    #
    # ncini['PO4'][tndx_ini, :, :, :] = po4
    #
    # ncbgcp1.close()

    #
    #
    # 3: Silicate
    #
    #

    esmnames1 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[2] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnames1)

    ncbgcs1 = netcdf(esmnames1, 'r')
    ncbgcs1o = ncbgcs1.variables

    print('Interpolate Silicate...')

    si = ece3cc.interp3d(ncbgcs1o, 'si', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['Si'][tndx_ini, :, :, :] = si

    ncbgcs1.close()

    #
    #
    # 3: Iron
    #
    #

    esmnamef1 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[3] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnamef1)

    ncbgcf1 = netcdf(esmnamef1, 'r')
    ncbgcf1o = ncbgcf1.variables

    print('Interpolate Iron...')

    fe = ece3cc.interp3d(ncbgcf1o, 'fe', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT)

    ncini['FER'][tndx_ini, :, :, :] = fe

    ncbgcf1.close()

    # #
    # #
    # # 3: Dissolved Oxygen
    # #
    # #
    #
    # esmnameo1 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[2] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
    #            + period + '.nc'
    # print(esmnameo1)
    #
    # ncbgco1 = netcdf(esmnameo1, 'r')
    # ncbgco1o = ncbgco1.variables
    #
    # print('Interpolate Dissolved Oxygen...')
    #
    # o2 = ece3cc.interp3d(ncbgco1o, 'o2', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
    #                      coefT, elemT)
    #
    # ncini['O2'][tndx_ini, :, :, :] = o2
    #
    # ncbgco1.close()

    #
    #
    # 3: Dissolved Inorganic Carbon
    #
    #

    esmnamed1 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[4] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnamed1)

    ncbgcd1 = netcdf(esmnamed1, 'r')
    ncbgcd1o = ncbgcd1.variables

    print('Interpolate Dissolved Inorganic Carbon...')

    dissic = ece3cc.interp3d(ncbgcd1o, 'dissic', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT,
                             LatT,
                             coefT, elemT)

    ncini['DIC'][tndx_ini, :, :, :] = dissic

    ncbgcd1.close()

    #
    #
    # 3: Dissolved Organic Carbon
    #
    #

    esmnamed2 = esmbgc_dir + '_' + vari_bgc_ssp_Omon[5] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
                + period + '.nc'
    print(esmnamed2)

    ncbgcd2 = netcdf(esmnamed2, 'r')
    ncbgcd2o = ncbgcd2.variables

    print('Interpolate Dissolved Organic Carbon...')

    dissoc = ece3cc.interp3d(ncbgcd1o, 'dissoc', tndx_glo, Nzgoodmin, depth, z_rho, iminT, imaxT, jminT, jmaxT, LonT,
                             LatT,
                             coefT, elemT)

    ncini['DOC'][tndx_ini, :, :, :] = dissoc

    ncbgcd2.close()

    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    esmnameu1 = esmphy_dir + '_' + variables_phy[2] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnameu1)

    ncphyu1 = netcdf(esmnameu1, 'r')
    ncphyu1o = ncphyu1.variables

    esmnamev1 = esmphy_dir + '_' + variables_phy[3] + '_' + frq + '_' + esm + '_' + case + '_' + exp + '_gn_' \
               + period + '.nc'
    print(esmnamev1)

    ncphyv1 = netcdf(esmnamev1, 'r')
    ncphyv1o = ncphyv1.variables

    cosa = np.cos(angle)
    sina = np.sin(angle)

    [u, v] = ece3cc.interp3d_uv(ncphyo, tndx_glo, Nzgoodmin, depth, z_rho, cosa, sina,
                                iminU, imaxU, jminU, jmaxU, LonU, LatU, coefU, elemU,
                                iminV, imaxV, jminV, jmaxV, LonV, LatV, coefV, elemV)

    ncini['u'][tndx_ini, :, :, :] = u
    ncini['v'][tndx_ini, :, :, :] = v

    ncphyu1.close()
    ncphyv1.close()

    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #

    (ubar, h0) = vgrd.vintegr(u, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h)
    (vbar, h0) = vgrd.vintegr(v, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h)

    ncini['ubar'][tndx_ini, :, :] = ubar
    ncini['vbar'][tndx_ini, :, :] = vbar

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
