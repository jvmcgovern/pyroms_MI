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
import cftime
import numpy as np
import matplotlib.pyplot as plt
import glob
from netCDF4 import Dataset as netcdf
from scipy.interpolate import griddata
from netCDF4 import date2index as d2i
from datetime import date, datetime
import PyCO2SYS as pyco2
import sys

# sys.path.insert(0,'')

from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from progressbar import *

if 1 == 1:
    # def main_func():

    #
    # #################### USERS DEFINED VARIABLES ########################
    #

    title = 'Initial file using GLORYS'

    NIVAvars = ['no3no2', 'po4', 'si', 'o2', 'dissic', 'talk']

    N = 20
    theta_s = 7.
    theta_b = 0.
    # hc = 20.
    hc = 50.
    vtransform = 2

    Yorig = 1990  # year origin of time : days since Yorig-01-01

    # Yini = 2017
    Yini = 1993
    Mini = 1
    Dini = 1

    # crocofiles_dir = my_home_dir + 'SWAG/Run_TEST/CROCO_FILES/'

    crocofiles_bcd_dir = '/media/dskfour/CS1KM_19932021_BCd/'
    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    ininame = crocofiles_bcd_dir + 'croco_ini_my_bc_MERCATOR_hc50m_MSL_ng' + '_Y' + str(Yini) + 'M' + str(Mini).zfill(
        2) + '.nc'
    # ininame = crocofiles_dir + 'croco_ini_MERCATOR_Y2017M11.nc'  # was hmin = 20m
    # ininame = crocofiles_dir + 'croco_ini_PHYBIO_CELTIC_h6.nc'  # was hmin = 6m
    # ininame = crocofiles_dir + 'croco_ini_PHYBIO_CELTIC_h8.nc'  # was hmin = 8m
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    # grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m
    grdname = crocofiles_dir + 'CELTIC_grd_v5_MSL_h10.nc'  # LAT to MSL included, INFOMAR block averaged hmin = 10m
    # grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m

    # Date in form YYYYMMDD
    # ini_date = '20050101'
    # ini_date = '20180101'
    # ini_date = '20171101'
    ini_date = str(Yini) + str(Mini).zfill(2) + str(Dini).zfill(2)

    # # CMEMS_IBI MY (Reanalysis data)
    glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    glorys_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    glorys_ending = '_R*_RE*.nc'

    # # CMEMS_GLO NRT (Near-real time data)
    # glorysfiles_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    # glorys_prefix = 'CMEMS_v5r1_IBI_PHY_NRT_PdE_01dav_'
    # glorys_ending = '_R*_AN*.nc'

    # glorys_bgc_prefix = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_'
    # glorysfiles_dir = '/media/dskthree/CMEMS_GLO/'
    # glorys_prefix = 'mercatorglorys12v1_gl12_mean_'
    # glorys_ending = '_R*_crop.nc'

    NORESMfiles_dir = '/media/dskone/CELTIC/NORESM/NORESMOCv1p2_C2CMI_1993_2020_NutO2CO2/'
    NORESM_prefix = 'NORESMOCv1p2_C2CMI_1993_2020_'
    NORESM_ending = '_reanal.nc'

    # lat
    # lon

    # Sample file name (note hyphen), and run stamp at end of name e.g. RYYYYMMDD.nc
    # daily: 'IBI36_cg_1d-m_20170605-20170605_3DT-bgc_hcst_R20170614.nc'
    # monthly: 'IBI36_cg_1m-m_201301_3DT-bgc_hcst.nc'

    # PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI_PISCESvars_012013_202112/monthly/'
    # PISCES24_prefix = 'IBI36_cg_1m-m_'
    # PISCES24_ending = '_3DT-bgc_hcst.nc'

    PISCES24files_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    PISCES24_prefix = 'CMEMS_v5r1_IBI_BIO_MY_PdE_01mav_'
    PISCES24_ending = '_R*_RE01.nc'

    broudir = '/media/dskone/VAL/BroullonClims/'
    brou_dic = 'TCO2_NNGv2LDEO_climatology_WOA18.nc'
    brou_tlk = 'AT_NNGv2_climatology_WOA18.nc'

    woadir = '/media/dskone/VAL/WOA18/'
    woa_sta = 'woa18_all_'
    woa_end = '_01.nc'

    # nav_lat
    # nav_lon

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation
    comp_delaunay = 1  # 1: compute delaunay triangulations - 0: use saved matrices (for debugging)

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)
    tndx_ini = 0  # time index in the CROCO initial file (should be only 0 here)

    #
    # ################## END USERS DEFINED VARIABLES ######################
    #

    #
    # Get the time in days since Yorig,1,1
    #

    date_str = (Yini, Mini, Dini)

    Tini = date.toordinal(date(Yini, Mini, Dini)) - date.toordinal(date(Yorig, 1, 1))
    # MERCATOR IBI DAILY MODEL TIMESTAMPED AT NOON, LEAVING INITIAL TIME TO MIDNIGHT TO ALLOW FOR
    # VALIDATION OF SST FROM HISTORY FILES
    # Tini = Tini + 0.5  # 12H
    Tini_str = "%06d" % Tini

    #
    # Get the GLORYS file name from the date (only one time step per file)
    #

    glorysname = glorysfiles_dir + glorys_prefix + ini_date + '_' + ini_date + glorys_ending

    glofile = sorted(glob.glob(glorysname))
    print(glorysname)

    NORESMnamebgc_ex = NORESMfiles_dir + NORESM_prefix + 'no3no2' + NORESM_ending

    print(NORESMnamebgc_ex)

    PISCESnamebgc = PISCES24files_dir + PISCES24_prefix + str(Yini) + str(Mini).zfill(2) + '*_*' + PISCES24_ending
    bgcfile = sorted(glob.glob(PISCESnamebgc))

    print(PISCESnamebgc)

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

    glor.create_ini_PHY_BGC(ininame, grdname, title, theta_s, theta_b, hc, N, Tini, vtransform)

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

    ncphy = netcdf(glofile[0], 'r')
    ncphyo = ncphy.variables
    depthg = np.array(ncphyo['depth'][:])

    ncpis = netcdf(bgcfile[0], 'r')
    ncpiso = ncpis.variables
    depthp = np.array(ncpiso['depth'][:])

    # WOA file sorting and opening
    woan_cur = netcdf(woadir + woa_sta + 'n' + str(Mini).zfill(2) + woa_end, 'r')
    woan_curo = woan_cur.variables
    depthnps = np.array(woan_curo['depth'][:])
    woap_cur = netcdf(woadir + woa_sta + 'p' + str(Mini).zfill(2) + woa_end, 'r')
    woap_curo = woap_cur.variables
    woas_cur = netcdf(woadir + woa_sta + 'i' + str(Mini).zfill(2) + woa_end, 'r')
    woas_curo = woas_cur.variables
    woao_cur = netcdf(woadir + woa_sta + 'o' + str(Mini).zfill(2) + woa_end, 'r')
    woao_curo = woao_cur.variables
    deptho = np.array(woao_curo['depth'][:])

    # Open Broullon file
    ncbroud = netcdf(broudir + brou_dic, 'r')
    ncbroudo = ncbroud.variables
    ncbrout = netcdf(broudir + brou_tlk, 'r')
    ncbrouto = ncbrout.variables
    depthb = np.array(ncbroudo['depth'][:])

    # Nitrate from NORESM
    ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
    ncnornio = ncnorni.variables
    depthn = np.array(ncnornio['depth'][:])
    Tinin = datetime(Yini, Mini, Dini)
    Tininx1 = d2i(Tinin, ncnornio['time'], select='before')
    Tininx2 = d2i(Tinin, ncnornio['time'], select='after')

    # Phosphate from NORESM
    ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
    ncnorpoo = ncnorpo.variables

    # Silicate from NORESM
    ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
    ncnorsio = ncnorsi.variables

    # Dissolved Oxygen from NORESM
    ncnorox = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[3] + NORESM_ending, 'r')
    ncnoroxo = ncnorox.variables

    # Dissolved Inorganic Carbon from NORESM
    ncnordc = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[4] + NORESM_ending, 'r')
    ncnordco = ncnordc.variables

    # Total Alkalinity from NORESM
    ncnoral = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[5] + NORESM_ending, 'r')
    ncnoralo = ncnoral.variables

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

    [Nz] = np.shape(depthg)

    #
    # get PISCES positions and indices
    #

    latP = np.array(ncpiso['latitude'][:])
    lonP = np.array(ncpiso['longitude'][:])

    iminP = glor.geo_idx(lonmin - 1, lonP)
    imaxP = glor.geo_idx(lonmax + 1, lonP)
    jminP = glor.geo_idx(latmin - 1, latP)
    jmaxP = glor.geo_idx(latmax + 1, latP)

    lonP = lonP[iminP:imaxP]
    latP = latP[jminP:jmaxP]
    (LonP, LatP) = np.meshgrid(lonP, latP)

    #
    # get WOA positions and indices
    #

    latW = np.array(woan_curo['lat'][:])
    lonW = np.array(woan_curo['lon'][:])

    iminW = glor.geo_idx(lonmin - 2, lonW)
    imaxW = glor.geo_idx(lonmax + 2, lonW)
    jminW = glor.geo_idx(latmin - 2, latW)
    jmaxW = glor.geo_idx(latmax + 2, latW)

    lonW = lonW[iminW:imaxW]
    latW = latW[jminW:jmaxW]
    (LonW, LatW) = np.meshgrid(lonW, latW)

    #
    # get BROULLON positions and indices
    #

    latB = np.array(ncbroudo['latitude'][:])
    lonB = np.array(ncbroudo['longitude'][:])

    iminB = glor.geo_idx(lonmin - 2, lonB)
    imaxB = glor.geo_idx(lonmax + 2, lonB)
    jminB = glor.geo_idx(latmin - 2, latB)
    jmaxB = glor.geo_idx(latmax + 2, latB)

    lonB = lonB[iminB:imaxB]
    latB = latB[jminB:jmaxB]
    (LonB, LatB) = np.meshgrid(lonB, latB)

    #
    # get NORESM positions and indices
    #

    LatN = np.array(ncnornio['lat'][:])
    LonN = np.array(ncnornio['lon'][:])

    iminN = 0
    imaxN = LonN.shape[1]
    jminN = 0
    jmaxN = LatN.shape[0]

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

        print('Compute Delaunay triangulation from PISCES points to CROCO rho_points...')
        [elemP, coefP] = glor.get_tri_coef(LonP, LatP, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefP, axis=2)
        coefP = coefP / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from WOA points to CROCO rho_points...')
        [elemW, coefW] = glor.get_tri_coef(LonW, LatW, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefP, axis=2)
        coefW = coefW / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from BROULLON points to CROCO rho_points...')
        [elemB, coefB] = glor.get_tri_coef(LonB, LatB, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefB, axis=2)
        coefB = coefB / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from NORESM points to CROCO rho_points...')
        [elemN, coefN] = glor.get_tri_coef(LonN, LatN, lon_rho, lat_rho, 1)
        coefnorm = np.sum(coefN, axis=2)
        coefN = coefN / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs.npz',
                 coefT=coefT, elemT=elemT,
                 coefU=coefU, elemU=elemU,
                 coefV=coefV, elemV=elemV,
                 coefP=coefP, elemP=elemP,
                 coefW=coefW, elemW=elemW,
                 coefB=coefB, elemB=elemB,
                 coefN=coefN, elemN=elemN)

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
        coefP = data['coefP']
        elemP = data['elemP']
        coefW = data['coefW']
        elemW = data['elemW']
        coefB = data['coefB']
        elemB = data['elemB']
        coefN = data['coefN']
        elemN = data['elemN']

    print('Delaunay triangulation done')

    #
    # 1: SSH
    #

    print('Interpolate SSH...')

    (zeta, NzGood) = glor.interp_tracers(ncphyo, 'zos', tndx_glo, -1, iminT, imaxT, jminT, jmaxT, LonT, LatT, coefT,
                                         elemT, '_FillValue')

    ncini['zeta'][tndx_ini, :, :] = zeta

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h, zeta, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h, zeta, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    # 2: Temperature
    #

    print('Interpolate Temperature...')

    temp = glor.interp3d(ncphyo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT,  '_FillValue')

    ncini['temp'][tndx_ini, :, :, :] = temp

    #
    # 3: Salinity
    #

    print('Interpolate Salinity...')

    salt = glor.interp3d(ncphyo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho, iminT, imaxT, jminT, jmaxT, LonT, LatT,
                         coefT, elemT, '_FillValue')

    ncini['salt'][tndx_ini, :, :, :] = salt

    #
    # 4a: DIC from BROULLON
    #
    print('Interpolate Dissolved Inorganic Carbon...')

    dissic = glor.interp3d(ncbroudo, 'TCO2_NNGv2LDEO', Mini - 1, Nzgoodmin, depthb, z_rho, iminB, imaxB, jminB, jmaxB,
                           LonB, LatB, coefB, elemB, [])

    ncini['DIC'][tndx_ini, :, :, :] = dissic*1.025

    #
    # 4a: TALK from BROULLON
    #
    print('Interpolate Total Alkalinity...')

    talk = glor.interp3d(ncbrouto, 'AT_NNGv2', Mini - 1, Nzgoodmin, depthb, z_rho, iminB, imaxB, jminB, jmaxB,
                         LonB, LatB, coefB, elemB, [])

    ncini['TALK'][tndx_ini, :, :, :] = talk*1.025

    # #
    # # 4b: Alk from CMEMS IBI
    # #
    #
    # print('Interpolate pH...')
    #
    # ph = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho, iminP, imaxP, jminP, jmaxP,
    #                    LonP, LatP, coefP, elemP, '_FillValue')

    #
    # 4c: O2 from CMEMS IBI
    #

    print('Interpolate Dissolved Oxygen...')

    o2 = glor.interp3d(woao_curo, 'o_an', tndx_glo, Nzgoodmin, deptho, z_rho, iminW, imaxW, jminW, jmaxW,
                       LonW, LatW, coefW, elemW, '_FillValue')

    ncini['O2'][tndx_ini, :, :, :] = o2

    #
    # 4d: PO4 from CMEMS IBI
    #

    print('Interpolate Dissolved Phosphate...')

    po4 = glor.interp3d(woap_curo, 'p_an', tndx_glo, Nzgoodmin, depthnps, z_rho, iminW, imaxW, jminW, jmaxW,
                        LonW, LatW, coefW, elemW, '_FillValue')

    ncini['PO4'][tndx_ini, :, :, :] = po4

    #
    # 4e: NO3 from CMEMS IBI
    #

    print('Interpolate Nitrate...')

    no3 = glor.interp3d(woan_curo, 'n_an', tndx_glo, Nzgoodmin, depthnps, z_rho, iminW, imaxW, jminW, jmaxW,
                       LonW, LatW, coefW, elemW, '_FillValue')

    ncini['NO3'][tndx_ini, :, :, :] = no3

    #
    # 4f: Si from CMEMS IBI
    #

    print('Interpolate Silicate...')

    si = glor.interp3d(woas_curo, 'i_an', tndx_glo, Nzgoodmin, depthnps, z_rho, iminW, imaxW, jminW, jmaxW,
                       LonW, LatW, coefW, elemW, '_FillValue')

    ncini['Si'][tndx_ini, :, :, :] = si

    #
    # 5i: Dissolved Iron from IBI PISCES
    #

    print('Interpolate Dissolved Iron...')

    fe = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho, iminP, imaxP, jminP, jmaxP,
                       LonP, LatP, coefP, elemP, '_FillValue')

    ncini['FER'][tndx_ini, :, :, :] = fe

    #
    # 5r: Ammonium from IBI PISCES
    #

    print('Interpolate Ammonium...')

    nh4 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho, iminP, imaxP, jminP, jmaxP,
                        LonP, LatP, coefP, elemP, '_FillValue')

    ncini['NH4'][tndx_ini, :, :, :] = nh4

    print('Interpolate Total Alkalinity...')

    # results = pyco2.sys(par1=ph, par1_type=3,
    #                     par2=dissic*1000, par2_type=2,
    #                     temperature=temp,
    #                     salinity=salt,
    #                     opt_pH_scale=1, opt_k_carbonic=4,
    #                     opt_k_bisulfate=1, opt_total_borate=1)
    # ncini['TALK'][tndx_ini, :, :, :] = results["alkalinity"]

    #
    # 6a: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #

    cosa = np.cos(angle)
    sina = np.sin(angle)

    [u, v] = glor.interp3d_uv(ncphyo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                              iminU, imaxU, jminU, jmaxU, LonU, LatU, coefU, elemU,
                              iminV, imaxV, jminV, jmaxV, LonV, LatV, coefV, elemV)

    ncini['u'][tndx_ini, :, :, :] = u
    ncini['v'][tndx_ini, :, :, :] = v

    #
    # 6b: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #

    (ubar, h0) = vgrd.vintegr(u, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h)
    (vbar, h0) = vgrd.vintegr(v, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h)

    ncini['ubar'][tndx_ini, :, :] = ubar
    ncini['vbar'][tndx_ini, :, :] = vbar

    ncphy.close()
    ncpis.close()
    ncnorni.close()
    ncnorpo.close()
    ncnorox.close()
    ncnorsi.close()
    ncnordc.close()
    ncnoral.close()
    ncin.close()

    # print('Figure')
    # plt.contourf(lon_rho, lat_rho, temp[N - 1, :, :])
    # plt.show()

#  return

#
# End
#
######################################################################

# main_func()
