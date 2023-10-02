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
import glob
import cftime
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import griddata
from scipy import interpolate
from netCDF4 import date2index as d2i
from netCDF4 import num2date as n2d
from datetime import date, datetime
from calendar import monthrange
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
    hc = 20.
    vtransform = 2

    Yorig = 1990  # year origin of time : days since Yorig-01-01

    Yini = 2017
    Mini = 1
    Dini = 1

    Ystart = 2017
    Mstart = 1
    Dstart = 1

    Yend = 2017
    Mend = 11
    Dend = 30

    # crocofiles_dir = my_home_dir + 'SWAG/Run_TEST/CROCO_FILES/'
    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    # grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m

    # Date in form YYYYMMDD
    # ini_date = '20050101'
    # ini_date = '20180101'
    # ini_date = '20171101'

    IBI_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    IBI_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    IBI_post = '_R20201201_RE01.nc'

    GLO_dir = '/media/dskthree/CMEMS_GLO/'
    GLO_prefix = 'mercatorglorys12v1_gl12_mean_'
    GLO_post = '_R*_crop.nc'

    hisdir_orig = '/media/dskone/CELTIC/CROCO_PHY_1p1_1719/'
    hisdir_J3 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_J3/'
    hisdir_FO = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_Fairall_Off/'
    hisdir_J5_GLS = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_J5_GLS/'
    hisdir_TSUV = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV/'
    hisdir_GLO = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_GLO/'
    hisdir_GLO_2 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_GLO_II/'
    hisdir_GLO_JS = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_GLO_JS/'
    hisdir_IBI_JS = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_JS/'

    his_prefix = 'croco_his_'

    sat_sstfiles_dir = '/media/dskone/VAL/TEMP/SST_L4/'
    # sat_sst_ending_core_ex = YYYYMMDD
    sat_sst_ending = '000000-IFR-L4_GHRSST-SSTfnd-ODYSSEA-ATL_005-v2.0-fv1.0.nc'

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation
    comp_delaunay = 1  # 1: compute delaunay triangulations - 0: use saved matrices (for debugging)

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    #
    # ################## END USERS DEFINED VARIABLES ######################
    #

    #
    # Get the time in days since Yorig,1,1
    #

    Tini = date.toordinal(date(Yini, Mini, Dini)) - date.toordinal(date(Yorig, 1, 1))
    Tini = Tini + 0.5  # 12H
    Tini_str = "%06d" % Tini

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
    # get a satellite SST subgrid
    #
    sat_sst_files = sorted(glob.glob(sat_sstfiles_dir + '????????' + sat_sst_ending))

    sat_sstex = sat_sst_files[0]

    ncphy = netcdf(sat_sstex, 'r')
    ncphyo = ncphy.variables

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncphyo['lat'][:])
    lonT = np.array(ncphyo['lon'][:])

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

    latU = np.array(ncphyo['lat'][:])
    lonU = np.array(ncphyo['lon'][:])

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

    latV = np.array(ncphyo['lat'][:])
    lonV = np.array(ncphyo['lon'][:])

    iminV = glor.geo_idx(lonmin - 1, lonV)
    imaxV = glor.geo_idx(lonmax + 1, lonV)
    jminV = glor.geo_idx(latmin - 1, latV)
    jmaxV = glor.geo_idx(latmax + 1, latV)

    lonV = lonV[iminV:imaxV]
    latV = latV[jminV:jmaxV]
    (LonV, LatV) = np.meshgrid(lonV, latV)

    ncphy.close()

    # Get IBI positions and indices
    ncibi = netcdf(IBI_dir + IBI_prefix + '20170101_20170101' + IBI_post, 'r')
    ncibio = ncibi.variables
    IBIlat = np.array(ncibio['latitude'][:])
    IBIlon = np.array(ncibio['longitude'][:])
    IBIdep = np.array(ncibio['depth'][:])
    ncibi.close()

    # Get GLO positions and indices
    GLOfiles = sorted(glob.glob(GLO_dir + GLO_prefix + '????????' + GLO_post))
    ncglo = netcdf(GLOfiles[0], 'r')
    ncgloo = ncglo.variables
    GLOlat = np.array(ncgloo['latitude'][:])
    GLOlon = np.array(ncgloo['longitude'][:])
    GLOdep = np.array(ncgloo['depth'][:])
    ncglo.close()

    ctdfil = '/media/dskone/CELTIC/IMI_CTD_34f3_19f6_97e1.nc'
    ncctf = netcdf(ctdfil, 'r')
    ncctfo = ncctf.variables

    timec = np.array(ncctfo['time'][:])
    deptc = np.array(ncctfo['depth'][:])
    lat_c = np.array(ncctfo['latitude'][:])
    lon_c = np.array(ncctfo['longitude'][:])
    tempc = np.array(ncctfo['temperature'][:])
    sal_c = np.array(ncctfo['salinity'][:])
    # tempm_or = np.zeros_like(tempc)
    # sal_m_or = np.zeros_like(sal_c)
    # tempm_J3 = np.zeros_like(tempc)
    # sal_m_J3 = np.zeros_like(sal_c)
    # tempm_FO = np.zeros_like(tempc)
    # sal_m_FO = np.zeros_like(sal_c)
    # tempm_J5_GLS = np.zeros_like(tempc)
    # sal_m_J5_GLS = np.zeros_like(sal_c)
    # tempm_TSUV = np.zeros_like(tempc)
    # sal_m_TSUV = np.zeros_like(sal_c)
    # tempm_GLO = np.zeros_like(tempc)
    # sal_m_GLO = np.zeros_like(sal_c)
    # tempm_GL2 = np.zeros_like(tempc)
    # sal_m_GL2 = np.zeros_like(sal_c)
    tempm_GJS = np.zeros_like(tempc)
    sal_m_GJS = np.zeros_like(sal_c)
    tempm_IJS = np.zeros_like(tempc)
    sal_m_IJS = np.zeros_like(sal_c)

    # Get year, month, day from timec
    timec_date = n2d(timec, units=ncctfo['time'].units, calendar=ncctfo['time'].calendar)
    timec_year = np.zeros(timec_date.shape[0], dtype=np.int32)
    timec_month = np.zeros(timec_date.shape[0], dtype=np.int32)
    timec_day = np.zeros(timec_date.shape[0], dtype=np.int32)
    ctd_latidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    ctd_lonidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    ibi_latidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    ibi_lonidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    glo_latidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    glo_lonidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    timec_ymdll = np.zeros([timec_date.shape[0], 9], dtype=np.int32)
    for td in range(0, len(timec_date)):
        timec_year[td] = timec_date[td].year
        timec_ymdll[td, 0] = timec_date[td].year
        timec_month[td] = timec_date[td].month
        timec_ymdll[td, 1] = timec_date[td].month
        timec_day[td] = timec_date[td].day
        timec_ymdll[td, 2] = timec_date[td].day
        ctd_latidx[td] = glor.geo_idx(lat_c[td], lat_rho[:, 0])
        timec_ymdll[td, 3] = glor.geo_idx(lat_c[td], lat_rho[:, 0])
        ctd_lonidx[td] = glor.geo_idx(lon_c[td], lon_rho[0, :])
        timec_ymdll[td, 4] = glor.geo_idx(lon_c[td], lon_rho[0, :])
        ibi_latidx[td] = glor.geo_idx(lat_c[td], IBIlat)
        timec_ymdll[td, 5] = glor.geo_idx(lat_c[td], IBIlat)
        ibi_lonidx[td] = glor.geo_idx(lon_c[td], IBIlon)
        timec_ymdll[td, 6] = glor.geo_idx(lon_c[td], IBIlon)
        glo_latidx[td] = glor.geo_idx(lat_c[td], GLOlat)
        timec_ymdll[td, 7] = glor.geo_idx(lat_c[td], GLOlat)
        glo_lonidx[td] = glor.geo_idx(lon_c[td], GLOlon)
        timec_ymdll[td, 8] = glor.geo_idx(lon_c[td], GLOlon)
    ncctf.close()

    # Get unique dates to cycle through the equivalent model output history files
    # Get unique grid coordinates for lat/lon combinations for each unique date
    timec_dates_lls = np.unique(timec_ymdll, axis=0)
    timec_dates = np.unique(timec_ymdll[:, 0:2], axis=0)
    depcount = np.zeros(timec_dates_lls.shape[0], dtype=np.int32)
    for dl in range(0, len(timec_dates_lls)):
        depths = deptc[(timec_year == timec_dates_lls[dl, 0]) &
                       (timec_month == timec_dates_lls[dl, 1]) &
                       (timec_day == timec_dates_lls[dl, 2]) &
                       (ctd_latidx == timec_dates_lls[dl, 3]) &
                       (ctd_lonidx == timec_dates_lls[dl, 4])] * (-1)
        depcount[dl] = len(depths)

        # his_orig = hisdir_orig + \
        #            his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchisor = netcdf(his_orig, 'r')
        # nchisoro = nchisor.variables
        # # Extract vertical section from history file
        # tempor = np.array(nchisoro['temp']
        #                   [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltor = np.array(nchisoro['salt']
        #                   [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_J3 = hisdir_J3 + \
        #          his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_J3 = netcdf(his_J3, 'r')
        # nchis_J3o = nchis_J3.variables
        # # Extract vertical section from history file
        # tempJ3 = np.array(nchis_J3o['temp']
        #                   [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltJ3 = np.array(nchis_J3o['salt']
        #                   [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_FO = hisdir_FO + \
        #          his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_FO = netcdf(his_FO, 'r')
        # nchis_FOo = nchis_FO.variables
        # # Extract vertical section from history file
        # tempFO = np.array(nchis_FOo['temp']
        #                   [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltFO = np.array(nchis_FOo['salt']
        #                   [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        # his_J5_GLS = hisdir_J5_GLS + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_J5_GLS = netcdf(his_J5_GLS, 'r')
        # nchis_J5_GLSo = nchis_J5_GLS.variables
        # # Extract vertical section from history file
        # tempJ5_GLS = np.array(nchis_J5_GLSo['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltJ5_GLS = np.array(nchis_J5_GLSo['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_TSUV = hisdir_TSUV + \
        #            his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_TSUV = netcdf(his_TSUV, 'r')
        # nchis_TSUVo = nchis_TSUV.variables
        # # Extract vertical section from history file
        # tempTSUV = np.array(nchis_TSUVo['temp']
        #                     [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltTSUV = np.array(nchis_TSUVo['salt']
        #                     [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        # his_GLO = hisdir_GLO + \
        #           his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_GLO = netcdf(his_GLO, 'r')
        # nchis_GLOo = nchis_GLO.variables
        # # Extract vertical section from history file
        # tempGLO = np.array(nchis_GLOo['temp']
        #                    [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltGLO = np.array(nchis_GLOo['salt']
        #                    [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_GLO2 = hisdir_GLO_2 + \
        #           his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_GLO2 = netcdf(his_GLO2, 'r')
        # nchis_GLO2o = nchis_GLO2.variables
        # # Extract vertical section from history file
        # tempGLO2 = np.array(nchis_GLO2o['temp']
        #                     [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltGLO2 = np.array(nchis_GLO2o['salt']
        #                     [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        his_GJS = hisdir_GLO_JS + \
                  his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        nchis_GJS = netcdf(his_GJS, 'r')
        nchis_GJSo = nchis_GJS.variables
        # Extract vertical section from history file
        tempGJS = np.array(nchis_GJSo['temp']
                           [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        saltGJS = np.array(nchis_GJSo['salt']
                           [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        his_IJS = hisdir_IBI_JS + \
                  his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        nchis_IJS = netcdf(his_IJS, 'r')
        nchis_IJSo = nchis_IJS.variables
        # Extract vertical section from history file
        tempIJS = np.array(nchis_IJSo['temp']
                           [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        saltIJS = np.array(nchis_IJSo['salt']
                           [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        # Interpolate section data onto depths from the CTD cast for that date
        z = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], 0, theta_s, theta_b, hc, N, 'w', vtransform)
        z[0] = np.floor(z[0])
        zmin = z.min()
        zmax = z.max()

        # ftor = interpolate.interp1d(z, np.append(tempor[0], tempor))
        # tmtoto = np.zeros_like(depths)
        # tmtoto[(depths >= zmin) & (depths <= zmax)] = ftor(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # ftJ3 = interpolate.interp1d(z, np.append(tempJ3[0], tempJ3))
        # tmtotj = np.zeros_like(depths)
        # tmtotj[(depths >= zmin) & (depths <= zmax)] = ftJ3(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # ftFO = interpolate.interp1d(z, np.append(tempFO[0], tempFO))
        # tmtotf = np.zeros_like(depths)
        # tmtotf[(depths >= zmin) & (depths <= zmax)] = ftFO(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # ftJ5_GLS = interpolate.interp1d(z, np.append(tempJ5_GLS[0], tempJ5_GLS))
        # tmtotg = np.zeros_like(depths)
        # tmtotg[(depths >= zmin) & (depths <= zmax)] = ftJ5_GLS(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # ftTSUV = interpolate.interp1d(z, np.append(tempTSUV[0], tempTSUV))
        # tmtott = np.zeros_like(depths)
        # tmtott[(depths >= zmin) & (depths <= zmax)] = ftTSUV(depths[(depths >= zmin) & (depths <= zmax)])

        # ftGLO = interpolate.interp1d(z, np.append(tempGLO[0], tempGLO))
        # tmtotglo = np.zeros_like(depths)
        # tmtotglo[(depths >= zmin) & (depths <= zmax)] = ftGLO(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # ftGLO2 = interpolate.interp1d(z, np.append(tempGLO2[0], tempGLO2))
        # tmtotglo2 = np.zeros_like(depths)
        # tmtotglo2[(depths >= zmin) & (depths <= zmax)] = ftGLO2(depths[(depths >= zmin) & (depths <= zmax)])
        #
        ftGJS = interpolate.interp1d(z, np.append(tempGJS[0], tempGJS))
        tmtotgjs = np.zeros_like(depths)
        tmtotgjs[(depths >= zmin) & (depths <= zmax)] = ftGJS(depths[(depths >= zmin) & (depths <= zmax)])

        ftIJS = interpolate.interp1d(z, np.append(tempIJS[0], tempIJS))
        tmtotijs = np.zeros_like(depths)
        tmtotijs[(depths >= zmin) & (depths <= zmax)] = ftIJS(depths[(depths >= zmin) & (depths <= zmax)])

        if depths.min() < zmin:
            # tmtoto[(depths < zmin)] = np.nan
            # tmtotj[(depths < zmin)] = np.nan
            # tmtotf[(depths < zmin)] = np.nan
            # tmtotg[(depths < zmin)] = np.nan
            # tmtott[(depths < zmin)] = np.nan
            # tmtotglo[(depths < zmin)] = np.nan
            # tmtotglo2[(depths < zmin)] = np.nan
            tmtotgjs[(depths < zmin)] = np.nan
            tmtotijs[(depths < zmin)] = np.nan
        if depths.max() > zmax:
            # tmtoto[(depths > zmax)] = np.nan
            # tmtotj[(depths > zmax)] = np.nan
            # tmtotf[(depths > zmax)] = np.nan
            # tmtotg[(depths > zmax)] = np.nan
            # tmtott[(depths > zmax)] = np.nan
            # tmtotglo[(depths > zmax)] = np.nan
            # tmtotglo2[(depths > zmax)] = np.nan
            tmtotgjs[(depths > zmax)] = np.nan
            tmtotijs[(depths > zmax)] = np.nan
        # tempm_or[(timec_year == timec_dates_lls[dl, 0]) &
        #          (timec_month == timec_dates_lls[dl, 1]) &
        #          (timec_day == timec_dates_lls[dl, 2]) &
        #          (ctd_latidx == timec_dates_lls[dl, 3]) &
        #          (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtoto
        # tempm_J3[(timec_year == timec_dates_lls[dl, 0]) &
        #          (timec_month == timec_dates_lls[dl, 1]) &
        #          (timec_day == timec_dates_lls[dl, 2]) &
        #          (ctd_latidx == timec_dates_lls[dl, 3]) &
        #          (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotj
        # tempm_FO[(timec_year == timec_dates_lls[dl, 0]) &
        #          (timec_month == timec_dates_lls[dl, 1]) &
        #          (timec_day == timec_dates_lls[dl, 2]) &
        #          (ctd_latidx == timec_dates_lls[dl, 3]) &
        #          (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotf
        # tempm_J5_GLS[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotg
        # tempm_TSUV[(timec_year == timec_dates_lls[dl, 0]) &
        #            (timec_month == timec_dates_lls[dl, 1]) &
        #            (timec_day == timec_dates_lls[dl, 2]) &
        #            (ctd_latidx == timec_dates_lls[dl, 3]) &
        #            (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtott
        # tempm_GLO[(timec_year == timec_dates_lls[dl, 0]) &
        #           (timec_month == timec_dates_lls[dl, 1]) &
        #           (timec_day == timec_dates_lls[dl, 2]) &
        #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        #           (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotglo
        # tempm_GL2[(timec_year == timec_dates_lls[dl, 0]) &
        #           (timec_month == timec_dates_lls[dl, 1]) &
        #           (timec_day == timec_dates_lls[dl, 2]) &
        #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        #           (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotglo2
        tempm_GJS[(timec_year == timec_dates_lls[dl, 0]) &
                  (timec_month == timec_dates_lls[dl, 1]) &
                  (timec_day == timec_dates_lls[dl, 2]) &
                  (ctd_latidx == timec_dates_lls[dl, 3]) &
                  (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotgjs
        tempm_IJS[(timec_year == timec_dates_lls[dl, 0]) &
                  (timec_month == timec_dates_lls[dl, 1]) &
                  (timec_day == timec_dates_lls[dl, 2]) &
                  (ctd_latidx == timec_dates_lls[dl, 3]) &
                  (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotijs

        # fsor = interpolate.interp1d(z, np.append(saltor[0], saltor))
        # smtoto = np.zeros_like(depths)
        # smtoto[(depths >= zmin) & (depths <= zmax)] = fsor(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fsJ3 = interpolate.interp1d(z, np.append(saltJ3[0], saltJ3))
        # smtotj = np.zeros_like(depths)
        # smtotj[(depths >= zmin) & (depths <= zmax)] = fsJ3(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fsFO = interpolate.interp1d(z, np.append(saltFO[0], saltFO))
        # smtotf = np.zeros_like(depths)
        # smtotf[(depths >= zmin) & (depths <= zmax)] = fsFO(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fsJ5_GLS = interpolate.interp1d(z, np.append(saltJ5_GLS[0], saltJ5_GLS))
        # smtotg = np.zeros_like(depths)
        # smtotg[(depths >= zmin) & (depths <= zmax)] = fsJ5_GLS(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fsTSUV = interpolate.interp1d(z, np.append(saltTSUV[0], saltTSUV))
        # smtott = np.zeros_like(depths)
        # smtott[(depths >= zmin) & (depths <= zmax)] = fsTSUV(depths[(depths >= zmin) & (depths <= zmax)])

        # fsGLO = interpolate.interp1d(z, np.append(saltGLO[0], saltGLO))
        # smtotglo = np.zeros_like(depths)
        # smtotglo[(depths >= zmin) & (depths <= zmax)] = fsGLO(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fsGLO2 = interpolate.interp1d(z, np.append(saltGLO2[0], saltGLO2))
        # smtotglo2 = np.zeros_like(depths)
        # smtotglo2[(depths >= zmin) & (depths <= zmax)] = fsGLO2(depths[(depths >= zmin) & (depths <= zmax)])
        #
        fsGJS = interpolate.interp1d(z, np.append(saltGJS[0], saltGJS))
        smtotgjs = np.zeros_like(depths)
        smtotgjs[(depths >= zmin) & (depths <= zmax)] = fsGJS(depths[(depths >= zmin) & (depths <= zmax)])

        fsIJS = interpolate.interp1d(z, np.append(saltIJS[0], saltIJS))
        smtotijs = np.zeros_like(depths)
        smtotijs[(depths >= zmin) & (depths <= zmax)] = fsIJS(depths[(depths >= zmin) & (depths <= zmax)])

        if depths.min() < zmin:
            # smtoto[(depths < zmin)] = np.nan
            # smtotj[(depths < zmin)] = np.nan
            # smtotf[(depths < zmin)] = np.nan
            # smtotg[(depths < zmin)] = np.nan
            # smtott[(depths < zmin)] = np.nan
            # smtotglo[(depths < zmin)] = np.nan
            # smtotglo2[(depths < zmin)] = np.nan
            smtotgjs[(depths < zmin)] = np.nan
            smtotijs[(depths < zmin)] = np.nan
        if depths.max() > zmax:
            # smtoto[(depths > zmax)] = np.nan
            # smtotj[(depths > zmax)] = np.nan
            # smtotf[(depths > zmax)] = np.nan
            # smtotg[(depths > zmax)] = np.nan
            # smtott[(depths > zmax)] = np.nan
            # smtotglo[(depths > zmax)] = np.nan
            # smtotglo2[(depths > zmax)] = np.nan
            smtotgjs[(depths > zmax)] = np.nan
            smtotijs[(depths > zmax)] = np.nan
        # sal_m_or[(timec_year == timec_dates_lls[dl, 0]) &
        #          (timec_month == timec_dates_lls[dl, 1]) &
        #          (timec_day == timec_dates_lls[dl, 2]) &
        #          (ctd_latidx == timec_dates_lls[dl, 3]) &
        #          (ctd_lonidx == timec_dates_lls[dl, 4])] = smtoto
        # sal_m_J3[(timec_year == timec_dates_lls[dl, 0]) &
        #          (timec_month == timec_dates_lls[dl, 1]) &
        #          (timec_day == timec_dates_lls[dl, 2]) &
        #          (ctd_latidx == timec_dates_lls[dl, 3]) &
        #          (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotj
        # sal_m_FO[(timec_year == timec_dates_lls[dl, 0]) &
        #          (timec_month == timec_dates_lls[dl, 1]) &
        #          (timec_day == timec_dates_lls[dl, 2]) &
        #          (ctd_latidx == timec_dates_lls[dl, 3]) &
        #          (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotf
        # sal_m_J5_GLS[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotg
        # sal_m_TSUV[(timec_year == timec_dates_lls[dl, 0]) &
        #            (timec_month == timec_dates_lls[dl, 1]) &
        #            (timec_day == timec_dates_lls[dl, 2]) &
        #            (ctd_latidx == timec_dates_lls[dl, 3]) &
        #            (ctd_lonidx == timec_dates_lls[dl, 4])] = smtott
        # sal_m_GLO[(timec_year == timec_dates_lls[dl, 0]) &
        #           (timec_month == timec_dates_lls[dl, 1]) &
        #           (timec_day == timec_dates_lls[dl, 2]) &
        #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        #           (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotglo
        # sal_m_GL2[(timec_year == timec_dates_lls[dl, 0]) &
        #           (timec_month == timec_dates_lls[dl, 1]) &
        #           (timec_day == timec_dates_lls[dl, 2]) &
        #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        #           (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotglo2
        sal_m_GJS[(timec_year == timec_dates_lls[dl, 0]) &
                  (timec_month == timec_dates_lls[dl, 1]) &
                  (timec_day == timec_dates_lls[dl, 2]) &
                  (ctd_latidx == timec_dates_lls[dl, 3]) &
                  (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotgjs
        sal_m_IJS[(timec_year == timec_dates_lls[dl, 0]) &
                  (timec_month == timec_dates_lls[dl, 1]) &
                  (timec_day == timec_dates_lls[dl, 2]) &
                  (ctd_latidx == timec_dates_lls[dl, 3]) &
                  (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotijs
        # nchisor.close()
        # nchis_J3.close()
        # nchis_FO.close()
        # nchis_J5_GLS.close()
        # nchis_TSUV.close()
        # nchis_GLO.close()
        # nchis_GLO2.close()
        nchis_GJS.close()
        nchis_IJS.close()

    trmnz = np.argwhere(tempm_GJS > 0).flatten()
    srmnz = np.argwhere(sal_m_GJS > 0).flatten()
    tsrmnz = np.intersect1d(trmnz, srmnz)

    timec = timec[tsrmnz]
    deptc = deptc[tsrmnz]
    lat_c = lat_c[tsrmnz]
    lon_c = lon_c[tsrmnz]
    tempc = tempc[tsrmnz]
    sal_c = sal_c[tsrmnz]
    # tempm_or = tempm_or[tsrmnz]
    # sal_m_or = sal_m_or[tsrmnz]
    # tempm_J3 = tempm_J3[tsrmnz]
    # sal_m_J3 = sal_m_J3[tsrmnz]
    # tempm_FO = tempm_FO[tsrmnz]
    # sal_m_FO = sal_m_FO[tsrmnz]
    # tempm_J5_GLS = tempm_J5_GLS[tsrmnz]
    # sal_m_J5_GLS = sal_m_J5_GLS[tsrmnz]
    # tempm_TSUV = tempm_TSUV[tsrmnz]
    # sal_m_TSUV = sal_m_TSUV[tsrmnz]
    # tempm_GLO = tempm_GLO[tsrmnz]
    # sal_m_GLO = sal_m_GLO[tsrmnz]
    # tempm_GL2 = tempm_GL2[tsrmnz]
    # sal_m_GL2 = sal_m_GL2[tsrmnz]
    tempm_GJS = tempm_GJS[tsrmnz]
    sal_m_GJS = sal_m_GJS[tsrmnz]
    tempm_IJS = tempm_IJS[tsrmnz]
    sal_m_IJS = sal_m_IJS[tsrmnz]

    # 2nd round of processing with timestamp revised
    ctdfil = '/media/dskone/CELTIC/IMI_CTD_34f3_19f6_97e1.nc'
    ncctf = netcdf(ctdfil, 'r')
    ncctfo = ncctf.variables

    timec_date = n2d(timec, units=ncctfo['time'].units, calendar=ncctfo['time'].calendar)
    timec_year = np.zeros(timec_date.shape[0], dtype=np.int32)
    timec_month = np.zeros(timec_date.shape[0], dtype=np.int32)
    timec_day = np.zeros(timec_date.shape[0], dtype=np.int32)
    ctd_latidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    ctd_lonidx = np.zeros(timec_date.shape[0], dtype=np.int32)
    timec_ymdll = np.zeros([timec_date.shape[0], 10], dtype=np.int32)
    for td in range(0, len(timec_date)):
        timec_year[td] = timec_date[td].year
        timec_ymdll[td, 0] = timec_date[td].year
        timec_month[td] = timec_date[td].month
        timec_ymdll[td, 1] = timec_date[td].month
        timec_day[td] = timec_date[td].day
        timec_ymdll[td, 2] = timec_date[td].day
        ctd_latidx[td] = glor.geo_idx(lat_c[td], lat_rho[:, 0])
        timec_ymdll[td, 3] = glor.geo_idx(lat_c[td], lat_rho[:, 0])
        ctd_lonidx[td] = glor.geo_idx(lon_c[td], lon_rho[0, :])
        timec_ymdll[td, 4] = glor.geo_idx(lon_c[td], lon_rho[0, :])
        ibi_latidx[td] = glor.geo_idx(lat_c[td], IBIlat)
        timec_ymdll[td, 5] = glor.geo_idx(lat_c[td], IBIlat)
        ibi_lonidx[td] = glor.geo_idx(lon_c[td], IBIlon)
        timec_ymdll[td, 6] = glor.geo_idx(lon_c[td], IBIlon)
        glo_latidx[td] = glor.geo_idx(lat_c[td], GLOlat)
        timec_ymdll[td, 7] = glor.geo_idx(lat_c[td], GLOlat)
        glo_lonidx[td] = glor.geo_idx(lon_c[td], GLOlon)
        timec_ymdll[td, 8] = glor.geo_idx(lon_c[td], GLOlon)
        timec_ymdll[td, 9] = timec_date[td].hour
    ncctf.close()

    # Get unique dates to cycle through the equivalent model output history files
    # Get unique grid coordinates for lat/lon combinations for each unique date
    timec_dates_lls = np.unique(timec_ymdll, axis=0)
    timec_dates = np.unique(timec_ymdll[:, 0:2], axis=0)
    timec_ll = np.unique(timec_ymdll[:, 3:5], axis=0)

    mapfname = hisdir_orig + 'croco_VAL' + '_CTD_positions_2017.jpg'
    fig = plt.figure(figsize=(7, 5))
    ax = plt.contourf(lon_rho, lat_rho, h, [20,  50])
    plt.scatter(lon_rho[0, timec_dates_lls[:, 4]], lat_rho[timec_dates_lls[:, 3], 0], s=0)
    ptlabel = np.zeros(timec_dates_lls.shape[0], dtype=np.int32)
    for p in range(0, timec_dates_lls.shape[0]):
        ptlabel[p] = p+1
    for i, txt in enumerate(ptlabel):
        plt.annotate(txt, (lon_rho[0, timec_dates_lls[i, 4]], lat_rho[timec_dates_lls[i, 3], 0]),
                     weight='bold', size=8)
    plt.savefig(mapfname, dpi='figure', format='jpg')

    for tdl in range(0, len(timec_dates_lls)):
        ibifile = IBI_dir + IBI_prefix + \
                  str(timec_dates_lls[tdl, 0]) + str(timec_dates_lls[tdl, 1]).zfill(2) + \
                  str(timec_dates_lls[tdl, 2]).zfill(2) + '_' + \
                  str(timec_dates_lls[tdl, 0]) + str(timec_dates_lls[tdl, 1]).zfill(2) + \
                  str(timec_dates_lls[tdl, 2]).zfill(2) + IBI_post
        glofiles = sorted(glob.glob(GLO_dir + GLO_prefix +
                                    str(timec_dates_lls[tdl, 0]) + str(timec_dates_lls[tdl, 1]).zfill(2) +
                                    str(timec_dates_lls[tdl, 2]).zfill(2) + GLO_post))
        glofile = glofiles[0]

        ncibi = netcdf(ibifile, 'r')
        ncibio = ncibi.variables
        # Extract vertical section from history file
        temp_ibi = np.array(ncibio['thetao']
                            [:, :, timec_dates_lls[tdl, 5], timec_dates_lls[tdl, 6]]).flatten()
        salt_ibi = np.array(ncibio['so']
                            [:, :, timec_dates_lls[tdl, 5], timec_dates_lls[tdl, 6]]).flatten()

        ncglo = netcdf(glofile, 'r')
        ncgloo = ncglo.variables
        # Extract vertical section from history file
        temp_glo = np.array(ncgloo['thetao']
                            [:, :, timec_dates_lls[tdl, 7], timec_dates_lls[tdl, 8]]).flatten()
        salt_glo = np.array(ncgloo['so']
                            [:, :, timec_dates_lls[tdl, 7], timec_dates_lls[tdl, 8]]).flatten()

        # Cycle through the individual dates to get repeated section plots for visual assessment of model performance
        # for d in range(0, len(timec_dates)):
        latsect = lat_rho[timec_dates_lls[tdl, 3], 0]
        latdegs = np.floor(latsect)
        latmins = np.round((latsect - latdegs) * 60)

        lonsect = lon_rho[0, timec_dates_lls[tdl, 4]]
        londegs = np.abs(np.ceil(lonsect))
        lonmins = np.round((abs(lonsect) - londegs) * 60)

        sectftname = hisdir_orig + 'croco_TEMP_sect' + '_' + str(int(tdl + 1)) + '_' + \
                     str(timec_dates_lls[tdl, 0]) + \
                     str(timec_dates_lls[tdl, 1]).zfill(2) + \
                     str(timec_dates_lls[tdl, 2]).zfill(2) + '_' + \
                     str(int(latdegs)).zfill(2) + 'd' + str(int(latmins)).zfill(2) + 'm' + 'N' + '_' + \
                     str(int(londegs)).zfill(2) + 'd' + str(int(lonmins)).zfill(2) + 'm' + 'W' + '.jpg'

        sectfsname = hisdir_orig + 'croco_SALT_sect' + '_' + str(int(tdl + 1)) + '_' + \
                     str(timec_dates_lls[tdl, 0]) + \
                     str(timec_dates_lls[tdl, 1]).zfill(2) + \
                     str(timec_dates_lls[tdl, 2]).zfill(2) + '_' + \
                     str(int(latdegs)).zfill(2) + 'd' + str(int(latmins)).zfill(2) + 'm' + 'N' + '_' + \
                     str(int(londegs)).zfill(2) + 'd' + str(int(lonmins)).zfill(2) + 'm' + 'W' + '.jpg'
        rwidx = np.argwhere(((timec_year == timec_dates_lls[tdl, 0]) &
                             (timec_month == timec_dates_lls[tdl, 1]) &
                             (timec_day == timec_dates_lls[tdl, 2]) &
                             (ctd_latidx == timec_dates_lls[tdl, 3]) &
                             (ctd_lonidx == timec_dates_lls[tdl, 4]))).flatten()

        plt.figure(figsize=(7, 5))
        sortidx = np.argsort(deptc[rwidx])
        dmin = np.min(deptc[rwidx][sortidx])
        dmax = np.max(deptc[rwidx][sortidx])
        ibididx = np.argwhere(IBIdep >= dmax)[0][0]-1
        glodidx = np.argwhere(GLOdep >= dmax)[0][0]-1
        # plt.plot(tempc[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
        #          tempm_or[rwidx][sortidx], deptc[rwidx][sortidx], 'b--',
        #          tempm_J3[rwidx][sortidx], deptc[rwidx][sortidx], 'b-.',
        #          tempm_FO[rwidx][sortidx], deptc[rwidx][sortidx], 'b:',
        #          tempm_J5_GLS[rwidx][sortidx], deptc[rwidx][sortidx], 'g',
        #          tempm_TSUV[rwidx][sortidx], deptc[rwidx][sortidx], 'g--',
        #          tempm_IJS[rwidx][sortidx], deptc[rwidx][sortidx], 'g:')
        #          # tempm_GLO[rwidx][sortidx], deptc[rwidx][sortidx], 'g:',
        #          # tempm_GL2[rwidx][sortidx], deptc[rwidx][sortidx], 'gh',
        #          # tempm_GJS[rwidx][sortidx], deptc[rwidx][sortidx], 'gs')
        # plt.legend(('CTD', 'M_Orig', 'M_J3', 'M_FO', 'M_J5_GLS', 'M_TSUV', 'M_IBI_JS'),
        #            loc='best', shadow=True)
        plt.plot(tempc[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
                 temp_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
                 tempm_IJS[rwidx][sortidx], deptc[rwidx][sortidx], 'b--',
                 temp_glo[0:glodidx], GLOdep[0:glodidx], 'g',
                 tempm_GJS[rwidx][sortidx], deptc[rwidx][sortidx], 'g--')
        plt.legend(('CTD', 'IBI_raw', 'IBI_croco', 'GLO_raw', 'GLO_croco'), loc='best', shadow=True)
        plt.title(str(timec_dates_lls[tdl, 9]) + 'o clock')
        plt.gca().invert_yaxis()
        plt.xlabel('Temp ($\degree$C)')
        plt.xlim(9.5, 16.5)
        plt.ylabel('Depth(m)')
        plt.savefig(sectftname, dpi='figure', format='jpg')

        plt.figure(figsize=(7, 5))
        # plt.plot(sal_c[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
        #          sal_m_or[rwidx][sortidx], deptc[rwidx][sortidx], 'b--',
        #          sal_m_J3[rwidx][sortidx], deptc[rwidx][sortidx], 'b-.',
        #          sal_m_FO[rwidx][sortidx], deptc[rwidx][sortidx], 'b:',
        #          sal_m_J5_GLS[rwidx][sortidx], deptc[rwidx][sortidx], 'g',
        #          sal_m_TSUV[rwidx][sortidx], deptc[rwidx][sortidx], 'g--',
        #          sal_m_IJS[rwidx][sortidx], deptc[rwidx][sortidx], 'g:')
        #          # sal_m_GLO[rwidx][sortidx], deptc[rwidx][sortidx], 'g:',
        #          # sal_m_GL2[rwidx][sortidx], deptc[rwidx][sortidx], 'gh',
        #          # sal_m_GJS[rwidx][sortidx], deptc[rwidx][sortidx], 'gs')
        # plt.legend(('CTD', 'M_Orig', 'M_J3', 'M_FO', 'M_J5_GLS', 'M_TSUV', 'M_IBI_JS'),
        #            loc='best', shadow=True)
        plt.plot(sal_c[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
                 salt_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
                 sal_m_IJS[rwidx][sortidx], deptc[rwidx][sortidx], 'b--',
                 salt_glo[0:glodidx], GLOdep[0:glodidx], 'g',
                 sal_m_GJS[rwidx][sortidx], deptc[rwidx][sortidx], 'g--')
        plt.legend(('CTD', 'IBI_raw', 'IBI_croco', 'GLO_raw', 'GLO_croco'), loc='best', shadow=True)
        plt.title(str(timec_dates_lls[tdl, 9]) + 'o clock')
        plt.gca().invert_yaxis()
        plt.xlabel('Salinity (psu)')
        # plt.xlim(32, 35.8)
        plt.ylabel('Depth(m)')

        plt.savefig(sectfsname, dpi='figure', format='jpg')
        ncibi.close()
        ncglo.close()
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

        np.savez('coeffs.npz',
                 coefT=coefT, elemT=elemT,
                 coefU=coefU, elemU=elemU,
                 coefV=coefV, elemV=elemV)

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
    # Get the GLORYS file name from the date (only one time step per file)
    #

    for iyear in range(Ystart, Yend + 1):
        for imonth in range(Mstart, Mend + 1):

            dpm = monthrange(iyear, imonth)

            sat_sst_files = sorted(glob.glob(sat_sstfiles_dir + str(iyear) + str(imonth).zfill(2) +
                                             '??' + sat_sst_ending))

            if len(sat_sst_files) == dpm[1]:

                # hisdir_stats = hisdir_TSUV
                # hisdir_stats = hisdir_GLO_JS
                hisdir_stats = hisdir_IBI_JS

                VALmname = hisdir_stats + 'croco_VALm' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                VALdname = hisdir_stats + 'croco_VALd' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                VALsmename = hisdir_stats + 'croco_VAL' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + 's_ME.jpg'

                print(' ')
                print(' Making validation file: ' + VALdname)
                print(' ')

                glor.create_VAL_PISCES_NORESM(VALdname, grdname, title, theta_s, theta_b, hc, N, Tini, vtransform)

                #
                # Open the CROCO validation data file for writing
                #

                ncvd = netcdf(VALdname, 'a')
                ncvdi = ncvd.variables

                for tday in range(0, dpm[1], 1):
                    sat_sstname = sat_sst_files[tday]

                    ncsatex = netcdf(sat_sstname, 'r')
                    ncsatexo = ncsatex.variables

                    print(sat_sstname)

                    #
                    # 1: SST
                    #

                    print('Interpolate Satellite SST...')

                    (sat_sst, NzGood) = glor.interp_sst(ncsatexo, 'analysed_sst', tndx_glo, -1, iminT, imaxT,
                                                        jminT, jmaxT, LonT, LatT, coefT, elemT)

                    ncvdi['temps'][tday, :, :] = sat_sst - 273.15

                    ncsatex.close()

                print(' ')
                print(' Making validation file: ' + VALmname)
                print(' ')

                glor.create_VAL_PISCES_NORESM(VALmname, grdname, title, theta_s, theta_b, hc, N, Tini, vtransform)

                hisfile = hisdir_stats + his_prefix + 'Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                nchis = netcdf(hisfile, 'r')
                nchiso = nchis.variables

                temps = np.array(nchiso['temp'][0:dpm[1], 0, :, :])
                salts = np.array(nchiso['salt'][0:dpm[1], 0, :, :])
                saltb = np.array(nchiso['salt'][0:dpm[1], len(nchiso['s_rho']) - 1, :, :])
                time = np.array(nchiso['time'][0:dpm[1]])

                temps_sat = np.array(ncvdi['temps'][:, :, :])

                for tslice in range(0, dpm[1], 1):
                    ncvdi['temps_t_bias'][tslice] = \
                        ((temps[tslice, :, :].flatten() -
                          temps_sat[tslice, :, :].flatten()).mean() / temps_sat[tslice, :, :].flatten().mean())
                    ncvdi['temps_t_ME'][tslice] = (temps[tslice, :, :].flatten() -
                                                   temps_sat[tslice, :, :].flatten()).mean()
                    ncvdi['temps_t_RMSD'][tslice] = (((temps[tslice, :, :].flatten() -
                                                       temps_sat[tslice, :, :].flatten()) ** 2).mean()) ** 0.5
                    ncvdi['temps_t_CC'][tslice] = np.corrcoef(temps[tslice, :, :].flatten(),
                                                              temps_sat[tslice, :, :].flatten())[0, 1]

                for la in range(0, lat_rho.shape[0], 1):
                    for lo in range(0, lon_rho.shape[1], 1):
                        ncvdi['temps_s_bias'][la, lo] = \
                            ((temps[:, la, lo].flatten() -
                              temps_sat[:, la, lo].flatten()).mean() / temps_sat[:, la, lo].flatten().mean())
                        ncvdi['temps_s_ME'][la, lo] = (temps[:, la, lo].flatten() -
                                                       temps_sat[:, la, lo].flatten()).mean()
                        ncvdi['temps_s_RMSD'][la, lo] = (((temps[:, la, lo].flatten() -
                                                           temps_sat[:, la, lo].flatten()) ** 2).mean()) ** 0.5
                        ncvdi['temps_s_CC'][la, lo] = np.corrcoef(temps[:, la, lo].flatten(),
                                                                  temps_sat[:, la, lo].flatten())[0, 1]

                ncmd = netcdf(VALmname, 'a')
                ncmdi = ncmd.variables
                ncmdi['temps_model-sat'][0:dpm[1], :, :] = temps - temps_sat
                ncmdi['salt_bottom-surf'][0:dpm[1], :, :] = saltb - salts
                ncmdi['temps'][0:dpm[1], :, :] = temps
                ncmdi['salts'][0:dpm[1], :, :] = salts
                ncmdi['saltb'][0:dpm[1], :, :] = saltb
                ncmdi['VAL_time'][0:dpm[1]] = time / 86400
                nchis.close()
                ncmd.close()

                ncvdi['VAL_time'][0:dpm[1]] = time / 86400
                ncvdi['temps_model-sat'][0:dpm[1], :, :] = temps - temps_sat
                ncvdi['salt_bottom-surf'][0:dpm[1], :, :] = saltb - salts

                spme = ncvdi['temps_s_ME'][:, :]
                spme2 = spme.filled(fill_value=0)

                # proj = ccrs.Mercator(central_longitude=-8.29, min_latitude=49, max_latitude=52.95)
                # ax = plt.axes(projection=proj)
                # ax = plt.axes(projection=ccrs.Mercator(central_longitude=-8.29, min_latitude=49, max_latitude=52.95))
                fig = plt.figure(figsize=(7, 5))
                ax = plt.contourf(lon_rho, lat_rho, spme2, [-1, -0.5, 0, 0.5, 1])
                plt.colorbar()
                # ax.set_global()
                # ax.coastlines()

                # coast_10m = cfeature.NaturalEarthFeature("physical", "land", "10m", edgecolor="k", facecolor="0.8")
                # ax.add_feature(coast_10m)
                plt.savefig(VALsmename, dpi='figure', format='jpg')
                # plt.savefig(VALsmename, dpi='figure', format='jpg', metadata=None,
                #             bbox_inches=None, pad_inches=0.1,
                #             facecolor='auto', edgecolor='auto',
                #             backend=None)
                ncvd.close()
