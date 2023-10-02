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

    Yini = 2013
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
    # grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m
    grdname = crocofiles_dir + 'CELTIC_grd_v5_MSL_h10.nc'  # LAT to MSL included, INFOMAR block averaged hmin = 10m
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

    # hisdir_orig = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_1719/'
    hisdir_orig = '/media/dskone/CELTIC/'
    hisdir_J3 = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_J3/'
    hisdir_FO = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_Fairall_Off/'
    hisdir_J5_GLS = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_J5_GLS/'
    hisdir_TSUV = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_TSUV/'
    hisdir_GLO = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_TSUV_GLO/'
    hisdir_GLO_2 = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_TSUV_GLO_II/'
    hisdir_GLO_JS = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_TSUV_GLO_JS/'
    hisdir_IBI_JS = '/media/dskthree/CROCO_PHY_1p1_Archive/CROCO_PHY_1p1_17_TSUV_IBI_JS/'
    hisdir_985849 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_985849/'
    hisdir_985962 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_985962/'
    hisdir_987878 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_987878/'
    hisdir_989174 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_989174/'
    hisdir_989667 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_989667/'
    hisdir_990003 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_990003/'
    hisdir_995440 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_995440/'
    hisdir_997146 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_997146/'
    hisdir_998468 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_998468/'
    hisdir_999199 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_999199/'
    hisdir_1001772 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1001772/'
    hisdir_1002301 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1002301/'
    hisdir_1002331 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1002331/'
    hisdir_1002634 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1002634/'
    hisdir_1002741 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1002741/'
    hisdir_1003228 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1003228/'
    hisdir_1004586 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_GLO_1004586/'
    hisdir_1004724 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_GLO_1004724/'
    hisdir_1004868 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_GLO_1004868/'
    hisdir_1005022 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_GLO_1005022/'
    hisdir_1005186 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_GLO_1005186/'
    hisdir_1005613 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_1319_TSUV_IBI_1005613/'
    hisdir_1006516 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_131519_TSUV_IBI_1006516/'
    hisdir_1006690 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_131519_TSUV_IBI_1006690/'
    hisdir_1007268 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_13_TSUV_IBI_1007268/'
    hisdir_1008287 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_1314_TSUV_IBI_1008287/'
    hisdir_1009495 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_1314_TSUV_IBI_1009495/'
    hisdir_1015382 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1015382/'
    hisdir_1015647 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1015647/'

    hisdir_1016234 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1016234/'
    hisdir_1016529 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1016529/'
    hisdir_1016656 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_1016656/'

    his_prefix = 'croco_his_'
    stn_prefix = 'croco_stn_'

    sat_sstfiles_dir = '/media/dskone/VAL/SAT_SST/L4/'
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

    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_34f3_19f6_97e1.nc'

    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2013.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2014.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2015.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2016.nc'
    ctdfil = '/media/dskone/CELTIC/IMI_CTD_2017.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2018.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2019.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2020.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2021.nc'
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

        # his_GJS = hisdir_GLO_JS + \
        #           his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_GJS = netcdf(his_GJS, 'r')
        # nchis_GJSo = nchis_GJS.variables
        # # Extract vertical section from history file
        # tempGJS = np.array(nchis_GJSo['temp']
        #                    [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltGJS = np.array(nchis_GJSo['salt']
        #                    [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zetaGJS = np.array(nchis_GJSo['zeta']
        #                    [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_IJS = hisdir_IBI_JS + \
        #           his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_IJS = netcdf(his_IJS, 'r')
        # nchis_IJSo = nchis_IJS.variables
        # # Extract vertical section from history file
        # tempIJS = np.array(nchis_IJSo['temp']
        #                    [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # saltIJS = np.array(nchis_IJSo['salt']
        #                    [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zetaIJS = np.array(nchis_IJSo['zeta']
        #                    [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        # Interpolate section data onto depths from the CTD cast for that date
        z = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], 0, theta_s, theta_b, hc, N, 'w', vtransform)
        z[0] = np.floor(z[0])
        zmin = z.min()
        zmax = z.max()

        # ftGJS = interpolate.interp1d(z, np.append(tempGJS[0], tempGJS))
        # tmtotgjs = np.zeros_like(depths)
        # tmtotgjs[(depths >= zmin) & (depths <= zmax)] = ftGJS(depths[(depths >= zmin) & (depths <= zmax)])

        # ftIJS = interpolate.interp1d(z, np.append(tempIJS[0], tempIJS))
        # tmtotijs = np.zeros_like(depths)
        # tmtotijs[(depths >= zmin) & (depths <= zmax)] = ftIJS(depths[(depths >= zmin) & (depths <= zmax)])

        # if depths.min() < zmin:
        #     # tmtotgjs[(depths < zmin)] = np.nan
        #     tmtotijs[(depths < zmin)] = np.nan
        # if depths.max() > zmax:
        #     # tmtotgjs[(depths > zmax)] = np.nan
        #     tmtotijs[(depths > zmax)] = np.nan
        # # tempm_GJS[(timec_year == timec_dates_lls[dl, 0]) &
        # #           (timec_month == timec_dates_lls[dl, 1]) &
        # #           (timec_day == timec_dates_lls[dl, 2]) &
        # #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        # #           (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotgjs
        # tempm_IJS[(timec_year == timec_dates_lls[dl, 0]) &
        #           (timec_month == timec_dates_lls[dl, 1]) &
        #           (timec_day == timec_dates_lls[dl, 2]) &
        #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        #           (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotijs

        # fsGJS = interpolate.interp1d(z, np.append(saltGJS[0], saltGJS))
        # smtotgjs = np.zeros_like(depths)
        # smtotgjs[(depths >= zmin) & (depths <= zmax)] = fsGJS(depths[(depths >= zmin) & (depths <= zmax)])

        # fsIJS = interpolate.interp1d(z, np.append(saltIJS[0], saltIJS))
        # smtotijs = np.zeros_like(depths)
        # smtotijs[(depths >= zmin) & (depths <= zmax)] = fsIJS(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # if depths.min() < zmin:
        #     # smtotgjs[(depths < zmin)] = np.nan
        #     smtotijs[(depths < zmin)] = np.nan
        # if depths.max() > zmax:
        #     # smtotgjs[(depths > zmax)] = np.nan
        #     smtotijs[(depths > zmax)] = np.nan
        # # sal_m_GJS[(timec_year == timec_dates_lls[dl, 0]) &
        # #           (timec_month == timec_dates_lls[dl, 1]) &
        # #           (timec_day == timec_dates_lls[dl, 2]) &
        # #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        # #           (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotgjs
        # sal_m_IJS[(timec_year == timec_dates_lls[dl, 0]) &
        #           (timec_month == timec_dates_lls[dl, 1]) &
        #           (timec_day == timec_dates_lls[dl, 2]) &
        #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        #           (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotijs
        # # nchis_GJS.close()
        # nchis_IJS.close()
    #
    # trmnz = np.argwhere(tempm_IJS > 0).flatten()
    # srmnz = np.argwhere(sal_m_IJS > 0).flatten()
    # tsrmnz = np.intersect1d(trmnz, srmnz)
    #
    # timec = timec[tsrmnz]
    # deptc = deptc[tsrmnz]
    # lat_c = lat_c[tsrmnz]
    # lon_c = lon_c[tsrmnz]
    # tempc = tempc[tsrmnz]
    # sal_c = sal_c[tsrmnz]
    # tempm_GJS = tempm_GJS[tsrmnz]
    # sal_m_GJS = sal_m_GJS[tsrmnz]
    # tempm_IJS = tempm_IJS[tsrmnz]
    # sal_m_IJS = sal_m_IJS[tsrmnz]

    # 2nd round of processing with timestamp revised
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

    ctd_surf_t = np.zeros(timec_dates_lls.shape[0])
    ctd_bott_t = np.zeros(timec_dates_lls.shape[0])
    ibi_surf_t = np.zeros(timec_dates_lls.shape[0])
    ibi_bott_t = np.zeros(timec_dates_lls.shape[0])
    glo_surf_t = np.zeros(timec_dates_lls.shape[0])
    glo_bott_t = np.zeros(timec_dates_lls.shape[0])
    cel_surf_t = np.zeros(timec_dates_lls.shape[0])
    cel_bott_t = np.zeros(timec_dates_lls.shape[0])

    ibi_t_ctd = np.zeros_like(tempc)
    glo_t_ctd = np.zeros_like(tempc)
    cel_t_ctd = np.zeros_like(tempc)

    ibi_s_ctd = np.zeros_like(sal_c)
    glo_s_ctd = np.zeros_like(sal_c)
    cel_s_ctd = np.zeros_like(sal_c)

    ctd_surf_s = np.zeros(timec_dates_lls.shape[0])
    ctd_bott_s = np.zeros(timec_dates_lls.shape[0])
    ibi_surf_s = np.zeros(timec_dates_lls.shape[0])
    ibi_bott_s = np.zeros(timec_dates_lls.shape[0])
    glo_surf_s = np.zeros(timec_dates_lls.shape[0])
    glo_bott_s = np.zeros(timec_dates_lls.shape[0])
    cel_surf_s = np.zeros(timec_dates_lls.shape[0])
    cel_bott_s = np.zeros(timec_dates_lls.shape[0])

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
        temp_ibi_fill = ncibio['thetao'].add_offset
        salt_ibi = np.array(ncibio['so']
                            [:, :, timec_dates_lls[tdl, 5], timec_dates_lls[tdl, 6]]).flatten()
        salt_ibi_fill = ncibio['so'].add_offset
        ncglo = netcdf(glofile, 'r')
        ncgloo = ncglo.variables
        # Extract vertical section from history file
        temp_glo = np.array(ncgloo['thetao']
                            [:, :, timec_dates_lls[tdl, 7], timec_dates_lls[tdl, 8]]).flatten()
        temp_glo_fill = ncgloo['thetao'].missing_value
        salt_glo = np.array(ncgloo['so']
                            [:, :, timec_dates_lls[tdl, 7], timec_dates_lls[tdl, 8]]).flatten()
        salt_glo_fill = ncgloo['so'].missing_value

        # hisdir_station = hisdir_985849
        # hisdir_station = hisdir_985962
        # hisdir_station = hisdir_987878
        # hisdir_station = hisdir_990003
        # hisdir_station = hisdir_995440
        # hisdir_station = hisdir_997146
        # hisdir_station = hisdir_998468
        # hisdir_station = hisdir_999199
        # hisdir_station = hisdir_1001772
        # hisdir_station = hisdir_1002301
        # hisdir_station = hisdir_1002331
        # hisdir_station = hisdir_1002634
        # hisdir_station = hisdir_1002741
        # hisdir_station = hisdir_1003228
        # hisdir_station = hisdir_1004586
        # hisdir_station = hisdir_1004724
        # hisdir_station = hisdir_1004868
        # hisdir_station = hisdir_1005022
        # hisdir_station = hisdir_1005186
        # hisdir_station = hisdir_1005613
        # hisdir_station = hisdir_1006516
        # hisdir_station = hisdir_1006690
        # hisdir_station = hisdir_1007268
        # hisdir_station = hisdir_1008287
        # hisdir_station = hisdir_1009495
        # hisdir_station = hisdir_1015382
        # hisdir_station = hisdir_1015647
        # hisdir_station = hisdir_1016234
        # hisdir_station = hisdir_1016529
        hisdir_station = hisdir_1016656

        crocostn_ibi = hisdir_station + stn_prefix + \
                       'Y' + str(timec_dates_lls[tdl, 0]) + \
                       'M' + str(timec_dates_lls[tdl, 1]).zfill(2) + '.nc'
        if tdl == 0 or tdl == 1:
            tdl2 = 0
        else:
            tdl2 = tdl - 1
        nccrocostn_ibi = netcdf(crocostn_ibi, 'r')
        nccrocostn_ibio = nccrocostn_ibi.variables
        latstn = np.array(nccrocostn_ibio['lat'][:])
        lonstn = np.array(nccrocostn_ibio['lon'][:])
        latgrd = lat_rho[timec_dates_lls[tdl, 3], timec_dates_lls[tdl, 4]]
        longrd = lon_rho[timec_dates_lls[tdl, 3], timec_dates_lls[tdl, 4]]
        latidx_stn = glor.geo_idx(latgrd, latstn)
        lonidx_stn = glor.geo_idx(longrd, lonstn)
        # tempcrocostn_ibi = np.array(nccrocostn_ibio['temp']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9]-12, tdl2, :])
        # saltcrocostn_ibi = np.array(nccrocostn_ibio['salt']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9]-12, tdl2, :])
        # zetacrocostn_ibi = np.array(nccrocostn_ibio['zeta']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9]-12, tdl2])
        # tempcrocostn_ibi = np.array(nccrocostn_ibio['temp']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9], tdl2, :])
        # saltcrocostn_ibi = np.array(nccrocostn_ibio['salt']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9], tdl2, :])
        # zetacrocostn_ibi = np.array(nccrocostn_ibio['zeta']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9], tdl2])
        tempcrocostn_ibi = np.array(nccrocostn_ibio['temp']
                                    [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9], latidx_stn, :])
        saltcrocostn_ibi = np.array(nccrocostn_ibio['salt']
                                    [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9], latidx_stn, :])
        zetacrocostn_ibi = np.array(nccrocostn_ibio['zeta']
                                    [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9], latidx_stn])
        # dibi = vgrd.zlevs(h[timec_dates_lls[tdl, 3], timec_dates_lls[tdl, 4]], zetacrocostn_ibi,
        #                   theta_s, theta_b, hc, N, 'r', vtransform)*-1
        dibi = vgrd.zlevs(h[timec_dates_lls[tdl, 3], timec_dates_lls[tdl, 4]], zetacrocostn_ibi,
                          theta_s, theta_b, 50, N, 'r', vtransform)*-1

        # crocostn_glo = hisdir_GLO_JS + stn_prefix + \
        #                'Y' + str(timec_dates_lls[tdl, 0]) + \
        #                'M' + str(timec_dates_lls[tdl, 1]).zfill(2) + '.nc'
        #
        # nccrocostn_glo = netcdf(crocostn_glo, 'r')
        # nccrocostn_gloo = nccrocostn_glo.variables
        # tempcrocostn_glo = np.array(nccrocostn_gloo['temp']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9]-12, tdl, :])
        # saltcrocostn_glo = np.array(nccrocostn_gloo['salt']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9]-12, tdl, :])
        # zetacrocostn_glo = np.array(nccrocostn_gloo['zeta']
        #                             [((timec_dates_lls[tdl, 2]-1)*24)+timec_dates_lls[tdl, 9]-12, tdl])
        # dglo = vgrd.zlevs(h[timec_dates_lls[tdl, 3], timec_dates_lls[tdl, 4]], zetacrocostn_glo,
        #                   theta_s, theta_b, hc, N, 'r', vtransform)*-1

        # Cycle through the individual dates to get repeated section plots for visual assessment of model performance
        # for d in range(0, len(timec_dates)):
        latsect = lat_rho[timec_dates_lls[tdl, 3], 0]
        latdegs = np.floor(latsect)
        latmins = np.round((latsect - latdegs) * 60)

        lonsect = lon_rho[0, timec_dates_lls[tdl, 4]]
        londegs = np.abs(np.ceil(lonsect))
        lonmins = np.round((abs(lonsect) - londegs) * 60)

        sectftname = hisdir_station + 'croco_TEMP_stn' + '_' + str(int(tdl + 1)) + '_' + \
                     str(timec_dates_lls[tdl, 0]) + \
                     str(timec_dates_lls[tdl, 1]).zfill(2) + \
                     str(timec_dates_lls[tdl, 2]).zfill(2) + '_' + \
                     str(int(latdegs)).zfill(2) + 'd' + str(int(latmins)).zfill(2) + 'm' + 'N' + '_' + \
                     str(int(londegs)).zfill(2) + 'd' + str(int(lonmins)).zfill(2) + 'm' + 'W' + '.jpg'

        sectfsname = hisdir_station + 'croco_SALT_stn' + '_' + str(int(tdl + 1)) + '_' + \
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
        idmin = np.min(IBIdep)
        gdmin = np.min(GLOdep)
        dmin = np.max((dmin, idmin, gdmin))
        dmax = np.max(deptc[rwidx][sortidx])
        ibididx = np.argwhere(IBIdep >= dmax)[0][0]-1
        glodidx = np.argwhere(GLOdep >= dmax)[0][0]-1
        glofidx = np.argwhere(temp_glo < 0)[0][0]-1
        glodidx = np.min((glodidx, glofidx))
        dmax_int = np.min((dmax, IBIdep[ibididx], GLOdep[glodidx]))
        ctd_dmax = np.argwhere(deptc[rwidx][sortidx] >= dmax_int)[0][0]-1

        ctdint_t = interpolate.interp1d(deptc[rwidx][sortidx], tempc[rwidx][sortidx], fill_value='extrapolate')
        ctd_surf_t[tdl] = ctdint_t(dmin)
        ctd_bott_t[tdl] = ctdint_t(dmax_int)

        if temp_ibi[0] != temp_ibi_fill:
            ibiint_t = interpolate.interp1d(IBIdep[temp_ibi != temp_ibi_fill], temp_ibi[temp_ibi != temp_ibi_fill],
                                            fill_value='extrapolate')
            ibi_surf_t[tdl] = ibiint_t(dmin)
            ibi_bott_t[tdl] = ibiint_t(dmax_int)
            ibi_tat_ctd = np.zeros_like(deptc[rwidx][sortidx])
            ibi_tat_ctd[0:ctd_dmax] = ibiint_t(deptc[rwidx][sortidx][0:ctd_dmax])
            ibi_tat_ctd[ibi_tat_ctd == 0] = np.nan
            ibi_t_ctd[rwidx[sortidx]] = ibi_tat_ctd
        else:
            ibi_surf_t[tdl] = np.nan
            ibi_bott_t[tdl] = np.nan
            ibi_t_ctd[rwidx[sortidx]] = np.nan

        if temp_glo[0] != temp_glo_fill:
            gloint_t = interpolate.interp1d(GLOdep[temp_glo != temp_glo_fill], temp_glo[temp_glo != temp_glo_fill],
                                            fill_value='extrapolate')
            glo_surf_t[tdl] = gloint_t(dmin)
            glo_bott_t[tdl] = gloint_t(dmax_int)
            glo_tat_ctd = np.zeros_like(deptc[rwidx][sortidx])
            glo_tat_ctd[0:ctd_dmax] = gloint_t(deptc[rwidx][sortidx][0:ctd_dmax])
            glo_tat_ctd[glo_tat_ctd == 0] = np.nan
            glo_t_ctd[rwidx[sortidx]] = glo_tat_ctd
        else:
            glo_surf_t[tdl] = np.nan
            glo_bott_t[tdl] = np.nan
            glo_t_ctd[rwidx[sortidx]] = np.nan

        if tempcrocostn_ibi[-1] != 0:
            celint_t = interpolate.interp1d(dibi, tempcrocostn_ibi, fill_value='extrapolate')
            cel_surf_t[tdl] = celint_t(dmin)
            cel_bott_t[tdl] = celint_t(dmax_int)
            cel_tat_ctd = np.zeros_like(deptc[rwidx][sortidx])
            cel_tat_ctd[0:ctd_dmax] = celint_t(deptc[rwidx][sortidx][0:ctd_dmax])
            cel_tat_ctd[cel_tat_ctd == 0] = np.nan
            cel_t_ctd[rwidx[sortidx]] = cel_tat_ctd
        else:
            cel_surf_t[tdl] = np.nan
            cel_bott_t[tdl] = np.nan
            cel_t_ctd[rwidx[sortidx]] = np.nan

        plt.plot(tempc[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
                 temp_ibi[temp_ibi != temp_ibi_fill], IBIdep[temp_ibi != temp_ibi_fill], 'b',
                 tempcrocostn_ibi, dibi, 'b--',
                 temp_glo[temp_glo != temp_glo_fill], GLOdep[temp_glo != temp_glo_fill], 'g')
                 # temp_glo[0:glodidx], GLOdep[0:glodidx], 'g',
                 # tempcrocostn_glo, dglo, 'g--')
        plt.legend(('CTD', 'IBI_raw', 'IBI_croco', 'GLO_raw'), loc='best', shadow=True)
        # plt.legend(('CTD', 'IBI_raw', 'IBI_croco', 'GLO_raw', 'GLO_croco'), loc='best', shadow=True)
        plt.title(str(timec_dates_lls[tdl, 9]) + 'o clock')
        plt.ylim(dmin, dmax)
        plt.gca().invert_yaxis()
        plt.xlabel('Temp ($\degree$C)')
        plt.xlim(9.5, 16.5)
        plt.ylabel('Depth(m)')
        plt.savefig(sectftname, dpi='figure', format='jpg')

        plt.figure(figsize=(7, 5))

        ctdint_s = interpolate.interp1d(deptc[rwidx][sortidx], sal_c[rwidx][sortidx], fill_value='extrapolate')
        ctd_surf_s[tdl] = ctdint_s(dmin)
        ctd_bott_s[tdl] = ctdint_s(dmax_int)

        if salt_ibi[0] != salt_ibi_fill:
            ibiint_s = interpolate.interp1d(IBIdep[salt_ibi != salt_ibi_fill], salt_ibi[salt_ibi != salt_ibi_fill],
                                            fill_value='extrapolate')
            ibi_surf_s[tdl] = ibiint_s(dmin)
            ibi_bott_s[tdl] = ibiint_s(dmax_int)
            ibi_sat_ctd = np.zeros_like(deptc[rwidx][sortidx])
            ibi_sat_ctd[0:ctd_dmax] = ibiint_s(deptc[rwidx][sortidx][0:ctd_dmax])
            ibi_sat_ctd[ibi_sat_ctd == 0] = np.nan
            ibi_s_ctd[rwidx[sortidx]] = ibi_sat_ctd
        else:
            ibi_surf_s[tdl] = np.nan
            ibi_bott_s[tdl] = np.nan
            ibi_s_ctd[rwidx[sortidx]] = np.nan

        if salt_glo[0] != salt_glo_fill:
            gloint_s = interpolate.interp1d(GLOdep[salt_glo != salt_glo_fill], salt_glo[salt_glo != salt_glo_fill],
                                            fill_value='extrapolate')
            glo_surf_s[tdl] = gloint_s(dmin)
            glo_bott_s[tdl] = gloint_s(dmax_int)
            glo_sat_ctd = np.zeros_like(deptc[rwidx][sortidx])
            glo_sat_ctd[0:ctd_dmax] = gloint_s(deptc[rwidx][sortidx][0:ctd_dmax])
            glo_sat_ctd[glo_sat_ctd == 0] = np.nan
            glo_s_ctd[rwidx[sortidx]] = glo_sat_ctd
        else:
            glo_surf_s[tdl] = np.nan
            glo_bott_s[tdl] = np.nan
            glo_s_ctd[rwidx[sortidx]] = np.nan

        if saltcrocostn_ibi[-1] != 0:
            celint_s = interpolate.interp1d(dibi, saltcrocostn_ibi, fill_value='extrapolate')
            cel_surf_s[tdl] = celint_s(dmin)
            cel_bott_s[tdl] = celint_s(dmax_int)
            cel_sat_ctd = np.zeros_like(deptc[rwidx][sortidx])
            cel_sat_ctd[0:ctd_dmax] = celint_s(deptc[rwidx][sortidx][0:ctd_dmax])
            cel_sat_ctd[cel_sat_ctd == 0] = np.nan
            cel_s_ctd[rwidx[sortidx]] = cel_sat_ctd
        else:
            cel_surf_s[tdl] = np.nan
            cel_bott_s[tdl] = np.nan
            cel_s_ctd[rwidx[sortidx]] = np.nan

        plt.plot(sal_c[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
                 salt_ibi[salt_ibi != salt_ibi_fill], IBIdep[salt_ibi != salt_ibi_fill], 'b',
                 saltcrocostn_ibi, dibi, 'b--',
                 salt_glo[salt_glo != salt_glo_fill], GLOdep[salt_glo != salt_glo_fill], 'g')
                 # salt_glo[0:glodidx], GLOdep[0:glodidx], 'g',
                 # saltcrocostn_glo, dglo, 'g--')
        plt.legend(('CTD', 'IBI_raw', 'IBI_croco', 'GLO_raw'), loc='best', shadow=True)
        # plt.legend(('CTD', 'IBI_raw', 'IBI_croco', 'GLO_raw', 'GLO_croco'), loc='best', shadow=True)
        plt.title(str(timec_dates_lls[tdl, 9]) + 'o clock')
        plt.ylim(dmin, dmax)
        plt.gca().invert_yaxis()
        plt.xlabel('Salinity (psu)')
        # plt.xlim(33.9, 35.7)
        plt.ylabel('Depth(m)')

        plt.savefig(sectfsname, dpi='figure', format='jpg')
        ncibi.close()
        ncglo.close()
        nccrocostn_ibi.close()
        # nccrocostn_glo.close()

        # ctd_surf_t, ctd_bott_t, ibi_surf_t, ibi_bott_t, glo_surf_t, glo_bott_t, cel_surf_t, cel_bott_t
        # ctd_surf_s, ctd_bott_s, ibi_surf_s, ibi_bott_s, glo_surf_s, glo_bott_s, cel_surf_s, cel_bott_s

    bias_yr = np.unique(timec_dates_lls[:, 0])
    bias_mths = np.unique(timec_dates_lls[:, 1])

    ibi_t_surfbias = ibi_surf_t - ctd_surf_t
    cel_t_surfbias = cel_surf_t - ctd_surf_t
    ibi_t_bottbias = ibi_bott_t - ctd_bott_t
    cel_t_bottbias = cel_bott_t - ctd_bott_t

    ibi_s_surfbias = ibi_surf_s - ctd_surf_s
    cel_s_surfbias = cel_surf_s - ctd_surf_s
    ibi_s_bottbias = ibi_bott_s - ctd_bott_s
    cel_s_bottbias = cel_bott_s - ctd_bott_s

    cm = plt.cm.get_cmap('RdYlBu_r')
    for mth in range(0, len(bias_mths)):
        maptname = hisdir_station + 'croco_tbias' + str(bias_yr[0]) + str(bias_mths[mth]).zfill(2) + '.jpg'
        mapsname = hisdir_station + 'croco_sbias' + str(bias_yr[0]) + str(bias_mths[mth]).zfill(2) + '.jpg'

        yrmthidx = np.argwhere((timec_dates_lls[:, 1] == bias_mths[mth]) &
                               (timec_dates_lls[:, 0] == bias_yr[0]) &
                               (ibi_t_surfbias != np.nan) & (cel_t_surfbias != np.nan) &
                               (ibi_t_bottbias != np.nan) & (cel_t_bottbias != np.nan) &
                               (ibi_s_surfbias != np.nan) & (cel_s_surfbias != np.nan) &
                               (ibi_s_bottbias != np.nan) & (cel_s_bottbias != np.nan)).flatten()

        t_vmin = np.nanmin((ibi_t_surfbias[yrmthidx], cel_t_surfbias[yrmthidx],
                            ibi_t_bottbias[yrmthidx], cel_t_bottbias[yrmthidx]))
        if abs(t_vmin) < 0.5:
            t_vmin = -0.5
        elif abs(t_vmin) < 1:
            t_vmin = -1
        else:
            t_vmin = np.ceil(abs(t_vmin)/0.5) * -0.5

        t_vmax = np.nanmax((ibi_t_surfbias[yrmthidx], cel_t_surfbias[yrmthidx],
                            ibi_t_bottbias[yrmthidx], cel_t_bottbias[yrmthidx]))
        if t_vmax < 0.5:
            t_vmax = 0.5
        elif t_vmax < 1:
            t_vmax = 1
        else:
            t_vmax = np.floor(t_vmax/0.5) * 0.5

        tlim = np.max((abs(t_vmin), t_vmax))
        t_vmin = tlim * -1
        t_vmax = tlim
        s_vmin = np.nanmin((ibi_s_surfbias[yrmthidx], cel_s_surfbias[yrmthidx],
                            ibi_s_bottbias[yrmthidx], cel_s_bottbias[yrmthidx]))
        if abs(s_vmin) < 0.25:
            s_vmin = -0.25
        elif abs(s_vmin) < 0.5:
            s_vmin = -0.5
        else:
            s_vmin = np.ceil(abs(s_vmin)/0.5) * -0.5
        s_vmax = np.nanmax((ibi_s_surfbias[yrmthidx], cel_s_surfbias[yrmthidx],
                            ibi_s_bottbias[yrmthidx], cel_s_bottbias[yrmthidx]))
        if s_vmax < 0.25:
            s_vmax = 0.25
        elif s_vmax < 0.5:
            s_vmax = 0.5
        else:
            s_vmax = np.ceil(abs(s_vmax)/0.5) * 0.5

        slim = np.max((abs(s_vmin), s_vmax))
        s_vmin = slim * -1
        s_vmax = slim

        fig, axs = plt.subplots(2, 2, figsize=(7, 5))
        fig.suptitle('Temp. bias($^\circ$C): ' +
                     str(bias_mths[mth]).zfill(2) + '/' + str(bias_yr[0]), weight='bold')
        axs[0, 0].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[0, 0].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[0, 0].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax1 = axs[0, 0].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=ibi_t_surfbias[yrmthidx], vmin=t_vmin, vmax=t_vmax, cmap=cm)
        axs[0, 0].set_title('IBI @ surf.')

        axs[0, 1].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[0, 1].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[0, 1].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax2 = axs[0, 1].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=cel_t_surfbias[yrmthidx], vmin=t_vmin, vmax=t_vmax, cmap=cm)
        axs[0, 1].set_title('CS1K @ surf.')

        axs[1, 0].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[1, 0].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[1, 0].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax3 = axs[1, 0].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=ibi_t_bottbias[yrmthidx], vmin=t_vmin, vmax=t_vmax, cmap=cm)
        axs[1, 0].set_title('IBI @ bott.')

        axs[1, 1].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[1, 1].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[1, 1].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax4 = axs[1, 1].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=cel_t_bottbias[yrmthidx], vmin=t_vmin, vmax=t_vmax, cmap=cm)
        axs[1, 1].set_title('CS1K @ bott.')

        fig.colorbar(ax4, ax=axs.ravel().tolist())
        plt.savefig(maptname, dpi='figure', format='jpg')

        fig, axs = plt.subplots(2, 2, figsize=(7, 5))
        fig.suptitle('Sal. bias (psu):' +
                     str(bias_mths[mth]).zfill(2) + '/' + str(bias_yr[0]), weight='bold')
        axs[0, 0].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[0, 0].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[0, 0].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax5 = axs[0, 0].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=ibi_s_surfbias[yrmthidx], vmin=s_vmin, vmax=s_vmax, cmap=cm)
        axs[0, 0].set_title('IBI @ surf.')

        axs[0, 1].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[0, 1].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[0, 1].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax6 = axs[0, 1].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=cel_s_surfbias[yrmthidx], vmin=s_vmin, vmax=s_vmax, cmap=cm)
        axs[0, 1].set_title('CS1K @ surf.')

        axs[1, 0].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[1, 0].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[1, 0].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax7 = axs[1, 0].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=ibi_s_bottbias[yrmthidx], vmin=s_vmin, vmax=s_vmax, cmap=cm)
        axs[1, 0].set_title('IBI @ bott.')

        axs[1, 1].contourf(lon_rho, lat_rho, h, [0, 10])
        axs0 = axs[1, 1].contour(lon_rho, lat_rho, h, [50, 100, 150, 200], colors='black',
                                 linewidths=[0.75, 0.5, 0.75, 0.5])
        axs[1, 1].clabel(axs0, levels=[50, 100, 150, 200], fmt='%1.0f', fontsize=6)
        ax8 = axs[1, 1].scatter(lon_rho[0, timec_dates_lls[yrmthidx, 4]], lat_rho[timec_dates_lls[yrmthidx, 3], 0], s=5,
                                c=cel_s_bottbias[yrmthidx], vmin=s_vmin, vmax=s_vmax, cmap=cm)
        axs[1, 1].set_title('CS1K @ bott.')

        fig.colorbar(ax8, ax=axs.ravel().tolist())
        plt.savefig(mapsname, dpi='figure', format='jpg')

    mrk = ["o", "v", "^", "<", ">", "s", "p", "P", "*", "X", "D", "8"]

    # T/S plot for July, October, November and December 2014
    Jul_ctdidx = np.argwhere((timec_ymdll[:, 1] == bias_mths[0]) &
                             (timec_ymdll[:, 0] == bias_yr[0]) &
                             (ibi_s_ctd != np.nan) & (ibi_t_ctd != np.nan) &
                             (cel_s_ctd != np.nan) & (cel_t_ctd != np.nan) &
                             (sal_c != np.nan) & (tempc != np.nan) &
                             (ibi_s_ctd != 0) & (ibi_t_ctd != 0) &
                             (cel_s_ctd != 0) & (cel_t_ctd != 0) &
                             (sal_c != 0) & (tempc != 0)).flatten()
    Oct_ctdidx = np.argwhere((timec_ymdll[:, 1] == bias_mths[1]) &
                             (timec_ymdll[:, 0] == bias_yr[0]) &
                             (ibi_s_ctd != np.nan) & (ibi_t_ctd != np.nan) &
                             (cel_s_ctd != np.nan) & (cel_t_ctd != np.nan) &
                             (sal_c != np.nan) & (tempc != np.nan) &
                             (ibi_s_ctd != 0) & (ibi_t_ctd != 0) &
                             (cel_s_ctd != 0) & (cel_t_ctd != 0) &
                             (sal_c != 0) & (tempc != 0)).flatten()
    Nov_ctdidx = np.argwhere((timec_ymdll[:, 1] == bias_mths[2]) &
                             (timec_ymdll[:, 0] == bias_yr[0]) &
                             (ibi_s_ctd != np.nan) & (ibi_t_ctd != np.nan) &
                             (cel_s_ctd != np.nan) & (cel_t_ctd != np.nan) &
                             (sal_c != np.nan) & (tempc != np.nan) &
                             (ibi_s_ctd != 0) & (ibi_t_ctd != 0) &
                             (cel_s_ctd != 0) & (cel_t_ctd != 0) &
                             (sal_c != 0) & (tempc != 0)).flatten()
    Dec_ctdidx = np.argwhere((timec_ymdll[:, 1] == bias_mths[3]) &
                             (timec_ymdll[:, 0] == bias_yr[0]) &
                             (ibi_s_ctd != np.nan) & (ibi_t_ctd != np.nan) &
                             (cel_s_ctd != np.nan) & (cel_t_ctd != np.nan) &
                             (sal_c != np.nan) & (tempc != np.nan) &
                             (ibi_s_ctd != 0) & (ibi_t_ctd != 0) &
                             (cel_s_ctd != 0) & (cel_t_ctd != 0) &
                             (sal_c != 0) & (tempc != 0)).flatten()

    s_vmin2 = np.nanmin([np.nanmin(sal_c[Jul_ctdidx]), np.nanmin(sal_c[Oct_ctdidx]), np.nanmin(sal_c[Nov_ctdidx]),
                         np.nanmin(sal_c[Dec_ctdidx]), np.nanmin(ibi_s_ctd[Jul_ctdidx]),
                         np.nanmin(ibi_s_ctd[Oct_ctdidx]), np.nanmin(ibi_s_ctd[Nov_ctdidx]),
                         np.nanmin(ibi_s_ctd[Dec_ctdidx])])
    s_vmax2 = np.nanmax([np.nanmax(sal_c[Jul_ctdidx]), np.nanmax(sal_c[Oct_ctdidx]), np.nanmax(sal_c[Nov_ctdidx]),
                         np.nanmax(sal_c[Dec_ctdidx]), np.nanmax(ibi_s_ctd[Jul_ctdidx]),
                         np.nanmax(ibi_s_ctd[Oct_ctdidx]), np.nanmax(ibi_s_ctd[Nov_ctdidx]),
                         np.nanmax(ibi_s_ctd[Dec_ctdidx])])
    s_vmin2 = np.floor(s_vmin2/0.1) * 0.1
    s_vmax2 = np.ceil(s_vmax2/0.1) * 0.1

    t_vmin2 = np.nanmin([np.nanmin(tempc[Jul_ctdidx]), np.nanmin(tempc[Oct_ctdidx]), np.nanmin(tempc[Nov_ctdidx]),
                         np.nanmin(tempc[Dec_ctdidx]), np.nanmin(ibi_t_ctd[Jul_ctdidx]),
                         np.nanmin(ibi_t_ctd[Oct_ctdidx]), np.nanmin(ibi_t_ctd[Nov_ctdidx]),
                         np.nanmin(ibi_t_ctd[Dec_ctdidx]), np.nanmin(cel_t_ctd[Jul_ctdidx]),
                         np.nanmin(cel_t_ctd[Oct_ctdidx]), np.nanmin(cel_t_ctd[Nov_ctdidx]),
                         np.nanmin(cel_t_ctd[Dec_ctdidx])])
    t_vmax2 = np.nanmax([np.nanmax(tempc[Jul_ctdidx]), np.nanmax(tempc[Oct_ctdidx]), np.nanmax(tempc[Nov_ctdidx]),
                         np.nanmax(tempc[Dec_ctdidx]), np.nanmax(ibi_t_ctd[Jul_ctdidx]),
                         np.nanmax(ibi_t_ctd[Oct_ctdidx]), np.nanmax(ibi_t_ctd[Nov_ctdidx]),
                         np.nanmax(ibi_t_ctd[Dec_ctdidx]), np.nanmax(cel_t_ctd[Jul_ctdidx]),
                         np.nanmax(cel_t_ctd[Oct_ctdidx]), np.nanmax(cel_t_ctd[Nov_ctdidx]),
                         np.nanmax(cel_t_ctd[Dec_ctdidx])])
    t_vmin2 = np.floor(t_vmin2/0.5) * 0.5
    t_vmax2 = np.ceil(t_vmax2/0.5) * 0.5

    maptsname = hisdir_station + 'croco_TS' + str(bias_yr[0]) + str(bias_mths[mth]).zfill(2) + '.jpg'
    fig, axs = plt.subplots(1, 3, figsize=(7, 5))
    axs[0].scatter(sal_c[Jul_ctdidx], tempc[Jul_ctdidx], s=5, marker=mrk[6])
    axs[0].scatter(sal_c[Oct_ctdidx], tempc[Oct_ctdidx], s=5, marker=mrk[9])
    axs[0].scatter(sal_c[Nov_ctdidx], tempc[Nov_ctdidx], s=5, marker=mrk[10])
    axs[0].scatter(sal_c[Dec_ctdidx], tempc[Dec_ctdidx], s=5, marker=mrk[11])
    axs[0].axis(xmin=s_vmin2, xmax=s_vmax2)
    axs[0].axis(ymin=t_vmin2, ymax=t_vmax2)
    axs[0].set_title('CTD')
    axs[1].scatter(ibi_s_ctd[Jul_ctdidx], ibi_t_ctd[Jul_ctdidx], s=5, marker=mrk[6])
    axs[1].scatter(ibi_s_ctd[Oct_ctdidx], ibi_t_ctd[Oct_ctdidx], s=5, marker=mrk[9])
    axs[1].scatter(ibi_s_ctd[Nov_ctdidx], ibi_t_ctd[Nov_ctdidx], s=5, marker=mrk[10])
    axs[1].scatter(ibi_s_ctd[Dec_ctdidx], ibi_t_ctd[Dec_ctdidx], s=5, marker=mrk[11])
    axs[1].axis(xmin=s_vmin2, xmax=s_vmax2)
    axs[1].axis(ymin=t_vmin2, ymax=t_vmax2)
    axs[1].set_title('IBI')
    july = axs[2].scatter(cel_s_ctd[Jul_ctdidx], cel_t_ctd[Jul_ctdidx], s=5, marker=mrk[6], label='Jul')
    octo = axs[2].scatter(cel_s_ctd[Oct_ctdidx], cel_t_ctd[Oct_ctdidx], s=5, marker=mrk[9], label='Oct')
    nove = axs[2].scatter(cel_s_ctd[Nov_ctdidx], cel_t_ctd[Nov_ctdidx], s=5, marker=mrk[10], label='Nov')
    dece = axs[2].scatter(cel_s_ctd[Dec_ctdidx], cel_t_ctd[Dec_ctdidx], s=5, marker=mrk[11], label='Dec')
    axs[2].axis(xmin=s_vmin2, xmax=s_vmax2)
    axs[2].axis(ymin=t_vmin2, ymax=t_vmax2)
    axs[2].set_title('CS1K')
    handles, labels = axs[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=4)
    plt.savefig(maptsname, dpi='figure', format='jpg')
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
                # hisdir_stats = hisdir_985849
                hisdir_stats = hisdir_985962

                # Labelling files by the fact model outputs now at midnight, not midday
                # Therefore, stats comparing like with like
                VALmname = hisdir_stats + 'croco_mn_VALm' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                VALdname = hisdir_stats + 'croco_mn_VALd' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                VALsmename = hisdir_stats + 'croco_mn_VAL' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + 's_ME.jpg'

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

                tempb = np.array(nchiso['temp'][0:dpm[1], 0, :, :])
                temps = np.array(nchiso['temp'][0:dpm[1], N-1, :, :])
                saltb = np.array(nchiso['salt'][0:dpm[1], 0, :, :])
                salts = np.array(nchiso['salt'][0:dpm[1], len(nchiso['s_rho']) - 1, :, :])
                time = np.array(nchiso['time'][0:dpm[1]])

                temps_sat = np.array(ncvdi['temps'][:, :, :])

                for tslice in range(0, dpm[1], 1):
                    temps_m_slice = temps[tslice, :, :].flatten()
                    temps_s_slice = temps_sat[tslice, :, :].flatten()
                    slice_idx = np.argwhere(temps_m_slice > 0)
                    ncvdi['temps_t_bias'][tslice] = \
                        ((temps_m_slice[slice_idx] - temps_s_slice[slice_idx]).mean() / temps_s_slice[slice_idx].mean())
                    ncvdi['temps_t_ME'][tslice] = (temps_m_slice[slice_idx] - temps_s_slice[slice_idx]).mean()
                    ncvdi['temps_t_RMSD'][tslice] = (((temps_m_slice[slice_idx] -
                                                       temps_s_slice[slice_idx]) ** 2).mean()) ** 0.5
                    ncvdi['temps_t_CC'][tslice] = np.corrcoef(temps_m_slice[slice_idx].flatten(),
                                                              temps_s_slice[slice_idx].flatten())[0, 1]

                for la in range(0, lat_rho.shape[0], 1):
                    for lo in range(0, lon_rho.shape[1], 1):
                        temps_m_slice = temps[:, la, lo].flatten()
                        temps_s_slice = temps_sat[:, la, lo].flatten()
                        if temps_m_slice.any() == 0:
                            ncvdi['temps_s_bias'][la, lo] = 0
                            ncvdi['temps_s_ME'][la, lo] = 0
                            ncvdi['temps_s_RMSD'][la, lo] = 0
                            ncvdi['temps_s_CC'][la, lo] = 1
                        else:
                            ncvdi['temps_s_bias'][la, lo] = \
                                ((temps_m_slice - temps_s_slice).mean() / temps_s_slice.mean())
                            ncvdi['temps_s_ME'][la, lo] = (temps_m_slice -
                                                           temps_s_slice).mean()
                            ncvdi['temps_s_RMSD'][la, lo] = (((temps_m_slice - temps_s_slice) ** 2).mean()) ** 0.5
                            ncvdi['temps_s_CC'][la, lo] = np.corrcoef(temps_m_slice.flatten(),
                                                                      temps_s_slice.flatten())[0, 1]

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
