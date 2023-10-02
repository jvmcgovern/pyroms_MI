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
    # ini_date = '20170101'

    IBI_dir = '/media/dskone/CELTIC/CMEMS_IBI/'
    IBI_prefix = 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_'
    IBI_post = '_R20201201_RE01.nc'

    # GLO_dir = '/media/dskthree/CMEMS_GLO/'
    # GLO_prefix = 'mercatorglorys12v1_gl12_mean_'
    # GLO_post = '_R*_crop.nc'
    hisdir_main = '/media/dskone/CELTIC/'
    hisdir_IBI_JS = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_JS/'
    hisdir_979519 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_979519/'
    hisdir_979685 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_979685/'
    hisdir_979807 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_979807/'
    hisdir_980116 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_980116/'
    hisdir_980558 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_980558/'
    hisdir_980771 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_980771/'

    hisdir_983014 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_983014/'
    hisdir_983076 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_983076/'
    hisdir_983204 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_983204/'
    hisdir_983286 = '/media/dskone/CELTIC/CROCO_PHY_1p1_17_TSUV_IBI_983286/'

    hisdir_980885 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_980885/'
    hisdir_981213 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_981213/'
    hisdir_981356 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_981356/'
    hisdir_981557 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_981557/'
    hisdir_981848 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_981848/'
    hisdir_982242 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_982242/'
    hisdir_982365 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_982365/'
    hisdir_982490 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_982490/'

    hisdir_983543 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_983543/'
    hisdir_983672 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_983672/'
    hisdir_983790 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_983790/'
    hisdir_983982 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_983982/'
    hisdir_985519 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_17_TSUV_IBI_985519/'
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

    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_34f3_19f6_97e1.nc'

    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2013.nc'
    ctdfil = '/media/dskone/CELTIC/IMI_CTD_2014.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2015.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2016.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2017.nc'
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
    tempm_IJS = np.zeros_like(tempc)
    sal_m_IJS = np.zeros_like(sal_c)
    tempm_979519 = np.zeros_like(tempc)
    sal_m_979519 = np.zeros_like(tempc)
    tempm_979685 = np.zeros_like(tempc)
    sal_m_979685 = np.zeros_like(tempc)
    tempm_979807 = np.zeros_like(tempc)
    sal_m_979807 = np.zeros_like(tempc)
    tempm_980116 = np.zeros_like(tempc)
    sal_m_980116 = np.zeros_like(tempc)
    tempm_980558 = np.zeros_like(tempc)
    sal_m_980558 = np.zeros_like(tempc)
    tempm_980771 = np.zeros_like(tempc)
    sal_m_980771 = np.zeros_like(tempc)

    tempm_983014 = np.zeros_like(tempc)
    sal_m_983014 = np.zeros_like(tempc)
    tempm_983076 = np.zeros_like(tempc)
    sal_m_983076 = np.zeros_like(tempc)
    tempm_983204 = np.zeros_like(tempc)
    sal_m_983204 = np.zeros_like(tempc)
    tempm_983286 = np.zeros_like(tempc)
    sal_m_983286 = np.zeros_like(tempc)

    tempm_980885 = np.zeros_like(tempc)
    sal_m_980885 = np.zeros_like(tempc)
    tempm_981213 = np.zeros_like(tempc)
    sal_m_981213 = np.zeros_like(tempc)
    tempm_981356 = np.zeros_like(tempc)
    sal_m_981356 = np.zeros_like(tempc)
    tempm_981557 = np.zeros_like(tempc)
    sal_m_981557 = np.zeros_like(tempc)
    tempm_981848 = np.zeros_like(tempc)
    sal_m_981848 = np.zeros_like(tempc)
    tempm_982242 = np.zeros_like(tempc)
    sal_m_982242 = np.zeros_like(tempc)
    tempm_982365 = np.zeros_like(tempc)
    sal_m_982365 = np.zeros_like(tempc)
    tempm_982490 = np.zeros_like(tempc)
    sal_m_982490 = np.zeros_like(tempc)

    tempm_983543 = np.zeros_like(tempc)
    sal_m_983543 = np.zeros_like(tempc)
    tempm_983672 = np.zeros_like(tempc)
    sal_m_983672 = np.zeros_like(tempc)
    tempm_983790 = np.zeros_like(tempc)
    sal_m_983790 = np.zeros_like(tempc)
    tempm_983982 = np.zeros_like(tempc)
    sal_m_983982 = np.zeros_like(tempc)

    tempm_990003 = np.zeros_like(tempc)
    sal_m_990003 = np.zeros_like(tempc)
    tempm_995440 = np.zeros_like(tempc)
    sal_m_995440 = np.zeros_like(tempc)

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

        # his_979519 = hisdir_979519 + \
        #           his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_979519 = netcdf(his_979519, 'r')
        # nchis_979519o = nchis_979519.variables
        # # Extract vertical section from history file
        # temp979519 = np.array(nchis_979519o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt979519 = np.array(nchis_979519o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta979519 = np.array(nchis_979519o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_979685 = hisdir_979685 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_979685 = netcdf(his_979685, 'r')
        # nchis_979685o = nchis_979685.variables
        # # Extract vertical section from history file
        # temp979685 = np.array(nchis_979685o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt979685 = np.array(nchis_979685o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta979685 = np.array(nchis_979685o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_979807 = hisdir_979807 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_979807 = netcdf(his_979807, 'r')
        # nchis_979807o = nchis_979807.variables
        # # Extract vertical section from history file
        # temp979807 = np.array(nchis_979807o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt979807 = np.array(nchis_979807o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta979807 = np.array(nchis_979807o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_980116 = hisdir_980116 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_980116 = netcdf(his_980116, 'r')
        # nchis_980116o = nchis_980116.variables
        # # Extract vertical section from history file
        # temp980116 = np.array(nchis_980116o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt980116 = np.array(nchis_980116o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta980116 = np.array(nchis_980116o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_980558 = hisdir_980558 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_980558 = netcdf(his_980558, 'r')
        # nchis_980558o = nchis_980558.variables
        # # Extract vertical section from history file
        # temp980558 = np.array(nchis_980558o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt980558 = np.array(nchis_980558o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta980558 = np.array(nchis_980558o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_980771 = hisdir_980771 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_980771 = netcdf(his_980771, 'r')
        # nchis_980771o = nchis_980771.variables
        # # Extract vertical section from history file
        # temp980771 = np.array(nchis_980771o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt980771 = np.array(nchis_980771o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta980771 = np.array(nchis_980771o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_980885 = hisdir_980885 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_980885 = netcdf(his_980885, 'r')
        # nchis_980885o = nchis_980885.variables
        # # Extract vertical section from history file
        # temp980885 = np.array(nchis_980885o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt980885 = np.array(nchis_980885o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta980885 = np.array(nchis_980885o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_981213 = hisdir_981213 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_981213 = netcdf(his_981213, 'r')
        # nchis_981213o = nchis_981213.variables
        # # Extract vertical section from history file
        # temp981213 = np.array(nchis_981213o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt981213 = np.array(nchis_981213o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta981213 = np.array(nchis_981213o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_981356 = hisdir_981356 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_981356 = netcdf(his_981356, 'r')
        # nchis_981356o = nchis_981356.variables
        # # Extract vertical section from history file
        # temp981356 = np.array(nchis_981356o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt981356 = np.array(nchis_981356o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta981356 = np.array(nchis_981356o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_981557 = hisdir_981557 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_981557 = netcdf(his_981557, 'r')
        # nchis_981557o = nchis_981557.variables
        # # Extract vertical section from history file
        # temp981557 = np.array(nchis_981557o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt981557 = np.array(nchis_981557o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta981557 = np.array(nchis_981557o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_981848 = hisdir_981848 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_981848 = netcdf(his_981848, 'r')
        # nchis_981848o = nchis_981848.variables
        # # Extract vertical section from history file
        # temp981848 = np.array(nchis_981848o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt981848 = np.array(nchis_981848o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta981848 = np.array(nchis_981848o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_982242 = hisdir_982242 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_982242 = netcdf(his_982242, 'r')
        # nchis_982242o = nchis_982242.variables
        # # Extract vertical section from history file
        # temp982242 = np.array(nchis_982242o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt982242 = np.array(nchis_982242o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta982242 = np.array(nchis_982242o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_982365 = hisdir_982365 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_982365 = netcdf(his_982365, 'r')
        # nchis_982365o = nchis_982365.variables
        # # Extract vertical section from history file
        # temp982365 = np.array(nchis_982365o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt982365 = np.array(nchis_982365o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta982365 = np.array(nchis_982365o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_982490 = hisdir_982490 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_982490 = netcdf(his_982490, 'r')
        # nchis_982490o = nchis_982490.variables
        # # Extract vertical section from history file
        # temp982490 = np.array(nchis_982490o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt982490 = np.array(nchis_982490o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta982490 = np.array(nchis_982490o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983014 = hisdir_983014 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983014 = netcdf(his_983014, 'r')
        # nchis_983014o = nchis_983014.variables
        # # Extract vertical section from history file
        # temp983014 = np.array(nchis_983014o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983014 = np.array(nchis_983014o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983014 = np.array(nchis_983014o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983076 = hisdir_983076 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983076 = netcdf(his_983076, 'r')
        # nchis_983076o = nchis_983076.variables
        # # Extract vertical section from history file
        # temp983076 = np.array(nchis_983076o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983076 = np.array(nchis_983076o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983076 = np.array(nchis_983076o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983204 = hisdir_983204 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983204 = netcdf(his_983204, 'r')
        # nchis_983204o = nchis_983204.variables
        # # Extract vertical section from history file
        # temp983204 = np.array(nchis_983204o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983204 = np.array(nchis_983204o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983204 = np.array(nchis_983204o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983286 = hisdir_983286 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983286 = netcdf(his_983286, 'r')
        # nchis_983286o = nchis_983286.variables
        # # Extract vertical section from history file
        # temp983286 = np.array(nchis_983286o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983286 = np.array(nchis_983286o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983286 = np.array(nchis_983286o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983543 = hisdir_983543 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983543 = netcdf(his_983543, 'r')
        # nchis_983543o = nchis_983543.variables
        # # Extract vertical section from history file
        # temp983543 = np.array(nchis_983543o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983543 = np.array(nchis_983543o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983543 = np.array(nchis_983543o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983672 = hisdir_983672 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983672 = netcdf(his_983672, 'r')
        # nchis_983672o = nchis_983672.variables
        # # Extract vertical section from history file
        # temp983672 = np.array(nchis_983672o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983672 = np.array(nchis_983672o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983672 = np.array(nchis_983672o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983790 = hisdir_983790 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983790 = netcdf(his_983790, 'r')
        # nchis_983790o = nchis_983790.variables
        # # Extract vertical section from history file
        # temp983790 = np.array(nchis_983790o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983790 = np.array(nchis_983790o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983790 = np.array(nchis_983790o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_983982 = hisdir_983982 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_983982 = netcdf(his_983982, 'r')
        # nchis_983982o = nchis_983982.variables
        # # Extract vertical section from history file
        # temp983982 = np.array(nchis_983982o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt983982 = np.array(nchis_983982o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta983982 = np.array(nchis_983982o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_990003 = hisdir_990003 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_990003 = netcdf(his_990003, 'r')
        # nchis_990003o = nchis_990003.variables
        # # Extract vertical section from history file
        # temp990003 = np.array(nchis_990003o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt990003 = np.array(nchis_990003o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta990003 = np.array(nchis_990003o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        #
        # his_995440 = hisdir_995440 + \
        #              his_prefix + 'Y' + str(timec_dates_lls[dl, 0]) + 'M' + str(timec_dates_lls[dl, 1]).zfill(2) + '.nc'
        # nchis_995440 = netcdf(his_995440, 'r')
        # nchis_995440o = nchis_995440.variables
        # # Extract vertical section from history file
        # temp995440 = np.array(nchis_995440o['temp']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # salt995440 = np.array(nchis_995440o['salt']
        #                       [timec_dates_lls[dl, 2], :, timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()
        # zeta995440 = np.array(nchis_995440o['zeta']
        #                       [timec_dates_lls[dl, 2], timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]]).flatten()

        # # Interpolate section data onto depths from the CTD cast for that date
        # # z = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], 0, theta_s, theta_b, hc, N, 'w', vtransform)
        # z = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], 0, theta_s, theta_b, hc, N, 'w', vtransform)
        # z1p1 = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], zeta983286,
        #                   theta_s, theta_b, hc, N, 'r', vtransform)
        # z1p1[0] = np.floor(z1p1[0])
        # z1p2p1 = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], zeta983982,
        #                     theta_s, theta_b, hc, N, 'r', vtransform)
        # z1p2p1[0] = np.floor(z1p2p1[0])
        # z1p2p1_hc50m = vgrd.zlevs(h[timec_dates_lls[dl, 3], timec_dates_lls[dl, 4]], zeta995440,
        #                     theta_s, theta_b, 50, N, 'r', vtransform)
        # z1p2p1_hc50m[0] = np.floor(z1p2p1_hc50m[0])
        #
        # z[0] = np.floor(z[0])
        # zmin = z.min()
        # zmax = z.max()

        # ftIJS = interpolate.interp1d(z1p1, tempIJS, fill_value="extrapolate")
        # # ftIJS = interpolate.interp1d(z, np.append(tempIJS[0], tempIJS))
        # tmtotijs = np.zeros_like(depths)
        # tmtotijs[(depths >= zmin) & (depths <= zmax)] = ftIJS(depths[(depths >= zmin) & (depths <= zmax)])

        # # ft979519 = interpolate.interp1d(z, np.append(temp979519[0], temp979519))
        # ft979519 = interpolate.interp1d(z1p1, temp979519, fill_value="extrapolate")
        # tmtot979519 = np.zeros_like(depths)
        # tmtot979519[(depths >= zmin) & (depths <= zmax)] = ft979519(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft979685 = interpolate.interp1d(z, np.append(temp979685[0], temp979685))
        # ft979685 = interpolate.interp1d(z1p1, temp979685, fill_value="extrapolate")
        # tmtot979685 = np.zeros_like(depths)
        # tmtot979685[(depths >= zmin) & (depths <= zmax)] = ft979685(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft979807 = interpolate.interp1d(z, np.append(temp979807[0], temp979807))
        # ft979807 = interpolate.interp1d(z1p1, temp979807, fill_value="extrapolate")
        # tmtot979807 = np.zeros_like(depths)
        # tmtot979807[(depths >= zmin) & (depths <= zmax)] = ft979807(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft980116 = interpolate.interp1d(z, np.append(temp980116[0], temp980116))
        # ft980116 = interpolate.interp1d(z1p1, temp980116, fill_value="extrapolate")
        # tmtot980116 = np.zeros_like(depths)
        # tmtot980116[(depths >= zmin) & (depths <= zmax)] = ft980116(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft980558 = interpolate.interp1d(z, np.append(temp980558[0], temp980558))
        # ft980558 = interpolate.interp1d(z1p1, temp980558, fill_value="extrapolate")
        # tmtot980558 = np.zeros_like(depths)
        # tmtot980558[(depths >= zmin) & (depths <= zmax)] = ft980558(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft980771 = interpolate.interp1d(z, np.append(temp980771[0], temp980771))
        # ft980771 = interpolate.interp1d(z1p1, temp980771, fill_value="extrapolate")
        # tmtot980771 = np.zeros_like(depths)
        # tmtot980771[(depths >= zmin) & (depths <= zmax)] = ft980771(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft980885 = interpolate.interp1d(z, np.append(temp980885[0], temp980885))
        # ft980885 = interpolate.interp1d(z1p2p1, temp980885, fill_value="extrapolate")
        # tmtot980885 = np.zeros_like(depths)
        # tmtot980885[(depths >= zmin) & (depths <= zmax)] = ft980885(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft981213 = interpolate.interp1d(z, np.append(temp981213[0], temp981213))
        # ft981213 = interpolate.interp1d(z1p2p1, temp981213, fill_value="extrapolate")
        # tmtot981213 = np.zeros_like(depths)
        # tmtot981213[(depths >= zmin) & (depths <= zmax)] = ft981213(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft981356 = interpolate.interp1d(z, np.append(temp981356[0], temp981356))
        # ft981356 = interpolate.interp1d(z1p2p1, temp981356, fill_value="extrapolate")
        # tmtot981356 = np.zeros_like(depths)
        # tmtot981356[(depths >= zmin) & (depths <= zmax)] = ft981356(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft981557 = interpolate.interp1d(z, np.append(temp981557[0], temp981557))
        # ft981557 = interpolate.interp1d(z1p2p1, temp981557, fill_value="extrapolate")
        # tmtot981557 = np.zeros_like(depths)
        # tmtot981557[(depths >= zmin) & (depths <= zmax)] = ft981557(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft981848 = interpolate.interp1d(z, np.append(temp981848[0], temp981848))
        # ft981848 = interpolate.interp1d(z1p2p1, temp981848, fill_value="extrapolate")
        # tmtot981848 = np.zeros_like(depths)
        # tmtot981848[(depths >= zmin) & (depths <= zmax)] = ft981848(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft982242 = interpolate.interp1d(z, np.append(temp982242[0], temp982242))
        # ft982242 = interpolate.interp1d(z1p2p1, temp982242, fill_value="extrapolate")
        # tmtot982242 = np.zeros_like(depths)
        # tmtot982242[(depths >= zmin) & (depths <= zmax)] = ft982242(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft982365 = interpolate.interp1d(z, np.append(temp982365[0], temp982365))
        # ft982365 = interpolate.interp1d(z1p2p1, temp982365, fill_value="extrapolate")
        # tmtot982365 = np.zeros_like(depths)
        # tmtot982365[(depths >= zmin) & (depths <= zmax)] = ft982365(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft982490 = interpolate.interp1d(z, np.append(temp982490[0], temp982490))
        # ft982490 = interpolate.interp1d(z1p2p1, temp982490, fill_value="extrapolate")
        # tmtot982490 = np.zeros_like(depths)
        # tmtot982490[(depths >= zmin) & (depths <= zmax)] = ft982490(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983014 = interpolate.interp1d(z, np.append(temp983014[0], temp983014))
        # ft983014 = interpolate.interp1d(z1p1, temp983014, fill_value="extrapolate")
        # tmtot983014 = np.zeros_like(depths)
        # tmtot983014[(depths >= zmin) & (depths <= zmax)] = ft983014(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983076 = interpolate.interp1d(z, np.append(temp983076[0], temp983076))
        # ft983076 = interpolate.interp1d(z1p1, temp983076, fill_value="extrapolate")
        # tmtot983076 = np.zeros_like(depths)
        # tmtot983076[(depths >= zmin) & (depths <= zmax)] = ft983076(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983204 = interpolate.interp1d(z, np.append(temp983204[0], temp983204))
        # ft983204 = interpolate.interp1d(z1p1, temp983204, fill_value="extrapolate")
        # tmtot983204 = np.zeros_like(depths)
        # tmtot983204[(depths >= zmin) & (depths <= zmax)] = ft983204(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983286 = interpolate.interp1d(z, np.append(temp983286[0], temp983286))
        # ft983286 = interpolate.interp1d(z1p1, temp983286, fill_value="extrapolate")
        # tmtot983286 = np.zeros_like(depths)
        # tmtot983286[(depths >= zmin) & (depths <= zmax)] = ft983286(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983543 = interpolate.interp1d(z, np.append(temp983543[0], temp983543))
        # ft983543 = interpolate.interp1d(z1p2p1, temp983543, fill_value="extrapolate")
        # tmtot983543 = np.zeros_like(depths)
        # tmtot983543[(depths >= zmin) & (depths <= zmax)] = ft983543(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983672 = interpolate.interp1d(z, np.append(temp983672[0], temp983672))
        # ft983672 = interpolate.interp1d(z1p2p1, temp983672, fill_value="extrapolate")
        # tmtot983672 = np.zeros_like(depths)
        # tmtot983672[(depths >= zmin) & (depths <= zmax)] = ft983672(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983790 = interpolate.interp1d(z, np.append(temp983790[0], temp983790))
        # ft983790 = interpolate.interp1d(z1p2p1, temp983790, fill_value="extrapolate")
        # tmtot983790 = np.zeros_like(depths)
        # tmtot983790[(depths >= zmin) & (depths <= zmax)] = ft983790(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft983982 = interpolate.interp1d(z, np.append(temp983982[0], temp983982))
        # ft983982 = interpolate.interp1d(z1p2p1, temp983982, fill_value="extrapolate")
        # tmtot983982 = np.zeros_like(depths)
        # tmtot983982[(depths >= zmin) & (depths <= zmax)] = ft983982(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft990003 = interpolate.interp1d(z, np.append(temp990003[0], temp990003))
        # ft990003 = interpolate.interp1d(z1p2p1, temp990003, fill_value="extrapolate")
        # tmtot990003 = np.zeros_like(depths)
        # tmtot990003[(depths >= zmin) & (depths <= zmax)] = ft990003(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # # ft995440 = interpolate.interp1d(z, np.append(temp995440[0], temp995440))
        # ft995440 = interpolate.interp1d(z1p2p1_hc50m, temp995440, fill_value="extrapolate")
        # tmtot995440 = np.zeros_like(depths)
        # tmtot995440[(depths >= zmin) & (depths <= zmax)] = ft995440(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # if depths.min() < zmin:
        #     # tmtotijs[(depths < zmin)] = np.nan
        #     tmtot979519[(depths < zmin)] = np.nan
        #     tmtot979685[(depths < zmin)] = np.nan
        #     tmtot979807[(depths < zmin)] = np.nan
        #     tmtot980116[(depths < zmin)] = np.nan
        #     tmtot980558[(depths < zmin)] = np.nan
        #     tmtot980771[(depths < zmin)] = np.nan
        #     tmtot980885[(depths < zmin)] = np.nan
        #     tmtot981213[(depths < zmin)] = np.nan
        #     tmtot981356[(depths < zmin)] = np.nan
        #     tmtot981557[(depths < zmin)] = np.nan
        #     tmtot981848[(depths < zmin)] = np.nan
        #     tmtot982242[(depths < zmin)] = np.nan
        #     tmtot982365[(depths < zmin)] = np.nan
        #     tmtot982490[(depths < zmin)] = np.nan
        #     tmtot983014[(depths < zmin)] = np.nan
        #     tmtot983076[(depths < zmin)] = np.nan
        #     tmtot983204[(depths < zmin)] = np.nan
        #     tmtot983286[(depths < zmin)] = np.nan
        #     tmtot983543[(depths < zmin)] = np.nan
        #     tmtot983672[(depths < zmin)] = np.nan
        #     tmtot983790[(depths < zmin)] = np.nan
        #     tmtot983982[(depths < zmin)] = np.nan
        #     tmtot990003[(depths < zmin)] = np.nan
        #     tmtot995440[(depths < zmin)] = np.nan
        # if depths.max() > zmax:
        #     # tmtotijs[(depths > zmax)] = np.nan
        #     tmtot979519[(depths > zmax)] = np.nan
        #     tmtot979685[(depths > zmax)] = np.nan
        #     tmtot979807[(depths > zmax)] = np.nan
        #     tmtot980116[(depths > zmax)] = np.nan
        #     tmtot980558[(depths > zmax)] = np.nan
        #     tmtot980771[(depths > zmax)] = np.nan
        #     tmtot980885[(depths > zmax)] = np.nan
        #     tmtot981213[(depths > zmax)] = np.nan
        #     tmtot981356[(depths > zmax)] = np.nan
        #     tmtot981557[(depths > zmax)] = np.nan
        #     tmtot981848[(depths > zmax)] = np.nan
        #     tmtot982242[(depths > zmax)] = np.nan
        #     tmtot982365[(depths > zmax)] = np.nan
        #     tmtot982490[(depths > zmax)] = np.nan
        #     tmtot983014[(depths > zmax)] = np.nan
        #     tmtot983076[(depths > zmax)] = np.nan
        #     tmtot983204[(depths > zmax)] = np.nan
        #     tmtot983286[(depths > zmax)] = np.nan
        #     tmtot983543[(depths > zmax)] = np.nan
        #     tmtot983672[(depths > zmax)] = np.nan
        #     tmtot983790[(depths > zmax)] = np.nan
        #     tmtot983982[(depths > zmax)] = np.nan
        #     tmtot990003[(depths > zmax)] = np.nan
        #     tmtot995440[(depths > zmax)] = np.nan
        # # tempm_IJS[(timec_year == timec_dates_lls[dl, 0]) &
        # #           (timec_month == timec_dates_lls[dl, 1]) &
        # #           (timec_day == timec_dates_lls[dl, 2]) &
        # #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        # #           (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtotijs
        # tempm_979519[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot979519
        # tempm_979685[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot979685
        # tempm_979807[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot979807
        # tempm_980116[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot980116
        # tempm_980558[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot980558
        # tempm_980771[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot980771
        # tempm_980885[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot980885
        # tempm_981213[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot981213
        # tempm_981356[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot981356
        # tempm_981557[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot981557
        # tempm_981848[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot981848
        # tempm_982242[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot982242
        # tempm_982365[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot982365
        # tempm_982490[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot982490
        # tempm_983014[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983014
        # tempm_983076[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983076
        # tempm_983204[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983204
        # tempm_983286[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983286
        # tempm_983543[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983543
        # tempm_983672[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983672
        # tempm_983790[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983790
        # tempm_983982[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot983982
        # tempm_990003[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot990003
        # tempm_995440[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = tmtot995440
        #
        # # fsIJS = interpolate.interp1d(z, np.append(saltIJS[0], saltIJS))
        # # smtotijs = np.zeros_like(depths)
        # # smtotijs[(depths >= zmin) & (depths <= zmax)] = fsIJS(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs979519 = interpolate.interp1d(z, np.append(salt979519[0], salt979519))
        # smtot979519 = np.zeros_like(depths)
        # smtot979519[(depths >= zmin) & (depths <= zmax)] = fs979519(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs979685 = interpolate.interp1d(z, np.append(salt979685[0], salt979685))
        # smtot979685 = np.zeros_like(depths)
        # smtot979685[(depths >= zmin) & (depths <= zmax)] = fs979685(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs979807 = interpolate.interp1d(z, np.append(salt979807[0], salt979807))
        # smtot979807 = np.zeros_like(depths)
        # smtot979807[(depths >= zmin) & (depths <= zmax)] = fs979807(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs980116 = interpolate.interp1d(z, np.append(salt980116[0], salt980116))
        # smtot980116 = np.zeros_like(depths)
        # smtot980116[(depths >= zmin) & (depths <= zmax)] = fs980116(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs980558 = interpolate.interp1d(z, np.append(salt980558[0], salt980558))
        # smtot980558 = np.zeros_like(depths)
        # smtot980558[(depths >= zmin) & (depths <= zmax)] = fs980558(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs980771 = interpolate.interp1d(z, np.append(salt980771[0], salt980771))
        # smtot980771 = np.zeros_like(depths)
        # smtot980771[(depths >= zmin) & (depths <= zmax)] = fs980771(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs980885 = interpolate.interp1d(z1p2p1, salt980885, fill_value="extrapolate")
        # smtot980885 = np.zeros_like(depths)
        # smtot980885[(depths >= zmin) & (depths <= zmax)] = fs980885(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs981213 = interpolate.interp1d(z1p2p1, salt981213, fill_value="extrapolate")
        # smtot981213 = np.zeros_like(depths)
        # smtot981213[(depths >= zmin) & (depths <= zmax)] = fs981213(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs981356 = interpolate.interp1d(z1p2p1, salt981356, fill_value="extrapolate")
        # smtot981356 = np.zeros_like(depths)
        # smtot981356[(depths >= zmin) & (depths <= zmax)] = fs981356(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs981557 = interpolate.interp1d(z1p2p1, salt981557, fill_value="extrapolate")
        # smtot981557 = np.zeros_like(depths)
        # smtot981557[(depths >= zmin) & (depths <= zmax)] = fs981557(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs981848 = interpolate.interp1d(z1p2p1, salt981848, fill_value="extrapolate")
        # smtot981848 = np.zeros_like(depths)
        # smtot981848[(depths >= zmin) & (depths <= zmax)] = fs981848(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs982242 = interpolate.interp1d(z1p2p1, salt982242, fill_value="extrapolate")
        # smtot982242 = np.zeros_like(depths)
        # smtot982242[(depths >= zmin) & (depths <= zmax)] = fs982242(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs982365 = interpolate.interp1d(z1p2p1, salt982365, fill_value="extrapolate")
        # smtot982365 = np.zeros_like(depths)
        # smtot982365[(depths >= zmin) & (depths <= zmax)] = fs982365(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs982490 = interpolate.interp1d(z1p2p1, salt982490, fill_value="extrapolate")
        # smtot982490 = np.zeros_like(depths)
        # smtot982490[(depths >= zmin) & (depths <= zmax)] = fs982490(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983014 = interpolate.interp1d(z, np.append(salt983014[0], salt983014))
        # smtot983014 = np.zeros_like(depths)
        # smtot983014[(depths >= zmin) & (depths <= zmax)] = fs983014(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983076 = interpolate.interp1d(z, np.append(salt983076[0], salt983076))
        # smtot983076 = np.zeros_like(depths)
        # smtot983076[(depths >= zmin) & (depths <= zmax)] = fs983076(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983204 = interpolate.interp1d(z, np.append(salt983204[0], salt983204))
        # smtot983204 = np.zeros_like(depths)
        # smtot983204[(depths >= zmin) & (depths <= zmax)] = fs983204(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983286 = interpolate.interp1d(z, np.append(salt983286[0], salt983286))
        # smtot983286 = np.zeros_like(depths)
        # smtot983286[(depths >= zmin) & (depths <= zmax)] = fs983286(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983543 = interpolate.interp1d(z1p2p1, salt983543, fill_value="extrapolate")
        # smtot983543 = np.zeros_like(depths)
        # smtot983543[(depths >= zmin) & (depths <= zmax)] = fs983543(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983672 = interpolate.interp1d(z1p2p1, salt983672, fill_value="extrapolate")
        # smtot983672 = np.zeros_like(depths)
        # smtot983672[(depths >= zmin) & (depths <= zmax)] = fs983672(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983790 = interpolate.interp1d(z1p2p1, salt983790, fill_value="extrapolate")
        # smtot983790 = np.zeros_like(depths)
        # smtot983790[(depths >= zmin) & (depths <= zmax)] = fs983790(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs983982 = interpolate.interp1d(z1p2p1, salt983982, fill_value="extrapolate")
        # smtot983982 = np.zeros_like(depths)
        # smtot983982[(depths >= zmin) & (depths <= zmax)] = fs983982(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs990003 = interpolate.interp1d(z1p2p1, salt990003, fill_value="extrapolate")
        # smtot990003 = np.zeros_like(depths)
        # smtot990003[(depths >= zmin) & (depths <= zmax)] = fs990003(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # fs995440 = interpolate.interp1d(z1p2p1_hc50m, salt995440, fill_value="extrapolate")
        # smtot995440 = np.zeros_like(depths)
        # smtot995440[(depths >= zmin) & (depths <= zmax)] = fs995440(depths[(depths >= zmin) & (depths <= zmax)])
        #
        # if depths.min() < zmin:
        #     # smtotijs[(depths < zmin)] = np.nan
        #     smtot979519[(depths < zmin)] = np.nan
        #     smtot979685[(depths < zmin)] = np.nan
        #     smtot979807[(depths < zmin)] = np.nan
        #     smtot980116[(depths < zmin)] = np.nan
        #     smtot980558[(depths < zmin)] = np.nan
        #     smtot980771[(depths < zmin)] = np.nan
        #     smtot980885[(depths < zmin)] = np.nan
        #     smtot981213[(depths < zmin)] = np.nan
        #     smtot981356[(depths < zmin)] = np.nan
        #     smtot981557[(depths < zmin)] = np.nan
        #     smtot981848[(depths < zmin)] = np.nan
        #     smtot982242[(depths < zmin)] = np.nan
        #     smtot982365[(depths < zmin)] = np.nan
        #     smtot982490[(depths < zmin)] = np.nan
        #     smtot983014[(depths < zmin)] = np.nan
        #     smtot983076[(depths < zmin)] = np.nan
        #     smtot983204[(depths < zmin)] = np.nan
        #     smtot983286[(depths < zmin)] = np.nan
        #     smtot983543[(depths < zmin)] = np.nan
        #     smtot983672[(depths < zmin)] = np.nan
        #     smtot983790[(depths < zmin)] = np.nan
        #     smtot983982[(depths < zmin)] = np.nan
        #     smtot990003[(depths < zmin)] = np.nan
        #     smtot995440[(depths < zmin)] = np.nan
        # if depths.max() > zmax:
        #     # smtotijs[(depths > zmax)] = np.nan
        #     smtot979519[(depths > zmax)] = np.nan
        #     smtot979685[(depths > zmax)] = np.nan
        #     smtot979807[(depths > zmax)] = np.nan
        #     smtot980116[(depths > zmax)] = np.nan
        #     smtot980558[(depths > zmax)] = np.nan
        #     smtot980771[(depths > zmax)] = np.nan
        #     smtot980885[(depths > zmax)] = np.nan
        #     smtot981213[(depths > zmax)] = np.nan
        #     smtot981356[(depths > zmax)] = np.nan
        #     smtot981557[(depths > zmax)] = np.nan
        #     smtot981848[(depths > zmax)] = np.nan
        #     smtot982242[(depths > zmax)] = np.nan
        #     smtot982365[(depths > zmax)] = np.nan
        #     smtot982490[(depths > zmax)] = np.nan
        #     smtot983014[(depths > zmax)] = np.nan
        #     smtot983076[(depths > zmax)] = np.nan
        #     smtot983204[(depths > zmax)] = np.nan
        #     smtot983286[(depths > zmax)] = np.nan
        #     smtot983543[(depths > zmax)] = np.nan
        #     smtot983672[(depths > zmax)] = np.nan
        #     smtot983790[(depths > zmax)] = np.nan
        #     smtot983982[(depths > zmax)] = np.nan
        #     smtot990003[(depths > zmax)] = np.nan
        #     smtot995440[(depths > zmax)] = np.nan
        # # sal_m_IJS[(timec_year == timec_dates_lls[dl, 0]) &
        # #           (timec_month == timec_dates_lls[dl, 1]) &
        # #           (timec_day == timec_dates_lls[dl, 2]) &
        # #           (ctd_latidx == timec_dates_lls[dl, 3]) &
        # #           (ctd_lonidx == timec_dates_lls[dl, 4])] = smtotijs
        # sal_m_979519[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot979519
        # sal_m_979685[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot979685
        # sal_m_979807[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot979807
        # sal_m_980116[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot980116
        # sal_m_980558[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot980558
        # sal_m_980771[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot980771
        # sal_m_980885[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot980885
        # sal_m_981213[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot981213
        # sal_m_981356[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot981356
        # sal_m_981557[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot981557
        # sal_m_981848[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot981848
        # sal_m_982242[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot982242
        # sal_m_982365[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot982365
        # sal_m_982490[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot982490
        # sal_m_983014[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983014
        # sal_m_983076[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983076
        # sal_m_983204[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983204
        # sal_m_983286[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983286
        # sal_m_983543[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983543
        # sal_m_983672[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983672
        # sal_m_983790[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983790
        # sal_m_983982[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot983982
        # sal_m_990003[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot990003
        # sal_m_995440[(timec_year == timec_dates_lls[dl, 0]) &
        #              (timec_month == timec_dates_lls[dl, 1]) &
        #              (timec_day == timec_dates_lls[dl, 2]) &
        #              (ctd_latidx == timec_dates_lls[dl, 3]) &
        #              (ctd_lonidx == timec_dates_lls[dl, 4])] = smtot995440
        # # nchis_IJS.close()
        # nchis_979519.close()
        # nchis_979685.close()
        # nchis_979807.close()
        # nchis_980116.close()
        # nchis_980558.close()
        # nchis_980771.close()
        # nchis_980885.close()
        # nchis_981213.close()
        # nchis_981356.close()
        # nchis_981557.close()
        # nchis_981848.close()
        # nchis_982242.close()
        # nchis_982365.close()
        # nchis_982490.close()
        # nchis_983014.close()
        # nchis_983076.close()
        # nchis_983204.close()
        # nchis_983286.close()
        # nchis_983543.close()
        # nchis_983672.close()
        # nchis_983790.close()
        # nchis_983982.close()
        # nchis_990003.close()
        # nchis_995440.close()

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
    # tempm_IJS = tempm_IJS[tsrmnz]
    # sal_m_IJS = sal_m_IJS[tsrmnz]
    # tempm_979519 = tempm_979519[tsrmnz]
    # sal_m_979519 = sal_m_979519[tsrmnz]
    # tempm_979685 = tempm_979685[tsrmnz]
    # sal_m_979685 = sal_m_979685[tsrmnz]
    # tempm_979807 = tempm_979807[tsrmnz]
    # sal_m_979807 = sal_m_979807[tsrmnz]
    # tempm_980116 = tempm_980116[tsrmnz]
    # sal_m_980116 = sal_m_980116[tsrmnz]
    # tempm_980558 = tempm_980558[tsrmnz]
    # sal_m_980558 = sal_m_980558[tsrmnz]
    # tempm_980771 = tempm_980771[tsrmnz]
    # sal_m_980771 = sal_m_980771[tsrmnz]
    # tempm_980885 = tempm_980885[tsrmnz]
    # sal_m_980885 = sal_m_980885[tsrmnz]
    # tempm_981213 = tempm_981213[tsrmnz]
    # sal_m_981213 = sal_m_981213[tsrmnz]
    # tempm_981356 = tempm_981356[tsrmnz]
    # sal_m_981356 = sal_m_981356[tsrmnz]
    # tempm_981557 = tempm_981557[tsrmnz]
    # sal_m_981557 = sal_m_981557[tsrmnz]
    # tempm_981848 = tempm_981848[tsrmnz]
    # sal_m_981848 = sal_m_981848[tsrmnz]
    # tempm_982242 = tempm_982242[tsrmnz]
    # sal_m_982242 = sal_m_982242[tsrmnz]
    # tempm_982365 = tempm_982365[tsrmnz]
    # sal_m_982365 = sal_m_982365[tsrmnz]
    # tempm_982490 = tempm_982490[tsrmnz]
    # sal_m_982490 = sal_m_982490[tsrmnz]
    #
    # tempm_983014 = tempm_983014[tsrmnz]
    # sal_m_983014 = sal_m_983014[tsrmnz]
    # tempm_983076 = tempm_983076[tsrmnz]
    # sal_m_983076 = sal_m_983076[tsrmnz]
    # tempm_983204 = tempm_983204[tsrmnz]
    # sal_m_983204 = sal_m_983204[tsrmnz]
    # tempm_983286 = tempm_983286[tsrmnz]
    # sal_m_983286 = sal_m_983286[tsrmnz]
    # tempm_983543 = tempm_983543[tsrmnz]
    # sal_m_983543 = sal_m_983543[tsrmnz]
    # tempm_983672 = tempm_983672[tsrmnz]
    # sal_m_983672 = sal_m_983672[tsrmnz]
    # tempm_983790 = tempm_983790[tsrmnz]
    # sal_m_983790 = sal_m_983790[tsrmnz]
    # tempm_983982 = tempm_983982[tsrmnz]
    # sal_m_983982 = sal_m_983982[tsrmnz]
    #
    # tempm_990003 = tempm_990003[tsrmnz]
    # sal_m_990003 = sal_m_990003[tsrmnz]
    # tempm_995440 = tempm_995440[tsrmnz]
    # sal_m_995440 = sal_m_995440[tsrmnz]

    # 2nd round of processing with timestamp revised
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_34f3_19f6_97e1.nc'

    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2013.nc'
    ctdfil = '/media/dskone/CELTIC/IMI_CTD_2014.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2015.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2016.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2017.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2018.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2019.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2020.nc'
    # ctdfil = '/media/dskone/CELTIC/IMI_CTD_2021.nc'

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
        # glo_latidx[td] = glor.geo_idx(lat_c[td], GLOlat)
        # timec_ymdll[td, 7] = glor.geo_idx(lat_c[td], GLOlat)
        # glo_lonidx[td] = glor.geo_idx(lon_c[td], GLOlon)
        # timec_ymdll[td, 8] = glor.geo_idx(lon_c[td], GLOlon)
        timec_ymdll[td, 9] = timec_date[td].hour
    ncctf.close()

    # Get unique dates to cycle through the equivalent model output history files
    # Get unique grid coordinates for lat/lon combinations for each unique date
    timec_dates_lls = np.unique(timec_ymdll, axis=0)
    timec_dates = np.unique(timec_ymdll[:, 0:2], axis=0)
    timec_ll = np.unique(timec_ymdll[:, 3:5], axis=0)

    mapfname = hisdir_main + 'croco_VAL' + '_CTD_positions_2019.jpg'
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

        ncibi = netcdf(ibifile, 'r')
        ncibio = ncibi.variables
        # Extract vertical section from history file
        temp_ibi = np.array(ncibio['thetao']
                            [:, :, timec_dates_lls[tdl, 5], timec_dates_lls[tdl, 6]]).flatten()
        salt_ibi = np.array(ncibio['so']
                            [:, :, timec_dates_lls[tdl, 5], timec_dates_lls[tdl, 6]]).flatten()

        # Cycle through the individual dates to get repeated section plots for visual assessment of model performance
        # for d in range(0, len(timec_dates)):
        latsect = lat_rho[timec_dates_lls[tdl, 3], 0]
        latdegs = np.floor(latsect)
        latmins = np.round((latsect - latdegs) * 60)

        lonsect = lon_rho[0, timec_dates_lls[tdl, 4]]
        londegs = np.abs(np.ceil(lonsect))
        lonmins = np.round((abs(lonsect) - londegs) * 60)

        sectftname = hisdir_main + 'croco1p2p1_TEMP_sect' + '_' + str(int(tdl + 1)) + '_' + \
                     str(timec_dates_lls[tdl, 0]) + \
                     str(timec_dates_lls[tdl, 1]).zfill(2) + \
                     str(timec_dates_lls[tdl, 2]).zfill(2) + '_' + \
                     str(int(latdegs)).zfill(2) + 'd' + str(int(latmins)).zfill(2) + 'm' + 'N' + '_' + \
                     str(int(londegs)).zfill(2) + 'd' + str(int(lonmins)).zfill(2) + 'm' + 'W' + '.jpg'
        # sectftname = hisdir_main + 'croco1p1_TEMP_sect' + '_' + str(int(tdl + 1)) + '_' + \

        sectfsname = hisdir_main + 'croco1p2p1_SALT_sect' + '_' + str(int(tdl + 1)) + '_' + \
                     str(timec_dates_lls[tdl, 0]) + \
                     str(timec_dates_lls[tdl, 1]).zfill(2) + \
                     str(timec_dates_lls[tdl, 2]).zfill(2) + '_' + \
                     str(int(latdegs)).zfill(2) + 'd' + str(int(latmins)).zfill(2) + 'm' + 'N' + '_' + \
                     str(int(londegs)).zfill(2) + 'd' + str(int(lonmins)).zfill(2) + 'm' + 'W' + '.jpg'
        # sectfsname = hisdir_main + 'croco1p1_SALT_sect' + '_' + str(int(tdl + 1)) + '_' + \

        rwidx = np.argwhere(((timec_year == timec_dates_lls[tdl, 0]) &
                             (timec_month == timec_dates_lls[tdl, 1]) &
                             (timec_day == timec_dates_lls[tdl, 2]) &
                             (ctd_latidx == timec_dates_lls[tdl, 3]) &
                             (ctd_lonidx == timec_dates_lls[tdl, 4]))).flatten()

        plt.figure(figsize=(7, 5))
        sortidx = np.argsort(deptc[rwidx])
        dmin = np.min(deptc[rwidx][sortidx])
        idmin = np.min(IBIdep)
        # dmin = np.max((dmin, idmin))
        dmin = idmin
        dmax = np.max(deptc[rwidx][sortidx])
        # ibididx = np.argwhere(IBIdep >= dmax)[0][0]-1
        # ibididx = np.argwhere(IBIdep >= dmax)[0][0]
        ibididx = np.argwhere(temp_ibi == 10)[0][0]-1
        # plt.plot(tempc[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
        #          temp_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
        #          tempm_IJS[rwidx][sortidx], deptc[rwidx][sortidx], 'b-',
        #          tempm_979519[rwidx][sortidx], deptc[rwidx][sortidx], 'b--',
        #          tempm_979685[rwidx][sortidx], deptc[rwidx][sortidx], 'b-.',
        #          tempm_979807[rwidx][sortidx], deptc[rwidx][sortidx], 'b:',
        #          tempm_980116[rwidx][sortidx], deptc[rwidx][sortidx], 'g-',
        #          tempm_980558[rwidx][sortidx], deptc[rwidx][sortidx], 'g--',
        #          tempm_980771[rwidx][sortidx], deptc[rwidx][sortidx], 'g-.',
        #          tempm_980885[rwidx][sortidx], deptc[rwidx][sortidx], 'g:')
        # # v1.1
        # plt.plot(tempc[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
        #          temp_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
        #          tempm_IJS[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # tempm_979519[rwidx][sortidx], deptc[rwidx][sortidx],
        #          tempm_979685[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # tempm_979807[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # tempm_980116[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # tempm_980558[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # tempm_980771[rwidx][sortidx], deptc[rwidx][sortidx],
        #          tempm_983014[rwidx][sortidx], deptc[rwidx][sortidx],
        #          tempm_983076[rwidx][sortidx], deptc[rwidx][sortidx],
        #          tempm_983204[rwidx][sortidx], deptc[rwidx][sortidx],
        #          tempm_983286[rwidx][sortidx], deptc[rwidx][sortidx])
        # # plt.legend(('CTD', 'IBI_raw', 'IBI_croco', '979519', '979685',
        # #             '979807', '980116', '980558', '980771', '983014',
        # #             '983076', '983204', '983286'), loc='best', shadow=True)
        # plt.legend(('CTD', 'IBI_raw', 'IBI_croco', '979685',
        #             '983014', '983076', '983204', '983286'), loc='best', shadow=True)
        # v1.2.1
        plt.plot(tempc[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
                 temp_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
                 # tempm_980885[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_981213[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_981356[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_981557[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_981848[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_982242[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_982365[rwidx][sortidx], deptc[rwidx][sortidx],
                 # tempm_982490[rwidx][sortidx], deptc[rwidx][sortidx])
                 # tempm_982490[rwidx][sortidx], deptc[rwidx][sortidx],
                 tempm_983543[rwidx][sortidx], deptc[rwidx][sortidx],
                 tempm_983672[rwidx][sortidx], deptc[rwidx][sortidx],
                 tempm_983790[rwidx][sortidx], deptc[rwidx][sortidx],
                 tempm_983982[rwidx][sortidx], deptc[rwidx][sortidx],
                 tempm_990003[rwidx][sortidx], deptc[rwidx][sortidx],
                 tempm_995440[rwidx][sortidx], deptc[rwidx][sortidx])
        # plt.legend(('CTD', 'IBI_raw', '980885', '981213', '981356',
        #             '981557', '981848', '982242', '982365', '982490',
        #             '983543', '983672', '983790', '983982'), loc='best', shadow=True)
        plt.legend(('CTD', 'IBI_raw', '983543', '983672', '983790', '983982',
                    '990003', '995440'), loc='best', shadow=True)
        plt.title(str(timec_dates_lls[tdl, 9]) + 'o clock')
        plt.ylim(dmin, dmax)
        plt.gca().invert_yaxis()
        plt.xlabel('Temp ($\degree$C)')
        plt.xlim(9.5, 16.5)
        plt.ylabel('Depth(m)')
        plt.savefig(sectftname, dpi='figure', format='jpg')

        plt.figure(figsize=(7, 5))
        # plt.plot(sal_c[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
        #          salt_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
        #          sal_m_IJS[rwidx][sortidx], deptc[rwidx][sortidx], 'b-',
        #          sal_m_979519[rwidx][sortidx], deptc[rwidx][sortidx], 'b--',
        #          sal_m_979685[rwidx][sortidx], deptc[rwidx][sortidx], 'b-.',
        #          sal_m_979807[rwidx][sortidx], deptc[rwidx][sortidx], 'b:',
        #          sal_m_980116[rwidx][sortidx], deptc[rwidx][sortidx], 'g-',
        #          sal_m_980558[rwidx][sortidx], deptc[rwidx][sortidx], 'g--',
        #          sal_m_980771[rwidx][sortidx], deptc[rwidx][sortidx], 'g-.',
        #          sal_m_980885[rwidx][sortidx], deptc[rwidx][sortidx], 'g:')
        # # v1.1
        # plt.plot(sal_c[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
        #          salt_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
        #          sal_m_IJS[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # sal_m_979519[rwidx][sortidx], deptc[rwidx][sortidx],
        #          sal_m_979685[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # sal_m_979807[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # sal_m_980116[rwidx][sortidx], deptc[rwidx][sortidx],
        #          # sal_m_980558[rwidx][sortidx], deptc[rwidx][sortidx],
        #          sal_m_983014[rwidx][sortidx], deptc[rwidx][sortidx],
        #          sal_m_983076[rwidx][sortidx], deptc[rwidx][sortidx],
        #          sal_m_983204[rwidx][sortidx], deptc[rwidx][sortidx],
        #          sal_m_983286[rwidx][sortidx], deptc[rwidx][sortidx])
        # # plt.legend(('CTD', 'IBI_raw', 'IBI_croco', '979519', '979685',
        # #             '979807', '980116', '980558', '980771', '983014',
        # #             '983076', '983204', '983286'), loc='best', shadow=True)
        # plt.legend(('CTD', 'IBI_raw', 'IBI_croco', '979685',
        #             '983014', '983076', '983204', '983286'), loc='best', shadow=True)
        # v1.2.1
        plt.plot(sal_c[rwidx][sortidx], deptc[rwidx][sortidx], 'k',
                 salt_ibi[0:ibididx], IBIdep[0:ibididx], 'b',
                 # sal_m_980885[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_981213[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_981356[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_981557[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_981848[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_982242[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_982365[rwidx][sortidx], deptc[rwidx][sortidx],
                 # sal_m_982490[rwidx][sortidx], deptc[rwidx][sortidx])
                 # sal_m_982490[rwidx][sortidx], deptc[rwidx][sortidx],
                 sal_m_983543[rwidx][sortidx], deptc[rwidx][sortidx],
                 sal_m_983672[rwidx][sortidx], deptc[rwidx][sortidx],
                 sal_m_983790[rwidx][sortidx], deptc[rwidx][sortidx],
                 sal_m_983982[rwidx][sortidx], deptc[rwidx][sortidx],
                 sal_m_990003[rwidx][sortidx], deptc[rwidx][sortidx],
                 sal_m_995440[rwidx][sortidx], deptc[rwidx][sortidx])
        # plt.legend(('CTD', 'IBI_raw', '980885', '981213', '981356',
        #             '981557', '981848', '982242', '982365', '982490',
        #             '983543', '983672', '983790', '983982'), loc='best', shadow=True)
        plt.legend(('CTD', 'IBI_raw', '983543', '983672', '983790', '983982',
                    '990003', '995440'), loc='best', shadow=True)
        plt.title(str(timec_dates_lls[tdl, 9]) + 'o clock')
        plt.ylim(dmin, dmax)
        plt.gca().invert_yaxis()
        plt.xlabel('Salinity (psu)')
        plt.xlim(34.4, 35.8)
        plt.ylabel('Depth(m)')

        plt.savefig(sectfsname, dpi='figure', format='jpg')
        ncibi.close()
        # nccrocostn_ibi.close()
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
                # hisdir_stats = hisdir_IBI_JS

                # hisdir_stats = [hisdir_979519, hisdir_979685, hisdir_979807, hisdir_980116, hisdir_980558,
                #                 hisdir_980771, hisdir_983014, hisdir_983076, hisdir_983204, hisdir_983286,
                #                 hisdir_980885, hisdir_981213, hisdir_981356, hisdir_981557, hisdir_981848,
                #                 hisdir_982242, hisdir_982365, hisdir_982490, hisdir_983543, hisdir_983672,
                #                 hisdir_983790, hisdir_983982, hisdir_985519, hisdir_985849, hisdir_985962]
                # hisdir_stats = [hisdir_985519, hisdir_985962]
                # hisdir_stats = [hisdir_987878, hisdir_989174, hisdir_989667, hisdir_990003, hisdir_995440]
                # hisdir_stats = [hisdir_981356, hisdir_982365, hisdir_983543, hisdir_983672,
                #                 hisdir_983790, hisdir_983982, hisdir_990003, hisdir_995440]
                # hisdir_stats = [hisdir_985519, hisdir_985962]
                # hisdir_stats = [hisdir_1005613, hisdir_1006516, hisdir_1006690,
                #                 hisdir_1007268, hisdir_1008287, hisdir_1009495]
                # hisdir_stats = [hisdir_1005613, hisdir_1006516, hisdir_1006690,
                #                 hisdir_1008287, hisdir_1009495]
                # hisdir_stats = [hisdir_1015382, hisdir_1015647]

                hisdir_stats = [hisdir_1016234, hisdir_1016529, hisdir_1016656]
                for hs in range(0, len(hisdir_stats)):
                    # Labelling files by the fact model outputs now at midnight, not midday
                    # Therefore, stats comparing like with like
                    VALmname = hisdir_stats[hs] + \
                               'croco_mn_VALm' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                    VALdname = hisdir_stats[hs] + \
                               'croco_mn_VALd' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
                    VALsmename = hisdir_stats[hs] + \
                                 'croco_mn_VAL_' + hisdir_stats[hs][-7:-1] + \
                                 '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + 's_ME.jpg'

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

                    hisfile = hisdir_stats[hs] + his_prefix + 'Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
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
                            ((temps_m_slice[slice_idx] - temps_s_slice[slice_idx]).mean()
                             / temps_s_slice[slice_idx].mean())
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
                    # ax = plt.axes(projection=ccrs.Mercator(central_longitude=-8.29,
                    #               min_latitude=49, max_latitude=52.95))
                    fig = plt.figure(figsize=(7, 5))
                    ax = plt.contourf(lon_rho, lat_rho, spme2, [-5, -2.5, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2.5, 5])
                    plt.colorbar()
                    # ax.set_global()
                    # ax.coastlines()

                    # coast_10m = cfeature.NaturalEarthFeature("physical", "land", "10m",
                    #                                           edgecolor="k", facecolor="0.8")
                    # ax.add_feature(coast_10m)
                    plt.savefig(VALsmename, dpi='figure', format='jpg')
                    # plt.savefig(VALsmename, dpi='figure', format='jpg', metadata=None,
                    #             bbox_inches=None, pad_inches=0.1,
                    #             facecolor='auto', edgecolor='auto',
                    #             backend=None)
                    ncvd.close()
