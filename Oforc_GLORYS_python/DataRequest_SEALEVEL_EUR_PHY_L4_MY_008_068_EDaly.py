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
import pandas

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

    GLO_prefix = 'dt_global_allsat_phy_l4_'
    GLO_post = '_*.nc'

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

    his_prefix = 'croco_his_'

    sat_sstfiles_dir = '/media/dskone/VAL/TEMP/SST_L4/'
    # sat_sst_ending_core_ex = YYYYMMDD
    sat_sst_ending = '000000-IFR-L4_GHRSST-SSTfnd-ODYSSEA-ATL_005-v2.0-fv1.0.nc'

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation
    comp_delaunay = 1  # 1: compute delaunay triangulations - 0: use saved matrices (for debugging)

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    chfile = '/media/dskone/EDaly_CMEMS_SEALEVEL_EUR_PHY_L4_MY_008_068/SurfVelsDataPoints2.xls'

    df = pandas.read_excel(chfile)

    ugos_out = np.zeros((df.shape[0], 1))

    vgos_out = np.zeros((df.shape[0], 1))

    glototal = sorted(glob.glob('/media/dskone/EDaly_CMEMS_SEALEVEL_EUR_PHY_L4_MY_008_068/' + GLO_prefix +
                                '????????' + GLO_post))
    glogen = glototal[0]
    ncg = netcdf(glogen, 'r')
    ncgo = ncg.variables
    latT = np.array(ncgo['latitude'][:])
    lonT = np.array(ncgo['longitude'][:])
    ncg.close()

    for rw in range(0, df.shape[0]):
        glofiles = sorted(glob.glob('/media/dskone/EDaly_CMEMS_SEALEVEL_EUR_PHY_L4_MY_008_068/' + GLO_prefix +
                                    str(df['Year'][rw]) + str(df['Month'][rw]).zfill(2) +
                                    str(df['Day'][rw]).zfill(2) + GLO_post))
        glofile = glofiles[0]
        ncf = netcdf(glofile, 'r')
        ncfo = ncf.variables

        #
        # get GLORYS positions and indices at T-points
        #

        lonidx = glor.geo_idx(df['Lon'][rw], lonT)
        latidx = glor.geo_idx(df['Lat'][rw], latT)

        ugos = ncfo['ugos'][0, latidx, lonidx]
        vgos = ncfo['vgos'][0, latidx, lonidx]

        ugos_out[rw] = ugos
        vgos_out[rw] = vgos

    df['U_gos'] = ugos_out[:, 0]
    df['V_gos'] = vgos_out[:, 0]

    df.to_excel(chfile[:-4] + '_out.xls')
