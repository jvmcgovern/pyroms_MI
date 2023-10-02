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

    GLO_prefix = 'mercatorglorys12v1_gl12_mean_'
    GLO_post = '_R*.nc'

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

    chfile = '/media/dskone/CHayes_CMEMS_GLO/latlongtimedata.xls'

    df = pandas.read_excel(chfile)

    MaxDepthGLO = np.zeros((df.shape[0], 1))

    Tsurf = np.zeros((df.shape[0], 1))
    T_20m = np.zeros((df.shape[0], 1))
    T_50m = np.zeros((df.shape[0], 1))
    Tbott = np.zeros((df.shape[0], 1))

    Ssurf = np.zeros((df.shape[0], 1))
    S_20m = np.zeros((df.shape[0], 1))
    S_50m = np.zeros((df.shape[0], 1))
    Sbott = np.zeros((df.shape[0], 1))

    glototal = sorted(glob.glob('/media/dskone/CHayes_CMEMS_GLO/' + GLO_prefix +
                                '????????' + GLO_post))
    glogen = glototal[0]
    ncg = netcdf(glogen, 'r')
    ncgo = ncg.variables
    depth = np.array(ncgo['depth'][:])
    latT = np.array(ncgo['latitude'][:])
    lonT = np.array(ncgo['longitude'][:])
    ncg.close()

    for rw in range(0, df.shape[0]-1):
        glofiles = sorted(glob.glob('/media/dskone/CHayes_CMEMS_GLO/' + GLO_prefix +
                                    str(df['Year'][rw]) + str(df['Month'][rw]).zfill(2) +
                                    str(df['Day'][rw]).zfill(2) + GLO_post))
        glofile = glofiles[0]
        ncf = netcdf(glofile, 'r')
        ncfo = ncf.variables

        #
        # get GLORYS positions and indices at T-points
        #

        lonidx = glor.geo_idx(df['StartLonDecStart'][rw], lonT)
        latidx = glor.geo_idx(df['StartLatDecStart'][rw], latT)

        salt = np.array(ncfo['so'][0, :, latidx, lonidx])
        temp = np.array(ncfo['thetao'][0, :, latidx, lonidx])

        maxdepthidx = np.argwhere(salt == min(salt))[0]-1
        Sbott[rw] = salt[maxdepthidx]
        Tbott[rw] = temp[maxdepthidx]
        MaxDepthGLO[rw] = depth[maxdepthidx]

        Ssurf[rw] = salt[0]
        Tsurf[rw] = temp[0]

        try:
            saltinterp = interpolate.interp1d(depth[0:maxdepthidx[0]+1], salt[0:maxdepthidx[0]+1])
            tempinterp = interpolate.interp1d(depth[0:maxdepthidx[0]+1], temp[0:maxdepthidx[0]+1])
            if MaxDepthGLO[rw] > 20:
                S_20m[rw] = saltinterp(20)
                T_20m[rw] = tempinterp(20)
            else:
                S_20m[rw] = Sbott[rw]
                T_20m[rw] = Tbott[rw]
            if MaxDepthGLO[rw] > 50:
                S_50m[rw] = saltinterp(50)
                T_50m[rw] = tempinterp(50)
            else:
                S_50m[rw] = Sbott[rw]
                T_50m[rw] = Tbott[rw]

        except:
            print('issues with ', str(rw))
        ncf.close()

    df['Sbott'] = Sbott[:, 0]
    df['Ssurf'] = Ssurf[:, 0]
    df['S_20m'] = S_20m[:, 0]
    df['S_50m'] = S_50m[:, 0]

    df['Tbott'] = Tbott[:, 0]
    df['Tsurf'] = Tsurf[:, 0]
    df['T_20m'] = T_20m[:, 0]
    df['T_50m'] = T_50m[:, 0]

    df['MaxDepthGLO'] = MaxDepthGLO[:, 0]

    df.to_excel(chfile[:-4] + '_out.xls')