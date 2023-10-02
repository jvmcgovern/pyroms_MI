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
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
from scipy.interpolate import griddata
from netCDF4 import date2index as d2i
from datetime import date, datetime
from calendar import monthrange
import sys

# sys.path.insert(0,'')

from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from progressbar import *

if 1 == 1:

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

    Yend = 2019
    Mend = 12
    Dend = 31

    # crocofiles_dir = my_home_dir + 'SWAG/Run_TEST/CROCO_FILES/'
    crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
    grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m
    # grdname = '/media/dskone/CELTIC/croco_grd_h6.nc'  # was hmin = 6m
    # grdname = '/media/dskone/CELTIC/croco_grd_h8.nc'  # was hmin = 8m

    hisdir1 = '/media/dskone/CELTIC/CROCO_PHY_1p1_1719/'
    hisdir2 = '/media/dskone/CELTIC/CROCO_PHY_1p2p1_1719/'
    his_prefix = 'croco_his_'

    sat_sstfiles_dir = '/media/dskone/VAL/TEMP/SST_L4/'
    # sat_sst_ending_core_ex = YYYYMMDD
    sat_sst_ending = '000000-IFR-L4_GHRSST-SSTfnd-ODYSSEA-ATL_005-v2.0-fv1.0.nc'

    Nzgoodmin = 4  # number minimum of good data points to use for the interpolation
    comp_delaunay = 1  # 1: compute delaunay triangulations - 0: use saved matrices (for debugging)

    tndx_glo = 0  # time index in the GLORYS file (should be only 0 here)

    for iyear in range(Ystart, Yend + 1):
        for imonth in range(Mstart, Mend + 1):

            VALmname = hisdir1 + 'croco_VALm' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'
            VALdname = hisdir1 + 'croco_VALd' + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'

            # open dataset
            ds = xr.open_dataset(VALmname)

            # Turn on chunking to activate dask and parallelize read/write.
            ds = ds.chunk({'VAL_time': 1})

            # Pick out some of the variables that will be included as coordinates
            ds = ds.set_coords(['Cs_r', 'Cs_w', 'hc', 'Vtransform'])

            # # Select a a subset of variables. Salt will be visualized, zeta is used to
            # # calculate the vertical coordinate
            # variables = ['salt', 'zeta']
            # ds[variables].isel(ocean_time=slice(47, None, 7 * 24),
            #                    xi_rho=slice(300, None)).to_netcdf('ROMS_example.nc', mode='w')

            # # load in the file
            # ds = xr.tutorial.open_dataset("ROMS_example.nc", chunks={"ocean_time": 1})

            # if ds.Vtransform == 1:
            #     Zo_rho = ds.hc * (ds.s_rho - ds.Cs_r) + ds.Cs_r * ds.h
            #     z_rho = Zo_rho + ds.zeta * (1 + Zo_rho / ds.h)
            # elif ds.Vtransform == 2:
            #     Zo_rho = (ds.hc * ds.s_rho + ds.Cs_r * ds.h) / (ds.hc + ds.h)
            #     z_rho = ds.zeta + (ds.zeta + ds.h) * Zo_rho
            #
            # ds.coords["z_rho"] = z_rho.transpose()  # needing transpose seems to be an xarray bug
            # ds.salt
            #
            # ds.saltb.isel(xi_rho=50, ocean_time=0).plot()
            # section = ds.salt.isel(xi_rho=50, eta_rho=slice(0, 167), ocean_time=0)
            # section.plot(x="lon_rho", y="z_rho", figsize=(15, 6), clim=(25, 35))
            # plt.ylim([-100, 1])
            #
            # ds.salt.isel(s_rho=-1, ocean_time=0).plot(x="lon_rho", y="lat_rho")
            # proj = ccrs.LambertConformal(central_longitude=-92, central_latitude=29)

            proj = ccrs.Mercator(central_longitude=-8.29, min_latitude=49, max_latitude=52.95)
            fig = plt.figure(figsize=(15, 5))
            ax = plt.axes(projection=proj)
            ds.temps_s_ME.isel().plot(x="xi_rho", y="eta_rho", transform=ccrs.PlateCarree())

            coast_10m = cfeature.NaturalEarthFeature("physical", "land", "10m", edgecolor="k", facecolor="0.8")
            ax.add_feature(coast_10m)
