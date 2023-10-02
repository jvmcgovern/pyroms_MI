# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Yellow_Sea/make_YELLOW_grd_v1.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2022
# Averaging of bathymetry instead of interpolation of high resolution data

import pandas as pd
from bathy_smoother import bathy_tools, bathy_smoothing
import numpy as np
from mpl_toolkits.basemap import Basemap  # need to execute: python -m pip install -U matplotlib==3.2
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import pyroms
from pyroms.grid import Gridgen
import matplotlib.pyplot as plt
import numpy.matlib
from osgeo import gdal
import time
import pickle
from netCDF4 import Dataset as netcdf

# import osr
# import os
# import matplotlib
# # from bathy_smoother import *
# from mpl_toolkits.basemap import shiftgrid  # need to execute: python -m pip install -U matplotlib==3.2
# import matplotlib.colors as colors
# from scipy.signal import medfilt2d
# import netCDF4
# filename = r'C:\Users\jmcgovern\PycharmProjects\pyroms_MI\CELTICSEA\Bathymetry_CelticSea_50m.tif'
# filename = r'C:\Users\jmcgovern\PycharmProjects\pyroms_MI\CELTICSEA\CelticSeaBathymetry_100m.tif'
# filename = r'C:\Users\jmcgovern\PycharmProjects\pyroms_MI\CELTICSEA\EMODnet_Bathy_100m.tif'
# filename = 'EMODnet_Bathy_100m.tif'
lat2mslname = '/media/dskone/CELTIC/BATHYMETRY/Final_MSL_LAT_separation_grid.xyz'
l2mtable = pd.read_csv(lat2mslname, names=['lmlon', 'lmlat', 'lat2msl'], header=None)
lmlon = l2mtable[['lmlon']].to_numpy().flatten()
lmlat = l2mtable[['lmlat']].to_numpy().flatten()
lat2msl = l2mtable[['lat2msl']].to_numpy().flatten()

# ulmlon = np.unique(lmlon)
# ulmlat = np.unique(lmlat)
# glat2msl = np.mean(lat2msl)*np.ones((ulmlon.shape[0], ulmlat.shape[0]))

filename = '/media/dskone/CELTIC/BATHYMETRY/EMODnet_Bathy_100m.tif'
bath = gdal.Open(filename)
bath_info = gdal.Info(filename)
band1 = bath.GetRasterBand(1)
nodataval = band1.GetNoDataValue()

# convert to a numpy array
data_array = bath.ReadAsArray().astype(np.float)

# replace missing values if necessary
if np.any(data_array == nodataval):
    data_array[data_array == nodataval] = np.nan

xoffset, px_w, rot1, yoffset, px_h, rot2 = bath.GetGeoTransform()

# supposing x and y are your pixel coordinate this
# is how to get the coordinate in space.


# x = 0  # all the way to 10560
# y = 0     # up to 5280

x = np.arange(data_array.shape[1])
x = np.matlib.repmat(x, data_array.shape[0], 1)

y = np.arange(data_array.shape[0])
y = np.matlib.repmat(y, data_array.shape[1], 1)
y = np.transpose(y)

posX = xoffset + (px_w * x) + (rot1 * y)
# Swapped px_h and rot2 between line above and below - gdal GeoTransForm apparently outputting them in reverse
posY = yoffset + (px_h * x) + (rot2 * y)
# shift to the center of the pixel
posX += px_w / 2.0
posY += px_h / 2.0

long = posX
lat = posY
lat = np.flipud(lat)
data_array = np.flipud(data_array)

# Grid dimension
Lm = 359  # Longitude
Mm = 439  # Latitude

lon0 = -10.75
lat0 = 52.95

lon1 = -10.75
lat1 = 49.

lon2 = -5.83
lat2 = 49.

lon3 = -5.83
lat3 = 52.95

map = Basemap(projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83, resolution='f')
map.drawcoastlines()
map.fillcontinents(color='coral', lake_color='aqua')
# draw parallels and meridians.
map.drawparallels(np.arange(49., 53., 1.))
map.drawmeridians(np.arange(-15., -5., 2.5))
map.drawmapboundary(fill_color='aqua')
plt.title("BIOCELTIC DOMAIN")
plt.show()

lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])
beta = np.array([1, 1, 1, 1])

# Generate the new grid
# Do this if you aren't going to move the grid corners interactively.
# hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
# AttributeError: module 'pyroms' has no attribute 'grid'

hgrd = Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)

# Do this if you are going to use the Boundary Interactor
# map.drawcoastlines()
# xp, yp = map(lonp, latp)
# bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mm+3,Lm+3), proj=map)
# hgrd=bry.grd

lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

lat_rho = hgrd.lat_rho
lon_rho = hgrd.lon_rho

latr = max(abs(lat_rho[1, 1] - lat_rho[2, 1]), abs(lat_rho[1, 1] - lat_rho[1, 2]))
lonr = max(abs(lon_rho[1, 1] - lon_rho[2, 1]), abs(lon_rho[1, 1] - lon_rho[1, 2]))

ltdelta = 0.5 * latr
ltmin = lat_rho - ltdelta
ltmax = lat_rho + ltdelta

lndelta = 0.5 * lonr
lnmin = lon_rho - lndelta
lnmax = lon_rho + lndelta

crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
grdname = crocofiles_dir + 'croco_grd.nc'  # was hmin = 20m
grd_pkl = crocofiles_dir + 'croco_crs_fin.pkl'

ncg = netcdf(grdname, 'r')
ncgrd = ncg.variables
mask_rho = np.array(ncgrd['mask_rho'][:])
ncg.close()

pickled = 1

if pickled == 0:
    crs_fin_grd = np.empty((lat_rho.shape[0], lon_rho.shape[1]), dtype=object)
    for lts in range(0, lat_rho.shape[0]):
        sta = time.time()
        for lns in range(0, lon_rho.shape[1]):
            if mask_rho[lts, lns] == 1.:
                crs_fin_grd[lts, lns] = np.argwhere(((lat >= ltmin[lts, lns]) &
                                                     (lat <= ltmax[lts, lns])) &
                                                    ((long >= lnmin[lts, lns]) &
                                                     (long <= lnmax[lts, lns])))
        dne = time.time()
        print(dne-sta)
    with open(grd_pkl, 'wb') as f:
        pickle.dump(crs_fin_grd, f)
else:
    with open(grd_pkl, 'rb') as f:
        crs_fin_grd = pickle.load(f)

# # Generate the mask
# for verts in map.coastsegs:
#     hgrd.mask_polygon(verts)
# # alternate version from johan.navarro.padron

for xx, yy in map.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy, np.float32)
    vv = np.zeros((xa.shape[0], 2))
    vv[:, 0] = xa
    vv[:, 1] = ya
    hgrd.mask_polygon(vv, mask_value=False)  # hgrd.mask_polygon(vv, mask_value=0)

# Edit the land mask interactively.
pyroms.grid.edit_mask_mesh(hgrd, proj=map)
# edit_mask_mesh_ij is a faster version using imshow  but no map projection.
coast = pyroms.utility.get_coast_from_map(map)
pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

# # Use the following to interpolate from etopo2 bathymetry.
# # generate the bathy
# # read in topo data (on a regular lat/lon grid)
# # this topo come with basemap so you should have it on your laptop.
# # just update datadir with the appropriate path
# # you can get this data from matplolib svn with
# # svn co https://matplotlib.svn.sourceforge.net/svnroot/matplotlib/trunk/htdocs/screenshots/data/"
#
# datadir = 'data/'
# topo = np.loadtxt(os.path.join(datadir, 'etopo20data.gz'))
# lons = np.loadtxt(os.path.join(datadir, 'etopo20lons.gz'))
# lats = np.loadtxt(os.path.join(datadir, 'etopo20lats.gz'))

# depth positive
topo = -data_array

# fix minimum depth
# hmin = 20
# hmin = 15
# hmin = 12
# hmin = 10
hmin = 10
topo = np.where(topo < hmin, hmin, topo)

# # interpolate new bathymetry (when gridded/regular)
# bath_int_fn = RegularGridInterpolator((lat[:, 0], long[0, :]), topo)
#
# h = bath_int_fn((hgrd.lat_rho.flat, hgrd.lon_rho.flat))
# h = np.reshape(h, hgrd.lon_rho.shape)

if pickled == 0:
    h = hmin * np.ones_like(mask_rho)
    for lts in range(0, lat_rho.shape[0]):
        # sta = time.time()
        for lns in range(0, lon_rho.shape[1]):
            if mask_rho[lts, lns] == 1.:
                h[lts, lns] = np.mean(topo[crs_fin_grd[lts, lns][:, 0], crs_fin_grd[lts, lns][:, 1]])
            else:
                h[lts, lns] = hmin/2.
        # dne = time.time()
        # print(dne-sta)
    with open(grd_pkl, 'wb') as f:
        pickle.dump(h, f)
else:
    with open(grd_pkl, 'rb') as f:
        h = pickle.load(f)

# interpolate irregularly gridded LAT to MSL data to model grid
lat2mslg = griddata((lmlat, lmlon), lat2msl, (hgrd.lat_rho.flat, hgrd.lon_rho.flat))
lat2mslg = np.reshape(lat2mslg, hgrd.lon_rho.shape)
lat2mslg = -lat2mslg

h = h + lat2mslg

# # interpolate new bathymetry (when ungridded/irregular)
# # lon, lat = np.meshgrid(long, lat) # use sparse matrix to limit computational overheads...
# lon, lat = np.meshgrid(long, lat, sparse=True)
# h = griddata((long.flat, lat.flat), topo.flat, (hgrd.lon_rho, hgrd.lat_rho), method='linear')

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# # set depth to hmin where masked
# # idx = np.where(hgrd.mask_rho == 0)
# idx = np.where(np.logical_or(np.equal(hgrd.mask_rho, 0), np.isnan(h)))
idx = np.where(np.isnan(h))
h[idx] = hmin

# shvI = np.where(np.greater_equal(h, 192))
# h[shvI] = h[shvI[0]-1, shvI[1]] + h[shvI[0]+1, shvI[1]] + h[shvI[0], shvI[1]-1] + h[shvI[0], shvI[1]+1]
# h[shvI] = h[shvI]/4
# shvII = np.where(np.greater_equal(h, 192))
# h[shvII] = h[shvII[0]-1, shvII[1]] + h[shvII[0]+1, shvII[1]] + h[shvII[0], shvII[1]-1] + h[shvII[0], shvII[1]+1]
# h[shvII] = h[shvII]/4
# shvIII = np.where(np.greater_equal(h, 187))
# h[shvIII] = h[shvIII[0]-1, shvIII[1]] + h[shvIII[0]+1, shvIII[1]] + h[shvIII[0], shvIII[1]-1] + h[shvIII[0], shvIII[1]+1]
# h[shvIII] = h[shvIII]/4

# save raw bathymetry
hraw = h.copy()

# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.2
# h = bathy_smoothing.smoothing_PlusMinus_rx0(hgrd.mask_rho, h, rx0_max, GridAreas)
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# vertical coordinate
theta_b = 0.0
theta_s = 7.0
# Tcline = hmin
Tcline = 50
N = 20

# This fn gives in ROMS nomenclature:
# vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=4
vgrd = pyroms.vgrid.s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=2
# vgrd = pyroms.vgrid.s_coordinate_5(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=5
# ROMS grid
# grd_name = 'CELTIC_II'

# grd_name = 'CELTIC_V'  # hmin = 20m
# grd_name = 'CELTIC_V_MSL_h8'  # hmin = 8m
grd_name = 'CELTIC_V_MSL_h10'  # hmin = 8m
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v1.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v2.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v3.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v5.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v5_h8.nc')
pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v5_MSL_h10.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v6_MSL_h8.nc')

# If there are two different coordinate systems, transformation needed
# **TransformPoint can only be used on single points, not arrays, loop needs to be introduced over array**
# # get CRS from dataset
# crs = osr.SpatialReference()
# crs.ImportFromWkt(bath.GetProjectionRef())
# # create lat/long crs with WGS84 datum
# crsGeo = osr.SpatialReference()
# crsGeo.ImportFromEPSG(4326)  # 4326 is the EPSG id of WGS84 lat/long crs
# t = osr.CoordinateTransformation(crs, crsGeo)
# (long, lat, z) = t.TransformPoint(posX, posY)
# print(long)
# print(lat)

# Too many points for mapping purposes, need to process to ROMS grid before plotting
# # Plot out data with Matplotlib's 'contour'
# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111)
# plt.contour(data_array, cmap="viridis",
#             levels=list(range(0, 5000, 100)))
# plt.title("Bathymetry data for Celtic Sea")
# cbar = plt.colorbar()
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()
