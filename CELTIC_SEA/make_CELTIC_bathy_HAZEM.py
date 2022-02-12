# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Yellow_Sea/make_YELLOW_grd_v1.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021
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
filename = 'EMODnet_Bathy_100m.tif'
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
# posY = yoffset + (rot2 * x) + (px_h * y)
# Swapped px_h and rot2 between line above and below - gdal GeoTransForm apparently outputting them in reverse
posY = yoffset + (px_h * x) + (rot2 * y)
# shift to the center of the pixel
posX += px_w / 2.0
posY += px_h / 2.0

long = posX
lat = posY
lat = np.flipud(lat)
data_array = np.flipud(data_array)

# # Grid dimension
# # Lm = 392  # Longitude
# # Mm = 350  # Latitude
# Lm = 359  # Longitude
# Mm = 439  # Latitude
#
# lon0 = -10.75
# # lat0 = 52.15
# lat0 = 52.95
#
# lon1 = -10.75
# lat1 = 49.
#
# # lon2 = -5.
# lon2 = -5.83
# lat2 = 49.
#
# # lon3 = -5.
# # lat3 = 52.15
# lon3 = -5.83
# lat3 = 52.95
#
# # map = Basemap(projection='lcc', lat_0=50.25, lat_1=49., lat_2=51, lon_0=-8.5,
# #               width=2000000, height=2000000, resolution='i')
#
# # map = Basemap(projection='lcc', lat_0=50.575, lat_1=50., lat_2=51., lon_0=-7.875,
# #               width=392000, height=350000, resolution='i')
#
# # map = Basemap(projection='merc', llcrnrlat=49, urcrnrlat=52.15, llcrnrlon=-10.75, urcrnrlon=-5, resolution='f')
# map = Basemap(projection='merc', llcrnrlat=49, urcrnrlat=52.95, llcrnrlon=-10.75, urcrnrlon=-5.83, resolution='f')
# map.drawcoastlines()
# map.fillcontinents(color='coral', lake_color='aqua')
# # draw parallels and meridians.
# map.drawparallels(np.arange(49., 53., 1.))
# map.drawmeridians(np.arange(-15., -5., 2.5))
# map.drawmapboundary(fill_color='aqua')
# plt.title("BIOCELTIC DOMAIN")
# plt.show()
#
# lonp = np.array([lon0, lon1, lon2, lon3])
# latp = np.array([lat0, lat1, lat2, lat3])
# beta = np.array([1, 1, 1, 1])

# Generate the new grid
# Do this if you aren't going to move the grid corners interactively.
# hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
# AttributeError: module 'pyroms' has no attribute 'grid'
#
# hgrd = Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
#
# # Do this if you are going to use the Boundary Interactor
# # map.drawcoastlines()
# # xp, yp = map(lonp, latp)
# # bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mm+3,Lm+3), proj=map)
# # hgrd=bry.grd
#
# lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
# hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)
#
# # # Generate the mask
# # for verts in map.coastsegs:
# #     hgrd.mask_polygon(verts)
# # # alternate version from johan.navarro.padron
#
# for xx, yy in map.coastpolygons:
#     xa = np.array(xx, np.float32)
#     ya = np.array(yy, np.float32)
#     vv = np.zeros((xa.shape[0], 2))
#     vv[:, 0] = xa
#     vv[:, 1] = ya
#     hgrd.mask_polygon(vv, mask_value=False)  # hgrd.mask_polygon(vv, mask_value=0)
#
# # Edit the land mask interactively.
# pyroms.grid.edit_mask_mesh(hgrd, proj=map)
# # edit_mask_mesh_ij is a faster version using imshow  but no map projection.
# coast = pyroms.utility.get_coast_from_map(map)
# pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

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
h = -data_array

# # fix minimum depth
# hmin = 20
# topo = np.where(topo < hmin, hmin, topo)
#
# # interpolate new bathymetry (when gridded/regular)
# bath_int_fn = RegularGridInterpolator((lat[:, 0], long[0, :]), topo)
#
# h = bath_int_fn((hgrd.lat_rho.flat, hgrd.lon_rho.flat))
# h = np.reshape(h, hgrd.lon_rho.shape)

# # interpolate new bathymetry (when ungridded/irregular)
# # lon, lat = np.meshgrid(long, lat) # use sparse matrix to limit computational overheads...
# lon, lat = np.meshgrid(long, lat, sparse=True)
# h = griddata((long.flat, lat.flat), topo.flat, (hgrd.lon_rho, hgrd.lat_rho), method='linear')

# insure that depth is always deeper than hmin
# h = np.where(h < hmin, hmin, h)

# # set depth to hmin where masked
# # idx = np.where(hgrd.mask_rho == 0)
# idx = np.where(np.logical_or(np.equal(hgrd.mask_rho, 0), np.isnan(h)))
# h[idx] = hmin

# shvI = np.where(np.greater_equal(h, 190))
# h[shvI] = h[shvI[0]-1, shvI[1]] + h[shvI[0]+1, shvI[1]] + h[shvI[0], shvI[1]-1] + h[shvI[0], shvI[1]+1]
# h[shvI] = h[shvI]/4
# shvII = np.where(np.greater_equal(h, 190))
# h[shvII] = h[shvII[0]-1, shvII[1]] + h[shvII[0]+1, shvII[1]] + h[shvII[0], shvII[1]-1] + h[shvII[0], shvII[1]+1]
# h[shvII] = h[shvII]/4
# shvIII = np.where(np.greater_equal(h, 185))
# h[shvIII] = h[shvIII[0]-1, shvIII[1]] + h[shvIII[0]+1, shvIII[1]] + h[shvIII[0], shvIII[1]-1] + h[shvIII[0],
# shvIII[1]+1]
# h[shvIII] = h[shvIII]/4

# save raw bathymetry
hraw = h.copy()

# # check bathymetry roughness
# RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
# print('Max Roughness value is: ', RoughMat.max())
#
# # smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
# rx0_max = 0.2
# # h = bathy_smoothing.smoothing_PlusMinus_rx0(hgrd.mask_rho, h, rx0_max, GridAreas)
# h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)
#
# # check bathymetry roughness again
# RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
# print('Max Roughness value is: ', RoughMat.max())

# vertical coordinate
# theta_b = 2
theta_b = 0.0
theta_s = 7.0
Tcline = 20
N = 20

# This fn gives in ROMS nomenclature:
# vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=4
vgrd = pyroms.vgrid.s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=2
# vgrd = pyroms.vgrid.s_coordinate_5(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=5
# ROMS grid
# grd_name = 'CELTIC_II'

grd_name = 'HAZEM_BATHY'
grd = pyroms.grid.ROMS_Grid(grd_name, h, vgrd)

# write grid to netcdf file
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v1.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v2.nc')
# pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v3.nc')
pyroms.grid.write_ROMS_grid(grd, filename='HAZEM_BATHY.nc')

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
