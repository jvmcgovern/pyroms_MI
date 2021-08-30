# To get lpsolve working, Try:
# 1: (didn't work)
#  conda install -c snorfalorpagus lpsolve
# 2:
# https://stackoverflow.com/questions/14906603/how-to-install-lpsolve-module-for-python-on-linux-ubuntu-10-04
# https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_Python_source.tar.gz/download
# https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_win64.zip/download
# 3:
# with intel mkl 2019 installed
# conda install -c conda-forge lpsolve55

# NB use the following: to ensure .pyd modules are recognised by pycharm
# https://intellij-support.jetbrains.com/hc/en-us/community/posts/360008700780-Pycharm-does-not-recognize-pyd-file-as-a-module
# https://intellij-support.jetbrains.com/hc/en-us/community/posts/360009758419-failed-to-load-pyd-package
# 'C:\\Program Files\\JetBrains\\PyCharm Community Edition 2020.2.3\\plugins\\python-ce\\helpers\\pydev'
# 'C:\\Program Files\\JetBrains\\PyCharm Community Edition 2020.2.3\\plugins\\python-ce\\helpers\\third_party\\thriftpy'
# 'C:\\Program Files\\JetBrains\\PyCharm Community Edition 2020.2.3\\plugins\\python-ce\\helpers\\pydev'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\python36.zip'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\DLLs'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\lib'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\lib\\site-packages'
# 'c:\\users\\jmcgovern\\pycharmprojects\\pyroms_mi\\pyroms'
# 'c:\\users\\jmcgovern\\pycharmprojects\\pyroms_mi\\pyroms_toolbox'
# 'c:\\users\\jmcgovern\\pycharmprojects\\pyroms_mi\\bathy_smoother'
# 'C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI'
# 'C:/Users/jmcgovern/PycharmProjects/pyroms_MI'
# import sys
#
# sys.path.append('C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI\\pyroms')
# sys.path.append('C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI\\pyroms_toolbox')
# sys.path.append('C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI\\bathy_smoother')

import osr
import os
from bathy_smoother import bathy_tools, bathy_smoothing
import numpy as np
from mpl_toolkits.basemap import Basemap  # need to execute: python -m pip install -U matplotlib==3.2
from scipy.interpolate import griddata
import pyroms

import numpy.matlib
import gdal

import matplotlib
import matplotlib.pyplot as plt
# from bathy_smoother import *
from mpl_toolkits.basemap import shiftgrid  # need to execute: python -m pip install -U matplotlib==3.2
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4

# filename = r'C:\Users\jmcgovern\PycharmProjects\pyroms_MI\CELTICSEA\Bathymetry_CelticSea_50m.tif'
# filename = r'C:\Users\jmcgovern\PycharmProjects\pyroms_MI\CELTICSEA\CelticSeaBathymetry_100m.tif'
filename = r'C:\Users\jmcgovern\PycharmProjects\pyroms_MI\CELTICSEA\EMODnet_Bathy_100m.tif'
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

posX = (px_w * x) + (rot1 * y) + xoffset
posY = (rot2 * x) + (px_h * y) + yoffset

# shift to the center of the pixel
posX += px_w / 2.0
posY += px_h / 2.0

long = posX
lat = posY

# Grid dimension
Lm = 473  # Longitude
Mm = 500  # Latitude

lon0 = -12.
lat0 = 52.5

lon1 = -12.
lat1 = 48.

lon2 = -5.
lat2 = 48.

lon3 = -5.
lat3 = 52.5

# map = Basemap(projection='lcc', lat_0=50.25, lat_1=49., lat_2=51, lon_0=-8.5,
#               width=2000000, height=2000000, resolution='i')

map = Basemap(projection='lcc', lat_0=50.25, lat_1=49., lat_2=51, lon_0=-8.5,
              width=473000, height=500000, resolution='i')

lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])
beta = np.array([1, 1, 1, 1])

# Generate the new grid
# Do this if you aren't going to move the grid corners interactively.
# hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)# AttributeError: module 'pyroms' has no attribute 'grid'
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)
# Do this if you are going to use the Boundary Interactor
# map.drawcoastlines()
# xp, yp = map(lonp, latp)
# bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mm+3,Lm+3), proj=map)
# hgrd=bry.grd

lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

# Generate the mask
# for verts in map.coastsegs:
#     hgrd.mask_polygon(verts)
# alternate version from johan.navarro.padron

for xx, yy in map.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy, np.float32)
    vv = np.zeros((xa.shape[0], 2))
    vv[:, 0] = xa
    vv[:, 1] = ya
    hgrd.mask_polygon(vv, mask_value=0)

# Edit the land mask interactively.
# pyroms.grid.edit_mask_mesh(hgrd, proj=map)
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
hmin = 5
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
lon, lat = np.meshgrid(long, lat)
h = griddata((lon.flat, lat.flat), topo.flat, (hgrd.lon_rho, hgrd.lat_rho), method='linear')

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

# save raw bathymetry
hraw = h.copy()

# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# vertical coordinate
theta_b = 2
theta_s = 7.0
Tcline = 50
N = 30
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)

# ROMS grid
grd_name = 'YELLOW'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename='YELLOW_grd_v1.nc')



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