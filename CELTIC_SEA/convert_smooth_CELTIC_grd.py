# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Yellow_Sea/make_YELLOW_grd_v1.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021
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
from netCDF4 import Dataset as netcdf

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

hgrd = Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=map)

lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

for xx, yy in map.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy, np.float32)
    vv = np.zeros((xa.shape[0], 2))
    vv[:, 0] = xa
    vv[:, 1] = ya
    hgrd.mask_polygon(vv, mask_value=False)

# Edit the land mask interactively.
pyroms.grid.edit_mask_mesh(hgrd, proj=map)
# edit_mask_mesh_ij is a faster version using imshow  but no map projection.
coast = pyroms.utility.get_coast_from_map(map)
pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

# depth positive
topo = -data_array

# fix minimum depth
# hmin = 20
# hmin = 15
# hmin = 12
# hmin = 10
hmin = 8
topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry (when gridded/regular)
bath_int_fn = RegularGridInterpolator((lat[:, 0], long[0, :]), topo)

h = bath_int_fn((hgrd.lat_rho.flat, hgrd.lon_rho.flat))
h = np.reshape(h, hgrd.lon_rho.shape)

# insure that depth is always deeper than hmin
h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
# idx = np.where(hgrd.mask_rho == 0)
idx = np.where(np.logical_or(np.equal(hgrd.mask_rho, 0), np.isnan(h)))
h[idx] = hmin

# # save raw bathymetry
# hraw = h.copy()

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

# Smoothed grid, hmin=20m, (datum LAT) smoothed using grid builder and providing stable results
ncgA = netcdf('/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/croco_grd.nc', 'r')
ncgrdA = ncgA.variables
hA = np.array(ncgrdA['h'][:])
ncgA.close()

# Unsmoothed grid, hmin=20m, (datum LAT); produced by PyROMS
ncgB = netcdf('/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/CELTIC_grd_v5_LAT_h20.nc', 'r')
ncgrdB = ncgB.variables
hB = np.array(ncgrdB['h'][:])
ncgB.close()

# Unsmoothed grid, hmin=20m, (datum MSL, using LAT to MSL from INFOMAR); produced by PyROMS
ncgC = netcdf('/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/CELTIC_grd_v5_MSL_h20.nc', 'r')
ncgrdC = ncgC.variables
hC = np.array(ncgrdC['h'][:])
ncgC.close()

h_smth = hA - hB
h_lat2msl = hC - hB

hraw = h + h_lat2msl
h = h + h_smth + h_lat2msl

# vertical coordinate
theta_b = 0.0
theta_s = 7.0
Tcline = hmin
N = 20

# This fn gives in ROMS nomenclature:
# vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=4
vgrd = pyroms.vgrid.s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=2
# vgrd = pyroms.vgrid.s_coordinate_5(h, theta_b, theta_s, Tcline, N, hraw=hraw)  # Vtransform=2 and Vstretching=5

# ROMS grid
grd_name = 'CELTIC_V_MSL_h8'  # hmin
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file

pyroms.grid.write_ROMS_grid(grd, filename='CELTIC_grd_v5_MSL_h8.nc')
