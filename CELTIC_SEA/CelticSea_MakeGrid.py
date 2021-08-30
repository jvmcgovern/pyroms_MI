import gdal
import osr
import numpy as np
import numpy.matlib
import matplotlib
import matplotlib.pyplot as plt

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
