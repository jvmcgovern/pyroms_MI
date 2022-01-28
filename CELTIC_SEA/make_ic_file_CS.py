# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Arctic_GLORYS/make_ic_file.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021
import matplotlib
import os
import subprocess
import pyroms
import pyroms_toolbox
from remap_CS import remap
from remap_uv_CS import remap_uv
matplotlib.use('Agg')

os.environ['PYROMS_GRIDID_FILE'] = '/home/pete/PycharmProjects/pyroms_MI/gridid_MI.txt'

# file = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/CMEMS_IBI' \
#        '/CMEMS_v5r1_IBI_BIO_MY_PdE_01dav_20180101_20180101_R20201201_RE01.nc'

file = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/CMEMS_IBI' \
       '/CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_20180101_20180101_R20201201_RE01.nc'

dst_dir = './'

print('Build IC file from the following file:')
print(file)
print(' ')

src_grd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_CMEMS_IBI(
    '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/IBI-MFC_005_002_mask_bathy.nc', name='IBI', area='global')
# Variables in IBI-MFC_005_002_mask_bathy.nc are:
# deptho      'sea floor depth below geoid' aka bathymetry (2D) (289*361)
# depth        is an array of the standard depths (1D)          (50 long)
# deptho_lev  'Model level number at sea floor' (2D)            (289*361)
# mask_vo     'sea_binary_mask' for v velocities
# mask_uo     'sea_binary_mask' for u velocities
# mask_thetao 'sea_binary_mask' for temperature
# latitude  1D array of values for grid definition             (361 lat)
# longitude 1D array of values for grid definition             (289 long)

dst_grd = pyroms.grid.get_ROMS_grid('CELTIC_V')

# remap
# To avoid Mercator Ocean Fillvalue of -32767 infiltrating the interpolation process, overwrite/modify the FillValue to
# 1e37 using the following commands in the ubuntu command line (use wildcard to apply to all files in folder
# ncatted -a _FillValue,thetao,m,d,1.e37 -a _FillValue,so,m,d,1.e37 -a _FillValue,zos,m,d,1.e37
# -a _FillValue,u0,m,d,1.e37 -a _FillValue,vo,m,d,1.e37 -a _FillValue,bottomT,m,d,1.e37
# CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_20180101_20180101_R20201201_RE01.nc

# Sometimes the FillValue attribute may disappear, meaning that you need to create it (c) instead of modifying (m) it
# Therefore, the repetitive combination
# -a _FillValue,thetao,m,d,1.e37 will need to be replaced with
#                      ^
#                      I
#                      v
# -a _FillValue,thetao,c,d,1.e37

zeta = remap(file, 'zos', src_grd, dst_grd, dst_dir=dst_dir)  # zos is sea surface height
dst_grd = pyroms.grid.get_ROMS_grid('CELTIC_V', zeta=zeta)
remap(file, 'thetao', src_grd, dst_grd, dst_dir=dst_dir)  # thetao is temperature degrees celsius
remap(file, 'so', src_grd, dst_grd, dst_dir=dst_dir)      # so is salinity psu
remap_uv(file, src_grd, dst_grd, dst_dir=dst_dir)         # velocities...

# merge file
ic_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_ic_' + dst_grd.name + '.nc'

out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_zos_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-O', out_file, ic_file)
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_thetao_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file)
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_so_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file)
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_u_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file)
subprocess.call(command)
os.remove(out_file)

out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_v_ic_' + dst_grd.name + '.nc'
command = ('ncks', '-a', '-A', out_file, ic_file)
subprocess.call(command)
os.remove(out_file)
