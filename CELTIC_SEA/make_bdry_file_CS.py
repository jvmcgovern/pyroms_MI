# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Palau_HYCOM/make_bdry_file.py and
# /home/pete/PycharmProjects/pyroms_MI/examples/CCS1_SODA3.3.1/make_bdry_file.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021
import matplotlib
import os
import sys
import subprocess
import numpy as np
import pyroms
import pyroms_toolbox
import glob
from remap_bdry_CS import remap_bdry
from remap_bdry_uv_CS import remap_bdry_uv
from nco import Nco
matplotlib.use('Agg')

os.environ['PYROMS_GRIDID_FILE'] = '/home/pete/PycharmProjects/pyroms_MI/gridid_MI.txt'

# IF running from python terminal outside of debug mode
# year = int(sys.argv[1])
# lst_year = [year]
year = 2018
lst_year = [2018]

data_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/CMEMS_IBI/'
dst_dir = './'

lst_file = []

for year in lst_year:
    year = np.str(year)
    lst = subprocess.getoutput('ls ' + data_dir + 'CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_' + year + '*')
    lst = lst.split()
    lst_file = lst_file + lst

print('Build OBC file from the following file list:')
print(lst_file)
print(' ')

src_grd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_CMEMS_IBI(
    '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/IBI-MFC_005_002_mask_bathy.nc', name='IBI', area='global')

dst_grd = pyroms.grid.get_ROMS_grid('CELTIC')

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

input_string = "/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/CMEMS_IBI/" \
               "CMEMS_v5r1_IBI_PHY_MY_PdE_01dav_201801*_201801*_R*_RE01.nc"
list_thetao = ['ncatted', '-a', '_FillValue,thetao,m,d,1.e37']
list_so = ['ncatted', '-a', '_FillValue,so,m,d,1.e37']
list_zos = ['ncatted', '-a', '_FillValue,zos,m,d,1.e37']
list_uo = ['ncatted', '-a', '_FillValue,uo,m,d,1.e37']
list_vo = ['ncatted', '-a', '_FillValue,vo,m,d,1.e37']
list2 = glob.glob(input_string)
for mf in range(len(list2)):
    cthetao = list_thetao.copy()
    cthetao = cthetao + [list2[mf]]
    subprocess.call(cthetao[:])
    cso = list_so.copy()
    cso = cso + [list2[mf]]
    subprocess.call(cso[:])
    czos = list_zos.copy()
    czos = czos + [list2[mf]]
    subprocess.call(czos[:])
    cuo = list_uo.copy()
    cuo = cuo + [list2[mf]]
    subprocess.call(cuo[:])
    cvo = list_vo.copy()
    cvo = cvo + [list2[mf]]
    subprocess.call(cvo[:])

for file in lst_file:
    zeta = remap_bdry(file, 'zos', src_grd, dst_grd, dst_dir=dst_dir)
    dst_grd = pyroms.grid.get_ROMS_grid('CELTIC', zeta=zeta)
    remap_bdry(file, 'thetao', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry(file, 'so', src_grd, dst_grd, dst_dir=dst_dir)
    remap_bdry_uv(file, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    bdry_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_bdry_' + dst_grd.name + '.nc'
    out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_zos_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-O', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_thetao_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_so_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_u_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)
    out_file = dst_dir + file.rsplit('/')[-1][11:22] + file.rsplit('/')[-1][32:40] + '_v_bdry_' + dst_grd.name + '.nc'
    command = ('ncks', '-a', '-A', out_file, bdry_file)
    subprocess.check_call(command)
    os.remove(out_file)

# # ncrcat -d ocean_time,0,100 IBI_PHY_MY_201801*_bdry_CELTIC.nc IBI_PHY_MY_201801_bdry_CELTIC.nc
# # The following combination may merge the daily boundary files into a single file
# input_string = "/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/" \
#                "IBI_PHY_MY_201801*_bdry_CELTIC.nc"
# list_inputs = ['ncrcat', '-d', 'ocean_time,0,100']
# list_fils = glob.glob(input_string)
# list_inputs = list_inputs + list_fils[:] + ['IBI_PHY_MY_201801_bdry_CELTIC.nc']
# subprocess.call(list_inputs[:])



