# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/make_tide/make_remap_weights_file.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021

import pyroms
import os
import pyroms_toolbox
from examples.make_tide.CGrid_TPXO8 import get_nc_CGrid_TPXO8, make_remap_grid_file
# import CGrid_TPXO8

# load the grid
# pth_tpxo = '/archive/u1/uaf/kate/tides/tpxo8/'
# srcgrd = CGrid_TPXO8.get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8atlas_30_v1.nc', \
#       xrange=(4400, 5100), yrange=(5200, 6500))

os.environ['PYROMS_GRIDID_FILE'] = '/home/pete/PycharmProjects/pyroms_MI/gridid_MI.txt'

pth_tpxo = '/home/pete/PycharmProjects/pyroms_MI/Tideforce/TPX/TPX_08/'
# srcgrd = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8atlas_30.nc', xrange=(4141, 4290), yrange=(121, 330))
# srcgrd_lr = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8_atlas6.nc', name='TPXO8atlas6',
#                                xrange=(829, 858), yrange=(25, 66))
# srcgrd = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8atlas_30.nc', xrange=(121, 330), yrange=(4141, 4290))
# srcgrd_lr = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8_atlas6.nc', name='TPXO8atlas6',
#                                xrange=(25, 66), yrange=(829, 858))
# srcgrd = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8atlas_30.nc', xrange=(4141, 4290), yrange=(121, 330))
# srcgrd_lr = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8_atlas6.nc', name='TPXO8atlas6',
#                                xrange=(829, 858), yrange=(25, 66))
srcgrd = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8atlas_30.nc', xrange=(4141, 4290), yrange=(5071, 5280))
srcgrd_lr = get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8_atlas6.nc', name='TPXO8atlas6',
                               xrange=(829, 858), yrange=(1009, 1056))
dstgrd = pyroms.grid.get_ROMS_grid('CELTIC')

# make remap grid file for scrip
make_remap_grid_file(srcgrd, Cpos='t')
make_remap_grid_file(srcgrd, Cpos='u')
make_remap_grid_file(srcgrd, Cpos='v')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_' + srcgrd.name + '_t.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_t.nc'
map1_name = srcgrd.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name = dstgrd.name + ' to ' + srcgrd.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method)


grid1_file = 'remap_grid_' + srcgrd.name + '_u.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_u_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_u.nc'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method)


grid1_file = 'remap_grid_' + srcgrd.name + '_v.nc'
grid2_file = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1 = 'remap_weights_' + srcgrd.name + '_to_' + dstgrd.name + '_bilinear_v_to_rho.nc'
interp_file2 = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd.name + '_bilinear_rho_to_v.nc'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method)

make_remap_grid_file(srcgrd_lr, Cpos='t')
make_remap_grid_file(srcgrd_lr, Cpos='u')
make_remap_grid_file(srcgrd_lr, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file_lr = 'remap_grid_' + srcgrd_lr.name + '_t.nc'
grid2_file_lr = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1_lr = 'remap_weights_' + srcgrd_lr.name + '_to_' + dstgrd.name + '_bilinear_t_to_rho.nc'
interp_file2_lr = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd_lr.name + '_bilinear_rho_to_t.nc'
map1_name_lr = srcgrd_lr.name + ' to ' + dstgrd.name + ' Bilinear Mapping'
map2_name_lr = dstgrd.name + ' to ' + srcgrd_lr.name + ' Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file_lr, grid2_file_lr,
                                       interp_file1_lr, interp_file2_lr,
                                       map1_name_lr, map2_name_lr,
                                       num_maps, map_method)


grid1_file_lr = 'remap_grid_' + srcgrd_lr.name + '_u.nc'
grid2_file_lr = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1_lr = 'remap_weights_' + srcgrd_lr.name + '_to_' + dstgrd.name + '_bilinear_u_to_rho.nc'
interp_file2_lr = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd_lr.name + '_bilinear_rho_to_u.nc'

pyroms.remapping.compute_remap_weights(grid1_file_lr, grid2_file_lr,
                                       interp_file1_lr, interp_file2_lr,
                                       map1_name_lr, map2_name_lr,
                                       num_maps, map_method)


grid1_file_lr = 'remap_grid_' + srcgrd_lr.name + '_v.nc'
grid2_file_lr = 'remap_grid_' + dstgrd.name + '_rho.nc'
interp_file1_lr = 'remap_weights_' + srcgrd_lr.name + '_to_' + dstgrd.name + '_bilinear_v_to_rho.nc'
interp_file2_lr = 'remap_weights_' + dstgrd.name + '_to_' + srcgrd_lr.name + '_bilinear_rho_to_v.nc'

pyroms.remapping.compute_remap_weights(grid1_file_lr, grid2_file_lr,
                                       interp_file1_lr, interp_file2_lr,
                                       map1_name_lr, map2_name_lr,
                                       num_maps, map_method)
