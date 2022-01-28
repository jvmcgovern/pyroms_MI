# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Arctic_GLORYS/make_remap_weights_file.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021

import pyroms
import pyroms_toolbox
import os
import matplotlib
matplotlib.use('Agg')
# import CGrid_GLORYS

# load the grid srcgrd = CGrid_GLORYS.get_nc_CGrid_GLORYS('/Volumes/P1/Data/GLORYS/data/GL2V1_mesh_mask_new.nc',
# name='GLORYS_CELTIC', area='npolar', ystart=690)
# srcgrd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_GLORYS(
# '/Volumes/P1/Data/GLORYS/data/GL2V1_mesh_mask_new.nc', name='GLORYS_CELTIC', area='npolar', ystart=690)

# srcgrd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_CMEMS_IBI(
#     '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/IBI-MFC_005_002_mask_bathy.nc', name='IBI', area='npolar')

srcgrd = pyroms_toolbox.CGrid_GLORYS.get_nc_CGrid_CMEMS_IBI(
    '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/IBI-MFC_005_002_mask_bathy.nc', name='IBI', area='global')

os.environ['PYROMS_GRIDID_FILE'] = '/home/pete/PycharmProjects/pyroms_MI/gridid_MI.txt'
dstgrd = pyroms.grid.get_ROMS_grid('CELTIC_V')

# make remap grid file for scrip
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='t')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='u')
pyroms_toolbox.CGrid_GLORYS.make_remap_grid_file(srcgrd, Cpos='v')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='rho')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='u')
pyroms.remapping.make_remap_grid_file(dstgrd, Cpos='v')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_IBI_t.nc'
grid2_file = 'remap_grid_CELTIC_V_rho.nc'
interp_file1 = 'remap_weights_IBI_to_CELTIC_V_bilinear_t_to_rho.nc'
interp_file2 = 'remap_weights_CELTIC_V_to_IBI_bilinear_rho_to_t.nc'
map1_name = 'IBI to CELTIC_V Bilinear Mapping'
map2_name = 'CELTIC_V to IBI Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method,
                                       grid1_periodic='.true.', grid2_periodic='.true.')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_IBI_u.nc'
grid2_file = 'remap_grid_CELTIC_V_rho.nc'
interp_file1 = 'remap_weights_IBI_to_CELTIC_V_bilinear_u_to_rho.nc'
interp_file2 = 'remap_weights_CELTIC_V_to_IBI_bilinear_rho_to_u.nc'
map1_name = 'IBI to CELTIC_V Bilinear Mapping'
map2_name = 'CELTIC_V to IBI Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method,
                                       grid1_periodic='.true.', grid2_periodic='.true.')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_IBI_v.nc'
grid2_file = 'remap_grid_CELTIC_V_rho.nc'
interp_file1 = 'remap_weights_IBI_to_CELTIC_V_bilinear_v_to_rho.nc'
interp_file2 = 'remap_weights_CELTIC_V_to_IBI_bilinear_rho_to_v.nc'
map1_name = 'IBI to CELTIC_V Bilinear Mapping'
map2_name = 'CELTIC_V to IBI Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method,
                                       grid1_periodic='.true.', grid2_periodic='.true.')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_IBI_t.nc'
grid2_file = 'remap_grid_CELTIC_V_u.nc'
interp_file1 = 'remap_weights_IBI_to_CELTIC_V_bilinear_t_to_u.nc'
interp_file2 = 'remap_weights_CELTIC_V_to_IBI_bilinear_u_to_t.nc'
map1_name = 'IBI to CELTIC_V Bilinear Mapping'
map2_name = 'CELTIC_V to IBI Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method,
                                       grid1_periodic='.true.', grid2_periodic='.true.')

# compute remap weights
# input namelist variables for bilinear remapping at rho points
grid1_file = 'remap_grid_IBI_t.nc'
grid2_file = 'remap_grid_CELTIC_V_v.nc'
interp_file1 = 'remap_weights_IBI_to_CELTIC_V_bilinear_t_to_v.nc'
interp_file2 = 'remap_weights_CELTIC_V_to_IBI_bilinear_v_to_t.nc'
map1_name = 'IBI to CELTIC_V Bilinear Mapping'
map2_name = 'CELTIC_V to IBI Bilinear Mapping'
num_maps = 1
map_method = 'bilinear'

pyroms.remapping.compute_remap_weights(grid1_file, grid2_file,
                                       interp_file1, interp_file2,
                                       map1_name, map2_name,
                                       num_maps, map_method,
                                       grid1_periodic='.true.', grid2_periodic='.true.')
