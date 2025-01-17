# File inspired by/modified from:
# /home/pete/PycharmProjects/pyroms_MI/examples/Arctic_GLORYS/remap_uv.py
# Joe McGovern, Marine Institute, Rinville West, Rinville, Oranmore, Co. Galway 2021

import numpy as np
import os
try:
    import netCDF4 as netCDF
except:
    import netCDF3 as netCDF
from datetime import datetime
# from matplotlib.dates import date2num, num2date
from matplotlib.dates import date2num
import pyroms
import pyroms_toolbox
# import _remapping
# import matplotlib.pyplot as plt
# import time


class nctime(object):
    pass


def remap_uv(src_file, src_grd, dst_grd, dmax=0, cdepth=0, kk=0, dst_dir='./'):
    # ystart = 690
    ystart = 0
    # get time
    nctime.long_name = 'time'
    # nctime.units = 'days since 1900-01-01 00:00:00'
    # # time reference "days since 1900-01-01 00:00:00"
    # ref = datetime(1900, 1, 1, 0, 0, 0)
    # ref = date2num(ref)
    nctime.units = 'days since 1968-05-23 00:00:00'
    # time reference "days since 1968-05-23 00:00:00"
    ref = datetime(1968, 5, 23, 0, 0, 0)
    ref = date2num(ref)
    # tag = src_file.rsplit('/')[-1].rsplit('_')[2]
    tag = src_file.rsplit('/')[-1].rsplit('_')[7]
    print(("tag:", tag))
    year = int(tag[:4])
    month = int(tag[4:6])
    day = int(tag[6:])
    time = datetime(year, month, day, 0, 0, 0)
    time = date2num(time)
    time = time - ref
    time = time + 0.5  # 1-day average

    # get dimensions
    Mp, Lp = dst_grd.hgrid.mask_rho.shape

    # create destination file
    dst_file = src_file.rsplit('/')[-1]
    dst_fileu = dst_dir + dst_file[11:22] + dst_file[32:40] + '_u_ic_' + dst_grd.name + '.nc'
    print('\nCreating destination file', dst_fileu)
    if os.path.exists(dst_fileu) is True:
        os.remove(dst_fileu)
    pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
    dst_filev = dst_dir + dst_file[11:22] + dst_file[32:40] + '_v_ic_' + dst_grd.name + '.nc'
    print('Creating destination file', dst_filev)
    if os.path.exists(dst_filev) is True:
        os.remove(dst_filev)
    pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)

    # open destination file
    ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    # load var
    cdf = netCDF.Dataset(src_file)
    # src_varu = cdf.variables['vozocrtx']
    # src_varv = cdf.variables['vomecrty']
    src_varu = cdf.variables['uo']
    src_varv = cdf.variables['vo']
    print("dims", src_varu.dimensions, src_varv.dimensions)

    # get missing value
    spval = src_varu._FillValue

    # ARCTIC grid sub-sample
    src_varu = src_varu[:]
    src_varv = src_varv[:]
    print("shape 1", src_varu.shape, src_varv.shape)
    src_varu = np.squeeze(src_varu)
    src_varv = np.squeeze(src_varv)
    print("shape 2", src_varu.shape, src_varv.shape)
    src_varu = src_varu[:, np.r_[ystart:np.size(src_varu, 1), -1], :]
    src_varv = src_varv[:, np.r_[ystart:np.size(src_varv, 1), -1], :]
    print("shape 3", src_varu.shape, src_varv.shape)

    # get weights file
    wts_file_a = 'remap_weights_IBI_to_CELTIC_bilinear_t_to_rho.nc'
    wts_file_u = 'remap_weights_IBI_to_CELTIC_bilinear_u_to_rho.nc'
    wts_file_v = 'remap_weights_IBI_to_CELTIC_bilinear_v_to_rho.nc'

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1, 0, 0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name + '_Z', dst_grd.hgrid, dst_zcoord)

    # create variable in destination file
    print('Creating variable u')
    ncu.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['u'].long_name = '3D u-momentum component'
    ncu.variables['u'].units = 'meter second-1'
    ncu.variables['u'].field = 'u-velocity, scalar, series'
    # create variable in destination file
    print('Creating variable ubar')
    ncu.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['ubar'].long_name = '2D u-momentum component'
    ncu.variables['ubar'].units = 'meter second-1'
    ncu.variables['ubar'].field = 'ubar-velocity, scalar, series'

    print('Creating variable v')
    ncv.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['v'].long_name = '3D v-momentum component'
    ncv.variables['v'].units = 'meter second-1'
    ncv.variables['v'].field = 'v-velocity, scalar, series'
    print('Creating variable vbar')
    ncv.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['vbar'].long_name = '2D v-momentum component'
    ncv.variables['vbar'].units = 'meter second-1'
    ncv.variables['vbar'].field = 'vbar-velocity, scalar, series'

    # remapping
    print('remapping and rotating u and v from', src_grd.name,
          'to', dst_grd.name)
    print('time =', time)

    # flood the grid
    print('flood the grid', src_varu.shape)
    src_uz = pyroms_toolbox.CGrid_GLORYS.flood(src_varu, src_grd, Cpos='u',
                                               spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
    src_vz = pyroms_toolbox.CGrid_GLORYS.flood(src_varv, src_grd, Cpos='v',
                                               spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
    # This IIIIIIII is how tracers are processed (for comparison with u/v processing above)
    #      VVVVVVVV
    # src_varz = pyroms_toolbox.CGrid_GLORYS.flood(src_var, src_grd, Cpos=['t', 'u' or 'v'],
    #                                              spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
    #
    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_uz = pyroms.remapping.remap(src_uz, wts_file_u, spval=spval)
    dst_vz = pyroms.remapping.remap(src_vz, wts_file_v, spval=spval)

    # This IIIIIIII is how tracers are processed (for comparison with u/v processing above)
    #      VVVVVVVV
    # dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)

    # vertical interpolation from standard z level to sigma
    print('vertical interpolation from standard z level to sigma')
    dst_u = pyroms.remapping.z2roms(dst_uz[::-1, :, :], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=False)
    dst_v = pyroms.remapping.z2roms(dst_vz[::-1, :, :], dst_grdz, dst_grd, Cpos='rho', spval=spval, flood=False)

    # This IIIIIIII is how tracers are processed (for comparison with u/v processing above)
    #      VVVVVVVV
    # dst_var = pyroms.remapping.z2roms(dst_varz[::-1, :, :], dst_grdz, dst_grd, Cpos=Cpos, spval=spval, flood=False)

    # rotate u,v fields
    src_angle = src_grd.angle
    src_angle = pyroms.remapping.remap(src_angle, wts_file_a)
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    U = dst_u + dst_v * 1j
    eitheta = np.exp(-1j * angle[:, :, :])
    U = U * eitheta
    dst_u = np.real(U)
    dst_v = np.imag(U)

    # move back to u,v points
    dst_u = 0.5 * (dst_u[:, :, :-1] + dst_u[:, :, 1:])
    dst_v = 0.5 * (dst_v[:, :-1, :] + dst_v[:, 1:, :])

    # spval
    idxu = np.where(dst_grd.hgrid.mask_u == 0)
    idxv = np.where(dst_grd.hgrid.mask_v == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u[n, idxu[0], idxu[1]] = spval
        dst_v[n, idxv[0], idxv[1]] = spval

    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u = 0.5 * (dst_grd.vgrid.z_w[0, :, :, :-1] + dst_grd.vgrid.z_w[0, :, :, 1:])
    z_v = 0.5 * (dst_grd.vgrid.z_w[0, :, :-1, :] + dst_grd.vgrid.z_w[0, :, 1:, :])

    dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
    dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))

    for i in range(dst_ubar.shape[1]):
        for j in range(dst_ubar.shape[0]):
            dst_ubar[j, i] = (dst_u[:, j, i] * np.diff(z_u[:, j, i])).sum() / -z_u[0, j, i]

    for i in range(dst_vbar.shape[1]):
        for j in range(dst_vbar.shape[0]):
            dst_vbar[j, i] = (dst_v[:, j, i] * np.diff(z_v[:, j, i])).sum() / -z_v[0, j, i]

    # spval
    dst_ubar[idxu[0], idxu[1]] = spval
    dst_vbar[idxv[0], idxv[1]] = spval

    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = time
    ncu.variables['u'][0] = dst_u
    ncu.variables['ubar'][0] = dst_ubar

    ncv.variables['ocean_time'][0] = time
    ncv.variables['v'][0] = dst_v
    ncv.variables['vbar'][0] = dst_vbar

    print(dst_u.shape)
    print(dst_ubar.shape)
    print(dst_v.shape)
    print(dst_vbar.shape)

    # close destination file
    ncu.close()
    ncv.close()
    cdf.close()
