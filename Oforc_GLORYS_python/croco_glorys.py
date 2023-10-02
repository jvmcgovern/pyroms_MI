#
#
#########################################################################
#
#
# croco_glorys.py
#
#
# A collection of functions for interpolation of GLORYS ocean reanalysis 
# data to create CROCO initial and boundary condition files
#
#
# Gustav Rautenbach, Steven Herbette, Pierrick Penven, 2021
#
#########################################################################
#
#
#  geo_idx = geo_idx(dd, dd_array)
#
#    Get the closest index of a particular coordinate
#
#
#  (elem,coef) = get_tri_coef(X, Y, newX, newY, verbose=0)
#
#    Get Delaunay linear interpolation pointers and coefficients
#
#
#  varnew =  horiz_interp_delaunay(lonold,latold,varold,lonnew,latnew,elem=0,coef=0)
#
#    Do an horizontal interpolation
#
#
#  Vout,NzGood = interp_tracers(nc,vname,l,k,imin,imax,jmin,jmax,Lon,Lat,coef,elem)
#
#    Remove the missing values from a gridded 2D field from a Netcdf file
#    and do an horizontal interpolation using Delaunay matrices (coef and elem)
#
#
#  vout=interp3d(nc,vname,tndx_glo,Nzgoodmin,depth,z_rho,imin,imax,jmin,jmax,Lon,Lat,coef,elem)
#
#    Do a full interpolation of a 3d variable from GLORYS to a CROCO sigma grid
#    1 - Horizontal interpolation on each GLORYS levels
#    2 - Vertical Interpolation from z to CROCO sigma levels
#
#
#  uout,vout=interp3d_uv(nc,tndx_glo,Nzgoodmin,depth,z_rho,cosa,sina,\
#                iminU,imaxU,jminU,jmaxU,LonU,LatU,coefU,elemU,\
#                iminV,imaxV,jminV,jmaxV,LonV,LatV,coefV,elemV)
#    Do a full interpolation of a 3d U and V from GLORYS to a CROCO sigma grid
#    Do a rotation to align with CROCO grid (if not rectangular)
#
#
#
#  create_inifile(ininame,grdname,title,theta_s,theta_b,hc,N,ini_time,vtransform)
#
#    This function create a Netcdf initial file
#
#
#  create_bryfile(bryname,grdname,title,obc...
#                          theta_s,theta_b,hc,N,...
#                          bry_time,cycle,clobber)
#
#    This function create a Netcdf boundary file
#
#
#  (LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,\
#   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,\
#   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry)\
#   = get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
#
#     Get Delaunay linear interpolation pointers and coefficients for each boundary
#
#  ncbry = interp_bry(obctype,ncglo,tndx_glo,ncbry,tndx_bry,\
#               h_bry,theta_s,theta_b,hc,N,vtransform,\
#	        Nzgoodmin,depth,angle_bry,\
#               LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,\
#               LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,\
#	        LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry)
#
#
#    Do a full interpolation of a 3d variable from GLORYS to a CROCO sigma grid
#    for each boundary
#
#
########################################################################
#
#
#  This is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 2 of the License,
#  or (at your option) any later version.
#
#  CROCOTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#
######################################################################
#
#

import numpy as np
from netCDF4 import Dataset as netcdf
from scipy.interpolate import griddata
import os
import sys
import time
# sys.path.insert(0,'/XXX/')
import croco_vgrid as vgrd
import croco_glorys as glor
from interp_Cgrid import *
from scipy.spatial import Delaunay
from progressbar import *
from netCDF4 import date2index as d2i
from datetime import date, datetime
import ttide as tt


#
#
#########################################################################
#
#
# get index of particular coordinate
#
#

def geo_idx(dd, dd_array):
    """
     - dd - the decimal degree (latitude or longitude)
     - dd_array - the list of decimal degrees to search.
     search for nearest decimal degree in an array of decimal degrees and return the index.
     np.argmin returns the indices of minium value along an axis.
     so subtract dd from all values in dd_array, take absolute value and find index of minium.
   """
    geo_idx = (np.abs(dd_array - dd)).argmin()
    return geo_idx


#
#
#########################################################################
#
#


#
#
#########################################################################
#
#

def get_tri_coef(X, Y, newX, newY, verbose=0):
    """
    Inputs:
        origin lon and lat 2d arrays (X,Y)
        child lon and lat 2d arrays (newX,newY)
    Ouputs:
        elem - pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
        coef - linear interpolation coefficients
    Use:
        To subsequently interpolate data from Fp to Fc, the following
        will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
        griddata calls. Since it avoids repeated triangulations and tsearches (that are done
        with every call to griddata) it should be much faster.
    """

    Xp = np.array([X.ravel(), Y.ravel()]).T
    Xc = np.array([newX.ravel(), newY.ravel()]).T

    # Compute Delaunay triangulation
    if verbose == 1:
        tstart = time.time()
    tri = Delaunay(Xp)
    if verbose == 1:
        print('Delaunay Triangulation', time.time() - tstart)

    # Compute enclosing simplex and barycentric coordinate (similar to tsearchn in MATLAB)
    npts = Xc.shape[0]
    p = np.zeros((npts, 3))

    points = tri.points[tri.vertices[tri.find_simplex(Xc)]]
    if verbose == 1:
        tstart = time.time()
    for i in progressbar(range(npts), '  Get_tri_coef: ', 40):

        if verbose == 1:
            print(np.float(i) / npts)

        if tri.find_simplex(Xc[i]) == -1:  # Point outside triangulation
            p[i, :] = p[i, :] * np.nan

        else:

            if verbose == 1:
                tstart = time.time()
            A = np.append(np.ones((3, 1)), points[i], axis=1)
            if verbose == 1:
                print('append A', time.time() - tstart)

            if verbose == 1:
                tstart = time.time()
            B = np.append(1., Xc[i])
            if verbose == 1:
                print('append B', time.time() - tstart)

            if verbose == 1:
                tstart = time.time()
            p[i, :] = np.linalg.lstsq(A.T, B.T)[0]
            if verbose == 1:
                print('solve', time.time() - tstart)

    if verbose == 1:
        print('Coef. computation 1', time.time() - tstart)

    if verbose == 1:
        tstart = time.time()
        elem = np.reshape(tri.vertices[tri.find_simplex(Xc)], (newX.shape[0], newY.shape[1], 3))
        coef = np.reshape(p, (newX.shape[0], newY.shape[1], 3))
    if verbose == 1:
        print('Coef. computation 2', time.time() - tstart)

    return elem, coef


#
#
#########################################################################
#
#


#
#
#########################################################################
#
#

def horiz_interp_delaunay(lonold, latold, varold, lonnew, latnew, elem=0, coef=0):
    #
    # horizontal interpolation
    #
    """
    lonold: 2D original longitude matrix
    latold: 2D original longitude matrix
    varold: 2D original datra matrix
    lonnew: 2D new longitude matrix
    latnew: 2D new longitude matrix
    print(len(args*))
    """

    # Horizontal Interpolation from croco variables' grid to new grid
    #: get interpolation coefficients
    if np.all(elem == 0):
        [elem, coef] = get_tri_coef(lonold, latold, lonnew, latnew)
        coefnorm = np.sum(coef, axis=2)
        coef = coef / coefnorm[:, :, np.newaxis]
        varnew = np.sum(coef * varold.ravel()[elem], 2)
        return elem, coef, varnew
    else:
        varnew = np.sum(coef * varold.ravel()[elem], 2)
        return varnew


#
#
#########################################################################
#
#


#
#
# #####################################################################
# ##### FUNCTION INTERP_TRACERS #######################################
# #####################################################################
#
#

def interp_tracers(nc, vname, l, k, imin, imax, jmin, jmax, Lon, Lat, coef, elem, bdi):
    #
    #
    #  Remove the missing values from a gridded 2D field
    #  and do an horizontal interpolation using Delaunay matrices (coef and elem)
    #
    #
    # #####################################################################
    #
    #
    #

    #
    # 1: Read data
    #
    if k == -1:
        Vin = np.array(nc[vname][l, jmin:jmax, imin:imax])
    else:
        Vin = np.array(nc[vname][l, k, jmin:jmax, imin:imax])

    #
    # 2: Remove bad values (using nearest values)
    #

    # igood = np.where(Vin < 1000.)
    # ibad = np.where(Vin >= 1000.)
    if bdi == 'add_offset':
        try:
            igood = np.where(Vin != nc[vname].add_offset)
            ibad = np.where(Vin == nc[vname].add_offset)
        except:
            igood = np.where(Vin != nc[vname]._FillValue)
            ibad = np.where(Vin == nc[vname]._FillValue)
    elif bdi == '_FillValue':
        try:
            igood = np.where(Vin != nc[vname]._FillValue)
            ibad = np.where(Vin == nc[vname]._FillValue)
        except:
            igood = np.where(Vin != nc[vname].add_offset)
            ibad = np.where(Vin == nc[vname].add_offset)
    else:
        try:
            igood = np.where(np.isfinite(Vin))
            ibad = np.where(np.isnan(Vin))
        except:
            igood = np.where(Vin >= 0)
            ibad = np.where(Vin < 0)

    NzGood = np.size(igood)
    Nbad = np.size(ibad)

    if NzGood == 0:

        print('Warning: no good data')
        Vin[:] = np.nan

    elif NzGood < 10:

        print('Warning: less than 10 good values')
        Vin[:] = np.mean(Vin[igood])

    elif Nbad > 0:
        # if k != -1 and k <= 3:
        if k != -1 and k > 0:
            # JVMCG 230223 had added the section below to prevent interpolation using fill vals (in case)
            Vin = np.array(nc[vname][l, :, jmin:jmax, imin:imax])

            # dchk = np.array(nc[vname][l, :, jmin:jmax, imin:imax])
            # maxdep = np.zeros_like(dchk[0, :, :])
            for chk in range(int(Nbad / 2)):
                Vin[k:, ibad[0][chk], ibad[1][chk]] = Vin[k - 1, ibad[0][chk], ibad[1][chk]]
            # for j in range(dchk.shape[1]):
            #     for i in range(dchk.shape[2]):
            #         try:
            #             maxdep[j, i] = np.argwhere(np.bitwise_or(dchk[:, j, i] < -100, dchk[:, j, i] > 10000))[0][0]
            #             Vin[int(maxdep[j, i]):, j, i] = Vin[int(maxdep[j, i] - 1), j, i]
            #         except:
            #             maxdep[j, i] = 0
            #             Vin[:, j, i] = Vin[:, j, i]
            Vin = Vin[k, :, :]

        Vin[ibad] = griddata((Lon[igood], Lat[igood]), Vin[igood], (Lon[ibad], Lat[ibad]), method='nearest')

    #
    # 3: 2D interpolation
    #

    Vout = np.sum(coef * Vin.ravel()[elem], 2)

    return Vout, NzGood


def interp_sst(nc, vname, l, k, imin, imax, jmin, jmax, Lon, Lat, coef, elem):
    #
    #
    #  Remove the missing values from a gridded 2D field
    #  and do an horizontal interpolation using Delaunay matrices (coef and elem)
    #
    #
    # #####################################################################
    #
    #
    #

    #
    # 1: Read data
    #

    if k == -1:
        Vin = np.array(nc[vname][l, jmin:jmax, imin:imax])
    else:
        Vin = np.array(nc[vname][l, k, jmin:jmax, imin:imax])

    #
    # 2: Remove bad values (using nearest values)
    #
    try:
        igood = np.where(Vin != nc[vname]._FillValue)
        ibad = np.where(Vin == nc[vname]._FillValue)
    except:
        igood = np.where(Vin != nc[vname].add_offset)
        ibad = np.where(Vin == nc[vname].add_offset)

    NzGood = np.size(igood)
    Nbad = np.size(ibad)

    if NzGood == 0:

        print('Warning: no good data')
        Vin[:] = np.nan

    elif NzGood < 10:

        print('Warning: less than 10 good values')
        Vin[:] = np.mean(Vin[igood])

    else:

        Vin[ibad] = griddata((Lon[igood], Lat[igood]), Vin[igood], (Lon[ibad], Lat[ibad]), method='nearest')

    #
    # 3: 2D interpolation
    #

    Vout = np.sum(coef * Vin.ravel()[elem], 2)

    return Vout, NzGood


#
#
#
# #####################################################################
# ##### END FUNCTION INTERP_TRACERS ###################################
# #####################################################################
#
#


#
#
# #####################################################################
# #### FUNCTION INTERP3D ##############################################
# #####################################################################
#
#

def interp3d(nc, vname, tndx_glo, Nzgoodmin, depth, z_rho, imin, imax, jmin, jmax, Lon, Lat, coef, elem, bdi):
    #
    #  Do a full interpolation of a 3d variable from GLORYS to a CROCO sigma grid
    #
    #  1 - Horizontal interpolation on each GLORYS levels
    #  2 - Vertical Interpolation from z to CROCO sigma levels
    #
    #
    ######################################################################
    #
    #
    #

    [N, M, L] = np.shape(z_rho)
    [Nz] = np.shape(depth)

    comp_horzinterp = 1  # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)

    if comp_horzinterp == 1:

        print('Horizontal interpolation over z levels')

        t3d = np.zeros((Nz, M, L))
        kgood = -1

        maxdep = -np.min(z_rho)
        mdidx = np.argwhere(depth > maxdep)[0][0]
        depth = depth[:mdidx]
        [Nz] = np.shape(depth)

        for k in progressbar(range(Nz), vname + ': ', 40):

            (t2d, Nzgood) = interp_tracers(nc, vname, tndx_glo, k, imin, imax, jmin, jmax,
                                           Lon, Lat, coef, elem, bdi)

            if Nzgood > Nzgoodmin:
                kgood = kgood + 1
                t3d[kgood, :, :] = t2d

        t3d = t3d[0:kgood, :, :]
        depth = depth[0:kgood]

        np.savez('t3d.npz', t3d=t3d, depth=depth)

    else:

        print('Load matrix...')
        data = np.load('t3d.npz')
        t3d = data['t3d']
        depth = data['depth']

    [Nz] = np.shape(depth)
    Z = -depth

    # croco_surfd = np.nanmax(z_rho)

    #
    # ----------------------------------------------------
    #  Vertical interpolation
    # ----------------------------------------------------
    #

    print('Vertical interpolation')

    #
    # Add a layer below the bottom and above the surface to avoid vertical extrapolations
    # and flip the matrices upside down (Z[Nz]=surface)
    #

    Z = np.flipud(np.concatenate(([100.], Z, [-10000.])))

    [Nz] = np.shape(Z)

    t3d = np.flipud(vgrd.add2layers(t3d))

    #
    # Do the vertical interpolations
    #

    vout = vgrd.ztosigma(t3d, Z, z_rho)

    return vout


#
#
# #####################################################################
# #### END FUNCTION INTERP3D ##########################################
# #####################################################################
#
#

#
#
# #####################################################################
# #### FUNCTION INTERP3D_UV ###########################################
# #####################################################################
#
#

def interp3d_uv(nc, tndx_glo, Nzgoodmin, depth, z_rho, cosa, sina,
                iminU, imaxU, jminU, jmaxU, LonU, LatU, coefU, elemU,
                iminV, imaxV, jminV, jmaxV, LonV, LatV, coefV, elemV):
    #
    #  Do a full interpolation of a 3d variable from GLORYS to a CROCO sigma grid
    #
    #  1 - Horizontal interpolation on each GLORYS levels
    #  2 - Vertical Interpolation from z to CROCO sigma levels
    #
    #
    # #####################################################################
    #
    #
    #

    [N, M, L] = np.shape(z_rho)
    [Nz] = np.shape(depth)

    comp_horzinterp = 1  # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)

    if comp_horzinterp == 1:

        print('Horizontal interpolation of u and v over z levels')

        u3d = np.zeros((Nz, M, L - 1))
        v3d = np.zeros((Nz, M - 1, L))
        kgood = -1

        maxdep = -np.min(z_rho)
        mdidx = np.argwhere(depth > maxdep)[0][0]
        depth = depth[:mdidx]
        [Nz] = np.shape(depth)

        for k in progressbar(range(Nz), ' uv : ', 40):

            (u2d, Nzgood_u) = interp_tracers(nc, 'uo', tndx_glo, k, iminU, imaxU,
                                             jminU, jmaxU, LonU, LatU, coefU, elemU, '_FillValue')
            (v2d, Nzgood_v) = interp_tracers(nc, 'vo', tndx_glo, k, iminV, imaxV,
                                             jminV, jmaxV, LonV, LatV, coefV, elemV, '_FillValue')

            Nzgood = np.min((Nzgood_u, Nzgood_v))

            if Nzgood > Nzgoodmin:
                kgood = kgood + 1

                #
                # Rotation and put to u-points and v-points
                #

                u3d[kgood, :, :] = rho2u_2d(u2d * cosa + v2d * sina)
                v3d[kgood, :, :] = rho2v_2d(v2d * cosa - u2d * sina)

        u3d = u3d[0:kgood, :, :]
        v3d = v3d[0:kgood, :, :]
        depth = depth[0:kgood]

        np.savez('u3d.npz', u3d=u3d, v3d=v3d, depth=depth)

    else:

        print('Load matrices...')
        data = np.load('u3d.npz')
        u3d = data['u3d']
        v3d = data['v3d']
        depth = data['depth']

    [Nz] = np.shape(depth)
    Z = -depth

    #
    # ----------------------------------------------------
    #  Vertical interpolation
    # ----------------------------------------------------
    #

    print('Vertical interpolation')

    #
    # Add a layer below the bottom and above the surface to avoid vertical extrapolations
    # and flip the matrices upside down (Z[Nz]=surface)
    #
    #
    Z = np.flipud(np.concatenate(([100.], Z, [-10000.])))

    [Nz] = np.shape(Z)

    u3d = np.flipud(vgrd.add2layers(u3d))
    v3d = np.flipud(vgrd.add2layers(v3d))

    #
    # Do the vertical interpolations
    #

    uout = vgrd.ztosigma(u3d, Z, rho2u_3d(z_rho))
    vout = vgrd.ztosigma(v3d, Z, rho2v_3d(z_rho))

    return uout, vout


def interp3d_uv_ESM(ncu, ncv, tndx_glo, Nzgoodmin, depth, z_rho, cosa, sina,
                    iminU, imaxU, jminU, jmaxU, LonU, LatU, coefU, elemU,
                    iminV, imaxV, jminV, jmaxV, LonV, LatV, coefV, elemV):
    #
    #  Do a full interpolation of a 3d variable from GLORYS to a CROCO sigma grid
    #
    #  1 - Horizontal interpolation on each GLORYS levels
    #  2 - Vertical Interpolation from z to CROCO sigma levels
    #
    #
    # #####################################################################
    #
    #
    #

    [N, M, L] = np.shape(z_rho)
    [Nz] = np.shape(depth)

    comp_horzinterp = 1  # if 1 compute horizontal interpolations - 0 use saved matrices (for debugging)

    if comp_horzinterp == 1:

        print('Horizontal interpolation of u and v over z levels')

        u3d = np.zeros((Nz, M, L - 1))
        v3d = np.zeros((Nz, M - 1, L))
        kgood = -1

        maxdep = -np.min(z_rho)
        mdidx = np.argwhere(depth > maxdep)[0][0]
        depth = depth[:mdidx]
        [Nz] = np.shape(depth)

        for k in progressbar(range(Nz), ' uv : ', 40):

            (u2d, Nzgood_u) = interp_tracers(ncu, 'uo_unrotated', tndx_glo, k, iminU, imaxU,
                                             jminU, jmaxU, LonU, LatU, coefU, elemU, '_FillValue')
            (v2d, Nzgood_v) = interp_tracers(ncv, 'vo_unrotated', tndx_glo, k, iminV, imaxV,
                                             jminV, jmaxV, LonV, LatV, coefV, elemV, '_FillValue')

            Nzgood = np.min((Nzgood_u, Nzgood_v))

            if Nzgood > Nzgoodmin:
                kgood = kgood + 1

                #
                # Rotation and put to u-points and v-points
                #

                u3d[kgood, :, :] = rho2u_2d(u2d * cosa + v2d * sina)
                v3d[kgood, :, :] = rho2v_2d(v2d * cosa - u2d * sina)

        u3d = u3d[0:kgood, :, :]
        v3d = v3d[0:kgood, :, :]
        depth = depth[0:kgood]

        np.savez('u3d.npz', u3d=u3d, v3d=v3d, depth=depth)

    else:

        print('Load matrices...')
        data = np.load('u3d.npz')
        u3d = data['u3d']
        v3d = data['v3d']
        depth = data['depth']

    [Nz] = np.shape(depth)
    Z = -depth

    #
    # ----------------------------------------------------
    #  Vertical interpolation
    # ----------------------------------------------------
    #

    print('Vertical interpolation')

    #
    # Add a layer below the bottom and above the surface to avoid vertical extrapolations
    # and flip the matrices upside down (Z[Nz]=surface)
    #
    #
    Z = np.flipud(np.concatenate(([100.], Z, [-10000.])))

    [Nz] = np.shape(Z)

    u3d = np.flipud(vgrd.add2layers(u3d))
    v3d = np.flipud(vgrd.add2layers(v3d))

    #
    # Do the vertical interpolations
    #

    uout = vgrd.ztosigma(u3d, Z, rho2u_3d(z_rho))
    vout = vgrd.ztosigma(v3d, Z, rho2v_3d(z_rho))

    return uout, vout


#
#
#
# #####################################################################
# ##### END FUNCTION INTERP3D_UV ######################################
# #####################################################################
#
#


#
#
# #####################################################################
# #### FUNCTION CREATE_INIFILE ########################################
# #####################################################################
#
#

def create_inifile(ininame, grdname, title, theta_s, theta_b, hc, N, ini_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    #
    #
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + ininame)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #

    type = 'INITIAL file'
    history = 'CROCO'
    if os.path.exists(ininame):
        os.remove(ininame)

    print('Create: ' + ininame)
    # nc = netcdf(ininame, 'w', format='NETCDF3_CLASSIC')
    nc = netcdf(ininame, 'w', format='NETCDF4_CLASSIC')

    #
    # Create global attributes
    #

    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = ininame
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_ini_time = nc.createDimension('ini_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_scrum_time = nc.createVariable('scrum_time', np.float64, ('ini_time',))
    nc_scrum_time.long_name = 'ini_time since initialization'
    nc_scrum_time.units = 'second'
    #
    nc_ocean_time = nc.createVariable('ocean_time', np.float64, ('ini_time',))
    nc_ocean_time.long_name = 'ini_time since initialization'
    nc_ocean_time.units = 'second'
    #
    nc_u = nc.createVariable('u', np.float64, ('ini_time', 's_rho', 'eta_u', 'xi_u',))
    nc_u.long_name = 'u-momentum component'
    nc_u.units = 'meter second-1'
    #
    nc_v = nc.createVariable('v', np.float64, ('ini_time', 's_rho', 'eta_v', 'xi_v',))
    nc_v.long_name = 'v-momentum component'
    nc_v.units = 'meter second-1'
    #
    nc_ubar = nc.createVariable('ubar', np.float64, ('ini_time', 'eta_u', 'xi_u',))
    nc_ubar.long_name = 'vertically integrated u-momentum component'
    nc_ubar.units = 'meter second-1'
    #
    nc_vbar = nc.createVariable('vbar', np.float64, ('ini_time', 'eta_v', 'xi_v',))
    nc_vbar.long_name = 'vertically integrated v-momentum component'
    nc_vbar.units = 'meter second-1'
    #
    nc_zeta = nc.createVariable('zeta', np.float64, ('ini_time', 'eta_rho', 'xi_rho',))
    nc_zeta.long_name = 'free-surface'
    nc_zeta.units = 'meter'
    #
    nc_temp = nc.createVariable('temp', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_temp.long_name = 'potential temperature'
    nc_temp.units = 'Celsius'
    #
    nc_salt = nc.createVariable('salt', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_salt.long_name = 'salinity'
    nc_salt.units = 'PSU'
    #
    nc_no3 = nc.createVariable('NO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_no3.long_name = 'nitrate'
    nc_no3.units = 'mmol/m3'
    #
    nc_po4 = nc.createVariable('PO4', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_po4.long_name = 'orthophosphate'
    nc_po4.units = 'mmol/m3'
    #
    nc_si = nc.createVariable('Si', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_si.long_name = 'silicate'
    nc_si.units = 'mmol/m3'
    #
    nc_o2 = nc.createVariable('O2', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_o2.long_name = 'dissolved oxygen'
    nc_o2.units = 'mmol/m3'
    #
    nc_dic = nc.createVariable('DIC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dic.long_name = 'dissolved inorganic carbon'
    nc_dic.units = 'mmol/m3'
    #
    nc_talk = nc.createVariable('TALK', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_talk.long_name = 'total alkalinity'
    nc_talk.units = 'mmol/m3'
    #
    # nc_nh4 = nc.createVariable('NH4', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nh4.long_name = 'ammonium'
    # nc_nh4.units = 'mmol/m3'
    #
    nc_fe = nc.createVariable('FER', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_fe.long_name = 'iron'
    nc_fe.units = 'mmol/m3'
    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = ini_time
    nc_tend[:] = ini_time
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_scrum_time[0] = ini_time * 24. * 3600.
    nc_ocean_time[0] = ini_time * 24. * 3600.
    nc_u[:] = 0
    nc_v[:] = 0
    nc_zeta[:] = 0
    nc_ubar[:] = 0
    nc_vbar[:] = 0
    nc_temp[:] = 0
    nc_salt[:] = 0
    nc_no3[:] = 0
    nc_po4[:] = 0
    nc_si[:] = 0
    nc_o2[:] = 0
    nc_dic[:] = 0
    nc_talk[:] = 0
    # nc_nh4[:] = 0
    nc_fe[:] = 0
    nc.close()
    #
    return


#
# #####################################################################
# ##### END FUNCTION CREATE_INIFILE ###################################
# #####################################################################
#


def create_ini_PHY_BGC(ininame, grdname, title, theta_s, theta_b, hc, N, ini_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    #
    #
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + ininame)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #

    type = 'INITIAL file'
    history = 'CROCO'
    if os.path.exists(ininame):
        os.remove(ininame)

    print('Create: ' + ininame)
    nc = netcdf(ininame, 'w', format='NETCDF3_CLASSIC')

    #
    # Create global attributes
    #

    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = ininame
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_ini_time = nc.createDimension('ini_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_scrum_time = nc.createVariable('scrum_time', np.float64, ('ini_time',))
    nc_scrum_time.long_name = 'ini_time since initialization'
    nc_scrum_time.units = 'second'
    #
    nc_ocean_time = nc.createVariable('ocean_time', np.float64, ('ini_time',))
    nc_ocean_time.long_name = 'ini_time since initialization'
    nc_ocean_time.units = 'second'
    #
    nc_u = nc.createVariable('u', np.float64, ('ini_time', 's_rho', 'eta_u', 'xi_u',))
    nc_u.long_name = 'u-momentum component'
    nc_u.units = 'meter second-1'
    #
    nc_v = nc.createVariable('v', np.float64, ('ini_time', 's_rho', 'eta_v', 'xi_v',))
    nc_v.long_name = 'v-momentum component'
    nc_v.units = 'meter second-1'
    #
    nc_ubar = nc.createVariable('ubar', np.float64, ('ini_time', 'eta_u', 'xi_u',))
    nc_ubar.long_name = 'vertically integrated u-momentum component'
    nc_ubar.units = 'meter second-1'
    #
    nc_vbar = nc.createVariable('vbar', np.float64, ('ini_time', 'eta_v', 'xi_v',))
    nc_vbar.long_name = 'vertically integrated v-momentum component'
    nc_vbar.units = 'meter second-1'
    #
    nc_zeta = nc.createVariable('zeta', np.float64, ('ini_time', 'eta_rho', 'xi_rho',))
    nc_zeta.long_name = 'free-surface'
    nc_zeta.units = 'meter'
    #
    nc_temp = nc.createVariable('temp', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_temp.long_name = 'potential temperature'
    nc_temp.units = 'Celsius'
    #
    nc_salt = nc.createVariable('salt', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_salt.long_name = 'salinity'
    nc_salt.units = 'PSU'
    #
    nc_no3 = nc.createVariable('NO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_no3.long_name = 'nitrate'
    nc_no3.units = 'mmol/m3'
    #
    nc_po4 = nc.createVariable('PO4', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_po4.long_name = 'orthophosphate'
    nc_po4.units = 'mmol/m3'
    #
    nc_si = nc.createVariable('Si', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_si.long_name = 'silicate'
    nc_si.units = 'mmol/m3'
    #
    nc_o2 = nc.createVariable('O2', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_o2.long_name = 'dissolved oxygen'
    nc_o2.units = 'mmol/m3'
    #
    nc_dic = nc.createVariable('DIC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dic.long_name = 'dissolved inorganic carbon'
    nc_dic.units = 'mmol/m3'
    #
    nc_nh4 = nc.createVariable('NH4', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nh4.long_name = 'ammonium'
    nc_nh4.units = 'mmol/m3'
    #
    nc_fe = nc.createVariable('FER', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_fe.long_name = 'iron'
    nc_fe.units = 'mmol/m3'

    nc_talk = nc.createVariable('TALK', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_talk.long_name = 'alkalinity'
    nc_talk.units = 'mmol/m3'

    # # nc_caco3 = nc.createVariable('CaCO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_caco3 = nc.createVariable('CACO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_caco3.long_name = 'calcite'
    # nc_caco3.units = 'mmol/m3'
    #
    # nc_poc = nc.createVariable('POC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_poc.long_name = 'particulate organic carbon'
    # nc_poc.units = 'mmol/m3'
    #
    # # nc_phy = nc.createVariable('PHY', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy = nc.createVariable('NANO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy.long_name = 'nanophytoplankton'
    # nc_phy.units = 'mmol/m3'
    #
    # nc_zoo = nc.createVariable('ZOO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo.long_name = 'microzooplankton'
    # nc_zoo.units = 'mmol/m3'
    #
    # nc_doc = nc.createVariable('DOC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_doc.long_name = 'dissolved organic carbon'
    # nc_doc.units = 'mmol/m3'
    #
    # # nc_phy2 = nc.createVariable('DIA', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy2 = nc.createVariable('DIA', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy2.long_name = 'diatom'
    # nc_phy2.units = 'mmol/m3'
    #
    # # nc_zoo2 = nc.createVariable('MESO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo2 = nc.createVariable('MESO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo2.long_name = 'mesozooplankton'
    # nc_zoo2.units = 'mmol/m3'
    #
    # # nc_bsi = nc.createVariable('BSi', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bsi = nc.createVariable('BSI', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bsi.long_name = 'biogenic silica'
    # nc_bsi.units = 'mmol/m3'
    #
    # # nc_bfe = nc.createVariable('BFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bfe = nc.createVariable('BFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bfe.long_name = 'big particle iron'
    # nc_bfe.units = 'mmol/m3'
    #
    # nc_goc = nc.createVariable('GOC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_goc.long_name = 'big particulate organic carbon'
    # nc_goc.units = 'mmol/m3'
    #
    # # nc_sfe = nc.createVariable('SFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_sfe = nc.createVariable('SFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_sfe.long_name = 'iron in the small particles'
    # nc_sfe.units = 'mmol/m3'
    #
    # # nc_dfe = nc.createVariable('DFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dfe = nc.createVariable('DFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dfe.long_name = 'iron content of the diatoms'
    # nc_dfe.units = 'mmol/m3'
    #
    # # nc_dsi = nc.createVariable('DSi', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dsi = nc.createVariable('DSI', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dsi.long_name = 'silicon content of the Diatoms'
    # nc_dsi.units = 'mmol/m3'
    #
    # # nc_nfe = nc.createVariable('NFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nfe = nc.createVariable('NFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nfe.long_name = 'nanophytoplankton iron'
    # nc_nfe.units = 'mmol/m3'
    #
    # nc_nchl = nc.createVariable('NCHL', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nchl.long_name = 'nanophytoplankton chlorophyll'
    # nc_nchl.units = 'mmol/m3'
    #
    # nc_dchl = nc.createVariable('DCHL', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dchl.long_name = 'diatom chlorophyll'
    # nc_dchl.units = 'mmol/m3'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = ini_time
    nc_tend[:] = ini_time
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_scrum_time[0] = ini_time * 24. * 3600.
    nc_ocean_time[0] = ini_time * 24. * 3600.
    nc_u[:] = 0
    nc_v[:] = 0
    nc_ubar[:] = 0
    nc_vbar[:] = 0
    nc_zeta[:] = 0
    nc_temp[:] = 0
    nc_salt[:] = 0
    nc_no3[:] = 0
    nc_po4[:] = 0
    nc_si[:] = 0
    nc_o2[:] = 0
    nc_dic[:] = 0
    nc_nh4[:] = 0
    nc_fe[:] = 0
    nc_talk[:] = 0
    # nc_caco3[:] = 0
    # nc_poc[:] = 0
    # nc_phy[:] = 0
    # nc_zoo[:] = 0
    # nc_doc[:] = 0
    # nc_phy2[:] = 0
    # nc_zoo2[:] = 0
    # nc_bsi[:] = 0
    # nc_bfe[:] = 0
    # nc_goc[:] = 0
    # nc_sfe[:] = 0
    # nc_dfe[:] = 0
    # nc_dsi[:] = 0
    # nc_nfe[:] = 0
    # nc_nchl[:] = 0
    # nc_dchl[:] = 0
    nc.close()
    #
    return

def create_ini_BFM(ininame, grdname, title, theta_s, theta_b, hc, N, ini_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    #
    #
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + ininame)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #

    type = 'INITIAL file'
    history = 'CROCO'
    if os.path.exists(ininame):
        os.remove(ininame)

    print('Create: ' + ininame)
    nc = netcdf(ininame, 'w', format='NETCDF3_CLASSIC')

    #
    # Create global attributes
    #

    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = ininame
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_ini_time = nc.createDimension('ini_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_scrum_time = nc.createVariable('scrum_time', np.float64, ('ini_time',))
    nc_scrum_time.long_name = 'ini_time since initialization'
    nc_scrum_time.units = 'second'
    #
    nc_ocean_time = nc.createVariable('ocean_time', np.float64, ('ini_time',))
    nc_ocean_time.long_name = 'ini_time since initialization'
    nc_ocean_time.units = 'second'
    #
    nc_u = nc.createVariable('u', np.float64, ('ini_time', 's_rho', 'eta_u', 'xi_u',))
    nc_u.long_name = 'u-momentum component'
    nc_u.units = 'meter second-1'
    #
    nc_v = nc.createVariable('v', np.float64, ('ini_time', 's_rho', 'eta_v', 'xi_v',))
    nc_v.long_name = 'v-momentum component'
    nc_v.units = 'meter second-1'
    #
    nc_ubar = nc.createVariable('ubar', np.float64, ('ini_time', 'eta_u', 'xi_u',))
    nc_ubar.long_name = 'vertically integrated u-momentum component'
    nc_ubar.units = 'meter second-1'
    #
    nc_vbar = nc.createVariable('vbar', np.float64, ('ini_time', 'eta_v', 'xi_v',))
    nc_vbar.long_name = 'vertically integrated v-momentum component'
    nc_vbar.units = 'meter second-1'
    #
    nc_zeta = nc.createVariable('zeta', np.float64, ('ini_time', 'eta_rho', 'xi_rho',))
    nc_zeta.long_name = 'free-surface'
    nc_zeta.units = 'meter'
    #
    nc_temp = nc.createVariable('temp', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_temp.long_name = 'potential temperature'
    nc_temp.units = 'Celsius'
    #
    nc_salt = nc.createVariable('salt', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_salt.long_name = 'salinity'
    nc_salt.units = 'PSU'
    #
    nc_no3 = nc.createVariable('N3n', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_no3.long_name = 'nitrate'
    nc_no3.units = 'mmol/m3'
    #
    nc_po4 = nc.createVariable('N1p', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_po4.long_name = 'orthophosphate'
    nc_po4.units = 'mmol/m3'
    #
    nc_si = nc.createVariable('N5s', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_si.long_name = 'silicate'
    nc_si.units = 'mmol/m3'
    #
    nc_o2 = nc.createVariable('O2o', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_o2.long_name = 'dissolved oxygen'
    nc_o2.units = 'mmol/m3'
    #
    nc_dic = nc.createVariable('O3c', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dic.long_name = 'dissolved inorganic carbon'
    nc_dic.units = 'mmol/m3'
    #
    nc_nh4 = nc.createVariable('N4n', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nh4.long_name = 'ammonium'
    nc_nh4.units = 'mmol/m3'
    #
    nc_fe = nc.createVariable('N7f', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_fe.long_name = 'iron'
    nc_fe.units = 'mmol/m3'

    nc_talk = nc.createVariable('O3h', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_talk.long_name = 'alkalinity'
    nc_talk.units = 'mmol/m3'

    # # nc_caco3 = nc.createVariable('CaCO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_caco3 = nc.createVariable('CACO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_caco3.long_name = 'calcite'
    # nc_caco3.units = 'mmol/m3'
    #
    # nc_poc = nc.createVariable('POC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_poc.long_name = 'particulate organic carbon'
    # nc_poc.units = 'mmol/m3'
    #
    # # nc_phy = nc.createVariable('PHY', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy = nc.createVariable('NANO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy.long_name = 'nanophytoplankton'
    # nc_phy.units = 'mmol/m3'
    #
    # nc_zoo = nc.createVariable('ZOO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo.long_name = 'microzooplankton'
    # nc_zoo.units = 'mmol/m3'
    #
    # nc_doc = nc.createVariable('DOC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_doc.long_name = 'dissolved organic carbon'
    # nc_doc.units = 'mmol/m3'
    #
    # # nc_phy2 = nc.createVariable('DIA', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy2 = nc.createVariable('DIA', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy2.long_name = 'diatom'
    # nc_phy2.units = 'mmol/m3'
    #
    # # nc_zoo2 = nc.createVariable('MESO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo2 = nc.createVariable('MESO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo2.long_name = 'mesozooplankton'
    # nc_zoo2.units = 'mmol/m3'
    #
    # # nc_bsi = nc.createVariable('BSi', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bsi = nc.createVariable('BSI', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bsi.long_name = 'biogenic silica'
    # nc_bsi.units = 'mmol/m3'
    #
    # # nc_bfe = nc.createVariable('BFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bfe = nc.createVariable('BFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bfe.long_name = 'big particle iron'
    # nc_bfe.units = 'mmol/m3'
    #
    # nc_goc = nc.createVariable('GOC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_goc.long_name = 'big particulate organic carbon'
    # nc_goc.units = 'mmol/m3'
    #
    # # nc_sfe = nc.createVariable('SFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_sfe = nc.createVariable('SFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_sfe.long_name = 'iron in the small particles'
    # nc_sfe.units = 'mmol/m3'
    #
    # # nc_dfe = nc.createVariable('DFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dfe = nc.createVariable('DFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dfe.long_name = 'iron content of the diatoms'
    # nc_dfe.units = 'mmol/m3'
    #
    # # nc_dsi = nc.createVariable('DSi', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dsi = nc.createVariable('DSI', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dsi.long_name = 'silicon content of the Diatoms'
    # nc_dsi.units = 'mmol/m3'
    #
    # # nc_nfe = nc.createVariable('NFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nfe = nc.createVariable('NFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nfe.long_name = 'nanophytoplankton iron'
    # nc_nfe.units = 'mmol/m3'
    #
    # nc_nchl = nc.createVariable('NCHL', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nchl.long_name = 'nanophytoplankton chlorophyll'
    # nc_nchl.units = 'mmol/m3'
    #
    # nc_dchl = nc.createVariable('DCHL', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dchl.long_name = 'diatom chlorophyll'
    # nc_dchl.units = 'mmol/m3'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = ini_time
    nc_tend[:] = ini_time
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_scrum_time[0] = ini_time * 24. * 3600.
    nc_ocean_time[0] = ini_time * 24. * 3600.
    nc_u[:] = 0
    nc_v[:] = 0
    nc_ubar[:] = 0
    nc_vbar[:] = 0
    nc_zeta[:] = 0
    nc_temp[:] = 0
    nc_salt[:] = 0
    nc_no3[:] = 0
    nc_po4[:] = 0
    nc_si[:] = 0
    nc_o2[:] = 0
    nc_dic[:] = 0
    nc_nh4[:] = 0
    nc_fe[:] = 0
    nc_talk[:] = 0
    # nc_caco3[:] = 0
    # nc_poc[:] = 0
    # nc_phy[:] = 0
    # nc_zoo[:] = 0
    # nc_doc[:] = 0
    # nc_phy2[:] = 0
    # nc_zoo2[:] = 0
    # nc_bsi[:] = 0
    # nc_bfe[:] = 0
    # nc_goc[:] = 0
    # nc_sfe[:] = 0
    # nc_dfe[:] = 0
    # nc_dsi[:] = 0
    # nc_nfe[:] = 0
    # nc_nchl[:] = 0
    # nc_dchl[:] = 0
    nc.close()
    #
    return


def create_ini_PISCES_NORESM(ininame, grdname, title, theta_s, theta_b, hc, N, ini_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    #
    #
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + ininame)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #

    type = 'INITIAL file'
    history = 'CROCO'
    if os.path.exists(ininame):
        os.remove(ininame)

    print('Create: ' + ininame)
    nc = netcdf(ininame, 'w', format='NETCDF3_CLASSIC')

    #
    # Create global attributes
    #

    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = ininame
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_ini_time = nc.createDimension('ini_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_scrum_time = nc.createVariable('scrum_time', np.float64, ('ini_time',))
    nc_scrum_time.long_name = 'ini_time since initialization'
    nc_scrum_time.units = 'second'
    #
    nc_ocean_time = nc.createVariable('ocean_time', np.float64, ('ini_time',))
    nc_ocean_time.long_name = 'ini_time since initialization'
    nc_ocean_time.units = 'second'
    #
    nc_u = nc.createVariable('u', np.float64, ('ini_time', 's_rho', 'eta_u', 'xi_u',))
    nc_u.long_name = 'u-momentum component'
    nc_u.units = 'meter second-1'
    #
    nc_v = nc.createVariable('v', np.float64, ('ini_time', 's_rho', 'eta_v', 'xi_v',))
    nc_v.long_name = 'v-momentum component'
    nc_v.units = 'meter second-1'
    #
    nc_ubar = nc.createVariable('ubar', np.float64, ('ini_time', 'eta_u', 'xi_u',))
    nc_ubar.long_name = 'vertically integrated u-momentum component'
    nc_ubar.units = 'meter second-1'
    #
    nc_vbar = nc.createVariable('vbar', np.float64, ('ini_time', 'eta_v', 'xi_v',))
    nc_vbar.long_name = 'vertically integrated v-momentum component'
    nc_vbar.units = 'meter second-1'
    #
    nc_zeta = nc.createVariable('zeta', np.float64, ('ini_time', 'eta_rho', 'xi_rho',))
    nc_zeta.long_name = 'free-surface'
    nc_zeta.units = 'meter'
    #
    nc_temp = nc.createVariable('temp', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_temp.long_name = 'potential temperature'
    nc_temp.units = 'Celsius'
    #
    nc_salt = nc.createVariable('salt', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_salt.long_name = 'salinity'
    nc_salt.units = 'PSU'
    #
    nc_no3 = nc.createVariable('NO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_no3.long_name = 'nitrate'
    nc_no3.units = 'mmol/m3'
    #
    nc_po4 = nc.createVariable('PO4', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_po4.long_name = 'orthophosphate'
    nc_po4.units = 'mmol/m3'
    #
    nc_si = nc.createVariable('Si', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_si.long_name = 'silicate'
    nc_si.units = 'mmol/m3'
    #
    nc_o2 = nc.createVariable('O2', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_o2.long_name = 'dissolved oxygen'
    nc_o2.units = 'mmol/m3'
    #
    nc_dic = nc.createVariable('DIC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dic.long_name = 'dissolved inorganic carbon'
    nc_dic.units = 'mmol/m3'
    #
    nc_nh4 = nc.createVariable('NH4', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nh4.long_name = 'ammonium'
    nc_nh4.units = 'mmol/m3'
    #
    # nc_fe = nc.createVariable('Fer', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_fe = nc.createVariable('FER', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_fe.long_name = 'iron'
    nc_fe.units = 'mmol/m3'

    # nc_talk = nc.createVariable('TALK', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_talk = nc.createVariable('TALK', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_talk.long_name = 'alkalinity'
    nc_talk.units = 'mmol/m3'

    # nc_caco3 = nc.createVariable('CaCO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_caco3 = nc.createVariable('CACO3', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_caco3.long_name = 'calcite'
    nc_caco3.units = 'mmol/m3'

    nc_poc = nc.createVariable('POC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_poc.long_name = 'particulate organic carbon'
    nc_poc.units = 'mmol/m3'

    # nc_phy = nc.createVariable('PHY', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_phy = nc.createVariable('NANO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_phy.long_name = 'nanophytoplankton'
    nc_phy.units = 'mmol/m3'

    nc_zoo = nc.createVariable('ZOO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_zoo.long_name = 'microzooplankton'
    nc_zoo.units = 'mmol/m3'

    nc_doc = nc.createVariable('DOC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_doc.long_name = 'dissolved organic carbon'
    nc_doc.units = 'mmol/m3'

    # nc_phy2 = nc.createVariable('DIA', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_phy2 = nc.createVariable('DIA', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_phy2.long_name = 'diatom'
    nc_phy2.units = 'mmol/m3'

    # nc_zoo2 = nc.createVariable('MESO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_zoo2 = nc.createVariable('MESO', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_zoo2.long_name = 'mesozooplankton'
    nc_zoo2.units = 'mmol/m3'

    # nc_bsi = nc.createVariable('BSi', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_bsi = nc.createVariable('BSI', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_bsi.long_name = 'biogenic silica'
    nc_bsi.units = 'mmol/m3'

    # nc_bfe = nc.createVariable('BFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_bfe = nc.createVariable('BFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_bfe.long_name = 'big particle iron'
    nc_bfe.units = 'mmol/m3'

    nc_goc = nc.createVariable('GOC', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_goc.long_name = 'big particulate organic carbon'
    nc_goc.units = 'mmol/m3'

    # nc_sfe = nc.createVariable('SFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_sfe = nc.createVariable('SFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_sfe.long_name = 'iron in the small particles'
    nc_sfe.units = 'mmol/m3'

    # nc_dfe = nc.createVariable('DFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dfe = nc.createVariable('DFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dfe.long_name = 'iron content of the diatoms'
    nc_dfe.units = 'mmol/m3'

    # nc_dsi = nc.createVariable('DSi', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dsi = nc.createVariable('DSI', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dsi.long_name = 'silicon content of the Diatoms'
    nc_dsi.units = 'mmol/m3'

    # nc_nfe = nc.createVariable('NFe', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nfe = nc.createVariable('NFE', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nfe.long_name = 'nanophytoplankton iron'
    nc_nfe.units = 'mmol/m3'

    nc_nchl = nc.createVariable('NCHL', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nchl.long_name = 'nanophytoplankton chlorophyll'
    nc_nchl.units = 'mmol/m3'

    nc_dchl = nc.createVariable('DCHL', np.float64, ('ini_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dchl.long_name = 'diatom chlorophyll'
    nc_dchl.units = 'mmol/m3'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = ini_time
    nc_tend[:] = ini_time
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_scrum_time[0] = ini_time * 24. * 3600.
    nc_ocean_time[0] = ini_time * 24. * 3600.
    nc_u[:] = 0
    nc_v[:] = 0
    nc_ubar[:] = 0
    nc_vbar[:] = 0
    nc_zeta[:] = 0
    nc_temp[:] = 0
    nc_salt[:] = 0
    nc_no3[:] = 0
    nc_po4[:] = 0
    nc_si[:] = 0
    nc_o2[:] = 0
    nc_dic[:] = 0
    nc_nh4[:] = 0
    nc_fe[:] = 0
    nc_talk[:] = 0
    nc_caco3[:] = 0
    nc_poc[:] = 0
    nc_phy[:] = 0
    nc_zoo[:] = 0
    nc_doc[:] = 0
    nc_phy2[:] = 0
    nc_zoo2[:] = 0
    nc_bsi[:] = 0
    nc_bfe[:] = 0
    nc_goc[:] = 0
    nc_sfe[:] = 0
    nc_dfe[:] = 0
    nc_dsi[:] = 0
    nc_nfe[:] = 0
    nc_nchl[:] = 0
    nc_dchl[:] = 0
    nc.close()
    #
    return


def create_ini_FENNEL_ROMS(ininame, grdname, title, theta_s, theta_b, hc, N, ini_time, vtransform, vstretching):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    #
    #
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + ininame)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #

    type = 'INITIAL file'
    history = 'ROMS Fennel'
    if os.path.exists(ininame):
        os.remove(ininame)

    print('Create: ' + ininame)
    nc = netcdf(ininame, 'w', format='NETCDF4_CLASSIC')

    #
    # Create global attributes
    #

    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = ininame
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    # nc_dim_ini_time = nc.createDimension('ini_time', 0)
    nc_dim_ini_time = nc.createDimension('ocean_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_scrum_time = nc.createVariable('scrum_time', np.float64, ('ocean_time',))
    nc_scrum_time.long_name = 'ocean_time since initialization'
    nc_scrum_time.units = 'second'
    #
    nc_ocean_time = nc.createVariable('ocean_time', np.float64, ('ocean_time',))
    nc_ocean_time.long_name = 'ocean_time since initialization'
    nc_ocean_time.units = 'second'
    #
    nc_u = nc.createVariable('u', np.float64, ('ocean_time', 's_rho', 'eta_u', 'xi_u',))
    nc_u.long_name = 'u-momentum component'
    nc_u.units = 'meter second-1'
    #
    nc_v = nc.createVariable('v', np.float64, ('ocean_time', 's_rho', 'eta_v', 'xi_v',))
    nc_v.long_name = 'v-momentum component'
    nc_v.units = 'meter second-1'
    #
    nc_ubar = nc.createVariable('ubar', np.float64, ('ocean_time', 'eta_u', 'xi_u',))
    nc_ubar.long_name = 'vertically integrated u-momentum component'
    nc_ubar.units = 'meter second-1'
    #
    nc_vbar = nc.createVariable('vbar', np.float64, ('ocean_time', 'eta_v', 'xi_v',))
    nc_vbar.long_name = 'vertically integrated v-momentum component'
    nc_vbar.units = 'meter second-1'
    #
    nc_zeta = nc.createVariable('zeta', np.float64, ('ocean_time', 'eta_rho', 'xi_rho',))
    nc_zeta.long_name = 'free-surface'
    nc_zeta.units = 'meter'
    #
    nc_temp = nc.createVariable('temp', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_temp.long_name = 'potential temperature'
    nc_temp.units = 'Celsius'
    #
    nc_salt = nc.createVariable('salt', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_salt.long_name = 'salinity'
    nc_salt.units = 'PSU'
    #
    nc_no3 = nc.createVariable('NO3', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_no3.long_name = 'nitrate'
    nc_no3.units = 'mmol/m3'
    #
    nc_o2 = nc.createVariable('oxygen', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_o2.long_name = 'dissolved oxygen'
    nc_o2.units = 'mmol/m3'
    #
    nc_dic = nc.createVariable('TIC', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_dic.long_name = 'dissolved inorganic carbon'
    nc_dic.units = 'mmol C/m3'
    #
    nc_nh4 = nc.createVariable('NH4', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_nh4.long_name = 'ammonium'
    nc_nh4.units = 'mmol/m3'

    nc_talk = nc.createVariable('alkalinity', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_talk.long_name = 'alkalinity'
    nc_talk.units = 'milliequivalents/m3'

    nc_phy = nc.createVariable('phytoplankton', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_phy.long_name = 'phytoplankton'
    nc_phy.units = 'mmol N/m3'

    nc_zoo = nc.createVariable('zooplankton', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_zoo.long_name = 'zooplankton'
    nc_zoo.units = 'mmol N/m3'

    nc_chl = nc.createVariable('chlorophyll', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_chl.long_name = 'chlorophyll'
    nc_chl.units = 'mg Chl/m3'
    #
    nc_ldn = nc.createVariable('LdetritusN', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_ldn.long_name = 'large fraction nitrogen detritus concentration'
    nc_ldn.units = 'mmol N/m3'
    #
    nc_sdn = nc.createVariable('SdetritusN', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_sdn.long_name = 'small fraction nitrogen detritus concentration'
    nc_sdn.units = 'mmol N/m3'
    #
    nc_ldc = nc.createVariable('LdetritusC', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_ldc.long_name = 'large fraction carbon detritus concentration'
    nc_ldc.units = 'mmol C/m3'
    #
    nc_sdc = nc.createVariable('SdetritusC', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',))
    nc_sdc.long_name = 'small fraction carbon detritus concentration'
    nc_sdc.units = 'mmol C/m3'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = vstretching
    nc_tstart[:] = ini_time
    nc_tend[:] = ini_time
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_scrum_time[0] = ini_time * 24. * 3600.
    nc_ocean_time[0] = ini_time * 24. * 3600.
    nc_u[:] = 0
    nc_v[:] = 0
    nc_ubar[:] = 0
    nc_vbar[:] = 0
    nc_zeta[:] = 0
    nc_temp[:] = 0
    nc_salt[:] = 0
    nc_no3[:] = 0
    nc_o2[:] = 0
    nc_dic[:] = 0
    nc_nh4[:] = 0
    nc_talk[:] = 0
    nc_phy[:] = 0
    nc_zoo[:] = 0.1
    nc_ldn[:] = 0.1
    nc_sdn[:] = 0.1
    nc_ldc[:] = 0.6625
    nc_sdc[:] = 0.6625

    nc.close()
    #
    return


def create_VAL_PISCES_NORESM(VALname, grdname, title, theta_s, theta_b, hc, N, VAL_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    #
    #
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + VALname)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #

    type = 'VALIDATION file'
    history = 'CROCO'
    if os.path.exists(VALname):
        os.remove(VALname)

    print('Create: ' + VALname)
    nc = netcdf(VALname, 'w', format='NETCDF3_CLASSIC')

    #
    # Create global attributes
    #

    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = VALname
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_VAL_time = nc.createDimension('VAL_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_VAL_time = nc.createVariable('VAL_time', np.float64, ('VAL_time',))
    nc_VAL_time.long_name = 'VAL_time since initialization'
    nc_VAL_time.units = 'days'

    # nc_u = nc.createVariable('u', np.float64, ('VAL_time', 's_rho', 'eta_u', 'xi_u',))
    # nc_u.long_name = 'u-momentum component'
    # nc_u.units = 'meter second-1'
    # #
    # nc_v = nc.createVariable('v', np.float64, ('VAL_time', 's_rho', 'eta_v', 'xi_v',))
    # nc_v.long_name = 'v-momentum component'
    # nc_v.units = 'meter second-1'
    # #
    # nc_ubar = nc.createVariable('ubar', np.float64, ('VAL_time', 'eta_u', 'xi_u',))
    # nc_ubar.long_name = 'vertically integrated u-momentum component'
    # nc_ubar.units = 'meter second-1'
    # #
    # nc_vbar = nc.createVariable('vbar', np.float64, ('VAL_time', 'eta_v', 'xi_v',))
    # nc_vbar.long_name = 'vertically integrated v-momentum component'
    # nc_vbar.units = 'meter second-1'
    #
    # nc_zeta = nc.createVariable('zeta', np.float64, ('VAL_time', 'eta_rho', 'xi_rho',))
    # nc_zeta.long_name = 'free-surface'
    # nc_zeta.units = 'meter'
    #
    # nc_temp = nc.createVariable('temp', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_temp.long_name = 'potential temperature'
    # nc_temp.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps', np.float64, ('VAL_time', 'eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_model-sat', np.float64, ('VAL_time', 'eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_t_RMSD', np.float64, ('VAL_time',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_t_ME', np.float64, ('VAL_time',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_t_CC', np.float64, ('VAL_time',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_t_bias', np.float64, ('VAL_time',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_s_RMSD', np.float64, ('eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_s_ME', np.float64, ('eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_s_CC', np.float64, ('eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_temps = nc.createVariable('temps_s_bias', np.float64, ('eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    #
    nc_salts = nc.createVariable('salts', np.float64, ('VAL_time', 'eta_rho', 'xi_rho',))
    nc_salts.long_name = 'sea surface salinity'
    nc_salts.units = 'PSU'
    #
    nc_saltb = nc.createVariable('saltb', np.float64, ('VAL_time', 'eta_rho', 'xi_rho',))
    nc_saltb.long_name = 'bottom salinity'
    nc_saltb.units = 'PSU'
    #
    nc_salts = nc.createVariable('salt_bottom-surf', np.float64, ('VAL_time', 'eta_rho', 'xi_rho',))
    nc_salts.long_name = 'sea surface salinity'
    nc_salts.units = 'PSU'
    #
    # nc_no3 = nc.createVariable('NO3', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_no3.long_name = 'nitrate'
    # nc_no3.units = 'mmol/m3'
    # #
    # nc_po4 = nc.createVariable('PO4', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_po4.long_name = 'orthophosphate'
    # nc_po4.units = 'mmol/m3'
    # #
    # nc_si = nc.createVariable('Si', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_si.long_name = 'silicate'
    # nc_si.units = 'mmol/m3'
    # #
    # nc_o2 = nc.createVariable('O2', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_o2.long_name = 'dissolved oxygen'
    # nc_o2.units = 'mmol/m3'
    # #
    # nc_dic = nc.createVariable('DIC', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dic.long_name = 'dissolved inorganic carbon'
    # nc_dic.units = 'mmol/m3'
    # #
    # nc_nh4 = nc.createVariable('NH4', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nh4.long_name = 'ammonium'
    # nc_nh4.units = 'mmol/m3'
    # #
    # nc_fe = nc.createVariable('Fer', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_fe.long_name = 'iron'
    # nc_fe.units = 'mmol/m3'
    #
    # nc_talk = nc.createVariable('TALK', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_talk.long_name = 'alkalinity'
    # nc_talk.units = 'mmol/m3'
    #
    # nc_caco3 = nc.createVariable('CACO3', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_caco3.long_name = 'calcite'
    # nc_caco3.units = 'mmol/m3'
    #
    # nc_poc = nc.createVariable('POC', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_poc.long_name = 'particulate organic carbon'
    # nc_poc.units = 'mmol/m3'
    #
    # nc_phy = nc.createVariable('PHY', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy.long_name = 'nanophytoplankton'
    # nc_phy.units = 'mmol/m3'
    #
    # nc_zoo = nc.createVariable('ZOO', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo.long_name = 'microzooplankton'
    # nc_zoo.units = 'mmol/m3'
    #
    # nc_doc = nc.createVariable('DOC', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_doc.long_name = 'dissolved organic carbon'
    # nc_doc.units = 'mmol/m3'
    #
    # nc_phy2 = nc.createVariable('DIA', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_phy2.long_name = 'diatom'
    # nc_phy2.units = 'mmol/m3'
    #
    # nc_zoo2 = nc.createVariable('MESO', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_zoo2.long_name = 'mesozooplankton'
    # nc_zoo2.units = 'mmol/m3'
    #
    # nc_bsi = nc.createVariable('BSi', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bsi.long_name = 'biogenic silica'
    # nc_bsi.units = 'mmol/m3'
    #
    # nc_bfe = nc.createVariable('BFe', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_bfe.long_name = 'big particle iron'
    # nc_bfe.units = 'mmol/m3'
    #
    # nc_goc = nc.createVariable('GOC', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_goc.long_name = 'big particulate organic carbon'
    # nc_goc.units = 'mmol/m3'
    #
    # nc_sfe = nc.createVariable('SFe', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_sfe.long_name = 'iron in the small particles'
    # nc_sfe.units = 'mmol/m3'
    #
    # nc_dfe = nc.createVariable('DFe', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dfe.long_name = 'iron content of the diatoms'
    # nc_dfe.units = 'mmol/m3'
    #
    # nc_dsi = nc.createVariable('DSi', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dsi.long_name = 'silicon content of the Diatoms'
    # nc_dsi.units = 'mmol/m3'
    #
    # nc_nfe = nc.createVariable('NFe', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nfe.long_name = 'nanophytoplankton iron'
    # nc_nfe.units = 'mmol/m3'
    #
    # nc_nchl = nc.createVariable('NCHL', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_nchl.long_name = 'nanophytoplankton chlorophyll'
    # nc_nchl.units = 'mmol/m3'
    #
    # nc_dchl = nc.createVariable('DCHL', np.float64, ('VAL_time', 's_rho', 'eta_rho', 'xi_rho',))
    # nc_dchl.long_name = 'diatom chlorophyll'
    # nc_dchl.units = 'mmol/m3'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(VAL_time)
    nc_tend[:] = np.max(VAL_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_VAL_time[:] = VAL_time
    # nc_u[:] = 0
    # nc_v[:] = 0
    # nc_ubar[:] = 0
    # nc_vbar[:] = 0
    # nc_zeta[:] = 0
    nc_temps[:] = 0
    nc_salts[:] = 0
    nc_saltb[:] = 0
    # nc_no3[:] = 0
    # nc_po4[:] = 0
    # nc_si[:] = 0
    # nc_o2[:] = 0
    # nc_dic[:] = 0
    # nc_nh4[:] = 0
    # nc_fe[:] = 0
    # nc_talk[:] = 0
    # nc_caco3[:] = 0
    # nc_poc[:] = 0
    # nc_phy[:] = 0
    # nc_zoo[:] = 0
    # nc_doc[:] = 0
    # nc_phy2[:] = 0
    # nc_zoo2[:] = 0
    # nc_bsi[:] = 0
    # nc_bfe[:] = 0
    # nc_goc[:] = 0
    # nc_sfe[:] = 0
    # nc_dfe[:] = 0
    # nc_dsi[:] = 0
    # nc_nfe[:] = 0
    # nc_nchl[:] = 0
    # nc_dchl[:] = 0
    nc.close()
    #
    return


def create_ERA5style_ESMatm(filname, varlab, lngnam, atmunits, atime, latlen, lonlen):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    ################################################################
    #
    #

    print(' ')
    print(' Creating ESM file  for var: ' + filname)

    #
    #  Create the file
    #

    if os.path.exists(filname):
        os.remove(filname)

    print('Create: ' + filname)
    nc = netcdf(filname, 'w', format='NETCDF4')

    #
    # Create global attributes
    #
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.type = 'ERA5-style atmospheric variable from ESM '

    #
    #  Create dimensions
    #
    nc_dim_xi_rho = nc.createDimension('lon', lonlen)
    nc_dim_eta_rho = nc.createDimension('lat', latlen)
    nc_dim_surf_time = nc.createDimension('time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    nc_surf_time = nc.createVariable('time', np.float64, ('time',))
    nc_surf_time.long_name = 'time'
    nc_surf_time.units = 'days since 1990-1-1 0:0:0'
    nc_surf_time.calendar = 'gregorian'
    #
    nc_lat = nc.createVariable('lat', np.float64, ('lat',))
    nc_lat.long_name = 'latitude'
    nc_lat.standard_name = 'latitude'
    nc_lat.units = 'degree_north'
    #
    nc_lon = nc.createVariable('lon', np.float64, ('lon',))
    nc_lon.long_name = 'longitude'
    nc_lon.standard_name = 'longitude'
    nc_lon.units = 'degree_east'
    #
    nc_varlab = nc.createVariable(varlab, np.float64, ('time', 'lat', 'lon',))
    nc_varlab.long_name = lngnam
    nc_varlab.units = atmunits
    nc_varlab.missing_value = 9999.
    nc_varlab.coordinates = 'lat lon time'

    #
    # Write variables
    #
    #
    nc_surf_time[:] = atime
    nc_varlab[:] = 0
    nc_lat[:] = 0
    nc_lon[:] = 0

    nc.close()
    #
    return

def create_surftemp_RCave(VALname, grdname, title, theta_s, theta_b, hc, N, surf_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + VALname)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #
    type = 'VALIDATION file'
    history = 'CROCO'
    if os.path.exists(VALname):
        os.remove(VALname)

    print('Create: ' + VALname)
    nc = netcdf(VALname, 'w', format='NETCDF3_CLASSIC')

    #
    # Create global attributes
    #
    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = VALname
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #
    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_surf_time = nc.createDimension('surf_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_surf_time = nc.createVariable('surf_time', np.float64, ('surf_time',))
    nc_surf_time.long_name = 'time'
    nc_surf_time.units = 'days since 1990-1-1 0:0:0'
    nc_surf_time.calendar = 'gregorian'
    #
    nc_lat_rho = nc.createVariable('lat_rho', np.float64, ('eta_rho', 'xi_rho',))
    nc_lat_rho.long_name = 'latitude'
    nc_lat_rho.standard_name = 'latitude'
    nc_lat_rho.units = 'degrees_north'
    nc_lat_rho.coordinates = 'xi_rho eta_rho'
    #
    nc_lon_rho = nc.createVariable('lon_rho', np.float64, ('eta_rho', 'xi_rho',))
    nc_lon_rho.long_name = 'longitude'
    nc_lon_rho.standard_name = 'longitude'
    nc_lon_rho.units = 'degrees_east'
    nc_lon_rho.coordinates = 'xi_rho eta_rho'
    #
    nc_temps = nc.createVariable('temps', np.float64, ('surf_time', 'eta_rho', 'xi_rho',))
    nc_temps.long_name = 'potential temperature'
    nc_temps.units = 'Celsius'
    nc_temps.coordinates = 'lat_rho lon_rho surf_time'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(surf_time)
    nc_tend[:] = np.max(surf_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_surf_time[:] = surf_time
    nc_temps[:] = 0
    nc_lat_rho[:] = 0
    nc_lon_rho[:] = 0

    nc.close()
    #
    return


def create_surfbgc_RCave(VALname, grdname, title, theta_s, theta_b, hc, N, bott_time, vtransform):
    #
    #
    ################################################################
    #
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   ininame      Netcdf initial file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   ini_time     Initial time.(Real)
    #   clobber      Switch to allow or not writing over an existing
    #                file.(character string)
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + VALname)
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()
    #
    #
    #
    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + 'm)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    #
    #  Create the initial file
    #
    type = 'VALIDATION file'
    history = 'CROCO'
    if os.path.exists(VALname):
        os.remove(VALname)

    print('Create: ' + VALname)
    nc = netcdf(VALname, 'w', format='NETCDF3_CLASSIC')

    #
    # Create global attributes
    #
    nc.title = title
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.ini_file = VALname
    nc.grd_file = grdname
    nc.type = 'CROCO initial file'

    #
    #  Create dimensions
    #
    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_bott_time = nc.createDimension('bott_time', 0)
    # nc_dim_surf_time = nc.createDimension('surf_time', 0)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    #
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_bott_time = nc.createVariable('bott_time', np.float64, ('bott_time',))
    nc_bott_time.long_name = 'time'
    nc_bott_time.units = 'days since 1990-1-1 0:0:0'
    nc_bott_time.calendar = 'gregorian'
    #
    # nc_surf_time = nc.createVariable('surf_time', np.float64, ('surf_time',))
    # nc_surf_time.long_name = 'time'
    # nc_surf_time.units = 'days since 1990-1-1 0:0:0'
    # nc_surf_time.calendar = 'gregorian'
    #
    nc_lat_rho = nc.createVariable('lat_rho', np.float64, ('eta_rho', 'xi_rho',))
    nc_lat_rho.long_name = 'latitude'
    nc_lat_rho.standard_name = 'latitude'
    nc_lat_rho.units = 'degrees_north'
    nc_lat_rho.coordinates = 'xi_rho eta_rho'
    #
    nc_lon_rho = nc.createVariable('lon_rho', np.float64, ('eta_rho', 'xi_rho',))
    nc_lon_rho.long_name = 'longitude'
    nc_lon_rho.standard_name = 'longitude'
    nc_lon_rho.units = 'degrees_east'
    nc_lon_rho.coordinates = 'xi_rho eta_rho'
    #
    # nc_dics = nc.createVariable('DIC', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    # nc_dics.long_name = 'bottace dissolved inorganic carbon'
    # nc_dics.units = 'umol C L-1'
    # nc_dics.coordinates = 'lat_rho lon_rho bott_time'
    # #
    # nc_talks = nc.createVariable('TALK', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    # nc_talks.long_name = 'total alkalinity at bottace'
    # nc_talks.units = 'umol C L-1'
    # nc_talks.coordinates = 'lat_rho lon_rho bott_time'
    #
    # nc_no3s = nc.createVariable('NO3', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    # nc_no3s.long_name = 'bottom nitrate'
    # nc_no3s.units = 'umol N L-1'
    # nc_no3s.coordinates = 'lat_rho lon_rho bott_time'
    # #
    # nc_o2s = nc.createVariable('O2', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    # nc_o2s.long_name = 'bottom dissolved oxygen'
    # nc_o2s.units = 'umol O L-1'
    # nc_o2s.coordinates = 'lat_rho lon_rho bott_time'
    # #
    # nc_doss = nc.createVariable('DO_SAT', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    # nc_doss.long_name = 'bottom dissolved oxygen saturation'
    # nc_doss.units = '%'
    # nc_doss.coordinates = 'lat_rho lon_rho bott_time'
    #
    # nc_no3s = nc.createVariable('NO3', np.float64, ('surf_time', 'eta_rho', 'xi_rho',))
    # nc_no3s.long_name = 'surface nitrate'
    # nc_no3s.units = 'umol N L-1'
    # nc_no3s.coordinates = 'lat_rho lon_rho surf_time'
    # #
    # nc_o2s = nc.createVariable('O2', np.float64, ('surf_time', 'eta_rho', 'xi_rho',))
    # nc_o2s.long_name = 'surface dissolved oxygen'
    # nc_o2s.units = 'umol O L-1'
    # nc_o2s.coordinates = 'lat_rho lon_rho surf_time'
    # #
    # nc_doss = nc.createVariable('DO_SAT', np.float64, ('surf_time', 'eta_rho', 'xi_rho',))
    # nc_doss.long_name = 'surface dissolved oxygen saturation'
    # nc_doss.units = '%'
    # nc_doss.coordinates = 'lat_rho lon_rho surf_time'
    #
    nc_temps = nc.createVariable('temp', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    nc_temps.long_name = 'bottom temperature'
    nc_temps.units = 'degC'
    nc_temps.coordinates = 'lat_rho lon_rho bott_time'
    #
    nc_salts = nc.createVariable('salt', np.float64, ('bott_time', 'eta_rho', 'xi_rho',))
    nc_salts.long_name = 'bottom salinity'
    nc_salts.units = 'psu'
    nc_salts.coordinates = 'lat_rho lon_rho bott_time'
    #
    # nc_temps = nc.createVariable('temp', np.float64, ('surf_time', 'eta_rho', 'xi_rho',))
    # nc_temps.long_name = 'surface temperature'
    # nc_temps.units = 'degC'
    # nc_temps.coordinates = 'lat_rho lon_rho surf_time'
    # #
    # nc_salts = nc.createVariable('salt', np.float64, ('surf_time', 'eta_rho', 'xi_rho',))
    # nc_salts.long_name = 'surface salinity'
    # nc_salts.units = 'psu'
    # nc_salts.coordinates = 'lat_rho lon_rho surf_time'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #

    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(bott_time)
    nc_tend[:] = np.max(bott_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_Cs_r[:] = Cs_r
    nc_bott_time[:] = bott_time
    # nc_surf_time[:] = bott_time
    # nc_no3s[:] = 0
    # nc_o2s[:] = 0
    # nc_doss[:] = 0
    nc_temps[:] = 0
    nc_salts[:] = 0
    nc_lat_rho[:] = 0
    nc_lon_rho[:] = 0

    nc.close()
    #
    return


#
#
# #####################################################################
# #### FUNCTION CREATE_BRYFILE ########################################
# #####################################################################
#
#

def create_bryfile(bryname, grdname, title, obc,
                   theta_s, theta_b, hc, N,
                   bry_time, cycle, vtransform):
    #
    #
    ################################################################
    #
    # function create_bryfile(bryname,grdname,title,obc...
    #                          theta_s,theta_b,hc,N,...
    #                          bry_time,cycle,clobber)
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   bryname      Netcdf climatology file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   obc          open boundaries flag (1=open , [S E N W]).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   bry_time     time.(vector)
    #   cycle        Length (days) for cycling the climatology.(Real)
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + bryname)
    print(' ')
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()

    #
    #
    #

    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + ' m)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    Nt = 0

    #
    #  Create the boundary file
    #
    type = 'BOUNDARY file'
    history = 'CROCO'

    if os.path.exists(bryname):
        os.remove(bryname)

    print('Create: ' + bryname)
    nc = netcdf(bryname, 'w', format='NETCDF4')

    #
    # set global attributes
    #

    nc.type = 'CROCO boundary file'
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.title = title
    nc.bry_file = bryname
    nc.grd_file = grdname

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)
    nc_dim_bry_time = nc.createDimension('bry_time', Nt)
    nc_dim_tclm_time = nc.createDimension('tclm_time', Nt)
    nc_dim_temp_time = nc.createDimension('temp_time', Nt)
    nc_dim_sclm_time = nc.createDimension('sclm_time', Nt)
    nc_dim_salt_time = nc.createDimension('salt_time', Nt)
    nc_dim_uclm_time = nc.createDimension('uclm_time', Nt)
    nc_dim_vclm_time = nc.createDimension('vclm_time', Nt)
    nc_dim_v2d_time = nc.createDimension('v2d_time', Nt)
    nc_dim_v3d_time = nc.createDimension('v3d_time', Nt)
    nc_dim_ssh_time = nc.createDimension('ssh_time', Nt)
    nc_dim_zeta_time = nc.createDimension('zeta_time', Nt)
    nc_dim_no3_time = nc.createDimension('no3_time', Nt)
    nc_dim_po4_time = nc.createDimension('po4_time', Nt)
    nc_dim_si_time = nc.createDimension('si_time', Nt)
    nc_dim_o2_time = nc.createDimension('o2_time', Nt)
    nc_dim_dic_time = nc.createDimension('dic_time', Nt)
    nc_dim_dic_time = nc.createDimension('talk_time', Nt)
    nc_dim_fe_time = nc.createDimension('fe_time', Nt)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #

    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    # s
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_bry_time = nc.createVariable('bry_time', np.float64, ('bry_time',))
    nc_bry_time.long_name = 'time for boundary climatology'
    nc_bry_time.units = 'day'
    nc_bry_time.calendar = 'XXX days in every year'
    nc_bry_time.cycle_length = cycle
    #
    nc_tclm_time = nc.createVariable('tclm_time', np.float64, ('tclm_time',))
    nc_tclm_time.long_name = 'time for temperature climatology'
    nc_tclm_time.units = 'day'
    nc_tclm_time.calendar = 'XXX days in every year'
    nc_tclm_time.cycle_length = cycle
    #
    nc_temp_time = nc.createVariable('temp_time', np.float64, ('temp_time',))
    nc_temp_time.long_name = 'time for temperature climatology'
    nc_temp_time.units = 'day'
    nc_temp_time.calendar = 'XXX days in every year'
    nc_temp_time.cycle_length = cycle
    #
    nc_sclm_time = nc.createVariable('sclm_time', np.float64, ('sclm_time',))
    nc_sclm_time.long_name = 'time for salinity climatology'
    nc_sclm_time.units = 'day'
    nc_sclm_time.calendar = 'XXX days in every year'
    nc_sclm_time.cycle_length = cycle
    #
    nc_salt_time = nc.createVariable('salt_time', np.float64, ('salt_time',))
    nc_salt_time.long_name = 'time for salinity climatology'
    nc_salt_time.units = 'day'
    nc_salt_time.calendar = 'XXX days in every year'
    nc_salt_time.cycle_length = cycle
    #
    nc_uclm_time = nc.createVariable('uclm_time', np.float64, ('uclm_time',))
    nc_uclm_time.long_name = 'time climatological u'
    nc_uclm_time.units = 'day'
    nc_uclm_time.calendar = 'XXX days in every year'
    nc_uclm_time.cycle_length = cycle
    #
    nc_vclm_time = nc.createVariable('vclm_time', np.float64, ('vclm_time',))
    nc_vclm_time.long_name = 'time climatological v'
    nc_vclm_time.units = 'day'
    nc_vclm_time.calendar = 'XXX days in every year'
    nc_vclm_time.cycle_length = cycle
    #
    nc_v2d_time = nc.createVariable('v2d_time', np.float64, ('v2d_time',))
    nc_v2d_time.long_name = 'time for 2D velocity climatology'
    nc_v2d_time.units = 'day'
    nc_v2d_time.calendar = 'XXX days in every year'
    nc_v2d_time.cycle_length = cycle
    #
    nc_v3d_time = nc.createVariable('v3d_time', np.float64, ('v3d_time',))
    nc_v3d_time.long_name = 'time for 3D velocity climatology'
    nc_v3d_time.units = 'day'
    nc_v3d_time.calendar = 'XXX days in every year'
    nc_v3d_time.cycle_length = cycle
    #
    nc_ssh_time = nc.createVariable('ssh_time', np.float64, ('ssh_time',))
    nc_ssh_time.long_name = 'time for sea surface height'
    nc_ssh_time.units = 'day'
    nc_ssh_time.calendar = 'XXX days in every year'
    nc_ssh_time.cycle_length = cycle
    #
    nc_zeta_time = nc.createVariable('zeta_time', np.float64, ('zeta_time',))
    nc_zeta_time.long_name = 'time for sea surface height'
    nc_zeta_time.units = 'day'
    nc_zeta_time.calendar = 'XXX days in every year'
    nc_zeta_time.cycle_length = cycle

    nc_no3_time = nc.createVariable('no3_time', np.float64, ('no3_time',))
    nc_no3_time.long_name = 'time for nitrate'
    nc_no3_time.units = 'day'
    nc_no3_time.calendar = 'XXX days in every year'
    nc_no3_time.cycle_length = cycle

    nc_po4_time = nc.createVariable('po4_time', np.float64, ('po4_time',))
    nc_po4_time.long_name = 'time for orthophosphate'
    nc_po4_time.units = 'day'
    nc_po4_time.calendar = 'XXX days in every year'
    nc_po4_time.cycle_length = cycle

    nc_si_time = nc.createVariable('si_time', np.float64, ('si_time',))
    nc_si_time.long_name = 'time for silicate'
    nc_si_time.units = 'day'
    nc_si_time.calendar = 'XXX days in every year'
    nc_si_time.cycle_length = cycle

    nc_o2_time = nc.createVariable('o2_time', np.float64, ('o2_time',))
    nc_o2_time.long_name = 'time for dissolved oxygen'
    nc_o2_time.units = 'day'
    nc_o2_time.calendar = 'XXX days in every year'
    nc_o2_time.cycle_length = cycle

    nc_dic_time = nc.createVariable('dic_time', np.float64, ('dic_time',))
    nc_dic_time.long_name = 'time for dissolved inorganic carbon'
    nc_dic_time.units = 'day'
    nc_dic_time.calendar = 'XXX days in every year'
    nc_dic_time.cycle_length = cycle

    nc_talk_time = nc.createVariable('talk_time', np.float64, ('talk_time',))
    nc_talk_time.long_name = 'time for total alkalinity'
    nc_talk_time.units = 'day'
    nc_talk_time.calendar = 'XXX days in every year'
    nc_talk_time.cycle_length = cycle

    nc_fe_time = nc.createVariable('fe_time', np.float64, ('fe_time',))
    nc_fe_time.long_name = 'time for iron'
    nc_fe_time.units = 'day'
    nc_fe_time.calendar = 'XXX days in every year'
    nc_fe_time.cycle_length = cycle

    #
    if obc[0] == 1:
        #
        #   Southern boundary
        #
        nc_temp_south = nc.createVariable('temp_south', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_south.long_name = 'southern boundary potential temperature'
        nc_temp_south.units = 'Celsius'
        nc_temp_south.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_south = nc.createVariable('salt_south', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_south.long_name = 'southern boundary salinity'
        nc_salt_south.units = 'PSU'
        nc_salt_south.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_south = nc.createVariable('u_south', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_south.long_name = 'southern boundary u-momentum component'
        nc_u_south.units = 'meter second-1'
        nc_u_south.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_south = nc.createVariable('v_south', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_south.long_name = 'southern boundary v-momentum component'
        nc_v_south.units = 'meter second-1'
        nc_v_south.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_south = nc.createVariable('ubar_south', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_south.long_name = 'southern boundary vertically integrated u-momentum component'
        nc_ubar_south.units = 'meter second-1'
        nc_ubar_south.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_south = nc.createVariable('vbar_south', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_south.long_name = 'southern boundary vertically integrated v-momentum component'
        nc_vbar_south.units = 'meter second-1'
        nc_vbar_south.coordinates = 'lon_v vclm_time'
        #
        nc_zeta_south = nc.createVariable('zeta_south', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_south.long_name = 'southern boundary sea surface height'
        nc_zeta_south.units = 'meter'
        nc_zeta_south.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_south = nc.createVariable('NO3_south', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_south.long_name = 'southern boundary nitrate'
        nc_no3_south.units = 'mmol/m3'
        nc_no3_south.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_south = nc.createVariable('PO4_south', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_south.long_name = 'southern boundary orthophosphate'
        nc_po4_south.units = 'mmol/m3'
        nc_po4_south.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_south = nc.createVariable('Si_south', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_south.long_name = 'southern boundary silicate'
        nc_si_south.units = 'mmol/m3'
        nc_si_south.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_south = nc.createVariable('O2_south', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_south.long_name = 'southern boundary dissolved oxygen'
        nc_o2_south.units = 'mmol/m3'
        nc_o2_south.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_south = nc.createVariable('DIC_south', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_south.long_name = 'southern boundary dissolved inorganic carbon'
        nc_dic_south.units = 'mmol/m3'
        nc_dic_south.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_talk_south = nc.createVariable('TALK_south', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_south.long_name = 'southern boundary total alkalinity'
        nc_talk_south.units = 'mmol/m3'
        nc_talk_south.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_fe_south = nc.createVariable('FER_south', np.float64, ('fe_time', 's_rho', 'xi_rho',))
        nc_fe_south.long_name = 'southern boundary iron'
        nc_fe_south.units = 'mmol/m3'
        #

    if obc[1] == 1:
        #
        #   Eastern boundary
        #
        nc_temp_east = nc.createVariable('temp_east', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_east.long_name = 'eastern boundary potential temperature'
        nc_temp_east.units = 'Celsius'
        nc_temp_east.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_east = nc.createVariable('salt_east', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_east.long_name = 'eastern boundary salinity'
        nc_salt_east.units = 'PSU'
        nc_salt_east.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_east = nc.createVariable('u_east', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_east.long_name = 'eastern boundary u-momentum component'
        nc_u_east.units = 'meter second-1'
        nc_u_east.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_east = nc.createVariable('v_east', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_east.long_name = 'eastern boundary v-momentum component'
        nc_v_east.units = 'meter second-1'
        nc_v_east.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_east = nc.createVariable('ubar_east', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_east.long_name = 'eastern boundary vertically integrated u-momentum component'
        nc_ubar_east.units = 'meter second-1'
        nc_ubar_east.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_east = nc.createVariable('vbar_east', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_east.long_name = 'eastern boundary vertically integrated v-momentum component'
        nc_vbar_east.units = 'meter second-1'
        nc_vbar_east.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_east = nc.createVariable('zeta_east', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_east.long_name = 'eastern boundary sea surface height'
        nc_zeta_east.units = 'meter'
        nc_zeta_east.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_east = nc.createVariable('NO3_east', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_east.long_name = 'eastern boundary nitrate'
        nc_no3_east.units = 'mmol/m3'
        nc_no3_east.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_east = nc.createVariable('PO4_east', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_east.long_name = 'eastern boundary orthophosphate'
        nc_po4_east.units = 'mmol/m3'
        nc_po4_east.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_east = nc.createVariable('Si_east', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_east.long_name = 'eastern boundary silicate'
        nc_si_east.units = 'mmol/m3'
        nc_si_east.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_east = nc.createVariable('O2_east', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_east.long_name = 'eastern boundary dissolved oxygen'
        nc_o2_east.units = 'mmol/m3'
        nc_o2_east.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_east = nc.createVariable('DIC_east', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_east.long_name = 'eastern boundary dissolved inorganic carbon'
        nc_dic_east.units = 'mmol/m3'
        nc_dic_east.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_talk_east = nc.createVariable('TALK_east', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_east.long_name = 'eastern boundary total alkalinity'
        nc_talk_east.units = 'mmol/m3'
        nc_talk_east.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_fe_east = nc.createVariable('FER_east', np.float64, ('fe_time', 's_rho', 'eta_rho',))
        nc_fe_east.long_name = 'eastern boundary iron'
        nc_fe_east.units = 'mmol/m3'
        #
        #

    if obc[2] == 1:
        #
        #   Northern boundary
        #
        nc_temp_north = nc.createVariable('temp_north', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_north.long_name = 'northern boundary potential temperature'
        nc_temp_north.units = 'Celsius'
        nc_temp_north.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_north = nc.createVariable('salt_north', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_north.long_name = 'northern boundary salinity'
        nc_salt_north.units = 'PSU'
        nc_salt_north.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_north = nc.createVariable('u_north', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_north.long_name = 'northern boundary u-momentum component'
        nc_u_north.units = 'meter second-1'
        nc_u_north.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_north = nc.createVariable('v_north', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_north.long_name = 'northern boundary v-momentum component'
        nc_v_north.units = 'meter second-1'
        nc_v_north.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_north = nc.createVariable('ubar_north', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_north.long_name = 'northern boundary vertically integrated u-momentum component'
        nc_ubar_north.units = 'meter second-1'
        nc_ubar_north.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_north = nc.createVariable('vbar_north', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_north.long_name = 'northern boundary vertically integrated v-momentum component'
        nc_vbar_north.units = 'meter second-1'
        nc_vbar_north.coordinates = 'lon_v vclm_time'

        nc_zeta_north = nc.createVariable('zeta_north', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_north.long_name = 'northern boundary sea surface height'
        nc_zeta_north.units = 'meter'
        nc_zeta_north.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_north = nc.createVariable('NO3_north', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_north.long_name = 'northern boundary nitrate'
        nc_no3_north.units = 'mmol/m3'
        nc_no3_north.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_north = nc.createVariable('PO4_north', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_north.long_name = 'northern boundary orthophosphate'
        nc_po4_north.units = 'mmol/m3'
        nc_po4_north.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_north = nc.createVariable('Si_north', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_north.long_name = 'northern boundary silicate'
        nc_si_north.units = 'mmol/m3'
        nc_si_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_north = nc.createVariable('O2_north', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_north.long_name = 'northern boundary dissolved oxygen'
        nc_o2_north.units = 'mmol/m3'
        nc_o2_north.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_north = nc.createVariable('DIC_north', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_north.long_name = 'northern boundary dissolved inorganic carbon'
        nc_dic_north.units = 'mmol/m3'
        nc_dic_north.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_talk_north = nc.createVariable('TALK_north', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_north.long_name = 'northern boundary total alkalinity'
        nc_talk_north.units = 'mmol/m3'
        nc_talk_north.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_fe_north = nc.createVariable('FER_north', np.float64, ('fe_time', 's_rho', 'xi_rho',))
        nc_fe_north.long_name = 'northern boundary iron'
        nc_fe_north.units = 'mmol/m3'
    #

    if obc[3] == 1:
        #
        #   Western boundary
        #
        nc_temp_west = nc.createVariable('temp_west', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_west.long_name = 'western boundary potential temperature'
        nc_temp_west.units = 'Celsius'
        nc_temp_west.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_west = nc.createVariable('salt_west', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_west.long_name = 'western boundary salinity'
        nc_salt_west.units = 'PSU'
        nc_salt_west.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_west = nc.createVariable('u_west', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_west.long_name = 'western boundary u-momentum component'
        nc_u_west.units = 'meter second-1'
        nc_u_west.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_west = nc.createVariable('v_west', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_west.long_name = 'western boundary v-momentum component'
        nc_v_west.units = 'meter second-1'
        nc_v_west.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_west = nc.createVariable('ubar_west', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_west.long_name = 'western boundary vertically integrated u-momentum component'
        nc_ubar_west.units = 'meter second-1'
        nc_ubar_west.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_west = nc.createVariable('vbar_west', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_west.long_name = 'western boundary vertically integrated v-momentum component'
        nc_vbar_west.units = 'meter second-1'
        nc_vbar_west.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_west = nc.createVariable('zeta_west', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_west.long_name = 'western boundary sea surface height'
        nc_zeta_west.units = 'meter'
        nc_zeta_west.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_west = nc.createVariable('NO3_west', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_west.long_name = 'western boundary nitrate'
        nc_no3_west.units = 'mmol/m3'
        nc_no3_west.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_west = nc.createVariable('PO4_west', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_west.long_name = 'western boundary orthophosphate'
        nc_po4_west.units = 'mmol/m3'
        nc_po4_west.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_west = nc.createVariable('Si_west', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_west.long_name = 'western boundary silicate'
        nc_si_west.units = 'mmol/m3'
        nc_si_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_west = nc.createVariable('O2_west', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_west.long_name = 'western boundary dissolved oxygen'
        nc_o2_west.units = 'mmol/m3'
        nc_o2_west.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_west = nc.createVariable('DIC_west', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_west.long_name = 'western boundary dissolved inorganic carbon'
        nc_dic_west.units = 'mmol/m3'
        nc_dic_west.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_talk_west = nc.createVariable('TALK_west', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_west.long_name = 'western boundary total alkalinity'
        nc_talk_west.units = 'mmol/m3'
        nc_talk_west.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_fe_west = nc.createVariable('FER_west', np.float64, ('fe_time', 's_rho', 'eta_rho',))
        nc_fe_west.long_name = 'western boundary iron'
        nc_fe_west.units = 'mmol/m3'
        #
        #

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #
    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(bry_time)
    nc_tend[:] = np.max(bry_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_sc_w[:] = sc_w
    nc_Cs_r[:] = Cs_r
    nc_Cs_w[:] = Cs_w
    nc_tclm_time[:] = bry_time
    nc_temp_time[:] = bry_time
    nc_sclm_time[:] = bry_time
    nc_salt_time[:] = bry_time
    nc_uclm_time[:] = bry_time
    nc_vclm_time[:] = bry_time
    nc_v2d_time[:] = bry_time
    nc_v3d_time[:] = bry_time
    nc_ssh_time[:] = bry_time
    nc_zeta_time[:] = bry_time
    nc_bry_time[:] = bry_time
    nc_no3_time[:] = bry_time
    nc_po4_time[:] = bry_time
    nc_si_time[:] = bry_time
    nc_o2_time[:] = bry_time
    nc_dic_time[:] = bry_time
    nc_talk_time[:] = bry_time
    # nc_nh4_time[:] = bry_time
    nc_fe_time[:] = bry_time

    if obc[0] == 1:
        nc_u_south[:] = 0.
        nc_v_south[:] = 0.
        nc_ubar_south[:] = 0.
        nc_vbar_south[:] = 0.
        nc_zeta_south[:] = 0.
        nc_temp_south[:] = 0.
        nc_salt_south[:] = 0.
        nc_no3_south[:] = 0.
        nc_po4_south[:] = 0.
        nc_si_south[:] = 0.
        nc_o2_south[:] = 0.
        nc_dic_south[:] = 0.
        nc_talk_south[:] = 0.
        nc_fe_south[:] = 0.

    if obc[1] == 1:
        nc_u_east[:] = 0.
        nc_v_east[:] = 0.
        nc_ubar_east[:] = 0.
        nc_vbar_east[:] = 0.
        nc_zeta_east[:] = 0.
        nc_temp_east[:] = 0.
        nc_salt_east[:] = 0.
        nc_no3_east[:] = 0.
        nc_po4_east[:] = 0.
        nc_si_east[:] = 0.
        nc_o2_east[:] = 0.
        nc_dic_east[:] = 0.
        nc_talk_east[:] = 0.
        nc_fe_east[:] = 0.

    if obc[2] == 1:
        nc_u_north[:] = 0.
        nc_v_north[:] = 0.
        nc_ubar_north[:] = 0.
        nc_vbar_north[:] = 0.
        nc_zeta_north[:] = 0.
        nc_temp_north[:] = 0.
        nc_salt_north[:] = 0.
        nc_no3_north[:] = 0.
        nc_po4_north[:] = 0.
        nc_si_north[:] = 0.
        nc_o2_north[:] = 0.
        nc_dic_north[:] = 0.
        nc_talk_north[:] = 0.
        nc_fe_north[:] = 0.

    if obc[3] == 1:
        nc_u_west[:] = 0.
        nc_v_west[:] = 0.
        nc_ubar_west[:] = 0.
        nc_vbar_west[:] = 0.
        nc_zeta_west[:] = 0.
        nc_temp_west[:] = 0.
        nc_salt_west[:] = 0.
        nc_no3_west[:] = 0.
        nc_po4_west[:] = 0.
        nc_si_west[:] = 0.
        nc_o2_west[:] = 0.
        nc_dic_west[:] = 0.
        nc_talk_west[:] = 0.
        nc_fe_west[:] = 0.
        #

    nc.close()

    return


def create_bryfile_PISCES_NORESM(bryname, grdname, title, obc,
                                 theta_s, theta_b, hc, N,
                                 bry_time, nor_time, pis_time, cycle, vtransform):
    #
    #
    ################################################################
    #
    # function create_bryfile(bryname,grdname,title,obc...
    #                          theta_s,theta_b,hc,N,...
    #                          bry_time,cycle,clobber)
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   bryname      Netcdf climatology file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   obc          open boundaries flag (1=open , [S E N W]).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   bry_time     time.(vector)
    #   cycle        Length (days) for cycling the climatology.(Real)
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + bryname)
    print(' ')
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()

    #
    #
    #

    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + ' m)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    Nt = 0

    #
    #  Create the boundary file
    #
    type = 'BOUNDARY file'
    history = 'CROCO'

    if os.path.exists(bryname):
        os.remove(bryname)

    print('Create: ' + bryname)
    nc = netcdf(bryname, 'w', format='NETCDF4')

    #
    # set global attributes
    #

    nc.type = 'CROCO boundary file'
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.title = title
    nc.bry_file = bryname
    nc.grd_file = grdname

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)

    # PHYSICS time dimension
    nc_dim_bry_time = nc.createDimension('bry_time', Nt)
    nc_dim_tclm_time = nc.createDimension('tclm_time', Nt)
    nc_dim_temp_time = nc.createDimension('temp_time', Nt)
    nc_dim_sclm_time = nc.createDimension('sclm_time', Nt)
    nc_dim_salt_time = nc.createDimension('salt_time', Nt)
    nc_dim_uclm_time = nc.createDimension('uclm_time', Nt)
    nc_dim_vclm_time = nc.createDimension('vclm_time', Nt)
    nc_dim_v2d_time = nc.createDimension('v2d_time', Nt)
    nc_dim_v3d_time = nc.createDimension('v3d_time', Nt)
    nc_dim_ssh_time = nc.createDimension('ssh_time', Nt)
    nc_dim_zeta_time = nc.createDimension('zeta_time', Nt)

    # # BGC time dimension
    nc_dim_no3_time = nc.createDimension('no3_time', Nt)
    nc_dim_po4_time = nc.createDimension('po4_time', Nt)
    nc_dim_si_time = nc.createDimension('si_time', Nt)
    nc_dim_o2_time = nc.createDimension('o2_time', Nt)
    nc_dim_dic_time = nc.createDimension('dic_time', Nt)
    # nc_dim_nh4_time = nc.createDimension('nh4_time', Nt)
    nc_dim_fer_time = nc.createDimension('fer_time', Nt)
    # nc_dim_doc_time = nc.createDimension('doc_time', Nt)
    nc_dim_talk_time = nc.createDimension('talk_time', Nt)
    # nc_dim_calc_time = nc.createDimension('calc_time', Nt)
    # nc_dim_poc_time = nc.createDimension('poc_time', Nt)
    # nc_dim_phy_time = nc.createDimension('phy_time', Nt)
    # nc_dim_zoo_time = nc.createDimension('zoo_time', Nt)
    # nc_dim_dia_time = nc.createDimension('dia_time', Nt)
    # nc_dim_mes_time = nc.createDimension('mes_time', Nt)
    # nc_dim_bsi_time = nc.createDimension('bsi_time', Nt)
    # nc_dim_bfe_time = nc.createDimension('bfe_time', Nt)
    # nc_dim_goc_time = nc.createDimension('goc_time', Nt)
    # nc_dim_sfe_time = nc.createDimension('sfe_time', Nt)
    # nc_dim_dfe_time = nc.createDimension('dfe_time', Nt)
    # nc_dim_dsi_time = nc.createDimension('dsi_time', Nt)
    # nc_dim_nfe_time = nc.createDimension('nfe_time', Nt)
    # nc_dim_nch_time = nc.createDimension('nch_time', Nt)
    # nc_dim_dch_time = nc.createDimension('dch_time', Nt)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #
    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    # s
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_bry_time = nc.createVariable('bry_time', np.float64, ('bry_time',))
    nc_bry_time.long_name = 'time for boundary climatology'
    nc_bry_time.units = 'day'
    nc_bry_time.calendar = 'XXX days in every year'
    nc_bry_time.cycle_length = cycle
    #
    nc_tclm_time = nc.createVariable('tclm_time', np.float64, ('tclm_time',))
    nc_tclm_time.long_name = 'time for temperature climatology'
    nc_tclm_time.units = 'day'
    nc_tclm_time.calendar = 'XXX days in every year'
    nc_tclm_time.cycle_length = cycle
    #
    nc_temp_time = nc.createVariable('temp_time', np.float64, ('temp_time',))
    nc_temp_time.long_name = 'time for temperature'
    nc_temp_time.units = 'day'
    nc_temp_time.calendar = 'XXX days in every year'
    nc_temp_time.cycle_length = cycle
    #
    nc_sclm_time = nc.createVariable('sclm_time', np.float64, ('sclm_time',))
    nc_sclm_time.long_name = 'time for salinity climatology'
    nc_sclm_time.units = 'day'
    nc_sclm_time.calendar = 'XXX days in every year'
    nc_sclm_time.cycle_length = cycle
    #
    nc_salt_time = nc.createVariable('salt_time', np.float64, ('salt_time',))
    nc_salt_time.long_name = 'time for salinity'
    nc_salt_time.units = 'day'
    nc_salt_time.calendar = 'XXX days in every year'
    nc_salt_time.cycle_length = cycle
    #
    nc_uclm_time = nc.createVariable('uclm_time', np.float64, ('uclm_time',))
    nc_uclm_time.long_name = 'time climatological u'
    nc_uclm_time.units = 'day'
    nc_uclm_time.calendar = 'XXX days in every year'
    nc_uclm_time.cycle_length = cycle
    #
    nc_vclm_time = nc.createVariable('vclm_time', np.float64, ('vclm_time',))
    nc_vclm_time.long_name = 'time climatological v'
    nc_vclm_time.units = 'day'
    nc_vclm_time.calendar = 'XXX days in every year'
    nc_vclm_time.cycle_length = cycle
    #
    nc_v2d_time = nc.createVariable('v2d_time', np.float64, ('v2d_time',))
    nc_v2d_time.long_name = 'time for 2D velocity'
    nc_v2d_time.units = 'day'
    nc_v2d_time.calendar = 'XXX days in every year'
    nc_v2d_time.cycle_length = cycle
    #
    nc_v3d_time = nc.createVariable('v3d_time', np.float64, ('v3d_time',))
    nc_v3d_time.long_name = 'time for 3D velocity'
    nc_v3d_time.units = 'day'
    nc_v3d_time.calendar = 'XXX days in every year'
    nc_v3d_time.cycle_length = cycle
    #
    nc_ssh_time = nc.createVariable('ssh_time', np.float64, ('ssh_time',))
    nc_ssh_time.long_name = 'time for sea surface height'
    nc_ssh_time.units = 'day'
    nc_ssh_time.calendar = 'XXX days in every year'
    nc_ssh_time.cycle_length = cycle
    #
    nc_zeta_time = nc.createVariable('zeta_time', np.float64, ('zeta_time',))
    nc_zeta_time.long_name = 'time for sea surface height'
    nc_zeta_time.units = 'day'
    nc_zeta_time.calendar = 'XXX days in every year'
    nc_zeta_time.cycle_length = cycle
    #
    nc_no3_time = nc.createVariable('no3_time', np.float64, ('no3_time',))
    nc_no3_time.long_name = 'time for nitrate'
    nc_no3_time.units = 'day'
    nc_no3_time.calendar = 'XXX days in every year'
    nc_no3_time.cycle_length = cycle
    #
    nc_po4_time = nc.createVariable('po4_time', np.float64, ('po4_time',))
    nc_po4_time.long_name = 'time for orthophosphate'
    nc_po4_time.units = 'day'
    nc_po4_time.calendar = 'XXX days in every year'
    nc_po4_time.cycle_length = cycle
    #
    nc_si_time = nc.createVariable('si_time', np.float64, ('si_time',))
    nc_si_time.long_name = 'time for silicate'
    nc_si_time.units = 'day'
    nc_si_time.calendar = 'XXX days in every year'
    nc_si_time.cycle_length = cycle
    #
    nc_o2_time = nc.createVariable('o2_time', np.float64, ('o2_time',))
    nc_o2_time.long_name = 'time for dissolved oxygen'
    nc_o2_time.units = 'day'
    nc_o2_time.calendar = 'XXX days in every year'
    nc_o2_time.cycle_length = cycle
    #
    nc_dic_time = nc.createVariable('dic_time', np.float64, ('dic_time',))
    nc_dic_time.long_name = 'time for dissolved inorganic carbon'
    nc_dic_time.units = 'day'
    nc_dic_time.calendar = 'XXX days in every year'
    nc_dic_time.cycle_length = cycle
    # #
    # nc_nh4_time = nc.createVariable('nh4_time', np.float64, ('nh4_time',))
    # nc_nh4_time.long_name = 'time for ammonium'
    # nc_nh4_time.units = 'day'
    # nc_nh4_time.calendar = 'XXX days in every year'
    # nc_nh4_time.cycle_length = cycle
    #
    nc_fer_time = nc.createVariable('fer_time', np.float64, ('fer_time',))
    nc_fer_time.long_name = 'time for iron'
    nc_fer_time.units = 'day'
    nc_fer_time.calendar = 'XXX days in every year'
    nc_fer_time.cycle_length = cycle
    #
    nc_talk_time = nc.createVariable('talk_time', np.float64, ('talk_time',))
    nc_talk_time.long_name = 'time for alkalinity'
    nc_talk_time.units = 'day'
    nc_talk_time.calendar = 'XXX days in every year'
    nc_talk_time.cycle_length = cycle
    # #
    # nc_calc_time = nc.createVariable('calc_time', np.float64, ('calc_time',))
    # nc_calc_time.long_name = 'time for calcite'
    # nc_calc_time.units = 'day'
    # nc_calc_time.calendar = 'XXX days in every year'
    # nc_calc_time.cycle_length = cycle
    # #
    # nc_poc_time = nc.createVariable('poc_time', np.float64, ('poc_time',))
    # nc_poc_time.long_name = 'time for particulate organic carbon'
    # nc_poc_time.units = 'day'
    # nc_poc_time.calendar = 'XXX days in every year'
    # nc_poc_time.cycle_length = cycle
    # #
    # nc_phy_time = nc.createVariable('phy_time', np.float64, ('phy_time',))
    # nc_phy_time.long_name = 'time for nanophytoplankton'
    # nc_phy_time.units = 'day'
    # nc_phy_time.calendar = 'XXX days in every year'
    # nc_phy_time.cycle_length = cycle
    # #
    # nc_zoo_time = nc.createVariable('zoo_time', np.float64, ('zoo_time',))
    # nc_zoo_time.long_name = 'time for microzooplankton'
    # nc_zoo_time.units = 'day'
    # nc_zoo_time.calendar = 'XXX days in every year'
    # nc_zoo_time.cycle_length = cycle
    # #
    # nc_doc_time = nc.createVariable('doc_time', np.float64, ('doc_time',))
    # nc_doc_time.long_name = 'time for dissolved organic carbon'
    # nc_doc_time.units = 'day'
    # nc_doc_time.calendar = 'XXX days in every year'
    # nc_doc_time.cycle_length = cycle
    # #
    # nc_dia_time = nc.createVariable('dia_time', np.float64, ('dia_time',))
    # nc_dia_time.long_name = 'time for diatom'
    # nc_dia_time.units = 'day'
    # nc_dia_time.calendar = 'XXX days in every year'
    # nc_dia_time.cycle_length = cycle
    # #
    # nc_mes_time = nc.createVariable('mes_time', np.float64, ('mes_time',))
    # nc_mes_time.long_name = 'time for mesozooplankton'
    # nc_mes_time.units = 'day'
    # nc_mes_time.calendar = 'XXX days in every year'
    # nc_mes_time.cycle_length = cycle
    # #
    # nc_bsi_time = nc.createVariable('bsi_time', np.float64, ('bsi_time',))
    # nc_bsi_time.long_name = 'time for biogenic silica'
    # nc_bsi_time.units = 'day'
    # nc_bsi_time.calendar = 'XXX days in every year'
    # nc_bsi_time.cycle_length = cycle
    # #
    # nc_bfe_time = nc.createVariable('bfe_time', np.float64, ('bfe_time',))
    # nc_bfe_time.long_name = 'time for big particle iron'
    # nc_bfe_time.units = 'day'
    # nc_bfe_time.calendar = 'XXX days in every year'
    # nc_bfe_time.cycle_length = cycle
    # #
    # nc_goc_time = nc.createVariable('goc_time', np.float64, ('goc_time',))
    # nc_goc_time.long_name = 'time for big particulate organic carbon'
    # nc_goc_time.units = 'day'
    # nc_goc_time.calendar = 'XXX days in every year'
    # nc_goc_time.cycle_length = cycle
    # #
    # nc_sfe_time = nc.createVariable('sfe_time', np.float64, ('sfe_time',))
    # nc_sfe_time.long_name = 'time for iron in the small particles'
    # nc_sfe_time.units = 'day'
    # nc_sfe_time.calendar = 'XXX days in every year'
    # nc_sfe_time.cycle_length = cycle
    # #
    # nc_dfe_time = nc.createVariable('dfe_time', np.float64, ('dfe_time',))
    # nc_dfe_time.long_name = 'time for iron content of the diatoms'
    # nc_dfe_time.units = 'day'
    # nc_dfe_time.calendar = 'XXX days in every year'
    # nc_dfe_time.cycle_length = cycle
    # #
    # nc_dsi_time = nc.createVariable('dsi_time', np.float64, ('dsi_time',))
    # nc_dsi_time.long_name = 'time for silicon content of the Diatoms'
    # nc_dsi_time.units = 'day'
    # nc_dsi_time.calendar = 'XXX days in every year'
    # nc_dsi_time.cycle_length = cycle
    # #
    # nc_nfe_time = nc.createVariable('nfe_time', np.float64, ('nfe_time',))
    # nc_nfe_time.long_name = 'time for nanophytoplankton iron'
    # nc_nfe_time.units = 'day'
    # nc_nfe_time.calendar = 'XXX days in every year'
    # nc_nfe_time.cycle_length = cycle
    # #
    # nc_nch_time = nc.createVariable('nch_time', np.float64, ('nch_time',))
    # nc_nch_time.long_name = 'time for nanophytoplankton chlorophyll'
    # nc_nch_time.units = 'day'
    # nc_nch_time.calendar = 'XXX days in every year'
    # nc_nch_time.cycle_length = cycle
    # #
    # nc_dch_time = nc.createVariable('dch_time', np.float64, ('dch_time',))
    # nc_dch_time.long_name = 'time for diatom chlorophyll'
    # nc_dch_time.units = 'day'
    # nc_dch_time.calendar = 'XXX days in every year'
    # nc_dch_time.cycle_length = cycle
    #

    if obc[0] == 1:
        #
        #   Southern boundary
        #
        nc_temp_south = nc.createVariable('temp_south', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_south.long_name = 'southern boundary potential temperature'
        nc_temp_south.units = 'Celsius'
        nc_temp_south.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_south = nc.createVariable('salt_south', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_south.long_name = 'southern boundary salinity'
        nc_salt_south.units = 'PSU'
        nc_salt_south.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_south = nc.createVariable('u_south', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_south.long_name = 'southern boundary u-momentum component'
        nc_u_south.units = 'meter second-1'
        nc_u_south.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_south = nc.createVariable('v_south', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_south.long_name = 'southern boundary v-momentum component'
        nc_v_south.units = 'meter second-1'
        nc_v_south.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_south = nc.createVariable('ubar_south', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_south.long_name = 'southern boundary vertically integrated u-momentum component'
        nc_ubar_south.units = 'meter second-1'
        nc_ubar_south.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_south = nc.createVariable('vbar_south', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_south.long_name = 'southern boundary vertically integrated v-momentum component'
        nc_vbar_south.units = 'meter second-1'
        nc_vbar_south.coordinates = 'lon_v vclm_time'
        #
        nc_zeta_south = nc.createVariable('zeta_south', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_south.long_name = 'southern boundary sea surface height'
        nc_zeta_south.units = 'meter'
        nc_zeta_south.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_south = nc.createVariable('NO3_south', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_south.long_name = 'southern boundary nitrate'
        nc_no3_south.units = 'mmol/m3'
        nc_no3_south.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_south = nc.createVariable('PO4_south', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_south.long_name = 'southern boundary orthophosphate'
        nc_po4_south.units = 'mmol/m3'
        nc_po4_south.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_south = nc.createVariable('Si_south', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_south.long_name = 'southern boundary silicate'
        nc_si_south.units = 'mmol/m3'
        nc_si_south.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_south = nc.createVariable('O2_south', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_south.long_name = 'southern boundary dissolved oxygen'
        nc_o2_south.units = 'mmol/m3'
        nc_o2_south.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_south = nc.createVariable('DIC_south', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_south.long_name = 'southern boundary dissolved inorganic carbon'
        nc_dic_south.units = 'mmol/m3'
        nc_dic_south.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_nh4_south = nc.createVariable('NH4_south', np.float64, ('nh4_time', 's_rho', 'xi_rho',))
        nc_nh4_south.long_name = 'southern boundary ammonium'
        nc_nh4_south.units = 'mmol/m3'
        nc_nh4_south.coordinates = 'lon_rho s_rho nh4_time'
        #
        nc_fer_south = nc.createVariable('FER_south', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_south.long_name = 'southern boundary iron'
        nc_fer_south.units = 'mmol/m3'
        nc_fer_south.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_doc_south = nc.createVariable('DOC_south', np.float64, ('doc_time', 's_rho', 'xi_rho',))
        nc_doc_south.long_name = 'southern boundary dissolved organic carbon'
        nc_doc_south.units = 'mmol/m3'
        nc_doc_south.coordinates = 'lon_rho s_rho doc_time'
        #
        nc_talk_south = nc.createVariable('TALK_south', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_south.long_name = 'southern boundary alkalinity'
        nc_talk_south.units = 'mmol/m3'
        nc_talk_south.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_caco3_south = nc.createVariable('CACO3_south', np.float64, ('calc_time', 's_rho', 'xi_rho',))
        nc_caco3_south.long_name = 'southern boundary calcite'
        nc_caco3_south.units = 'mmol/m3'
        nc_caco3_south.coordinates = 'lon_rho s_rho calc_time'
        #
        nc_poc_south = nc.createVariable('POC_south', np.float64, ('poc_time', 's_rho', 'xi_rho',))
        nc_poc_south.long_name = 'southern boundary particulate organic carbon'
        nc_poc_south.units = 'mmol/m3'
        nc_poc_south.coordinates = 'lon_rho s_rho poc_time'
        #
        nc_phy_south = nc.createVariable('NANO_south', np.float64, ('phy_time', 's_rho', 'xi_rho',))
        nc_phy_south.long_name = 'southern boundary nanophytoplankton'
        nc_phy_south.units = 'mmol/m3'
        nc_phy_south.coordinates = 'lon_rho s_rho phy_time'
        #
        nc_zoo_south = nc.createVariable('ZOO_south', np.float64, ('zoo_time', 's_rho', 'xi_rho',))
        nc_zoo_south.long_name = 'southern boundary microzooplankton'
        nc_zoo_south.units = 'mmol/m3'
        nc_zoo_south.coordinates = 'lon_rho s_rho zoo_time'
        #
        nc_phy2_south = nc.createVariable('DIA_south', np.float64, ('dia_time', 's_rho', 'xi_rho',))
        nc_phy2_south.long_name = 'southern boundary diatom'
        nc_phy2_south.units = 'mmol/m3'
        nc_phy2_south.coordinates = 'lon_rho s_rho dia_time'
        #
        nc_zoo2_south = nc.createVariable('MESO_south', np.float64, ('mes_time', 's_rho', 'xi_rho',))
        nc_zoo2_south.long_name = 'southern boundary mesozooplankton'
        nc_zoo2_south.units = 'mmol/m3'
        nc_zoo2_south.coordinates = 'lon_rho s_rho mes_time'
        #
        nc_bsi_south = nc.createVariable('BSI_south', np.float64, ('bsi_time', 's_rho', 'xi_rho',))
        nc_bsi_south.long_name = 'southern boundary biogenic silica'
        nc_bsi_south.units = 'mmol/m3'
        nc_bsi_south.coordinates = 'lon_rho s_rho bsi_time'
        #
        nc_bfe_south = nc.createVariable('BFE_south', np.float64, ('bfe_time', 's_rho', 'xi_rho',))
        nc_bfe_south.long_name = 'southern boundary big particle iron'
        nc_bfe_south.units = 'mmol/m3'
        nc_bfe_south.coordinates = 'lon_rho s_rho bfe_time'
        #
        nc_goc_south = nc.createVariable('GOC_south', np.float64, ('goc_time', 's_rho', 'xi_rho',))
        nc_goc_south.long_name = 'southern boundary big particulate organic carbon'
        nc_goc_south.units = 'mmol/m3'
        nc_goc_south.coordinates = 'lon_rho s_rho goc_time'
        #
        nc_sfe_south = nc.createVariable('SFE_south', np.float64, ('sfe_time', 's_rho', 'xi_rho',))
        nc_sfe_south.long_name = 'southern boundary iron in the small particles'
        nc_sfe_south.units = 'mmol/m3'
        nc_sfe_south.coordinates = 'lon_rho s_rho sfe_time'
        #
        nc_dfe_south = nc.createVariable('DFE_south', np.float64, ('dfe_time', 's_rho', 'xi_rho',))
        nc_dfe_south.long_name = 'southern boundary iron content of the diatoms'
        nc_dfe_south.units = 'mmol/m3'
        nc_dfe_south.coordinates = 'lon_rho s_rho dfe_time'
        #
        nc_dsi_south = nc.createVariable('DSI_south', np.float64, ('dsi_time', 's_rho', 'xi_rho',))
        nc_dsi_south.long_name = 'southern boundary silicon content of the Diatoms'
        nc_dsi_south.units = 'mmol/m3'
        nc_dsi_south.coordinates = 'lon_rho s_rho dsi_time'
        #
        nc_nfe_south = nc.createVariable('NFE_south', np.float64, ('nfe_time', 's_rho', 'xi_rho',))
        nc_nfe_south.long_name = 'southern boundary nanophytoplankton iron'
        nc_nfe_south.units = 'mmol/m3'
        nc_nfe_south.coordinates = 'lon_rho s_rho nfe_time'
        #
        nc_nchl_south = nc.createVariable('NCHL_south', np.float64, ('nch_time', 's_rho', 'xi_rho',))
        nc_nchl_south.long_name = 'southern boundary nanophytoplankton chlorophyll'
        nc_nchl_south.units = 'mmol/m3'
        nc_nchl_south.coordinates = 'lon_rho s_rho nch_time'
        #
        nc_dchl_south = nc.createVariable('DCHL_south', np.float64, ('dch_time', 's_rho', 'xi_rho',))
        nc_dchl_south.long_name = 'southern boundary diatom chlorophyll'
        nc_dchl_south.units = 'mmol/m3'
        nc_dchl_south.coordinates = 'lon_rho s_rho dch_time'
        #
    if obc[1] == 1:
        #
        #   Eastern boundary
        #
        nc_temp_east = nc.createVariable('temp_east', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_east.long_name = 'eastern boundary potential temperature'
        nc_temp_east.units = 'Celsius'
        nc_temp_east.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_east = nc.createVariable('salt_east', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_east.long_name = 'eastern boundary salinity'
        nc_salt_east.units = 'PSU'
        nc_salt_east.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_east = nc.createVariable('u_east', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_east.long_name = 'eastern boundary u-momentum component'
        nc_u_east.units = 'meter second-1'
        nc_u_east.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_east = nc.createVariable('v_east', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_east.long_name = 'eastern boundary v-momentum component'
        nc_v_east.units = 'meter second-1'
        nc_v_east.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_east = nc.createVariable('ubar_east', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_east.long_name = 'eastern boundary vertically integrated u-momentum component'
        nc_ubar_east.units = 'meter second-1'
        nc_ubar_east.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_east = nc.createVariable('vbar_east', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_east.long_name = 'eastern boundary vertically integrated v-momentum component'
        nc_vbar_east.units = 'meter second-1'
        nc_vbar_east.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_east = nc.createVariable('zeta_east', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_east.long_name = 'eastern boundary sea surface height'
        nc_zeta_east.units = 'meter'
        nc_zeta_east.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_east = nc.createVariable('NO3_east', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_east.long_name = 'eastern boundary nitrate'
        nc_no3_east.units = 'mmol/m3'
        nc_no3_east.coordinates = 'lat_rho s_rho no3_time'
        #
        nc_po4_east = nc.createVariable('PO4_east', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_east.long_name = 'eastern boundary orthophosphate'
        nc_po4_east.units = 'mmol/m3'
        nc_po4_east.coordinates = 'lat_rho s_rho po4_time'
        #
        nc_si_east = nc.createVariable('Si_east', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_east.long_name = 'eastern boundary silicate'
        nc_si_east.units = 'mmol/m3'
        nc_si_east.coordinates = 'lat_rho s_rho si_time'
        #
        nc_o2_east = nc.createVariable('O2_east', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_east.long_name = 'eastern boundary dissolved oxygen'
        nc_o2_east.units = 'mmol/m3'
        nc_o2_east.coordinates = 'lat_rho s_rho o2_time'
        #
        nc_dic_east = nc.createVariable('DIC_east', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_east.long_name = 'eastern boundary dissolved inorganic carbon'
        nc_dic_east.units = 'mmol/m3'
        nc_dic_east.coordinates = 'lat_rho s_rho dic_time'
        #
        nc_nh4_east = nc.createVariable('NH4_east', np.float64, ('nh4_time', 's_rho', 'eta_rho',))
        nc_nh4_east.long_name = 'eastern boundary total alkalinity'
        nc_nh4_east.units = 'mmol/m3'
        nc_nh4_east.coordinates = 'lat_rho s_rho nh4_time'
        #
        nc_fer_east = nc.createVariable('FER_east', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_east.long_name = 'eastern boundary iron'
        nc_fer_east.units = 'mmol/m3'
        nc_fer_east.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_doc_east = nc.createVariable('DOC_east', np.float64, ('doc_time', 's_rho', 'eta_rho',))
        nc_doc_east.long_name = 'eastern boundary dissolved organic carbon'
        nc_doc_east.units = 'mmol/m3'
        nc_doc_east.coordinates = 'lon_rho s_rho doc_time'
        #
        nc_talk_east = nc.createVariable('TALK_east', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_east.long_name = 'eastern boundary alkalinity'
        nc_talk_east.units = 'mmol/m3'
        nc_talk_east.coordinates = 'lat_rho s_rho talk_time'
        #
        nc_caco3_east = nc.createVariable('CACO3_east', np.float64, ('calc_time', 's_rho', 'eta_rho',))
        nc_caco3_east.long_name = 'eastern boundary calcite'
        nc_caco3_east.units = 'mmol/m3'
        nc_caco3_east.coordinates = 'lat_rho s_rho calc_time'
        #
        nc_poc_east = nc.createVariable('POC_east', np.float64, ('poc_time', 's_rho', 'eta_rho',))
        nc_poc_east.long_name = 'eastern boundary particulate organic carbon'
        nc_poc_east.units = 'mmol/m3'
        nc_poc_east.coordinates = 'lat_rho s_rho poc_time'
        #
        nc_phy_east = nc.createVariable('NANO_east', np.float64, ('phy_time', 's_rho', 'eta_rho',))
        nc_phy_east.long_name = 'eastern boundary nanophytoplankton'
        nc_phy_east.units = 'mmol/m3'
        nc_phy_east.coordinates = 'lat_rho s_rho phy_time'
        #
        nc_zoo_east = nc.createVariable('ZOO_east', np.float64, ('zoo_time', 's_rho', 'eta_rho',))
        nc_zoo_east.long_name = 'eastern boundary microzooplankton'
        nc_zoo_east.units = 'mmol/m3'
        nc_zoo_east.coordinates = 'lat_rho s_rho zoo_time'
        #
        nc_phy2_east = nc.createVariable('DIA_east', np.float64, ('dia_time', 's_rho', 'eta_rho',))
        nc_phy2_east.long_name = 'eastern boundary diatom'
        nc_phy2_east.units = 'mmol/m3'
        nc_phy2_east.coordinates = 'lat_rho s_rho dia_time'
        #
        nc_zoo2_east = nc.createVariable('MESO_east', np.float64, ('mes_time', 's_rho', 'eta_rho',))
        nc_zoo2_east.long_name = 'eastern boundary mesozooplankton'
        nc_zoo2_east.units = 'mmol/m3'
        nc_zoo2_east.coordinates = 'lat_rho s_rho mes_time'
        #
        nc_bsi_east = nc.createVariable('BSI_east', np.float64, ('bsi_time', 's_rho', 'eta_rho',))
        nc_bsi_east.long_name = 'eastern boundary biogenic silica'
        nc_bsi_east.units = 'mmol/m3'
        nc_bsi_east.coordinates = 'lat_rho s_rho bsi_time'
        #
        nc_bfe_east = nc.createVariable('BFE_east', np.float64, ('bfe_time', 's_rho', 'eta_rho',))
        nc_bfe_east.long_name = 'eastern boundary big particle iron'
        nc_bfe_east.units = 'mmol/m3'
        nc_bfe_east.coordinates = 'lat_rho s_rho bfe_time'
        #
        nc_goc_east = nc.createVariable('GOC_east', np.float64, ('goc_time', 's_rho', 'eta_rho',))
        nc_goc_east.long_name = 'eastern boundary big particulate organic carbon'
        nc_goc_east.units = 'mmol/m3'
        nc_goc_east.coordinates = 'lat_rho s_rho goc_time'
        #
        nc_sfe_east = nc.createVariable('SFE_east', np.float64, ('sfe_time', 's_rho', 'eta_rho',))
        nc_sfe_east.long_name = 'eastern boundary iron in the small particles'
        nc_sfe_east.units = 'mmol/m3'
        nc_sfe_east.coordinates = 'lat_rho s_rho sfe_time'
        #
        nc_dfe_east = nc.createVariable('DFE_east', np.float64, ('dfe_time', 's_rho', 'eta_rho',))
        nc_dfe_east.long_name = 'eastern boundary iron content of the diatoms'
        nc_dfe_east.units = 'mmol/m3'
        nc_dfe_east.coordinates = 'lat_rho s_rho dfe_time'
        #
        nc_dsi_east = nc.createVariable('DSI_east', np.float64, ('dsi_time', 's_rho', 'eta_rho',))
        nc_dsi_east.long_name = 'eastern boundary silicon content of the Diatoms'
        nc_dsi_east.units = 'mmol/m3'
        nc_dsi_east.coordinates = 'lat_rho s_rho dsi_time'
        #
        nc_nfe_east = nc.createVariable('NFE_east', np.float64, ('nfe_time', 's_rho', 'eta_rho',))
        nc_nfe_east.long_name = 'eastern boundary nanophytoplankton iron'
        nc_nfe_east.units = 'mmol/m3'
        nc_nfe_east.coordinates = 'lat_rho s_rho nfe_time'
        #
        nc_nchl_east = nc.createVariable('NCHL_east', np.float64, ('nch_time', 's_rho', 'eta_rho',))
        nc_nchl_east.long_name = 'eastern boundary nanophytoplankton chlorophyll'
        nc_nchl_east.units = 'mmol/m3'
        nc_nchl_east.coordinates = 'lat_rho s_rho nch_time'
        #
        nc_dchl_east = nc.createVariable('DCHL_east', np.float64, ('dch_time', 's_rho', 'eta_rho',))
        nc_dchl_east.long_name = 'eastern boundary diatom chlorophyll'
        nc_dchl_east.units = 'mmol/m3'
        nc_dchl_east.coordinates = 'lat_rho s_rho dch_time'

    if obc[2] == 1:
        #
        #   Northern boundary
        #
        nc_temp_north = nc.createVariable('temp_north', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_north.long_name = 'northern boundary potential temperature'
        nc_temp_north.units = 'Celsius'
        nc_temp_north.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_north = nc.createVariable('salt_north', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_north.long_name = 'northern boundary salinity'
        nc_salt_north.units = 'PSU'
        nc_salt_north.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_north = nc.createVariable('u_north', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_north.long_name = 'northern boundary u-momentum component'
        nc_u_north.units = 'meter second-1'
        nc_u_north.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_north = nc.createVariable('v_north', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_north.long_name = 'northern boundary v-momentum component'
        nc_v_north.units = 'meter second-1'
        nc_v_north.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_north = nc.createVariable('ubar_north', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_north.long_name = 'northern boundary vertically integrated u-momentum component'
        nc_ubar_north.units = 'meter second-1'
        nc_ubar_north.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_north = nc.createVariable('vbar_north', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_north.long_name = 'northern boundary vertically integrated v-momentum component'
        nc_vbar_north.units = 'meter second-1'
        nc_vbar_north.coordinates = 'lon_v vclm_time'

        nc_zeta_north = nc.createVariable('zeta_north', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_north.long_name = 'northern boundary sea surface height'
        nc_zeta_north.units = 'meter'
        nc_zeta_north.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_north = nc.createVariable('NO3_north', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_north.long_name = 'northern boundary nitrate'
        nc_no3_north.units = 'mmol/m3'
        nc_no3_north.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_north = nc.createVariable('PO4_north', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_north.long_name = 'northern boundary orthophosphate'
        nc_po4_north.units = 'mmol/m3'
        nc_po4_north.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_north = nc.createVariable('Si_north', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_north.long_name = 'northern boundary silicate'
        nc_si_north.units = 'mmol/m3'
        nc_si_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_north = nc.createVariable('O2_north', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_north.long_name = 'northern boundary dissolved oxygen'
        nc_o2_north.units = 'mmol/m3'
        nc_o2_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_north = nc.createVariable('DIC_north', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_north.long_name = 'northern boundary dissolved inorganic carbon'
        nc_dic_north.units = 'mmol/m3'
        nc_dic_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_nh4_north = nc.createVariable('NH4_north', np.float64, ('nh4_time', 's_rho', 'xi_rho',))
        nc_nh4_north.long_name = 'northern boundary total alkalinity'
        nc_nh4_north.units = 'mmol/m3'
        nc_nh4_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_north = nc.createVariable('FER_north', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_north.long_name = 'northern boundary iron'
        nc_fer_north.units = 'mmol/m3'
        nc_fer_north.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_doc_north = nc.createVariable('DOC_north', np.float64, ('doc_time', 's_rho', 'xi_rho',))
        nc_doc_north.long_name = 'northern boundary dissolved organic carbon'
        nc_doc_north.units = 'mmol/m3'
        nc_doc_north.coordinates = 'lon_rho s_rho doc_time'
        #
        nc_talk_north = nc.createVariable('TALK_north', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_north.long_name = 'northern boundary alkalinity'
        nc_talk_north.units = 'mmol/m3'
        nc_talk_north.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_caco3_north = nc.createVariable('CACO3_north', np.float64, ('calc_time', 's_rho', 'xi_rho',))
        nc_caco3_north.long_name = 'northern boundary calcite'
        nc_caco3_north.units = 'mmol/m3'
        nc_caco3_north.coordinates = 'lon_rho s_rho calc_time'
        #
        nc_poc_north = nc.createVariable('POC_north', np.float64, ('poc_time', 's_rho', 'xi_rho',))
        nc_poc_north.long_name = 'northern boundary particulate organic carbon'
        nc_poc_north.units = 'mmol/m3'
        nc_poc_north.coordinates = 'lon_rho s_rho poc_time'
        #
        nc_phy_north = nc.createVariable('NANO_north', np.float64, ('phy_time', 's_rho', 'xi_rho',))
        nc_phy_north.long_name = 'northern boundary nanophytoplankton'
        nc_phy_north.units = 'mmol/m3'
        nc_phy_north.coordinates = 'lon_rho s_rho phy_time'
        #
        nc_zoo_north = nc.createVariable('ZOO_north', np.float64, ('zoo_time', 's_rho', 'xi_rho',))
        nc_zoo_north.long_name = 'northern boundary microzooplankton'
        nc_zoo_north.units = 'mmol/m3'
        nc_zoo_north.coordinates = 'lon_rho s_rho zoo_time'
        #
        nc_phy2_north = nc.createVariable('DIA_north', np.float64, ('dia_time', 's_rho', 'xi_rho',))
        nc_phy2_north.long_name = 'northern boundary diatom'
        nc_phy2_north.units = 'mmol/m3'
        nc_phy2_north.coordinates = 'lon_rho s_rho dia_time'
        #
        nc_zoo2_north = nc.createVariable('MESO_north', np.float64, ('mes_time', 's_rho', 'xi_rho',))
        nc_zoo2_north.long_name = 'northern boundary mesozooplankton'
        nc_zoo2_north.units = 'mmol/m3'
        nc_zoo2_north.coordinates = 'lon_rho s_rho mes_time'
        #
        nc_bsi_north = nc.createVariable('BSI_north', np.float64, ('bsi_time', 's_rho', 'xi_rho',))
        nc_bsi_north.long_name = 'northern boundary biogenic silica'
        nc_bsi_north.units = 'mmol/m3'
        nc_bsi_north.coordinates = 'lon_rho s_rho bsi_time'
        #
        nc_bfe_north = nc.createVariable('BFE_north', np.float64, ('bfe_time', 's_rho', 'xi_rho',))
        nc_bfe_north.long_name = 'northern boundary big particle iron'
        nc_bfe_north.units = 'mmol/m3'
        nc_bfe_north.coordinates = 'lon_rho s_rho bfe_time'
        #
        nc_goc_north = nc.createVariable('GOC_north', np.float64, ('goc_time', 's_rho', 'xi_rho',))
        nc_goc_north.long_name = 'northern boundary big particulate organic carbon'
        nc_goc_north.units = 'mmol/m3'
        nc_goc_north.coordinates = 'lon_rho s_rho goc_time'
        #
        nc_sfe_north = nc.createVariable('SFE_north', np.float64, ('sfe_time', 's_rho', 'xi_rho',))
        nc_sfe_north.long_name = 'northern boundary iron in the small particles'
        nc_sfe_north.units = 'mmol/m3'
        nc_sfe_north.coordinates = 'lon_rho s_rho sfe_time'
        #
        nc_dfe_north = nc.createVariable('DFE_north', np.float64, ('dfe_time', 's_rho', 'xi_rho',))
        nc_dfe_north.long_name = 'northern boundary iron content of the diatoms'
        nc_dfe_north.units = 'mmol/m3'
        nc_dfe_north.coordinates = 'lon_rho s_rho dfe_time'
        #
        nc_dsi_north = nc.createVariable('DSI_north', np.float64, ('dsi_time', 's_rho', 'xi_rho',))
        nc_dsi_north.long_name = 'northern boundary silicon content of the Diatoms'
        nc_dsi_north.units = 'mmol/m3'
        nc_dsi_north.coordinates = 'lon_rho s_rho dsi_time'
        #
        nc_nfe_north = nc.createVariable('NFE_north', np.float64, ('nfe_time', 's_rho', 'xi_rho',))
        nc_nfe_north.long_name = 'northern boundary nanophytoplankton iron'
        nc_nfe_north.units = 'mmol/m3'
        nc_nfe_north.coordinates = 'lon_rho s_rho nfe_time'
        #
        nc_nchl_north = nc.createVariable('NCHL_north', np.float64, ('nch_time', 's_rho', 'xi_rho',))
        nc_nchl_north.long_name = 'northern boundary nanophytoplankton chlorophyll'
        nc_nchl_north.units = 'mmol/m3'
        nc_nchl_north.coordinates = 'lon_rho s_rho nch_time'
        #
        nc_dchl_north = nc.createVariable('DCHL_north', np.float64, ('dch_time', 's_rho', 'xi_rho',))
        nc_dchl_north.long_name = 'northern boundary diatom chlorophyll'
        nc_dchl_north.units = 'mmol/m3'
        nc_dchl_north.coordinates = 'lon_rho s_rho dch_time'
        #

    if obc[3] == 1:
        #
        #   Western boundary
        #
        nc_temp_west = nc.createVariable('temp_west', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_west.long_name = 'western boundary potential temperature'
        nc_temp_west.units = 'Celsius'
        nc_temp_west.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_west = nc.createVariable('salt_west', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_west.long_name = 'western boundary salinity'
        nc_salt_west.units = 'PSU'
        nc_salt_west.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_west = nc.createVariable('u_west', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_west.long_name = 'western boundary u-momentum component'
        nc_u_west.units = 'meter second-1'
        nc_u_west.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_west = nc.createVariable('v_west', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_west.long_name = 'western boundary v-momentum component'
        nc_v_west.units = 'meter second-1'
        nc_v_west.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_west = nc.createVariable('ubar_west', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_west.long_name = 'western boundary vertically integrated u-momentum component'
        nc_ubar_west.units = 'meter second-1'
        nc_ubar_west.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_west = nc.createVariable('vbar_west', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_west.long_name = 'western boundary vertically integrated v-momentum component'
        nc_vbar_west.units = 'meter second-1'
        nc_vbar_west.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_west = nc.createVariable('zeta_west', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_west.long_name = 'western boundary sea surface height'
        nc_zeta_west.units = 'meter'
        nc_zeta_west.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_west = nc.createVariable('NO3_west', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_west.long_name = 'western boundary nitrate'
        nc_no3_west.units = 'mmol/m3'
        nc_no3_west.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_west = nc.createVariable('PO4_west', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_west.long_name = 'western boundary orthophosphate'
        nc_po4_west.units = 'mmol/m3'
        nc_po4_west.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_west = nc.createVariable('Si_west', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_west.long_name = 'western boundary silicate'
        nc_si_west.units = 'mmol/m3'
        nc_si_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_west = nc.createVariable('O2_west', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_west.long_name = 'western boundary dissolved oxygen'
        nc_o2_west.units = 'mmol/m3'
        nc_o2_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_west = nc.createVariable('DIC_west', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_west.long_name = 'western boundary dissolved inorganic carbon'
        nc_dic_west.units = 'mmol/m3'
        nc_dic_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_nh4_west = nc.createVariable('NH4_west', np.float64, ('nh4_time', 's_rho', 'eta_rho',))
        nc_nh4_west.long_name = 'western boundary total alkalinity'
        nc_nh4_west.units = 'mmol/m3'
        nc_nh4_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_west = nc.createVariable('FER_west', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_west.long_name = 'western boundary iron'
        nc_fer_west.units = 'mmol/m3'
        nc_fer_west.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_doc_west = nc.createVariable('DOC_west', np.float64, ('doc_time', 's_rho', 'eta_rho',))
        nc_doc_west.long_name = 'western boundary dissolved organic carbon'
        nc_doc_west.units = 'mmol/m3'
        nc_doc_west.coordinates = 'lon_rho s_rho doc_time'
        #
        nc_talk_west = nc.createVariable('TALK_west', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_west.long_name = 'western boundary alkalinity'
        nc_talk_west.units = 'mmol/m3'
        nc_talk_west.coordinates = 'lat_rho s_rho talk_time'
        #
        nc_caco3_west = nc.createVariable('CACO3_west', np.float64, ('calc_time', 's_rho', 'eta_rho',))
        nc_caco3_west.long_name = 'western boundary calcite'
        nc_caco3_west.units = 'mmol/m3'
        nc_caco3_west.coordinates = 'lat_rho s_rho calc_time'
        #
        nc_poc_west = nc.createVariable('POC_west', np.float64, ('poc_time', 's_rho', 'eta_rho',))
        nc_poc_west.long_name = 'western boundary particulate organic carbon'
        nc_poc_west.units = 'mmol/m3'
        nc_poc_west.coordinates = 'lat_rho s_rho poc_time'
        #
        nc_phy_west = nc.createVariable('NANO_west', np.float64, ('phy_time', 's_rho', 'eta_rho',))
        nc_phy_west.long_name = 'western boundary nanophytoplankton'
        nc_phy_west.units = 'mmol/m3'
        nc_phy_west.coordinates = 'lat_rho s_rho phy_time'
        #
        nc_zoo_west = nc.createVariable('ZOO_west', np.float64, ('zoo_time', 's_rho', 'eta_rho',))
        nc_zoo_west.long_name = 'western boundary microzooplankton'
        nc_zoo_west.units = 'mmol/m3'
        nc_zoo_west.coordinates = 'lat_rho s_rho zoo_time'
        #
        nc_phy2_west = nc.createVariable('DIA_west', np.float64, ('dia_time', 's_rho', 'eta_rho',))
        nc_phy2_west.long_name = 'western boundary diatom'
        nc_phy2_west.units = 'mmol/m3'
        nc_phy2_west.coordinates = 'lat_rho s_rho dia_time'
        #
        nc_zoo2_west = nc.createVariable('MESO_west', np.float64, ('mes_time', 's_rho', 'eta_rho',))
        nc_zoo2_west.long_name = 'western boundary mesozooplankton'
        nc_zoo2_west.units = 'mmol/m3'
        nc_zoo2_west.coordinates = 'lat_rho s_rho mes_time'
        #
        nc_bsi_west = nc.createVariable('BSI_west', np.float64, ('bsi_time', 's_rho', 'eta_rho',))
        nc_bsi_west.long_name = 'western boundary biogenic silica'
        nc_bsi_west.units = 'mmol/m3'
        nc_bsi_west.coordinates = 'lat_rho s_rho bsi_time'
        #
        nc_bfe_west = nc.createVariable('BFE_west', np.float64, ('bfe_time', 's_rho', 'eta_rho',))
        nc_bfe_west.long_name = 'western boundary big particle iron'
        nc_bfe_west.units = 'mmol/m3'
        nc_bfe_west.coordinates = 'lat_rho s_rho bfe_time'
        #
        nc_goc_west = nc.createVariable('GOC_west', np.float64, ('goc_time', 's_rho', 'eta_rho',))
        nc_goc_west.long_name = 'western boundary big particulate organic carbon'
        nc_goc_west.units = 'mmol/m3'
        nc_goc_west.coordinates = 'lat_rho s_rho goc_time'
        #
        nc_sfe_west = nc.createVariable('SFE_west', np.float64, ('sfe_time', 's_rho', 'eta_rho',))
        nc_sfe_west.long_name = 'western boundary iron in the small particles'
        nc_sfe_west.units = 'mmol/m3'
        nc_sfe_west.coordinates = 'lat_rho s_rho sfe_time'
        #
        nc_dfe_west = nc.createVariable('DFE_west', np.float64, ('dfe_time', 's_rho', 'eta_rho',))
        nc_dfe_west.long_name = 'western boundary iron content of the diatoms'
        nc_dfe_west.units = 'mmol/m3'
        nc_dfe_west.coordinates = 'lat_rho s_rho dfe_time'
        #
        nc_dsi_west = nc.createVariable('DSI_west', np.float64, ('dsi_time', 's_rho', 'eta_rho',))
        nc_dsi_west.long_name = 'western boundary silicon content of the Diatoms'
        nc_dsi_west.units = 'mmol/m3'
        nc_dsi_west.coordinates = 'lat_rho s_rho dsi_time'
        #
        nc_nfe_west = nc.createVariable('NFE_west', np.float64, ('nfe_time', 's_rho', 'eta_rho',))
        nc_nfe_west.long_name = 'western boundary nanophytoplankton iron'
        nc_nfe_west.units = 'mmol/m3'
        nc_nfe_west.coordinates = 'lat_rho s_rho nfe_time'
        #
        nc_nchl_west = nc.createVariable('NCHL_west', np.float64, ('nch_time', 's_rho', 'eta_rho',))
        nc_nchl_west.long_name = 'western boundary nanophytoplankton chlorophyll'
        nc_nchl_west.units = 'mmol/m3'
        nc_nchl_west.coordinates = 'lat_rho s_rho nch_time'
        #
        nc_dchl_west = nc.createVariable('DCHL_west', np.float64, ('dch_time', 's_rho', 'eta_rho',))
        nc_dchl_west.long_name = 'western boundary diatom chlorophyll'
        nc_dchl_west.units = 'mmol/m3'
        nc_dchl_west.coordinates = 'lat_rho s_rho dch_time'
        #
    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #
    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(bry_time)
    nc_tend[:] = np.max(bry_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_sc_w[:] = sc_w
    nc_Cs_r[:] = Cs_r
    nc_Cs_w[:] = Cs_w
    nc_tclm_time[:] = bry_time
    nc_temp_time[:] = bry_time
    nc_sclm_time[:] = bry_time
    nc_salt_time[:] = bry_time
    nc_uclm_time[:] = bry_time
    nc_vclm_time[:] = bry_time
    nc_v2d_time[:] = bry_time
    nc_v3d_time[:] = bry_time
    nc_ssh_time[:] = bry_time
    nc_zeta_time[:] = bry_time
    nc_bry_time[:] = bry_time

    nc_no3_time[:] = nor_time
    nc_po4_time[:] = nor_time
    nc_si_time[:] = nor_time
    nc_o2_time[:] = nor_time
    nc_dic_time[:] = nor_time
    nc_talk_time[:] = nor_time

    # nc_nh4_time[:] = pis_time
    nc_fer_time[:] = pis_time
    # nc_calc_time[:] = pis_time
    # nc_poc_time[:] = pis_time
    # nc_phy_time[:] = pis_time
    # nc_zoo_time[:] = pis_time
    # nc_doc_time[:] = pis_time
    # nc_dia_time[:] = pis_time
    # nc_mes_time[:] = pis_time
    # nc_bsi_time[:] = pis_time
    # nc_bfe_time[:] = pis_time
    # nc_goc_time[:] = pis_time
    # nc_sfe_time[:] = pis_time
    # nc_dfe_time[:] = pis_time
    # nc_dsi_time[:] = pis_time
    # nc_nfe_time[:] = pis_time
    # nc_nch_time[:] = pis_time
    # nc_dch_time[:] = pis_time

    if obc[0] == 1:
        nc_u_south[:] = 0.
        nc_v_south[:] = 0.
        nc_ubar_south[:] = 0.
        nc_vbar_south[:] = 0.
        nc_zeta_south[:] = 0.
        nc_temp_south[:] = 0.
        nc_salt_south[:] = 0.
        nc_no3_south[:] = 0.
        nc_po4_south[:] = 0.
        nc_si_south[:] = 0.
        nc_o2_south[:] = 0.
        nc_dic_south[:] = 0.
        nc_nh4_south[:] = 0.
        nc_fer_south[:] = 0.
        nc_doc_south[:] = 0.
        nc_talk_south[:] = 0.
        nc_caco3_south[:] = 0.
        nc_poc_south[:] = 0.
        nc_phy_south[:] = 0.
        nc_zoo_south[:] = 0.
        nc_phy2_south[:] = 0.
        nc_zoo2_south[:] = 0.
        nc_bsi_south[:] = 0.
        nc_bfe_south[:] = 0.
        nc_goc_south[:] = 0.
        nc_sfe_south[:] = 0.
        nc_dfe_south[:] = 0.
        nc_dsi_south[:] = 0.
        nc_nfe_south[:] = 0.
        nc_nchl_south[:] = 0.
        nc_dchl_south[:] = 0.

    if obc[1] == 1:
        nc_u_east[:] = 0.
        nc_v_east[:] = 0.
        nc_ubar_east[:] = 0.
        nc_vbar_east[:] = 0.
        nc_zeta_east[:] = 0.
        nc_temp_east[:] = 0.
        nc_salt_east[:] = 0.
        nc_no3_east[:] = 0.
        nc_po4_east[:] = 0.
        nc_si_east[:] = 0.
        nc_o2_east[:] = 0.
        nc_dic_east[:] = 0.
        nc_nh4_east[:] = 0.
        nc_fer_east[:] = 0.
        nc_doc_east[:] = 0.
        nc_talk_east[:] = 0.
        nc_caco3_east[:] = 0.
        nc_poc_east[:] = 0.
        nc_phy_east[:] = 0.
        nc_zoo_east[:] = 0.
        nc_phy2_east[:] = 0.
        nc_zoo2_east[:] = 0.
        nc_bsi_east[:] = 0.
        nc_bfe_east[:] = 0.
        nc_goc_east[:] = 0.
        nc_sfe_east[:] = 0.
        nc_dfe_east[:] = 0.
        nc_dsi_east[:] = 0.
        nc_nfe_east[:] = 0.
        nc_nchl_east[:] = 0.
        nc_dchl_east[:] = 0.

    if obc[2] == 1:
        nc_u_north[:] = 0.
        nc_v_north[:] = 0.
        nc_ubar_north[:] = 0.
        nc_vbar_north[:] = 0.
        nc_zeta_north[:] = 0.
        nc_temp_north[:] = 0.
        nc_salt_north[:] = 0.
        nc_no3_north[:] = 0.
        nc_po4_north[:] = 0.
        nc_si_north[:] = 0.
        nc_o2_north[:] = 0.
        nc_dic_north[:] = 0.
        nc_nh4_north[:] = 0.
        nc_fer_north[:] = 0.
        nc_doc_north[:] = 0.
        nc_talk_north[:] = 0.
        nc_caco3_north[:] = 0.
        nc_poc_north[:] = 0.
        nc_phy_north[:] = 0.
        nc_zoo_north[:] = 0.
        nc_phy2_north[:] = 0.
        nc_zoo2_north[:] = 0.
        nc_bsi_north[:] = 0.
        nc_bfe_north[:] = 0.
        nc_goc_north[:] = 0.
        nc_sfe_north[:] = 0.
        nc_dfe_north[:] = 0.
        nc_dsi_north[:] = 0.
        nc_nfe_north[:] = 0.
        nc_nchl_north[:] = 0.
        nc_dchl_north[:] = 0.

    if obc[3] == 1:
        nc_u_west[:] = 0.
        nc_v_west[:] = 0.
        nc_ubar_west[:] = 0.
        nc_vbar_west[:] = 0.
        nc_zeta_west[:] = 0.
        nc_temp_west[:] = 0.
        nc_salt_west[:] = 0.
        nc_no3_west[:] = 0.
        nc_po4_west[:] = 0.
        nc_si_west[:] = 0.
        nc_o2_west[:] = 0.
        nc_dic_west[:] = 0.
        nc_nh4_west[:] = 0.
        nc_fer_west[:] = 0.
        nc_doc_west[:] = 0.
        nc_talk_west[:] = 0.
        nc_caco3_west[:] = 0.
        nc_poc_west[:] = 0.
        nc_phy_west[:] = 0.
        nc_zoo_west[:] = 0.
        nc_phy2_west[:] = 0.
        nc_zoo2_west[:] = 0.
        nc_bsi_west[:] = 0.
        nc_bfe_west[:] = 0.
        nc_goc_west[:] = 0.
        nc_sfe_west[:] = 0.
        nc_dfe_west[:] = 0.
        nc_dsi_west[:] = 0.
        nc_nfe_west[:] = 0.
        nc_nchl_west[:] = 0.
        nc_dchl_west[:] = 0.

    nc.close()

    return


def create_bryfile_BGC_PHY(bryname, grdname, title, obc,
                           theta_s, theta_b, hc, N,
                           bry_time, nor_time, pis_time, cycle, vtransform):
    #
    #
    ################################################################
    #
    # function create_bryfile(bryname,grdname,title,obc...
    #                          theta_s,theta_b,hc,N,...
    #                          bry_time,cycle,clobber)
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   bryname      Netcdf climatology file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   obc          open boundaries flag (1=open , [S E N W]).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   bry_time     time.(vector)
    #   cycle        Length (days) for cycling the climatology.(Real)
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + bryname)
    print(' ')
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()

    #
    #
    #

    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + ' m)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    Nt = 0

    #
    #  Create the boundary file
    #
    type = 'BOUNDARY file'
    history = 'CROCO'

    if os.path.exists(bryname):
        os.remove(bryname)

    print('Create: ' + bryname)
    nc = netcdf(bryname, 'w', format='NETCDF4')

    #
    # set global attributes
    #

    nc.type = 'CROCO boundary file'
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.title = title
    nc.bry_file = bryname
    nc.grd_file = grdname

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)

    # PHYSICS time dimension
    nc_dim_bry_time = nc.createDimension('bry_time', Nt)
    nc_dim_tmpm_time = nc.createDimension('tmpm_time', Nt)
    nc_dim_temp_time = nc.createDimension('temp_time', Nt)
    nc_dim_salm_time = nc.createDimension('salm_time', Nt)
    nc_dim_salt_time = nc.createDimension('salt_time', Nt)
    nc_dim_uclm_time = nc.createDimension('uclm_time', Nt)
    nc_dim_vclm_time = nc.createDimension('vclm_time', Nt)
    nc_dim_v2d_time = nc.createDimension('v2d_time', Nt)
    nc_dim_v3d_time = nc.createDimension('v3d_time', Nt)
    nc_dim_ssh_time = nc.createDimension('ssh_time', Nt)
    nc_dim_zeta_time = nc.createDimension('zeta_time', Nt)

    # # BGC time dimension
    nc_dim_no3_time = nc.createDimension('no3_time', Nt)
    nc_dim_po4_time = nc.createDimension('po4_time', Nt)
    nc_dim_si_time = nc.createDimension('si_time', Nt)
    nc_dim_o2_time = nc.createDimension('o2_time', Nt)
    nc_dim_dic_time = nc.createDimension('dic_time', Nt)
    nc_dim_nh4_time = nc.createDimension('nh4_time', Nt)
    nc_dim_fer_time = nc.createDimension('fer_time', Nt)
    nc_dim_talk_time = nc.createDimension('talk_time', Nt)
    nc_dim_ph_time = nc.createDimension('ph_time', Nt)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #
    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    # s
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_bry_time = nc.createVariable('bry_time', np.float64, ('bry_time',))
    nc_bry_time.long_name = 'time for boundary climatology'
    nc_bry_time.units = 'day'
    nc_bry_time.calendar = 'XXX days in every year'
    nc_bry_time.cycle_length = cycle
    #
    nc_tmpm_time = nc.createVariable('tmpm_time', np.float64, ('tmpm_time',))
    nc_tmpm_time.long_name = 'time for temperature monthly'
    nc_tmpm_time.units = 'day'
    nc_tmpm_time.calendar = 'XXX days in every year'
    nc_tmpm_time.cycle_length = cycle
    #
    nc_temp_time = nc.createVariable('temp_time', np.float64, ('temp_time',))
    nc_temp_time.long_name = 'time for temperature'
    nc_temp_time.units = 'day'
    nc_temp_time.calendar = 'XXX days in every year'
    nc_temp_time.cycle_length = cycle
    #
    nc_salm_time = nc.createVariable('salm_time', np.float64, ('salm_time',))
    nc_salm_time.long_name = 'time for salinity monthly'
    nc_salm_time.units = 'day'
    nc_salm_time.calendar = 'XXX days in every year'
    nc_salm_time.cycle_length = cycle
    #
    nc_salt_time = nc.createVariable('salt_time', np.float64, ('salt_time',))
    nc_salt_time.long_name = 'time for salinity'
    nc_salt_time.units = 'day'
    nc_salt_time.calendar = 'XXX days in every year'
    nc_salt_time.cycle_length = cycle
    #
    nc_uclm_time = nc.createVariable('uclm_time', np.float64, ('uclm_time',))
    nc_uclm_time.long_name = 'time climatological u'
    nc_uclm_time.units = 'day'
    nc_uclm_time.calendar = 'XXX days in every year'
    nc_uclm_time.cycle_length = cycle
    #
    nc_vclm_time = nc.createVariable('vclm_time', np.float64, ('vclm_time',))
    nc_vclm_time.long_name = 'time climatological v'
    nc_vclm_time.units = 'day'
    nc_vclm_time.calendar = 'XXX days in every year'
    nc_vclm_time.cycle_length = cycle
    #
    nc_v2d_time = nc.createVariable('v2d_time', np.float64, ('v2d_time',))
    nc_v2d_time.long_name = 'time for 2D velocity'
    nc_v2d_time.units = 'day'
    nc_v2d_time.calendar = 'XXX days in every year'
    nc_v2d_time.cycle_length = cycle
    #
    nc_v3d_time = nc.createVariable('v3d_time', np.float64, ('v3d_time',))
    nc_v3d_time.long_name = 'time for 3D velocity'
    nc_v3d_time.units = 'day'
    nc_v3d_time.calendar = 'XXX days in every year'
    nc_v3d_time.cycle_length = cycle
    #
    nc_ssh_time = nc.createVariable('ssh_time', np.float64, ('ssh_time',))
    nc_ssh_time.long_name = 'time for sea surface height'
    nc_ssh_time.units = 'day'
    nc_ssh_time.calendar = 'XXX days in every year'
    nc_ssh_time.cycle_length = cycle
    #
    nc_zeta_time = nc.createVariable('zeta_time', np.float64, ('zeta_time',))
    nc_zeta_time.long_name = 'time for sea surface height'
    nc_zeta_time.units = 'day'
    nc_zeta_time.calendar = 'XXX days in every year'
    nc_zeta_time.cycle_length = cycle
    #
    nc_no3_time = nc.createVariable('no3_time', np.float64, ('no3_time',))
    nc_no3_time.long_name = 'time for nitrate'
    nc_no3_time.units = 'day'
    nc_no3_time.calendar = 'XXX days in every year'
    nc_no3_time.cycle_length = cycle
    #
    nc_po4_time = nc.createVariable('po4_time', np.float64, ('po4_time',))
    nc_po4_time.long_name = 'time for orthophosphate'
    nc_po4_time.units = 'day'
    nc_po4_time.calendar = 'XXX days in every year'
    nc_po4_time.cycle_length = cycle
    #
    nc_si_time = nc.createVariable('si_time', np.float64, ('si_time',))
    nc_si_time.long_name = 'time for silicate'
    nc_si_time.units = 'day'
    nc_si_time.calendar = 'XXX days in every year'
    nc_si_time.cycle_length = cycle
    #
    nc_o2_time = nc.createVariable('o2_time', np.float64, ('o2_time',))
    nc_o2_time.long_name = 'time for dissolved oxygen'
    nc_o2_time.units = 'day'
    nc_o2_time.calendar = 'XXX days in every year'
    nc_o2_time.cycle_length = cycle
    #
    nc_dic_time = nc.createVariable('dic_time', np.float64, ('dic_time',))
    nc_dic_time.long_name = 'time for dissolved inorganic carbon'
    nc_dic_time.units = 'day'
    nc_dic_time.calendar = 'XXX days in every year'
    nc_dic_time.cycle_length = cycle
    #
    nc_nh4_time = nc.createVariable('nh4_time', np.float64, ('nh4_time',))
    nc_nh4_time.long_name = 'time for ammonium'
    nc_nh4_time.units = 'day'
    nc_nh4_time.calendar = 'XXX days in every year'
    nc_nh4_time.cycle_length = cycle
    #
    nc_fer_time = nc.createVariable('fer_time', np.float64, ('fer_time',))
    nc_fer_time.long_name = 'time for iron'
    nc_fer_time.units = 'day'
    nc_fer_time.calendar = 'XXX days in every year'
    nc_fer_time.cycle_length = cycle
    #
    nc_talk_time = nc.createVariable('talk_time', np.float64, ('talk_time',))
    nc_talk_time.long_name = 'time for alkalinity'
    nc_talk_time.units = 'day'
    nc_talk_time.calendar = 'XXX days in every year'
    nc_talk_time.cycle_length = cycle
    #
    nc_ph_time = nc.createVariable('ph_time', np.float64, ('ph_time',))
    nc_ph_time.long_name = 'time for ph'
    nc_ph_time.units = 'day'
    nc_ph_time.calendar = 'XXX days in every year'
    nc_ph_time.cycle_length = cycle
    #

    if obc[0] == 1:
        #
        #   Southern boundary
        #
        nc_temp_south = nc.createVariable('temp_south', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_south.long_name = 'southern boundary potential temperature'
        nc_temp_south.units = 'Celsius'
        nc_temp_south.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_south = nc.createVariable('salt_south', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_south.long_name = 'southern boundary salinity'
        nc_salt_south.units = 'PSU'
        nc_salt_south.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_south = nc.createVariable('u_south', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_south.long_name = 'southern boundary u-momentum component'
        nc_u_south.units = 'meter second-1'
        nc_u_south.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_south = nc.createVariable('v_south', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_south.long_name = 'southern boundary v-momentum component'
        nc_v_south.units = 'meter second-1'
        nc_v_south.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_south = nc.createVariable('ubar_south', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_south.long_name = 'southern boundary vertically integrated u-momentum component'
        nc_ubar_south.units = 'meter second-1'
        nc_ubar_south.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_south = nc.createVariable('vbar_south', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_south.long_name = 'southern boundary vertically integrated v-momentum component'
        nc_vbar_south.units = 'meter second-1'
        nc_vbar_south.coordinates = 'lon_v vclm_time'
        #
        nc_zeta_south = nc.createVariable('zeta_south', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_south.long_name = 'southern boundary sea surface height'
        nc_zeta_south.units = 'meter'
        nc_zeta_south.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_south = nc.createVariable('NO3_south', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_south.long_name = 'southern boundary nitrate'
        nc_no3_south.units = 'mmol/m3'
        nc_no3_south.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_south = nc.createVariable('PO4_south', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_south.long_name = 'southern boundary orthophosphate'
        nc_po4_south.units = 'mmol/m3'
        nc_po4_south.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_south = nc.createVariable('Si_south', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_south.long_name = 'southern boundary silicate'
        nc_si_south.units = 'mmol/m3'
        nc_si_south.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_south = nc.createVariable('O2_south', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_south.long_name = 'southern boundary dissolved oxygen'
        nc_o2_south.units = 'mmol/m3'
        nc_o2_south.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_south = nc.createVariable('DIC_south', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_south.long_name = 'southern boundary dissolved inorganic carbon'
        nc_dic_south.units = 'mmol/m3'
        nc_dic_south.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_nh4_south = nc.createVariable('NH4_south', np.float64, ('nh4_time', 's_rho', 'xi_rho',))
        nc_nh4_south.long_name = 'southern boundary ammonium'
        nc_nh4_south.units = 'mmol/m3'
        nc_nh4_south.coordinates = 'lon_rho s_rho nh4_time'
        #
        nc_fer_south = nc.createVariable('FER_south', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_south.long_name = 'southern boundary iron'
        nc_fer_south.units = 'mmol/m3'
        nc_fer_south.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_talk_south = nc.createVariable('TALK_south', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_south.long_name = 'southern boundary alkalinity'
        nc_talk_south.units = 'mmol/m3'
        nc_talk_south.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_ph_south = nc.createVariable('PH_south', np.float64, ('ph_time', 's_rho', 'xi_rho',))
        nc_ph_south.long_name = 'southern boundary pH'
        nc_ph_south.units = '-'
        nc_ph_south.coordinates = 'lon_rho s_rho ph_time'
        #
        nc_tmpm_south = nc.createVariable('tmpm_south', np.float64, ('tmpm_time', 's_rho', 'xi_rho',))
        nc_tmpm_south.long_name = 'southern boundary potential temperature monthly'
        nc_tmpm_south.units = 'Celsius'
        nc_tmpm_south.coordinates = 'lon_rho s_rho tmpm_time'
        #
        nc_salm_south = nc.createVariable('salm_south', np.float64, ('salm_time', 's_rho', 'xi_rho',))
        nc_salm_south.long_name = 'southern boundary salinity monthly'
        nc_salm_south.units = 'PSU'
        nc_salm_south.coordinates = 'lon_rho s_rho salm_time'
        #
    if obc[1] == 1:
        #
        #   Eastern boundary
        #
        nc_temp_east = nc.createVariable('temp_east', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_east.long_name = 'eastern boundary potential temperature'
        nc_temp_east.units = 'Celsius'
        nc_temp_east.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_east = nc.createVariable('salt_east', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_east.long_name = 'eastern boundary salinity'
        nc_salt_east.units = 'PSU'
        nc_salt_east.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_east = nc.createVariable('u_east', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_east.long_name = 'eastern boundary u-momentum component'
        nc_u_east.units = 'meter second-1'
        nc_u_east.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_east = nc.createVariable('v_east', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_east.long_name = 'eastern boundary v-momentum component'
        nc_v_east.units = 'meter second-1'
        nc_v_east.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_east = nc.createVariable('ubar_east', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_east.long_name = 'eastern boundary vertically integrated u-momentum component'
        nc_ubar_east.units = 'meter second-1'
        nc_ubar_east.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_east = nc.createVariable('vbar_east', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_east.long_name = 'eastern boundary vertically integrated v-momentum component'
        nc_vbar_east.units = 'meter second-1'
        nc_vbar_east.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_east = nc.createVariable('zeta_east', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_east.long_name = 'eastern boundary sea surface height'
        nc_zeta_east.units = 'meter'
        nc_zeta_east.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_east = nc.createVariable('NO3_east', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_east.long_name = 'eastern boundary nitrate'
        nc_no3_east.units = 'mmol/m3'
        nc_no3_east.coordinates = 'lat_rho s_rho no3_time'
        #
        nc_po4_east = nc.createVariable('PO4_east', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_east.long_name = 'eastern boundary orthophosphate'
        nc_po4_east.units = 'mmol/m3'
        nc_po4_east.coordinates = 'lat_rho s_rho po4_time'
        #
        nc_si_east = nc.createVariable('Si_east', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_east.long_name = 'eastern boundary silicate'
        nc_si_east.units = 'mmol/m3'
        nc_si_east.coordinates = 'lat_rho s_rho si_time'
        #
        nc_o2_east = nc.createVariable('O2_east', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_east.long_name = 'eastern boundary dissolved oxygen'
        nc_o2_east.units = 'mmol/m3'
        nc_o2_east.coordinates = 'lat_rho s_rho o2_time'
        #
        nc_dic_east = nc.createVariable('DIC_east', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_east.long_name = 'eastern boundary dissolved inorganic carbon'
        nc_dic_east.units = 'mmol/m3'
        nc_dic_east.coordinates = 'lat_rho s_rho dic_time'
        #
        nc_nh4_east = nc.createVariable('NH4_east', np.float64, ('nh4_time', 's_rho', 'eta_rho',))
        nc_nh4_east.long_name = 'eastern boundary total alkalinity'
        nc_nh4_east.units = 'mmol/m3'
        nc_nh4_east.coordinates = 'lat_rho s_rho nh4_time'
        #
        nc_fer_east = nc.createVariable('FER_east', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_east.long_name = 'eastern boundary iron'
        nc_fer_east.units = 'mmol/m3'
        nc_fer_east.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_talk_east = nc.createVariable('TALK_east', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_east.long_name = 'eastern boundary alkalinity'
        nc_talk_east.units = 'mmol/m3'
        nc_talk_east.coordinates = 'lat_rho s_rho talk_time'
        #
        nc_ph_east = nc.createVariable('PH_east', np.float64, ('ph_time', 's_rho', 'eta_rho',))
        nc_ph_east.long_name = 'eastern boundary pH'
        nc_ph_east.units = '-'
        nc_ph_east.coordinates = 'lat_rho s_rho ph_time'
        #
        nc_tmpm_east = nc.createVariable('tmpm_east', np.float64, ('tmpm_time', 's_rho', 'eta_rho',))
        nc_tmpm_east.long_name = 'eastern boundary potential temperature monthly'
        nc_tmpm_east.units = 'Celsius'
        nc_tmpm_east.coordinates = 'lat_rho s_rho tmpm_time'
        #
        nc_salm_east = nc.createVariable('salm_east', np.float64, ('salm_time', 's_rho', 'eta_rho',))
        nc_salm_east.long_name = 'eastern boundary salinity monthly'
        nc_salm_east.units = 'PSU'
        nc_salm_east.coordinates = 'lat_rho s_rho salm_time'
        #

    if obc[2] == 1:
        #
        #   Northern boundary
        #
        nc_temp_north = nc.createVariable('temp_north', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_north.long_name = 'northern boundary potential temperature'
        nc_temp_north.units = 'Celsius'
        nc_temp_north.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_north = nc.createVariable('salt_north', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_north.long_name = 'northern boundary salinity'
        nc_salt_north.units = 'PSU'
        nc_salt_north.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_north = nc.createVariable('u_north', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_north.long_name = 'northern boundary u-momentum component'
        nc_u_north.units = 'meter second-1'
        nc_u_north.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_north = nc.createVariable('v_north', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_north.long_name = 'northern boundary v-momentum component'
        nc_v_north.units = 'meter second-1'
        nc_v_north.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_north = nc.createVariable('ubar_north', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_north.long_name = 'northern boundary vertically integrated u-momentum component'
        nc_ubar_north.units = 'meter second-1'
        nc_ubar_north.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_north = nc.createVariable('vbar_north', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_north.long_name = 'northern boundary vertically integrated v-momentum component'
        nc_vbar_north.units = 'meter second-1'
        nc_vbar_north.coordinates = 'lon_v vclm_time'

        nc_zeta_north = nc.createVariable('zeta_north', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_north.long_name = 'northern boundary sea surface height'
        nc_zeta_north.units = 'meter'
        nc_zeta_north.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_north = nc.createVariable('NO3_north', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_north.long_name = 'northern boundary nitrate'
        nc_no3_north.units = 'mmol/m3'
        nc_no3_north.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_north = nc.createVariable('PO4_north', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_north.long_name = 'northern boundary orthophosphate'
        nc_po4_north.units = 'mmol/m3'
        nc_po4_north.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_north = nc.createVariable('Si_north', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_north.long_name = 'northern boundary silicate'
        nc_si_north.units = 'mmol/m3'
        nc_si_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_north = nc.createVariable('O2_north', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_north.long_name = 'northern boundary dissolved oxygen'
        nc_o2_north.units = 'mmol/m3'
        nc_o2_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_north = nc.createVariable('DIC_north', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_north.long_name = 'northern boundary dissolved inorganic carbon'
        nc_dic_north.units = 'mmol/m3'
        nc_dic_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_nh4_north = nc.createVariable('NH4_north', np.float64, ('nh4_time', 's_rho', 'xi_rho',))
        nc_nh4_north.long_name = 'northern boundary total alkalinity'
        nc_nh4_north.units = 'mmol/m3'
        nc_nh4_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_north = nc.createVariable('FER_north', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_north.long_name = 'northern boundary iron'
        nc_fer_north.units = 'mmol/m3'
        nc_fer_north.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_talk_north = nc.createVariable('TALK_north', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_north.long_name = 'northern boundary alkalinity'
        nc_talk_north.units = 'mmol/m3'
        nc_talk_north.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_ph_north = nc.createVariable('PH_north', np.float64, ('ph_time', 's_rho', 'xi_rho',))
        nc_ph_north.long_name = 'northern boundary pH'
        nc_ph_north.units = '-'
        nc_ph_north.coordinates = 'lon_rho s_rho ph_time'
        #
        nc_tmpm_north = nc.createVariable('tmpm_north', np.float64, ('tmpm_time', 's_rho', 'xi_rho',))
        nc_tmpm_north.long_name = 'northern boundary potential temperature monthly'
        nc_tmpm_north.units = 'Celsius'
        nc_tmpm_north.coordinates = 'lon_rho s_rho tmpm_time'
        #
        nc_salm_north = nc.createVariable('salm_north', np.float64, ('salm_time', 's_rho', 'xi_rho',))
        nc_salm_north.long_name = 'northern boundary salinity monthly'
        nc_salm_north.units = 'PSU'
        nc_salm_north.coordinates = 'lon_rho s_rho salm_time'

    if obc[3] == 1:
        #
        #   Western boundary
        #
        nc_temp_west = nc.createVariable('temp_west', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_west.long_name = 'western boundary potential temperature'
        nc_temp_west.units = 'Celsius'
        nc_temp_west.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_west = nc.createVariable('salt_west', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_west.long_name = 'western boundary salinity'
        nc_salt_west.units = 'PSU'
        nc_salt_west.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_west = nc.createVariable('u_west', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_west.long_name = 'western boundary u-momentum component'
        nc_u_west.units = 'meter second-1'
        nc_u_west.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_west = nc.createVariable('v_west', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_west.long_name = 'western boundary v-momentum component'
        nc_v_west.units = 'meter second-1'
        nc_v_west.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_west = nc.createVariable('ubar_west', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_west.long_name = 'western boundary vertically integrated u-momentum component'
        nc_ubar_west.units = 'meter second-1'
        nc_ubar_west.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_west = nc.createVariable('vbar_west', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_west.long_name = 'western boundary vertically integrated v-momentum component'
        nc_vbar_west.units = 'meter second-1'
        nc_vbar_west.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_west = nc.createVariable('zeta_west', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_west.long_name = 'western boundary sea surface height'
        nc_zeta_west.units = 'meter'
        nc_zeta_west.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_west = nc.createVariable('NO3_west', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_west.long_name = 'western boundary nitrate'
        nc_no3_west.units = 'mmol/m3'
        nc_no3_west.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_west = nc.createVariable('PO4_west', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_west.long_name = 'western boundary orthophosphate'
        nc_po4_west.units = 'mmol/m3'
        nc_po4_west.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_west = nc.createVariable('Si_west', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_west.long_name = 'western boundary silicate'
        nc_si_west.units = 'mmol/m3'
        nc_si_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_west = nc.createVariable('O2_west', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_west.long_name = 'western boundary dissolved oxygen'
        nc_o2_west.units = 'mmol/m3'
        nc_o2_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_west = nc.createVariable('DIC_west', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_west.long_name = 'western boundary dissolved inorganic carbon'
        nc_dic_west.units = 'mmol/m3'
        nc_dic_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_nh4_west = nc.createVariable('NH4_west', np.float64, ('nh4_time', 's_rho', 'eta_rho',))
        nc_nh4_west.long_name = 'western boundary total alkalinity'
        nc_nh4_west.units = 'mmol/m3'
        nc_nh4_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_west = nc.createVariable('FER_west', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_west.long_name = 'western boundary iron'
        nc_fer_west.units = 'mmol/m3'
        nc_fer_west.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_talk_west = nc.createVariable('TALK_west', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_west.long_name = 'western boundary alkalinity'
        nc_talk_west.units = 'mmol/m3'
        nc_talk_west.coordinates = 'lat_rho s_rho talk_time'
        #
        nc_ph_west = nc.createVariable('PH_west', np.float64, ('ph_time', 's_rho', 'eta_rho',))
        nc_ph_west.long_name = 'western boundary pH'
        nc_ph_west.units = '-'
        nc_ph_west.coordinates = 'lat_rho s_rho ph_time'
        #
        nc_tmpm_west = nc.createVariable('tmpm_west', np.float64, ('tmpm_time', 's_rho', 'eta_rho',))
        nc_tmpm_west.long_name = 'western boundary potential temperature monthly'
        nc_tmpm_west.units = 'Celsius'
        nc_tmpm_west.coordinates = 'lat_rho s_rho tmpm_time'
        #
        nc_salm_west = nc.createVariable('salm_west', np.float64, ('salm_time', 's_rho', 'eta_rho',))
        nc_salm_west.long_name = 'western boundary salinity monthly'
        nc_salm_west.units = 'PSU'
        nc_salm_west.coordinates = 'lat_rho s_rho salm_time'
    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #
    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(bry_time)
    nc_tend[:] = np.max(bry_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_sc_w[:] = sc_w
    nc_Cs_r[:] = Cs_r
    nc_Cs_w[:] = Cs_w
    nc_tmpm_time[:] = nor_time
    nc_temp_time[:] = bry_time
    nc_salm_time[:] = nor_time
    nc_salt_time[:] = bry_time
    nc_uclm_time[:] = bry_time
    nc_vclm_time[:] = bry_time
    nc_v2d_time[:] = bry_time
    nc_v3d_time[:] = bry_time
    nc_ssh_time[:] = bry_time
    nc_zeta_time[:] = bry_time
    nc_bry_time[:] = bry_time

    nc_no3_time[:] = nor_time
    nc_po4_time[:] = nor_time
    nc_si_time[:] = nor_time
    nc_o2_time[:] = nor_time
    nc_dic_time[:] = nor_time
    nc_talk_time[:] = nor_time
    nc_nh4_time[:] = pis_time
    nc_fer_time[:] = pis_time
    nc_ph_time[:] = pis_time

    if obc[0] == 1:
        nc_u_south[:] = 0.
        nc_v_south[:] = 0.
        nc_ubar_south[:] = 0.
        nc_vbar_south[:] = 0.
        nc_zeta_south[:] = 0.
        nc_temp_south[:] = 0.
        nc_salt_south[:] = 0.
        nc_no3_south[:] = 0.
        nc_po4_south[:] = 0.
        nc_si_south[:] = 0.
        nc_o2_south[:] = 0.
        nc_dic_south[:] = 0.
        nc_nh4_south[:] = 0.
        nc_fer_south[:] = 0.
        nc_talk_south[:] = 0.
        nc_ph_south[:] = 0.
        nc_tmpm_south[:] = 0.
        nc_salm_south[:] = 0.

    if obc[1] == 1:
        nc_u_east[:] = 0.
        nc_v_east[:] = 0.
        nc_ubar_east[:] = 0.
        nc_vbar_east[:] = 0.
        nc_zeta_east[:] = 0.
        nc_temp_east[:] = 0.
        nc_salt_east[:] = 0.
        nc_no3_east[:] = 0.
        nc_po4_east[:] = 0.
        nc_si_east[:] = 0.
        nc_o2_east[:] = 0.
        nc_dic_east[:] = 0.
        nc_nh4_east[:] = 0.
        nc_fer_east[:] = 0.
        nc_talk_east[:] = 0.
        nc_ph_east[:] = 0.
        nc_tmpm_east[:] = 0.
        nc_salm_east[:] = 0.

    if obc[2] == 1:
        nc_u_north[:] = 0.
        nc_v_north[:] = 0.
        nc_ubar_north[:] = 0.
        nc_vbar_north[:] = 0.
        nc_zeta_north[:] = 0.
        nc_temp_north[:] = 0.
        nc_salt_north[:] = 0.
        nc_no3_north[:] = 0.
        nc_po4_north[:] = 0.
        nc_si_north[:] = 0.
        nc_o2_north[:] = 0.
        nc_dic_north[:] = 0.
        nc_nh4_north[:] = 0.
        nc_fer_north[:] = 0.
        nc_talk_north[:] = 0.
        nc_ph_north[:] = 0.
        nc_tmpm_north[:] = 0.
        nc_salm_north[:] = 0.

    if obc[3] == 1:
        nc_u_west[:] = 0.
        nc_v_west[:] = 0.
        nc_ubar_west[:] = 0.
        nc_vbar_west[:] = 0.
        nc_zeta_west[:] = 0.
        nc_temp_west[:] = 0.
        nc_salt_west[:] = 0.
        nc_no3_west[:] = 0.
        nc_po4_west[:] = 0.
        nc_si_west[:] = 0.
        nc_o2_west[:] = 0.
        nc_dic_west[:] = 0.
        nc_nh4_west[:] = 0.
        nc_fer_west[:] = 0.
        nc_talk_west[:] = 0.
        nc_ph_west[:] = 0.
        nc_tmpm_west[:] = 0.
        nc_salm_west[:] = 0.

    nc.close()

    return

def create_bryfile_BFM(bryname, grdname, title, obc,
                       theta_s, theta_b, hc, N,
                       bry_time, nor_time, pis_time, cycle, vtransform):
    #
    #
    ################################################################
    #
    # function create_bryfile(bryname,grdname,title,obc...
    #                          theta_s,theta_b,hc,N,...
    #                          bry_time,cycle,clobber)
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   bryname      Netcdf climatology file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   obc          open boundaries flag (1=open , [S E N W]).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   bry_time     time.(vector)
    #   cycle        Length (days) for cycling the climatology.(Real)
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + bryname)
    print(' ')
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()

    #
    #
    #

    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + ' m)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    Nt = 0

    #
    #  Create the boundary file
    #
    type = 'BOUNDARY file'
    history = 'CROCO'

    if os.path.exists(bryname):
        os.remove(bryname)

    print('Create: ' + bryname)
    nc = netcdf(bryname, 'w', format='NETCDF4')

    #
    # set global attributes
    #

    nc.type = 'CROCO boundary file'
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.title = title
    nc.bry_file = bryname
    nc.grd_file = grdname

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)

    # PHYSICS time dimension
    nc_dim_bry_time = nc.createDimension('bry_time', Nt)
    nc_dim_tmpm_time = nc.createDimension('tmpm_time', Nt)
    nc_dim_temp_time = nc.createDimension('temp_time', Nt)
    nc_dim_salm_time = nc.createDimension('salm_time', Nt)
    nc_dim_salt_time = nc.createDimension('salt_time', Nt)
    nc_dim_uclm_time = nc.createDimension('uclm_time', Nt)
    nc_dim_vclm_time = nc.createDimension('vclm_time', Nt)
    nc_dim_v2d_time = nc.createDimension('v2d_time', Nt)
    nc_dim_v3d_time = nc.createDimension('v3d_time', Nt)
    nc_dim_ssh_time = nc.createDimension('ssh_time', Nt)
    nc_dim_zeta_time = nc.createDimension('zeta_time', Nt)

    # # BGC time dimension
    nc_dim_no3_time = nc.createDimension('no3_time', Nt)
    nc_dim_po4_time = nc.createDimension('po4_time', Nt)
    nc_dim_si_time = nc.createDimension('si_time', Nt)
    nc_dim_o2_time = nc.createDimension('o2_time', Nt)
    nc_dim_dic_time = nc.createDimension('dic_time', Nt)
    nc_dim_nh4_time = nc.createDimension('nh4_time', Nt)
    nc_dim_fer_time = nc.createDimension('fer_time', Nt)
    nc_dim_talk_time = nc.createDimension('talk_time', Nt)
    nc_dim_ph_time = nc.createDimension('ph_time', Nt)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #
    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    # s
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_bry_time = nc.createVariable('bry_time', np.float64, ('bry_time',))
    nc_bry_time.long_name = 'time for boundary climatology'
    nc_bry_time.units = 'day'
    nc_bry_time.calendar = 'XXX days in every year'
    nc_bry_time.cycle_length = cycle
    #
    nc_tmpm_time = nc.createVariable('tmpm_time', np.float64, ('tmpm_time',))
    nc_tmpm_time.long_name = 'time for temperature monthly'
    nc_tmpm_time.units = 'day'
    nc_tmpm_time.calendar = 'XXX days in every year'
    nc_tmpm_time.cycle_length = cycle
    #
    nc_temp_time = nc.createVariable('temp_time', np.float64, ('temp_time',))
    nc_temp_time.long_name = 'time for temperature'
    nc_temp_time.units = 'day'
    nc_temp_time.calendar = 'XXX days in every year'
    nc_temp_time.cycle_length = cycle
    #
    nc_salm_time = nc.createVariable('salm_time', np.float64, ('salm_time',))
    nc_salm_time.long_name = 'time for salinity monthly'
    nc_salm_time.units = 'day'
    nc_salm_time.calendar = 'XXX days in every year'
    nc_salm_time.cycle_length = cycle
    #
    nc_salt_time = nc.createVariable('salt_time', np.float64, ('salt_time',))
    nc_salt_time.long_name = 'time for salinity'
    nc_salt_time.units = 'day'
    nc_salt_time.calendar = 'XXX days in every year'
    nc_salt_time.cycle_length = cycle
    #
    nc_uclm_time = nc.createVariable('uclm_time', np.float64, ('uclm_time',))
    nc_uclm_time.long_name = 'time climatological u'
    nc_uclm_time.units = 'day'
    nc_uclm_time.calendar = 'XXX days in every year'
    nc_uclm_time.cycle_length = cycle
    #
    nc_vclm_time = nc.createVariable('vclm_time', np.float64, ('vclm_time',))
    nc_vclm_time.long_name = 'time climatological v'
    nc_vclm_time.units = 'day'
    nc_vclm_time.calendar = 'XXX days in every year'
    nc_vclm_time.cycle_length = cycle
    #
    nc_v2d_time = nc.createVariable('v2d_time', np.float64, ('v2d_time',))
    nc_v2d_time.long_name = 'time for 2D velocity'
    nc_v2d_time.units = 'day'
    nc_v2d_time.calendar = 'XXX days in every year'
    nc_v2d_time.cycle_length = cycle
    #
    nc_v3d_time = nc.createVariable('v3d_time', np.float64, ('v3d_time',))
    nc_v3d_time.long_name = 'time for 3D velocity'
    nc_v3d_time.units = 'day'
    nc_v3d_time.calendar = 'XXX days in every year'
    nc_v3d_time.cycle_length = cycle
    #
    nc_ssh_time = nc.createVariable('ssh_time', np.float64, ('ssh_time',))
    nc_ssh_time.long_name = 'time for sea surface height'
    nc_ssh_time.units = 'day'
    nc_ssh_time.calendar = 'XXX days in every year'
    nc_ssh_time.cycle_length = cycle
    #
    nc_zeta_time = nc.createVariable('zeta_time', np.float64, ('zeta_time',))
    nc_zeta_time.long_name = 'time for sea surface height'
    nc_zeta_time.units = 'day'
    nc_zeta_time.calendar = 'XXX days in every year'
    nc_zeta_time.cycle_length = cycle
    #
    nc_no3_time = nc.createVariable('no3_time', np.float64, ('no3_time',))
    nc_no3_time.long_name = 'time for nitrate'
    nc_no3_time.units = 'day'
    nc_no3_time.calendar = 'XXX days in every year'
    nc_no3_time.cycle_length = cycle
    #
    nc_po4_time = nc.createVariable('po4_time', np.float64, ('po4_time',))
    nc_po4_time.long_name = 'time for orthophosphate'
    nc_po4_time.units = 'day'
    nc_po4_time.calendar = 'XXX days in every year'
    nc_po4_time.cycle_length = cycle
    #
    nc_si_time = nc.createVariable('si_time', np.float64, ('si_time',))
    nc_si_time.long_name = 'time for silicate'
    nc_si_time.units = 'day'
    nc_si_time.calendar = 'XXX days in every year'
    nc_si_time.cycle_length = cycle
    #
    nc_o2_time = nc.createVariable('o2_time', np.float64, ('o2_time',))
    nc_o2_time.long_name = 'time for dissolved oxygen'
    nc_o2_time.units = 'day'
    nc_o2_time.calendar = 'XXX days in every year'
    nc_o2_time.cycle_length = cycle
    #
    nc_dic_time = nc.createVariable('dic_time', np.float64, ('dic_time',))
    nc_dic_time.long_name = 'time for dissolved inorganic carbon'
    nc_dic_time.units = 'day'
    nc_dic_time.calendar = 'XXX days in every year'
    nc_dic_time.cycle_length = cycle
    #
    nc_nh4_time = nc.createVariable('nh4_time', np.float64, ('nh4_time',))
    nc_nh4_time.long_name = 'time for ammonium'
    nc_nh4_time.units = 'day'
    nc_nh4_time.calendar = 'XXX days in every year'
    nc_nh4_time.cycle_length = cycle
    #
    nc_fer_time = nc.createVariable('fer_time', np.float64, ('fer_time',))
    nc_fer_time.long_name = 'time for iron'
    nc_fer_time.units = 'day'
    nc_fer_time.calendar = 'XXX days in every year'
    nc_fer_time.cycle_length = cycle
    #
    nc_talk_time = nc.createVariable('talk_time', np.float64, ('talk_time',))
    nc_talk_time.long_name = 'time for alkalinity'
    nc_talk_time.units = 'day'
    nc_talk_time.calendar = 'XXX days in every year'
    nc_talk_time.cycle_length = cycle
    #
    nc_ph_time = nc.createVariable('ph_time', np.float64, ('ph_time',))
    nc_ph_time.long_name = 'time for ph'
    nc_ph_time.units = 'day'
    nc_ph_time.calendar = 'XXX days in every year'
    nc_ph_time.cycle_length = cycle
    #

    if obc[0] == 1:
        #
        #   Southern boundary
        #
        nc_temp_south = nc.createVariable('temp_south', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_south.long_name = 'southern boundary potential temperature'
        nc_temp_south.units = 'Celsius'
        nc_temp_south.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_south = nc.createVariable('salt_south', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_south.long_name = 'southern boundary salinity'
        nc_salt_south.units = 'PSU'
        nc_salt_south.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_south = nc.createVariable('u_south', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_south.long_name = 'southern boundary u-momentum component'
        nc_u_south.units = 'meter second-1'
        nc_u_south.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_south = nc.createVariable('v_south', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_south.long_name = 'southern boundary v-momentum component'
        nc_v_south.units = 'meter second-1'
        nc_v_south.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_south = nc.createVariable('ubar_south', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_south.long_name = 'southern boundary vertically integrated u-momentum component'
        nc_ubar_south.units = 'meter second-1'
        nc_ubar_south.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_south = nc.createVariable('vbar_south', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_south.long_name = 'southern boundary vertically integrated v-momentum component'
        nc_vbar_south.units = 'meter second-1'
        nc_vbar_south.coordinates = 'lon_v vclm_time'
        #
        nc_zeta_south = nc.createVariable('zeta_south', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_south.long_name = 'southern boundary sea surface height'
        nc_zeta_south.units = 'meter'
        nc_zeta_south.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_south = nc.createVariable('N3n_south', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_south.long_name = 'southern boundary nitrate'
        nc_no3_south.units = 'mmol/m3'
        nc_no3_south.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_south = nc.createVariable('N1p_south', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_south.long_name = 'southern boundary orthophosphate'
        nc_po4_south.units = 'mmol/m3'
        nc_po4_south.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_south = nc.createVariable('N5s_south', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_south.long_name = 'southern boundary silicate'
        nc_si_south.units = 'mmol/m3'
        nc_si_south.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_south = nc.createVariable('O2o_south', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_south.long_name = 'southern boundary dissolved oxygen'
        nc_o2_south.units = 'mmol/m3'
        nc_o2_south.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_south = nc.createVariable('O3c_south', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_south.long_name = 'southern boundary dissolved inorganic carbon'
        nc_dic_south.units = 'mmol/m3'
        nc_dic_south.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_nh4_south = nc.createVariable('N4n_south', np.float64, ('nh4_time', 's_rho', 'xi_rho',))
        nc_nh4_south.long_name = 'southern boundary ammonium'
        nc_nh4_south.units = 'mmol/m3'
        nc_nh4_south.coordinates = 'lon_rho s_rho nh4_time'
        #
        nc_fer_south = nc.createVariable('N7f_south', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_south.long_name = 'southern boundary iron'
        nc_fer_south.units = 'mmol/m3'
        nc_fer_south.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_talk_south = nc.createVariable('O3h_south', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_south.long_name = 'southern boundary alkalinity'
        nc_talk_south.units = 'mmol/m3'
        nc_talk_south.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_ph_south = nc.createVariable('PH_south', np.float64, ('ph_time', 's_rho', 'xi_rho',))
        nc_ph_south.long_name = 'southern boundary pH'
        nc_ph_south.units = '-'
        nc_ph_south.coordinates = 'lon_rho s_rho ph_time'
        #
        nc_tmpm_south = nc.createVariable('tmpm_south', np.float64, ('tmpm_time', 's_rho', 'xi_rho',))
        nc_tmpm_south.long_name = 'southern boundary potential temperature monthly'
        nc_tmpm_south.units = 'Celsius'
        nc_tmpm_south.coordinates = 'lon_rho s_rho tmpm_time'
        #
        nc_salm_south = nc.createVariable('salm_south', np.float64, ('salm_time', 's_rho', 'xi_rho',))
        nc_salm_south.long_name = 'southern boundary salinity monthly'
        nc_salm_south.units = 'PSU'
        nc_salm_south.coordinates = 'lon_rho s_rho salm_time'
        #
    if obc[1] == 1:
        #
        #   Eastern boundary
        #
        nc_temp_east = nc.createVariable('temp_east', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_east.long_name = 'eastern boundary potential temperature'
        nc_temp_east.units = 'Celsius'
        nc_temp_east.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_east = nc.createVariable('salt_east', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_east.long_name = 'eastern boundary salinity'
        nc_salt_east.units = 'PSU'
        nc_salt_east.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_east = nc.createVariable('u_east', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_east.long_name = 'eastern boundary u-momentum component'
        nc_u_east.units = 'meter second-1'
        nc_u_east.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_east = nc.createVariable('v_east', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_east.long_name = 'eastern boundary v-momentum component'
        nc_v_east.units = 'meter second-1'
        nc_v_east.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_east = nc.createVariable('ubar_east', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_east.long_name = 'eastern boundary vertically integrated u-momentum component'
        nc_ubar_east.units = 'meter second-1'
        nc_ubar_east.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_east = nc.createVariable('vbar_east', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_east.long_name = 'eastern boundary vertically integrated v-momentum component'
        nc_vbar_east.units = 'meter second-1'
        nc_vbar_east.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_east = nc.createVariable('zeta_east', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_east.long_name = 'eastern boundary sea surface height'
        nc_zeta_east.units = 'meter'
        nc_zeta_east.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_east = nc.createVariable('N3n_east', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_east.long_name = 'eastern boundary nitrate'
        nc_no3_east.units = 'mmol/m3'
        nc_no3_east.coordinates = 'lat_rho s_rho no3_time'
        #
        nc_po4_east = nc.createVariable('N1p_east', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_east.long_name = 'eastern boundary orthophosphate'
        nc_po4_east.units = 'mmol/m3'
        nc_po4_east.coordinates = 'lat_rho s_rho po4_time'
        #
        nc_si_east = nc.createVariable('N5s_east', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_east.long_name = 'eastern boundary silicate'
        nc_si_east.units = 'mmol/m3'
        nc_si_east.coordinates = 'lat_rho s_rho si_time'
        #
        nc_o2_east = nc.createVariable('O2o_east', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_east.long_name = 'eastern boundary dissolved oxygen'
        nc_o2_east.units = 'mmol/m3'
        nc_o2_east.coordinates = 'lat_rho s_rho o2_time'
        #
        nc_dic_east = nc.createVariable('O3c_east', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_east.long_name = 'eastern boundary dissolved inorganic carbon'
        nc_dic_east.units = 'mmol/m3'
        nc_dic_east.coordinates = 'lat_rho s_rho dic_time'
        #
        nc_nh4_east = nc.createVariable('N4n_east', np.float64, ('nh4_time', 's_rho', 'eta_rho',))
        nc_nh4_east.long_name = 'eastern boundary total alkalinity'
        nc_nh4_east.units = 'mmol/m3'
        nc_nh4_east.coordinates = 'lat_rho s_rho nh4_time'
        #
        nc_fer_east = nc.createVariable('N7f_east', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_east.long_name = 'eastern boundary iron'
        nc_fer_east.units = 'mmol/m3'
        nc_fer_east.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_talk_east = nc.createVariable('O3h_east', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_east.long_name = 'eastern boundary alkalinity'
        nc_talk_east.units = 'mmol/m3'
        nc_talk_east.coordinates = 'lat_rho s_rho talk_time'
        #
        nc_ph_east = nc.createVariable('PH_east', np.float64, ('ph_time', 's_rho', 'eta_rho',))
        nc_ph_east.long_name = 'eastern boundary pH'
        nc_ph_east.units = '-'
        nc_ph_east.coordinates = 'lat_rho s_rho ph_time'
        #
        nc_tmpm_east = nc.createVariable('tmpm_east', np.float64, ('tmpm_time', 's_rho', 'eta_rho',))
        nc_tmpm_east.long_name = 'eastern boundary potential temperature monthly'
        nc_tmpm_east.units = 'Celsius'
        nc_tmpm_east.coordinates = 'lat_rho s_rho tmpm_time'
        #
        nc_salm_east = nc.createVariable('salm_east', np.float64, ('salm_time', 's_rho', 'eta_rho',))
        nc_salm_east.long_name = 'eastern boundary salinity monthly'
        nc_salm_east.units = 'PSU'
        nc_salm_east.coordinates = 'lat_rho s_rho salm_time'
        #

    if obc[2] == 1:
        #
        #   Northern boundary
        #
        nc_temp_north = nc.createVariable('temp_north', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_north.long_name = 'northern boundary potential temperature'
        nc_temp_north.units = 'Celsius'
        nc_temp_north.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_north = nc.createVariable('salt_north', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_north.long_name = 'northern boundary salinity'
        nc_salt_north.units = 'PSU'
        nc_salt_north.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_north = nc.createVariable('u_north', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_north.long_name = 'northern boundary u-momentum component'
        nc_u_north.units = 'meter second-1'
        nc_u_north.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_north = nc.createVariable('v_north', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_north.long_name = 'northern boundary v-momentum component'
        nc_v_north.units = 'meter second-1'
        nc_v_north.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_north = nc.createVariable('ubar_north', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_north.long_name = 'northern boundary vertically integrated u-momentum component'
        nc_ubar_north.units = 'meter second-1'
        nc_ubar_north.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_north = nc.createVariable('vbar_north', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_north.long_name = 'northern boundary vertically integrated v-momentum component'
        nc_vbar_north.units = 'meter second-1'
        nc_vbar_north.coordinates = 'lon_v vclm_time'

        nc_zeta_north = nc.createVariable('zeta_north', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_north.long_name = 'northern boundary sea surface height'
        nc_zeta_north.units = 'meter'
        nc_zeta_north.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_north = nc.createVariable('N3n_north', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_north.long_name = 'northern boundary nitrate'
        nc_no3_north.units = 'mmol/m3'
        nc_no3_north.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_north = nc.createVariable('N1p_north', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_north.long_name = 'northern boundary orthophosphate'
        nc_po4_north.units = 'mmol/m3'
        nc_po4_north.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_north = nc.createVariable('N5s_north', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_north.long_name = 'northern boundary silicate'
        nc_si_north.units = 'mmol/m3'
        nc_si_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_north = nc.createVariable('O2o_north', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_north.long_name = 'northern boundary dissolved oxygen'
        nc_o2_north.units = 'mmol/m3'
        nc_o2_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_north = nc.createVariable('O3c_north', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_north.long_name = 'northern boundary dissolved inorganic carbon'
        nc_dic_north.units = 'mmol/m3'
        nc_dic_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_nh4_north = nc.createVariable('N4n_north', np.float64, ('nh4_time', 's_rho', 'xi_rho',))
        nc_nh4_north.long_name = 'northern boundary total alkalinity'
        nc_nh4_north.units = 'mmol/m3'
        nc_nh4_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_north = nc.createVariable('N7f_north', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_north.long_name = 'northern boundary iron'
        nc_fer_north.units = 'mmol/m3'
        nc_fer_north.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_talk_north = nc.createVariable('O3h_north', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_north.long_name = 'northern boundary alkalinity'
        nc_talk_north.units = 'mmol/m3'
        nc_talk_north.coordinates = 'lon_rho s_rho talk_time'
        #
        nc_ph_north = nc.createVariable('PH_north', np.float64, ('ph_time', 's_rho', 'xi_rho',))
        nc_ph_north.long_name = 'northern boundary pH'
        nc_ph_north.units = '-'
        nc_ph_north.coordinates = 'lon_rho s_rho ph_time'
        #
        nc_tmpm_north = nc.createVariable('tmpm_north', np.float64, ('tmpm_time', 's_rho', 'xi_rho',))
        nc_tmpm_north.long_name = 'northern boundary potential temperature monthly'
        nc_tmpm_north.units = 'Celsius'
        nc_tmpm_north.coordinates = 'lon_rho s_rho tmpm_time'
        #
        nc_salm_north = nc.createVariable('salm_north', np.float64, ('salm_time', 's_rho', 'xi_rho',))
        nc_salm_north.long_name = 'northern boundary salinity monthly'
        nc_salm_north.units = 'PSU'
        nc_salm_north.coordinates = 'lon_rho s_rho salm_time'

    if obc[3] == 1:
        #
        #   Western boundary
        #
        nc_temp_west = nc.createVariable('temp_west', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_west.long_name = 'western boundary potential temperature'
        nc_temp_west.units = 'Celsius'
        nc_temp_west.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_west = nc.createVariable('salt_west', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_west.long_name = 'western boundary salinity'
        nc_salt_west.units = 'PSU'
        nc_salt_west.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_west = nc.createVariable('u_west', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_west.long_name = 'western boundary u-momentum component'
        nc_u_west.units = 'meter second-1'
        nc_u_west.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_west = nc.createVariable('v_west', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_west.long_name = 'western boundary v-momentum component'
        nc_v_west.units = 'meter second-1'
        nc_v_west.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_west = nc.createVariable('ubar_west', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_west.long_name = 'western boundary vertically integrated u-momentum component'
        nc_ubar_west.units = 'meter second-1'
        nc_ubar_west.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_west = nc.createVariable('vbar_west', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_west.long_name = 'western boundary vertically integrated v-momentum component'
        nc_vbar_west.units = 'meter second-1'
        nc_vbar_west.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_west = nc.createVariable('zeta_west', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_west.long_name = 'western boundary sea surface height'
        nc_zeta_west.units = 'meter'
        nc_zeta_west.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_west = nc.createVariable('N3n_west', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_west.long_name = 'western boundary nitrate'
        nc_no3_west.units = 'mmol/m3'
        nc_no3_west.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_west = nc.createVariable('N1p_west', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_west.long_name = 'western boundary orthophosphate'
        nc_po4_west.units = 'mmol/m3'
        nc_po4_west.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_west = nc.createVariable('N5s_west', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_west.long_name = 'western boundary silicate'
        nc_si_west.units = 'mmol/m3'
        nc_si_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_west = nc.createVariable('O2o_west', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_west.long_name = 'western boundary dissolved oxygen'
        nc_o2_west.units = 'mmol/m3'
        nc_o2_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_west = nc.createVariable('O3c_west', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_west.long_name = 'western boundary dissolved inorganic carbon'
        nc_dic_west.units = 'mmol/m3'
        nc_dic_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_nh4_west = nc.createVariable('N4n_west', np.float64, ('nh4_time', 's_rho', 'eta_rho',))
        nc_nh4_west.long_name = 'western boundary total alkalinity'
        nc_nh4_west.units = 'mmol/m3'
        nc_nh4_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_west = nc.createVariable('N7f_west', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_west.long_name = 'western boundary iron'
        nc_fer_west.units = 'mmol/m3'
        nc_fer_west.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_talk_west = nc.createVariable('O3h_west', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_west.long_name = 'western boundary alkalinity'
        nc_talk_west.units = 'mmol/m3'
        nc_talk_west.coordinates = 'lat_rho s_rho talk_time'
        #
        nc_ph_west = nc.createVariable('PH_west', np.float64, ('ph_time', 's_rho', 'eta_rho',))
        nc_ph_west.long_name = 'western boundary pH'
        nc_ph_west.units = '-'
        nc_ph_west.coordinates = 'lat_rho s_rho ph_time'
        #
        nc_tmpm_west = nc.createVariable('tmpm_west', np.float64, ('tmpm_time', 's_rho', 'eta_rho',))
        nc_tmpm_west.long_name = 'western boundary potential temperature monthly'
        nc_tmpm_west.units = 'Celsius'
        nc_tmpm_west.coordinates = 'lat_rho s_rho tmpm_time'
        #
        nc_salm_west = nc.createVariable('salm_west', np.float64, ('salm_time', 's_rho', 'eta_rho',))
        nc_salm_west.long_name = 'western boundary salinity monthly'
        nc_salm_west.units = 'PSU'
        nc_salm_west.coordinates = 'lat_rho s_rho salm_time'
    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #
    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(bry_time)
    nc_tend[:] = np.max(bry_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_sc_w[:] = sc_w
    nc_Cs_r[:] = Cs_r
    nc_Cs_w[:] = Cs_w
    nc_tmpm_time[:] = nor_time
    nc_temp_time[:] = bry_time
    nc_salm_time[:] = nor_time
    nc_salt_time[:] = bry_time
    nc_uclm_time[:] = bry_time
    nc_vclm_time[:] = bry_time
    nc_v2d_time[:] = bry_time
    nc_v3d_time[:] = bry_time
    nc_ssh_time[:] = bry_time
    nc_zeta_time[:] = bry_time
    nc_bry_time[:] = bry_time

    nc_no3_time[:] = nor_time
    nc_po4_time[:] = nor_time
    nc_si_time[:] = nor_time
    nc_o2_time[:] = nor_time
    nc_dic_time[:] = nor_time
    nc_talk_time[:] = nor_time
    nc_nh4_time[:] = pis_time
    nc_fer_time[:] = pis_time
    nc_ph_time[:] = pis_time

    if obc[0] == 1:
        nc_u_south[:] = 0.
        nc_v_south[:] = 0.
        nc_ubar_south[:] = 0.
        nc_vbar_south[:] = 0.
        nc_zeta_south[:] = 0.
        nc_temp_south[:] = 0.
        nc_salt_south[:] = 0.
        nc_no3_south[:] = 0.
        nc_po4_south[:] = 0.
        nc_si_south[:] = 0.
        nc_o2_south[:] = 0.
        nc_dic_south[:] = 0.
        nc_nh4_south[:] = 0.
        nc_fer_south[:] = 0.
        nc_talk_south[:] = 0.
        nc_ph_south[:] = 0.
        nc_tmpm_south[:] = 0.
        nc_salm_south[:] = 0.

    if obc[1] == 1:
        nc_u_east[:] = 0.
        nc_v_east[:] = 0.
        nc_ubar_east[:] = 0.
        nc_vbar_east[:] = 0.
        nc_zeta_east[:] = 0.
        nc_temp_east[:] = 0.
        nc_salt_east[:] = 0.
        nc_no3_east[:] = 0.
        nc_po4_east[:] = 0.
        nc_si_east[:] = 0.
        nc_o2_east[:] = 0.
        nc_dic_east[:] = 0.
        nc_nh4_east[:] = 0.
        nc_fer_east[:] = 0.
        nc_talk_east[:] = 0.
        nc_ph_east[:] = 0.
        nc_tmpm_east[:] = 0.
        nc_salm_east[:] = 0.

    if obc[2] == 1:
        nc_u_north[:] = 0.
        nc_v_north[:] = 0.
        nc_ubar_north[:] = 0.
        nc_vbar_north[:] = 0.
        nc_zeta_north[:] = 0.
        nc_temp_north[:] = 0.
        nc_salt_north[:] = 0.
        nc_no3_north[:] = 0.
        nc_po4_north[:] = 0.
        nc_si_north[:] = 0.
        nc_o2_north[:] = 0.
        nc_dic_north[:] = 0.
        nc_nh4_north[:] = 0.
        nc_fer_north[:] = 0.
        nc_talk_north[:] = 0.
        nc_ph_north[:] = 0.
        nc_tmpm_north[:] = 0.
        nc_salm_north[:] = 0.

    if obc[3] == 1:
        nc_u_west[:] = 0.
        nc_v_west[:] = 0.
        nc_ubar_west[:] = 0.
        nc_vbar_west[:] = 0.
        nc_zeta_west[:] = 0.
        nc_temp_west[:] = 0.
        nc_salt_west[:] = 0.
        nc_no3_west[:] = 0.
        nc_po4_west[:] = 0.
        nc_si_west[:] = 0.
        nc_o2_west[:] = 0.
        nc_dic_west[:] = 0.
        nc_nh4_west[:] = 0.
        nc_fer_west[:] = 0.
        nc_talk_west[:] = 0.
        nc_ph_west[:] = 0.
        nc_tmpm_west[:] = 0.
        nc_salm_west[:] = 0.

    nc.close()

    return


def create_bryfile_CMIP6(bryname, grdname, title, obc,
                         theta_s, theta_b, hc, N,
                         bry_time, nor_time, pis_time, cycle, vtransform):
    #
    #
    ################################################################
    #
    # function create_bryfile(bryname,grdname,title,obc...
    #                          theta_s,theta_b,hc,N,...
    #                          bry_time,cycle,clobber)
    #
    #   This function create the header of a Netcdf climatology
    #   file.
    #
    #   Input:
    #
    #   bryname      Netcdf climatology file name (character string).
    #   grdname      Netcdf grid file name (character string).
    #   obc          open boundaries flag (1=open , [S E N W]).
    #   theta_s      S-coordinate surface control parameter.(Real)
    #   theta_b      S-coordinate bottom control parameter.(Real)
    #   hc           Width (m) of surface or bottom boundary layer
    #                where higher vertical resolution is required
    #                during stretching.(Real)
    #   N            Number of vertical levels.(Integer)
    #   bry_time     time.(vector)
    #   cycle        Length (days) for cycling the climatology.(Real)
    #
    #
    ################################################################
    #
    #

    print(' ')
    print(' Creating the file : ' + bryname)
    print(' ')
    print(' VTRANSFORM = ' + str(vtransform))

    #
    #  Read the grid file and check the topography
    #

    nc = netcdf(grdname, 'r')
    h = nc.variables['h'][:]
    maskr = nc.variables['mask_rho'][:]
    [Mp, Lp] = np.shape(h)
    nc.close()

    #
    #
    #

    hmin = np.min(h[np.where(maskr == 1)])
    if vtransform == 1:
        if hc > hmin:
            print('Error: hc (' + hc + ' m) > hmin (' + hmin + ' m)')

    L = Lp - 1
    M = Mp - 1
    Np = N + 1

    Nt = 0

    #
    #  Create the boundary file
    #
    type = 'BOUNDARY file'
    history = 'CROCO'

    if os.path.exists(bryname):
        os.remove(bryname)

    print('Create: ' + bryname)
    nc = netcdf(bryname, 'w', format='NETCDF4')

    #
    # set global attributes
    #

    nc.type = 'CROCO boundary file'
    nc.history = 'Created ' + str(time.ctime(time.time()))
    nc.title = title
    nc.bry_file = bryname
    nc.grd_file = grdname

    #
    #  Create dimensions
    #

    nc_dim_xi_u = nc.createDimension('xi_u', L)
    nc_dim_xi_v = nc.createDimension('xi_v', Lp)
    nc_dim_xi_rho = nc.createDimension('xi_rho', Lp)
    nc_dim_eta_u = nc.createDimension('eta_u', Mp)
    nc_dim_eta_v = nc.createDimension('eta_v', M)
    nc_dim_eta_rho = nc.createDimension('eta_rho', Mp)
    nc_dim_s_rho = nc.createDimension('s_rho', N)
    nc_dim_s_w = nc.createDimension('s_w', Np)
    nc_dim_tracer = nc.createDimension('tracer', 2)

    # PHYSICS time dimension
    nc_dim_bry_time = nc.createDimension('bry_time', Nt)
    nc_dim_temp_time = nc.createDimension('temp_time', Nt)
    nc_dim_salt_time = nc.createDimension('salt_time', Nt)
    nc_dim_uclm_time = nc.createDimension('uclm_time', Nt)
    nc_dim_vclm_time = nc.createDimension('vclm_time', Nt)
    nc_dim_v2d_time = nc.createDimension('v2d_time', Nt)
    nc_dim_v3d_time = nc.createDimension('v3d_time', Nt)
    nc_dim_ssh_time = nc.createDimension('ssh_time', Nt)
    nc_dim_zeta_time = nc.createDimension('zeta_time', Nt)

    # # BGC time dimension
    nc_dim_no3_time = nc.createDimension('no3_time', Nt)
    nc_dim_po4_time = nc.createDimension('po4_time', Nt)
    nc_dim_si_time = nc.createDimension('si_time', Nt)
    nc_dim_o2_time = nc.createDimension('o2_time', Nt)
    nc_dim_dic_time = nc.createDimension('dic_time', Nt)
    nc_dim_nh4_time = nc.createDimension('nh4_time', Nt)
    nc_dim_fer_time = nc.createDimension('fer_time', Nt)
    nc_dim_talk_time = nc.createDimension('talk_time', Nt)
    nc_dim_one = nc.createDimension('one', 1)

    #
    #  Create variables and attributes
    #
    nc_spherical = nc.createVariable('spherical', 'S1', ('one',))
    nc_spherical.long_name = 'grid type logical switch'
    nc_spherical.flag_values = 'T, F'
    nc_spherical.flag_meanings = 'spherical Cartesian'
    #
    nc_Vtransform = nc.createVariable('Vtransform', 'i4', ('one',))
    nc_Vtransform.long_name = 'vertical terrain-following transformation equation'
    #
    nc_Vstretching = nc.createVariable('Vstretching', 'i4', ('one',))
    nc_Vstretching.long_name = 'vertical terrain-following stretching function'
    #
    nc_tstart = nc.createVariable('tstart', np.float64, ('one',))
    nc_tstart.long_name = 'start processing day'
    nc_tstart.units = 'day'
    #
    nc_tend = nc.createVariable('tend', np.float64, ('one',))
    nc_tend.long_name = 'end processing day'
    nc_tend.units = 'day'
    #
    nc_theta_s = nc.createVariable('theta_s', np.float64, ('one',))
    nc_theta_s.long_name = 'S-coordinate surface control parameter'
    nc_theta_s.units = 'nondimensional'
    #
    nc_theta_b = nc.createVariable('theta_b', np.float64, ('one',))
    nc_theta_b.long_name = 'S-coordinate bottom control parameter'
    nc_theta_b.units = 'nondimensional'
    #
    nc_Tcline = nc.createVariable('Tcline', np.float64, ('one',))
    nc_Tcline.long_name = 'S-coordinate surface/bottom layer width'
    nc_Tcline.units = 'meter'
    #
    nc_hc = nc.createVariable('hc', np.float64, ('one',))
    nc_hc.long_name = 'S-coordinate parameter, critical depth'
    nc_hc.units = 'meter'
    # s
    nc_sc_r = nc.createVariable('sc_r', np.float64, ('s_rho',))
    nc_sc_r.long_name = 'S-coordinate at RHO-points'
    nc_sc_r.valid_min = -1.
    nc_sc_r.valid_max = 0.
    nc_sc_r.positive = 'up'

    if vtransform == 1:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_r.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_r.formula_terms = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
    #
    nc_sc_w = nc.createVariable('sc_w', np.float64, ('s_w',))
    nc_sc_w.long_name = 'S-coordinate at W-points'
    nc_sc_w.valid_min = -1.
    nc_sc_w.valid_max = 0.
    nc_sc_w.positive = 'up'

    if vtransform == 1:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g1'
    elif vtransform == 2:
        nc_sc_w.standard_name = 'ocean_s_coordinate_g2'

    nc_sc_w.formula_terms = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
    #
    nc_Cs_r = nc.createVariable('Cs_r', np.float64, ('s_rho',))
    nc_Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
    nc_Cs_r.units = 'nondimensional'
    nc_Cs_r.valid_min = -1
    nc_Cs_r.valid_max = 0
    #
    nc_Cs_w = nc.createVariable('Cs_w', np.float64, ('s_w',))
    nc_Cs_w.long_name = 'S-coordinate stretching curves at W-points'
    nc_Cs_w.units = 'nondimensional'
    nc_Cs_w.valid_min = -1
    nc_Cs_w.valid_max = 0
    #
    nc_bry_time = nc.createVariable('bry_time', np.float64, ('bry_time',))
    nc_bry_time.long_name = 'time for boundary climatology'
    nc_bry_time.units = 'day'
    nc_bry_time.calendar = 'XXX days in every year'
    nc_bry_time.cycle_length = cycle
    #
    nc_temp_time = nc.createVariable('temp_time', np.float64, ('temp_time',))
    nc_temp_time.long_name = 'time for temperature'
    nc_temp_time.units = 'day'
    nc_temp_time.calendar = 'XXX days in every year'
    nc_temp_time.cycle_length = cycle
    #
    nc_salt_time = nc.createVariable('salt_time', np.float64, ('salt_time',))
    nc_salt_time.long_name = 'time for salinity'
    nc_salt_time.units = 'day'
    nc_salt_time.calendar = 'XXX days in every year'
    nc_salt_time.cycle_length = cycle
    #
    nc_uclm_time = nc.createVariable('uclm_time', np.float64, ('uclm_time',))
    nc_uclm_time.long_name = 'time climatological u'
    nc_uclm_time.units = 'day'
    nc_uclm_time.calendar = 'XXX days in every year'
    nc_uclm_time.cycle_length = cycle
    #
    nc_vclm_time = nc.createVariable('vclm_time', np.float64, ('vclm_time',))
    nc_vclm_time.long_name = 'time climatological v'
    nc_vclm_time.units = 'day'
    nc_vclm_time.calendar = 'XXX days in every year'
    nc_vclm_time.cycle_length = cycle
    #
    nc_v2d_time = nc.createVariable('v2d_time', np.float64, ('v2d_time',))
    nc_v2d_time.long_name = 'time for 2D velocity'
    nc_v2d_time.units = 'day'
    nc_v2d_time.calendar = 'XXX days in every year'
    nc_v2d_time.cycle_length = cycle
    #
    nc_v3d_time = nc.createVariable('v3d_time', np.float64, ('v3d_time',))
    nc_v3d_time.long_name = 'time for 3D velocity'
    nc_v3d_time.units = 'day'
    nc_v3d_time.calendar = 'XXX days in every year'
    nc_v3d_time.cycle_length = cycle
    #
    nc_ssh_time = nc.createVariable('ssh_time', np.float64, ('ssh_time',))
    nc_ssh_time.long_name = 'time for sea surface height'
    nc_ssh_time.units = 'day'
    nc_ssh_time.calendar = 'XXX days in every year'
    nc_ssh_time.cycle_length = cycle
    #
    nc_zeta_time = nc.createVariable('zeta_time', np.float64, ('zeta_time',))
    nc_zeta_time.long_name = 'time for sea surface height'
    nc_zeta_time.units = 'day'
    nc_zeta_time.calendar = 'XXX days in every year'
    nc_zeta_time.cycle_length = cycle
    #
    nc_no3_time = nc.createVariable('no3_time', np.float64, ('no3_time',))
    nc_no3_time.long_name = 'time for nitrate'
    nc_no3_time.units = 'day'
    nc_no3_time.calendar = 'XXX days in every year'
    nc_no3_time.cycle_length = cycle
    #
    nc_po4_time = nc.createVariable('po4_time', np.float64, ('po4_time',))
    nc_po4_time.long_name = 'time for orthophosphate'
    nc_po4_time.units = 'day'
    nc_po4_time.calendar = 'XXX days in every year'
    nc_po4_time.cycle_length = cycle
    #
    nc_si_time = nc.createVariable('si_time', np.float64, ('si_time',))
    nc_si_time.long_name = 'time for silicate'
    nc_si_time.units = 'day'
    nc_si_time.calendar = 'XXX days in every year'
    nc_si_time.cycle_length = cycle
    #
    nc_o2_time = nc.createVariable('o2_time', np.float64, ('o2_time',))
    nc_o2_time.long_name = 'time for dissolved oxygen'
    nc_o2_time.units = 'day'
    nc_o2_time.calendar = 'XXX days in every year'
    nc_o2_time.cycle_length = cycle
    #
    nc_dic_time = nc.createVariable('dic_time', np.float64, ('dic_time',))
    nc_dic_time.long_name = 'time for dissolved inorganic carbon'
    nc_dic_time.units = 'day'
    nc_dic_time.calendar = 'XXX days in every year'
    nc_dic_time.cycle_length = cycle
    #
    nc_fer_time = nc.createVariable('fer_time', np.float64, ('fer_time',))
    nc_fer_time.long_name = 'time for iron'
    nc_fer_time.units = 'day'
    nc_fer_time.calendar = 'XXX days in every year'
    nc_fer_time.cycle_length = cycle
    #
    nc_talk_time = nc.createVariable('talk_time', np.float64, ('talk_time',))
    nc_talk_time.long_name = 'time for alkalinity'
    nc_talk_time.units = 'day'
    nc_talk_time.calendar = 'XXX days in every year'
    nc_talk_time.cycle_length = cycle
    #

    if obc[0] == 1:
        #
        #   Southern boundary
        #
        nc_temp_south = nc.createVariable('temp_south', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_south.long_name = 'southern boundary potential temperature'
        nc_temp_south.units = 'Celsius'
        nc_temp_south.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_south = nc.createVariable('salt_south', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_south.long_name = 'southern boundary salinity'
        nc_salt_south.units = 'PSU'
        nc_salt_south.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_south = nc.createVariable('u_south', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_south.long_name = 'southern boundary u-momentum component'
        nc_u_south.units = 'meter second-1'
        nc_u_south.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_south = nc.createVariable('v_south', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_south.long_name = 'southern boundary v-momentum component'
        nc_v_south.units = 'meter second-1'
        nc_v_south.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_south = nc.createVariable('ubar_south', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_south.long_name = 'southern boundary vertically integrated u-momentum component'
        nc_ubar_south.units = 'meter second-1'
        nc_ubar_south.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_south = nc.createVariable('vbar_south', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_south.long_name = 'southern boundary vertically integrated v-momentum component'
        nc_vbar_south.units = 'meter second-1'
        nc_vbar_south.coordinates = 'lon_v vclm_time'
        #
        nc_zeta_south = nc.createVariable('zeta_south', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_south.long_name = 'southern boundary sea surface height'
        nc_zeta_south.units = 'meter'
        nc_zeta_south.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_south = nc.createVariable('NO3_south', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_south.long_name = 'southern boundary nitrate'
        nc_no3_south.units = 'mmol/m3'
        nc_no3_south.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_south = nc.createVariable('PO4_south', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_south.long_name = 'southern boundary orthophosphate'
        nc_po4_south.units = 'mmol/m3'
        nc_po4_south.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_south = nc.createVariable('Si_south', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_south.long_name = 'southern boundary silicate'
        nc_si_south.units = 'mmol/m3'
        nc_si_south.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_south = nc.createVariable('O2_south', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_south.long_name = 'southern boundary dissolved oxygen'
        nc_o2_south.units = 'mmol/m3'
        nc_o2_south.coordinates = 'lon_rho s_rho o2_time'
        #
        nc_dic_south = nc.createVariable('DIC_south', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_south.long_name = 'southern boundary dissolved inorganic carbon'
        nc_dic_south.units = 'mmol/m3'
        nc_dic_south.coordinates = 'lon_rho s_rho dic_time'
        #
        nc_fer_south = nc.createVariable('FER_south', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_south.long_name = 'southern boundary iron'
        nc_fer_south.units = 'mmol/m3'
        nc_fer_south.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_talk_south = nc.createVariable('TALK_south', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_south.long_name = 'southern boundary alkalinity'
        nc_talk_south.units = 'mmol/m3'
        nc_talk_south.coordinates = 'lon_rho s_rho talk_time'
        #
    if obc[1] == 1:
        #
        #   Eastern boundary
        #
        nc_temp_east = nc.createVariable('temp_east', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_east.long_name = 'eastern boundary potential temperature'
        nc_temp_east.units = 'Celsius'
        nc_temp_east.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_east = nc.createVariable('salt_east', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_east.long_name = 'eastern boundary salinity'
        nc_salt_east.units = 'PSU'
        nc_salt_east.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_east = nc.createVariable('u_east', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_east.long_name = 'eastern boundary u-momentum component'
        nc_u_east.units = 'meter second-1'
        nc_u_east.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_east = nc.createVariable('v_east', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_east.long_name = 'eastern boundary v-momentum component'
        nc_v_east.units = 'meter second-1'
        nc_v_east.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_east = nc.createVariable('ubar_east', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_east.long_name = 'eastern boundary vertically integrated u-momentum component'
        nc_ubar_east.units = 'meter second-1'
        nc_ubar_east.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_east = nc.createVariable('vbar_east', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_east.long_name = 'eastern boundary vertically integrated v-momentum component'
        nc_vbar_east.units = 'meter second-1'
        nc_vbar_east.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_east = nc.createVariable('zeta_east', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_east.long_name = 'eastern boundary sea surface height'
        nc_zeta_east.units = 'meter'
        nc_zeta_east.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_east = nc.createVariable('NO3_east', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_east.long_name = 'eastern boundary nitrate'
        nc_no3_east.units = 'mmol/m3'
        nc_no3_east.coordinates = 'lat_rho s_rho no3_time'
        #
        nc_po4_east = nc.createVariable('PO4_east', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_east.long_name = 'eastern boundary orthophosphate'
        nc_po4_east.units = 'mmol/m3'
        nc_po4_east.coordinates = 'lat_rho s_rho po4_time'
        #
        nc_si_east = nc.createVariable('Si_east', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_east.long_name = 'eastern boundary silicate'
        nc_si_east.units = 'mmol/m3'
        nc_si_east.coordinates = 'lat_rho s_rho si_time'
        #
        nc_o2_east = nc.createVariable('O2_east', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_east.long_name = 'eastern boundary dissolved oxygen'
        nc_o2_east.units = 'mmol/m3'
        nc_o2_east.coordinates = 'lat_rho s_rho o2_time'
        #
        nc_dic_east = nc.createVariable('DIC_east', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_east.long_name = 'eastern boundary dissolved inorganic carbon'
        nc_dic_east.units = 'mmol/m3'
        nc_dic_east.coordinates = 'lat_rho s_rho dic_time'
        #
        nc_fer_east = nc.createVariable('FER_east', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_east.long_name = 'eastern boundary iron'
        nc_fer_east.units = 'mmol/m3'
        nc_fer_east.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_talk_east = nc.createVariable('TALK_east', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_east.long_name = 'eastern boundary alkalinity'
        nc_talk_east.units = 'mmol/m3'
        nc_talk_east.coordinates = 'lat_rho s_rho talk_time'
        #

    if obc[2] == 1:
        #
        #   Northern boundary
        #
        nc_temp_north = nc.createVariable('temp_north', np.float64, ('temp_time', 's_rho', 'xi_rho',))
        nc_temp_north.long_name = 'northern boundary potential temperature'
        nc_temp_north.units = 'Celsius'
        nc_temp_north.coordinates = 'lon_rho s_rho temp_time'
        #
        nc_salt_north = nc.createVariable('salt_north', np.float64, ('salt_time', 's_rho', 'xi_rho',))
        nc_salt_north.long_name = 'northern boundary salinity'
        nc_salt_north.units = 'PSU'
        nc_salt_north.coordinates = 'lon_rho s_rho salt_time'
        #
        nc_u_north = nc.createVariable('u_north', np.float64, ('v3d_time', 's_rho', 'xi_u',))
        nc_u_north.long_name = 'northern boundary u-momentum component'
        nc_u_north.units = 'meter second-1'
        nc_u_north.coordinates = 'lon_u s_rho u_bry_time'
        #
        nc_v_north = nc.createVariable('v_north', np.float64, ('v3d_time', 's_rho', 'xi_rho',))
        nc_v_north.long_name = 'northern boundary v-momentum component'
        nc_v_north.units = 'meter second-1'
        nc_v_north.coordinates = 'lon_v s_rho vclm_time'
        #
        nc_ubar_north = nc.createVariable('ubar_north', np.float64, ('v2d_time', 'xi_u',))
        nc_ubar_north.long_name = 'northern boundary vertically integrated u-momentum component'
        nc_ubar_north.units = 'meter second-1'
        nc_ubar_north.coordinates = 'lon_u uclm_time'
        #
        nc_vbar_north = nc.createVariable('vbar_north', np.float64, ('v2d_time', 'xi_rho',))
        nc_vbar_north.long_name = 'northern boundary vertically integrated v-momentum component'
        nc_vbar_north.units = 'meter second-1'
        nc_vbar_north.coordinates = 'lon_v vclm_time'

        nc_zeta_north = nc.createVariable('zeta_north', np.float64, ('zeta_time', 'xi_rho',))
        nc_zeta_north.long_name = 'northern boundary sea surface height'
        nc_zeta_north.units = 'meter'
        nc_zeta_north.coordinates = 'lon_rho zeta_time'
        #
        nc_no3_north = nc.createVariable('NO3_north', np.float64, ('no3_time', 's_rho', 'xi_rho',))
        nc_no3_north.long_name = 'northern boundary nitrate'
        nc_no3_north.units = 'mmol/m3'
        nc_no3_north.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_north = nc.createVariable('PO4_north', np.float64, ('po4_time', 's_rho', 'xi_rho',))
        nc_po4_north.long_name = 'northern boundary orthophosphate'
        nc_po4_north.units = 'mmol/m3'
        nc_po4_north.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_north = nc.createVariable('Si_north', np.float64, ('si_time', 's_rho', 'xi_rho',))
        nc_si_north.long_name = 'northern boundary silicate'
        nc_si_north.units = 'mmol/m3'
        nc_si_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_north = nc.createVariable('O2_north', np.float64, ('o2_time', 's_rho', 'xi_rho',))
        nc_o2_north.long_name = 'northern boundary dissolved oxygen'
        nc_o2_north.units = 'mmol/m3'
        nc_o2_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_north = nc.createVariable('DIC_north', np.float64, ('dic_time', 's_rho', 'xi_rho',))
        nc_dic_north.long_name = 'northern boundary dissolved inorganic carbon'
        nc_dic_north.units = 'mmol/m3'
        nc_dic_north.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_north = nc.createVariable('FER_north', np.float64, ('fer_time', 's_rho', 'xi_rho',))
        nc_fer_north.long_name = 'northern boundary iron'
        nc_fer_north.units = 'mmol/m3'
        nc_fer_north.coordinates = 'lon_rho s_rho fer_time'
        #
        nc_talk_north = nc.createVariable('TALK_north', np.float64, ('talk_time', 's_rho', 'xi_rho',))
        nc_talk_north.long_name = 'northern boundary alkalinity'
        nc_talk_north.units = 'mmol/m3'
        nc_talk_north.coordinates = 'lon_rho s_rho talk_time'

    if obc[3] == 1:
        #
        #   Western boundary
        #
        nc_temp_west = nc.createVariable('temp_west', np.float64, ('temp_time', 's_rho', 'eta_rho',))
        nc_temp_west.long_name = 'western boundary potential temperature'
        nc_temp_west.units = 'Celsius'
        nc_temp_west.coordinates = 'lat_rho s_rho temp_time'
        #
        nc_salt_west = nc.createVariable('salt_west', np.float64, ('salt_time', 's_rho', 'eta_rho',))
        nc_salt_west.long_name = 'western boundary salinity'
        nc_salt_west.units = 'PSU'
        nc_salt_west.coordinates = 'lat_rho s_rho salt_time'
        #
        nc_u_west = nc.createVariable('u_west', np.float64, ('v3d_time', 's_rho', 'eta_rho',))
        nc_u_west.long_name = 'western boundary u-momentum component'
        nc_u_west.units = 'meter second-1'
        nc_u_west.coordinates = 'lat_u s_rho u_bry_time'
        #
        nc_v_west = nc.createVariable('v_west', np.float64, ('v3d_time', 's_rho', 'eta_v',))
        nc_v_west.long_name = 'western boundary v-momentum component'
        nc_v_west.units = 'meter second-1'
        nc_v_west.coordinates = 'lat_v s_rho vclm_time'
        #
        nc_ubar_west = nc.createVariable('ubar_west', np.float64, ('v2d_time', 'eta_rho',))
        nc_ubar_west.long_name = 'western boundary vertically integrated u-momentum component'
        nc_ubar_west.units = 'meter second-1'
        nc_ubar_west.coordinates = 'lat_u uclm_time'
        #
        nc_vbar_west = nc.createVariable('vbar_west', np.float64, ('v2d_time', 'eta_v',))
        nc_vbar_west.long_name = 'western boundary vertically integrated v-momentum component'
        nc_vbar_west.units = 'meter second-1'
        nc_vbar_west.coordinates = 'lat_v vclm_time'
        #
        nc_zeta_west = nc.createVariable('zeta_west', np.float64, ('zeta_time', 'eta_rho',))
        nc_zeta_west.long_name = 'western boundary sea surface height'
        nc_zeta_west.units = 'meter'
        nc_zeta_west.coordinates = 'lat_rho zeta_time'
        #
        nc_no3_west = nc.createVariable('NO3_west', np.float64, ('no3_time', 's_rho', 'eta_rho',))
        nc_no3_west.long_name = 'western boundary nitrate'
        nc_no3_west.units = 'mmol/m3'
        nc_no3_west.coordinates = 'lon_rho s_rho no3_time'
        #
        nc_po4_west = nc.createVariable('PO4_west', np.float64, ('po4_time', 's_rho', 'eta_rho',))
        nc_po4_west.long_name = 'western boundary orthophosphate'
        nc_po4_west.units = 'mmol/m3'
        nc_po4_west.coordinates = 'lon_rho s_rho po4_time'
        #
        nc_si_west = nc.createVariable('Si_west', np.float64, ('si_time', 's_rho', 'eta_rho',))
        nc_si_west.long_name = 'western boundary silicate'
        nc_si_west.units = 'mmol/m3'
        nc_si_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_o2_west = nc.createVariable('O2_west', np.float64, ('o2_time', 's_rho', 'eta_rho',))
        nc_o2_west.long_name = 'western boundary dissolved oxygen'
        nc_o2_west.units = 'mmol/m3'
        nc_o2_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_dic_west = nc.createVariable('DIC_west', np.float64, ('dic_time', 's_rho', 'eta_rho',))
        nc_dic_west.long_name = 'western boundary dissolved inorganic carbon'
        nc_dic_west.units = 'mmol/m3'
        nc_dic_west.coordinates = 'lon_rho s_rho si_time'
        #
        nc_fer_west = nc.createVariable('FER_west', np.float64, ('fer_time', 's_rho', 'eta_rho',))
        nc_fer_west.long_name = 'western boundary iron'
        nc_fer_west.units = 'mmol/m3'
        nc_fer_west.coordinates = 'lat_rho s_rho fer_time'
        #
        nc_talk_west = nc.createVariable('TALK_west', np.float64, ('talk_time', 's_rho', 'eta_rho',))
        nc_talk_west.long_name = 'western boundary alkalinity'
        nc_talk_west.units = 'mmol/m3'
        nc_talk_west.coordinates = 'lat_rho s_rho talk_time'

    #
    # Compute S coordinates
    #

    (sc_r, Cs_r, sc_w, Cs_w) = vgrd.scoordinate(theta_s, theta_b, N, hc, vtransform)
    print('vtransform = ' + str(vtransform))

    #
    # Write variables
    #
    nc_spherical[:] = 'T'
    nc_Vtransform[:] = vtransform
    nc_Vstretching[:] = 1
    nc_tstart[:] = np.min(bry_time)
    nc_tend[:] = np.max(bry_time)
    nc_theta_s[:] = theta_s
    nc_theta_b[:] = theta_b
    nc_Tcline[:] = hc
    nc_hc[:] = hc
    nc_sc_r[:] = sc_r
    nc_sc_w[:] = sc_w
    nc_Cs_r[:] = Cs_r
    nc_Cs_w[:] = Cs_w
    nc_temp_time[:] = bry_time
    nc_salt_time[:] = bry_time
    nc_uclm_time[:] = bry_time
    nc_vclm_time[:] = bry_time
    nc_v2d_time[:] = bry_time
    nc_v3d_time[:] = bry_time
    nc_ssh_time[:] = bry_time
    nc_zeta_time[:] = bry_time
    nc_bry_time[:] = bry_time

    nc_no3_time[:] = nor_time
    nc_po4_time[:] = nor_time
    nc_si_time[:] = nor_time
    nc_o2_time[:] = nor_time
    nc_dic_time[:] = nor_time
    nc_talk_time[:] = nor_time
    nc_fer_time[:] = pis_time

    if obc[0] == 1:
        nc_u_south[:] = 0.
        nc_v_south[:] = 0.
        nc_ubar_south[:] = 0.
        nc_vbar_south[:] = 0.
        nc_zeta_south[:] = 0.
        nc_temp_south[:] = 0.
        nc_salt_south[:] = 0.
        nc_no3_south[:] = 0.
        nc_po4_south[:] = 0.
        nc_si_south[:] = 0.
        nc_o2_south[:] = 0.
        nc_dic_south[:] = 0.
        nc_fer_south[:] = 0.
        nc_talk_south[:] = 0.

    if obc[1] == 1:
        nc_u_east[:] = 0.
        nc_v_east[:] = 0.
        nc_ubar_east[:] = 0.
        nc_vbar_east[:] = 0.
        nc_zeta_east[:] = 0.
        nc_temp_east[:] = 0.
        nc_salt_east[:] = 0.
        nc_no3_east[:] = 0.
        nc_po4_east[:] = 0.
        nc_si_east[:] = 0.
        nc_o2_east[:] = 0.
        nc_dic_east[:] = 0.
        nc_fer_east[:] = 0.
        nc_talk_east[:] = 0.

    if obc[2] == 1:
        nc_u_north[:] = 0.
        nc_v_north[:] = 0.
        nc_ubar_north[:] = 0.
        nc_vbar_north[:] = 0.
        nc_zeta_north[:] = 0.
        nc_temp_north[:] = 0.
        nc_salt_north[:] = 0.
        nc_no3_north[:] = 0.
        nc_po4_north[:] = 0.
        nc_si_north[:] = 0.
        nc_o2_north[:] = 0.
        nc_dic_north[:] = 0.
        nc_fer_north[:] = 0.
        nc_talk_north[:] = 0.

    if obc[3] == 1:
        nc_u_west[:] = 0.
        nc_v_west[:] = 0.
        nc_ubar_west[:] = 0.
        nc_vbar_west[:] = 0.
        nc_zeta_west[:] = 0.
        nc_temp_west[:] = 0.
        nc_salt_west[:] = 0.
        nc_no3_west[:] = 0.
        nc_po4_west[:] = 0.
        nc_si_west[:] = 0.
        nc_o2_west[:] = 0.
        nc_dic_west[:] = 0.
        nc_fer_west[:] = 0.
        nc_talk_west[:] = 0.

    nc.close()

    return


#
#
#
# #####################################################################
# ##### END FUNCTION CREATE_BRYFILE ###################################
# #####################################################################
#
#


#
#
# #####################################################################
# #### FUNCTION GET_DELAUNAY_BRY ######################################
# #####################################################################
#
#

def get_delaunay_bry(lon_bry, lat_bry, dl, ncglo):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonmin = np.min(lon_bry) - dl
    lonmax = np.max(lon_bry) + dl
    latmin = np.min(lat_bry) - dl
    latmax = np.max(lat_bry) + dl

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncglo['latitude'][:])
    lonT = np.array(ncglo['longitude'][:])

    iminT_bry = glor.geo_idx(lonmin, lonT)
    imaxT_bry = glor.geo_idx(lonmax, lonT)
    jminT_bry = glor.geo_idx(latmin, latT)
    jmaxT_bry = glor.geo_idx(latmax, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncglo['latitude'][:])
    lonU = np.array(ncglo['longitude'][:])

    iminU_bry = glor.geo_idx(lonmin, lonU)
    imaxU_bry = glor.geo_idx(lonmax, lonU)
    jminU_bry = glor.geo_idx(latmin, latU)
    jmaxU_bry = glor.geo_idx(latmax, latU)

    lonU = lonU[iminU_bry:imaxU_bry]
    latU = latU[jminU_bry:jmaxU_bry]
    (LonU_bry, LatU_bry) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncglo['latitude'][:])
    lonV = np.array(ncglo['longitude'][:])

    iminV_bry = glor.geo_idx(lonmin, lonV)
    imaxV_bry = glor.geo_idx(lonmax, lonV)
    jminV_bry = glor.geo_idx(latmin, latV)
    jmaxV_bry = glor.geo_idx(latmax, latV)

    lonV = lonV[iminV_bry:imaxV_bry]
    latV = latV[jminV_bry:jmaxV_bry]
    (LonV_bry, LatV_bry) = np.meshgrid(lonV, latV)

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU_bry, coefU_bry] = glor.get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV_bry, coefV_bry] = glor.get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry,
                 coefU_bry=coefU_bry, elemU_bry=elemU_bry,
                 coefV_bry=coefV_bry, elemV_bry=elemV_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']
        coefU_bry = data['coefU_bry']
        elemU_bry = data['elemU_bry']
        coefV_bry = data['coefV_bry']
        elemV_bry = data['elemV_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
            LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
            LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry)


def get_delaunay_bry_PISCES_NORESM(lon_bry, lat_bry, dl, ncglo, ncpiso, ncnoro):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonmin = np.min(lon_bry) - dl
    lonmax = np.max(lon_bry) + dl
    latmin = np.min(lat_bry) - dl
    latmax = np.max(lat_bry) + dl

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncglo['latitude'][:])
    lonT = np.array(ncglo['longitude'][:])

    iminT_bry = glor.geo_idx(lonmin, lonT)
    imaxT_bry = glor.geo_idx(lonmax, lonT)
    jminT_bry = glor.geo_idx(latmin, latT)
    jmaxT_bry = glor.geo_idx(latmax, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncglo['latitude'][:])
    lonU = np.array(ncglo['longitude'][:])

    iminU_bry = glor.geo_idx(lonmin, lonU)
    imaxU_bry = glor.geo_idx(lonmax, lonU)
    jminU_bry = glor.geo_idx(latmin, latU)
    jmaxU_bry = glor.geo_idx(latmax, latU)

    lonU = lonU[iminU_bry:imaxU_bry]
    latU = latU[jminU_bry:jmaxU_bry]
    (LonU_bry, LatU_bry) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncglo['latitude'][:])
    lonV = np.array(ncglo['longitude'][:])

    iminV_bry = glor.geo_idx(lonmin, lonV)
    imaxV_bry = glor.geo_idx(lonmax, lonV)
    jminV_bry = glor.geo_idx(latmin, latV)
    jmaxV_bry = glor.geo_idx(latmax, latV)

    lonV = lonV[iminV_bry:imaxV_bry]
    latV = latV[jminV_bry:jmaxV_bry]
    (LonV_bry, LatV_bry) = np.meshgrid(lonV, latV)

    #
    # get PISCES positions and indices
    #
    LatP_bry = np.array(ncpiso['nav_lat'][:])
    LonP_bry = np.array(ncpiso['nav_lon'][:])

    iminP_bry = 0
    imaxP_bry = LonP_bry.shape[1]
    jminP_bry = 0
    jmaxP_bry = LatP_bry.shape[0]

    #
    # get NORESM positions and indices
    #
    LatN_bry = np.array(ncnoro['lat'][:])
    LonN_bry = np.array(ncnoro['lon'][:])

    iminN_bry = 0
    imaxN_bry = LonN_bry.shape[1]
    jminN_bry = 0
    jmaxN_bry = LatN_bry.shape[0]

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU_bry, coefU_bry] = glor.get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV_bry, coefV_bry] = glor.get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from PISCES points to CROCO rho_points...')
        [elemP_bry, coefP_bry] = glor.get_tri_coef(LonP_bry, LatP_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefP_bry, axis=2)
        coefP_bry = coefP_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from NORESM points to CROCO rho_points...')
        [elemN_bry, coefN_bry] = glor.get_tri_coef(LonN_bry, LatN_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefN_bry, axis=2)
        coefN_bry = coefN_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry,
                 coefU_bry=coefU_bry, elemU_bry=elemU_bry,
                 coefV_bry=coefV_bry, elemV_bry=elemV_bry,
                 coefN_bry=coefN_bry, elemN_bry=elemN_bry,
                 coefP_bry=coefP_bry, elemP_bry=elemP_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']
        coefU_bry = data['coefU_bry']
        elemU_bry = data['elemU_bry']
        coefV_bry = data['coefV_bry']
        elemV_bry = data['elemV_bry']
        coefN_bry = data['coefN_bry']
        elemN_bry = data['elemN_bry']
        coefP_bry = data['coefP_bry']
        elemP_bry = data['elemP_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
            LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
            LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
            LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
            LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)


def get_delaunay_bry_BGC_PHY(lon_bry, lat_bry, dl, ncglo, ncpiso, ncnoro):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonmin = np.min(lon_bry) - dl
    lonmax = np.max(lon_bry) + dl
    latmin = np.min(lat_bry) - dl
    latmax = np.max(lat_bry) + dl

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncglo['latitude'][:])
    lonT = np.array(ncglo['longitude'][:])

    iminT_bry = glor.geo_idx(lonmin, lonT)
    imaxT_bry = glor.geo_idx(lonmax, lonT)
    jminT_bry = glor.geo_idx(latmin, latT)
    jmaxT_bry = glor.geo_idx(latmax, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncglo['latitude'][:])
    lonU = np.array(ncglo['longitude'][:])

    iminU_bry = glor.geo_idx(lonmin, lonU)
    imaxU_bry = glor.geo_idx(lonmax, lonU)
    jminU_bry = glor.geo_idx(latmin, latU)
    jmaxU_bry = glor.geo_idx(latmax, latU)

    lonU = lonU[iminU_bry:imaxU_bry]
    latU = latU[jminU_bry:jmaxU_bry]
    (LonU_bry, LatU_bry) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncglo['latitude'][:])
    lonV = np.array(ncglo['longitude'][:])

    iminV_bry = glor.geo_idx(lonmin, lonV)
    imaxV_bry = glor.geo_idx(lonmax, lonV)
    jminV_bry = glor.geo_idx(latmin, latV)
    jmaxV_bry = glor.geo_idx(latmax, latV)

    lonV = lonV[iminV_bry:imaxV_bry]
    latV = latV[jminV_bry:jmaxV_bry]
    (LonV_bry, LatV_bry) = np.meshgrid(lonV, latV)

    #
    # get PISCES positions and indices
    #
    latP = np.array(ncpiso['latitude'][:])
    lonP = np.array(ncpiso['longitude'][:])

    iminP_bry = glor.geo_idx(lonmin, lonP)
    imaxP_bry = glor.geo_idx(lonmax, lonP)
    jminP_bry = glor.geo_idx(latmin, latP)
    jmaxP_bry = glor.geo_idx(latmax, latP)

    lonP = lonP[iminP_bry:imaxP_bry]
    latP = latP[jminP_bry:jmaxP_bry]
    (LonP_bry, LatP_bry) = np.meshgrid(lonP, latP)

    #
    # get NORESM positions and indices
    #
    LatN_bry = np.array(ncnoro['lat'][:])
    LonN_bry = np.array(ncnoro['lon'][:])

    iminN_bry = 0
    imaxN_bry = LonN_bry.shape[1]
    jminN_bry = 0
    jmaxN_bry = LatN_bry.shape[0]

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU_bry, coefU_bry] = glor.get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV_bry, coefV_bry] = glor.get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from PISCES points to CROCO rho_points...')
        [elemP_bry, coefP_bry] = glor.get_tri_coef(LonP_bry, LatP_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefP_bry, axis=2)
        coefP_bry = coefP_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from NORESM points to CROCO rho_points...')
        [elemN_bry, coefN_bry] = glor.get_tri_coef(LonN_bry, LatN_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefN_bry, axis=2)
        coefN_bry = coefN_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry,
                 coefU_bry=coefU_bry, elemU_bry=elemU_bry,
                 coefV_bry=coefV_bry, elemV_bry=elemV_bry,
                 coefN_bry=coefN_bry, elemN_bry=elemN_bry,
                 coefP_bry=coefP_bry, elemP_bry=elemP_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']
        coefU_bry = data['coefU_bry']
        elemU_bry = data['elemU_bry']
        coefV_bry = data['coefV_bry']
        elemV_bry = data['elemV_bry']
        coefN_bry = data['coefN_bry']
        elemN_bry = data['elemN_bry']
        coefP_bry = data['coefP_bry']
        elemP_bry = data['elemP_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
            LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
            LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
            LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
            LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)

def get_delaunay_bry_IBIclim(lon_bry, lat_bry, dlg, dlp, ncglo, ncpiso):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonmin = np.min(lon_bry) - dlg
    lonmax = np.max(lon_bry) + dlg
    latmin = np.min(lat_bry) - dlg
    latmax = np.max(lat_bry) + dlg

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncglo['latitude'][:])
    lonT = np.array(ncglo['longitude'][:])

    iminT_bry = glor.geo_idx(lonmin, lonT)
    imaxT_bry = glor.geo_idx(lonmax, lonT)
    jminT_bry = glor.geo_idx(latmin, latT)
    jmaxT_bry = glor.geo_idx(latmax, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncglo['latitude'][:])
    lonU = np.array(ncglo['longitude'][:])

    iminU_bry = glor.geo_idx(lonmin, lonU)
    imaxU_bry = glor.geo_idx(lonmax, lonU)
    jminU_bry = glor.geo_idx(latmin, latU)
    jmaxU_bry = glor.geo_idx(latmax, latU)

    lonU = lonU[iminU_bry:imaxU_bry]
    latU = latU[jminU_bry:jmaxU_bry]
    (LonU_bry, LatU_bry) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncglo['latitude'][:])
    lonV = np.array(ncglo['longitude'][:])

    iminV_bry = glor.geo_idx(lonmin, lonV)
    imaxV_bry = glor.geo_idx(lonmax, lonV)
    jminV_bry = glor.geo_idx(latmin, latV)
    jmaxV_bry = glor.geo_idx(latmax, latV)

    lonV = lonV[iminV_bry:imaxV_bry]
    latV = latV[jminV_bry:jmaxV_bry]
    (LonV_bry, LatV_bry) = np.meshgrid(lonV, latV)

    #
    # get PISCES positions and indices
    #

    lonminp = np.min(lon_bry) - dlp
    lonmaxp = np.max(lon_bry) + dlp
    latminp = np.min(lat_bry) - dlp
    latmaxp = np.max(lat_bry) + dlp

    latP = np.array(ncpiso['lat'][:])
    lonP = np.array(ncpiso['lon'][:])

    iminP_bry = glor.geo_idx(lonminp, lonP)
    imaxP_bry = glor.geo_idx(lonmaxp, lonP)
    jminP_bry = glor.geo_idx(latminp, latP)
    jmaxP_bry = glor.geo_idx(latmaxp, latP)

    lonP = lonP[iminP_bry:imaxP_bry]
    latP = latP[jminP_bry:jmaxP_bry]
    (LonP_bry, LatP_bry) = np.meshgrid(lonP, latP)

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU_bry, coefU_bry] = glor.get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV_bry, coefV_bry] = glor.get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from PISCES points to CROCO rho_points...')
        [elemP_bry, coefP_bry] = glor.get_tri_coef(LonP_bry, LatP_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefP_bry, axis=2)
        coefP_bry = coefP_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry,
                 coefU_bry=coefU_bry, elemU_bry=elemU_bry,
                 coefV_bry=coefV_bry, elemV_bry=elemV_bry,
                 coefP_bry=coefP_bry, elemP_bry=elemP_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']
        coefU_bry = data['coefU_bry']
        elemU_bry = data['elemU_bry']
        coefV_bry = data['coefV_bry']
        elemV_bry = data['elemV_bry']
        coefP_bry = data['coefP_bry']
        elemP_bry = data['elemP_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
            LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
            LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
            LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)

def get_delaunay_bry_WOA_BROU(lon_bry, lat_bry, dlg, dlp, dlw, dlb, ncglo, ncpiso, ncwoao, ncbrouo):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonming = np.min(lon_bry) - dlg
    lonmaxg = np.max(lon_bry) + dlg
    latming = np.min(lat_bry) - dlg
    latmaxg = np.max(lat_bry) + dlg

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncglo['latitude'][:])
    lonT = np.array(ncglo['longitude'][:])

    iminT_bry = glor.geo_idx(lonming, lonT)
    imaxT_bry = glor.geo_idx(lonmaxg, lonT)
    jminT_bry = glor.geo_idx(latming, latT)
    jmaxT_bry = glor.geo_idx(latmaxg, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncglo['latitude'][:])
    lonU = np.array(ncglo['longitude'][:])

    iminU_bry = glor.geo_idx(lonming, lonU)
    imaxU_bry = glor.geo_idx(lonmaxg, lonU)
    jminU_bry = glor.geo_idx(latming, latU)
    jmaxU_bry = glor.geo_idx(latmaxg, latU)

    lonU = lonU[iminU_bry:imaxU_bry]
    latU = latU[jminU_bry:jmaxU_bry]
    (LonU_bry, LatU_bry) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncglo['latitude'][:])
    lonV = np.array(ncglo['longitude'][:])

    iminV_bry = glor.geo_idx(lonming, lonV)
    imaxV_bry = glor.geo_idx(lonmaxg, lonV)
    jminV_bry = glor.geo_idx(latming, latV)
    jmaxV_bry = glor.geo_idx(latmaxg, latV)

    lonV = lonV[iminV_bry:imaxV_bry]
    latV = latV[jminV_bry:jmaxV_bry]
    (LonV_bry, LatV_bry) = np.meshgrid(lonV, latV)

    #
    # get PISCES positions and indices
    #

    lonminp = np.min(lon_bry) - dlp
    lonmaxp = np.max(lon_bry) + dlp
    latminp = np.min(lat_bry) - dlp
    latmaxp = np.max(lat_bry) + dlp

    latP = np.array(ncpiso['latitude'][:])
    lonP = np.array(ncpiso['longitude'][:])

    iminP_bry = glor.geo_idx(lonminp, lonP)
    imaxP_bry = glor.geo_idx(lonmaxp, lonP)
    jminP_bry = glor.geo_idx(latminp, latP)
    jmaxP_bry = glor.geo_idx(latmaxp, latP)

    lonP = lonP[iminP_bry:imaxP_bry]
    latP = latP[jminP_bry:jmaxP_bry]
    (LonP_bry, LatP_bry) = np.meshgrid(lonP, latP)

    #
    # get WOA positions and indices
    #

    lonminw = np.min(lon_bry) - dlw
    lonmaxw = np.max(lon_bry) + dlw
    latminw = np.min(lat_bry) - dlw
    latmaxw = np.max(lat_bry) + dlw

    latW = np.array(ncwoao['lat'][:])
    lonW = np.array(ncwoao['lon'][:])

    iminW_bry = glor.geo_idx(lonminw, lonW)
    imaxW_bry = glor.geo_idx(lonmaxw, lonW)
    jminW_bry = glor.geo_idx(latminw, latW)
    jmaxW_bry = glor.geo_idx(latmaxw, latW)

    lonW = lonW[iminW_bry:imaxW_bry]
    latW = latW[jminW_bry:jmaxW_bry]
    (LonW_bry, LatW_bry) = np.meshgrid(lonW, latW)

    #
    # get BROULLON positions and indices
    #

    lonminb = np.min(lon_bry) - dlb
    lonmaxb = np.max(lon_bry) + dlb
    latminb = np.min(lat_bry) - dlb
    latmaxb = np.max(lat_bry) + dlb

    latB = np.array(ncbrouo['latitude'][:])
    lonB = np.array(ncbrouo['longitude'][:])

    iminB_bry = glor.geo_idx(lonminb, lonB)
    imaxB_bry = glor.geo_idx(lonmaxb, lonB)
    jminB_bry = glor.geo_idx(latminb, latB)
    jmaxB_bry = glor.geo_idx(latmaxb, latB)

    lonB = lonB[iminB_bry:imaxB_bry]
    latB = latB[jminB_bry:jmaxB_bry]
    (LonB_bry, LatB_bry) = np.meshgrid(lonB, latB)

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU_bry, coefU_bry] = glor.get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV_bry, coefV_bry] = glor.get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from PISCES points to CROCO rho_points...')
        [elemP_bry, coefP_bry] = glor.get_tri_coef(LonP_bry, LatP_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefP_bry, axis=2)
        coefP_bry = coefP_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from WOA points to CROCO rho_points...')
        [elemW_bry, coefW_bry] = glor.get_tri_coef(LonW_bry, LatW_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefW_bry, axis=2)
        coefW_bry = coefW_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from BROULLON points to CROCO rho_points...')
        [elemB_bry, coefB_bry] = glor.get_tri_coef(LonB_bry, LatB_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefB_bry, axis=2)
        coefB_bry = coefB_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry,
                 coefU_bry=coefU_bry, elemU_bry=elemU_bry,
                 coefV_bry=coefV_bry, elemV_bry=elemV_bry,
                 coefP_bry=coefP_bry, elemP_bry=elemP_bry,
                 coefW_bry=coefW_bry, elemW_bry=elemW_bry,
                 coefB_bry=coefB_bry, elemB_bry=elemB_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']
        coefU_bry = data['coefU_bry']
        elemU_bry = data['elemU_bry']
        coefV_bry = data['coefV_bry']
        elemV_bry = data['elemV_bry']
        coefP_bry = data['coefP_bry']
        elemP_bry = data['elemP_bry']
        coefW_bry = data['coefW_bry']
        elemW_bry = data['elemW_bry']
        coefB_bry = data['coefB_bry']
        elemB_bry = data['elemB_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
            LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
            LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
            LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry,
            LonW_bry, LatW_bry, iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry, elemW_bry, coefW_bry,
            LonB_bry, LatB_bry, iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry, elemB_bry, coefB_bry)


def get_delaunay_bry_ESM(lon_bry, lat_bry, dl, lon_esm, lat_esm):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonmin = np.min(lon_bry) - dl
    lonmax = np.max(lon_bry) + dl
    latmin = np.min(lat_bry) - dl
    latmax = np.max(lat_bry) + dl

    #
    # get GLORYS positions and indices at T-points
    #

    latT = lat_esm
    lonT = lon_esm

    iminT_bry = glor.geo_idx(lonmin, lonT)
    imaxT_bry = glor.geo_idx(lonmax, lonT)
    jminT_bry = glor.geo_idx(latmin, latT)
    jmaxT_bry = glor.geo_idx(latmax, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry)

def get_delaunay_bry_GLO_BGC_PHY(lon_bry, lat_bry, dl, ncglo, ncpiso):
    #
    #
    # ###############################################################
    #
    # function get_delaunay_bry(lon_bry,lat_bry,dl,ncglo)
    #
    #   This function computes the delaunay matrices for the interpolations
    #   at the boundaies
    #
    #   Input:
    #
    #   lon_bry      Longitudes of the boundary (vector).
    #   lat_bry      Latitudes of the boundary (vector).
    #   dl           extra extension in deg (real).
    #   ncglo        netcdf structure poiting to the GLORYS file
    #
    #   Output:
    #
    #   LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    #   LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    #   LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    # ###############################################################
    #
    #

    comp_delaunay = 1

    lonmin = np.min(lon_bry) - dl
    lonmax = np.max(lon_bry) + dl
    latmin = np.min(lat_bry) - dl
    latmax = np.max(lat_bry) + dl

    #
    # get GLORYS positions and indices at T-points
    #

    latT = np.array(ncglo['latitude'][:])
    lonT = np.array(ncglo['longitude'][:])

    iminT_bry = glor.geo_idx(lonmin, lonT)
    imaxT_bry = glor.geo_idx(lonmax, lonT)
    jminT_bry = glor.geo_idx(latmin, latT)
    jmaxT_bry = glor.geo_idx(latmax, latT)

    lonT = lonT[iminT_bry:imaxT_bry]
    latT = latT[jminT_bry:jmaxT_bry]
    (LonT_bry, LatT_bry) = np.meshgrid(lonT, latT)

    #
    # get GLORYS positions and indices at U-points
    #

    latU = np.array(ncglo['latitude'][:])
    lonU = np.array(ncglo['longitude'][:])

    iminU_bry = glor.geo_idx(lonmin, lonU)
    imaxU_bry = glor.geo_idx(lonmax, lonU)
    jminU_bry = glor.geo_idx(latmin, latU)
    jmaxU_bry = glor.geo_idx(latmax, latU)

    lonU = lonU[iminU_bry:imaxU_bry]
    latU = latU[jminU_bry:jmaxU_bry]
    (LonU_bry, LatU_bry) = np.meshgrid(lonU, latU)

    #
    # get GLORYS positions and indices at V-points
    #

    latV = np.array(ncglo['latitude'][:])
    lonV = np.array(ncglo['longitude'][:])

    iminV_bry = glor.geo_idx(lonmin, lonV)
    imaxV_bry = glor.geo_idx(lonmax, lonV)
    jminV_bry = glor.geo_idx(latmin, latV)
    jmaxV_bry = glor.geo_idx(latmax, latV)

    lonV = lonV[iminV_bry:imaxV_bry]
    latV = latV[jminV_bry:jmaxV_bry]
    (LonV_bry, LatV_bry) = np.meshgrid(lonV, latV)

    #
    # get PISCES positions and indices
    #
    latP = np.array(ncpiso['latitude'][:])
    lonP = np.array(ncpiso['longitude'][:])

    iminP_bry = glor.geo_idx(lonmin, lonP)
    imaxP_bry = glor.geo_idx(lonmax, lonP)
    jminP_bry = glor.geo_idx(latmin, latP)
    jmaxP_bry = glor.geo_idx(latmax, latP)

    lonP = lonP[iminP_bry:imaxP_bry]
    latP = latP[jminP_bry:jmaxP_bry]
    (LonP_bry, LatP_bry) = np.meshgrid(lonP, latP)

    #
    # Horizontal and vertical interp/extrapolations
    #

    print(' ')
    print(' Interpolations / extrapolations')
    print(' ')

    #
    # Get the 2D interpolation coefficients
    #

    if comp_delaunay == 1:

        print('Compute Delaunay triangulation from GLORYS T-points to CROCO rho_points...')
        [elemT_bry, coefT_bry] = glor.get_tri_coef(LonT_bry, LatT_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefT_bry, axis=2)
        coefT_bry = coefT_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS U-points to CROCO rho_points...')
        [elemU_bry, coefU_bry] = glor.get_tri_coef(LonU_bry, LatU_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefU_bry, axis=2)
        coefU_bry = coefU_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from GLORYS V-points to CROCO rho_points...')
        [elemV_bry, coefV_bry] = glor.get_tri_coef(LonV_bry, LatV_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefV_bry, axis=2)
        coefV_bry = coefV_bry / coefnorm[:, :, np.newaxis]

        print('Compute Delaunay triangulation from PISCES points to CROCO rho_points...')
        [elemP_bry, coefP_bry] = glor.get_tri_coef(LonP_bry, LatP_bry, lon_bry, lat_bry, 1)
        coefnorm = np.sum(coefP_bry, axis=2)
        coefP_bry = coefP_bry / coefnorm[:, :, np.newaxis]

        #
        # Save the Delaunay triangulation matrices
        #

        np.savez('coeffs_bry_GLO_PHY_BGC.npz',
                 coefT_bry=coefT_bry, elemT_bry=elemT_bry,
                 coefU_bry=coefU_bry, elemU_bry=elemU_bry,
                 coefV_bry=coefV_bry, elemV_bry=elemV_bry,
                 coefP_bry=coefP_bry, elemP_bry=elemP_bry)

    else:

        #
        # Load the Delaunay triangulation matrices
        #

        print('Load Delaunay triangulation...')
        data = np.load('coeffs_bry.npz')
        coefT_bry = data['coefT_bry']
        elemT_bry = data['elemT_bry']
        coefU_bry = data['coefU_bry']
        elemU_bry = data['elemU_bry']
        coefV_bry = data['coefV_bry']
        elemV_bry = data['elemV_bry']
        coefP_bry = data['coefP_bry']
        elemP_bry = data['elemP_bry']

    print('Delaunay triangulation done')

    return (LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
            LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
            LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
            LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)


#
#
#
# #####################################################################
# ##### END FUNCTION GET_DELAUNAY_BRY #################################
# #####################################################################
#
#


#
#
# #####################################################################
# #### FUNCTION INTERP_BRY ############################################
# #####################################################################
#
#

def interp_bry(obctype, ncglo, ncbgco, tndx_glo, ncbry, tndx_bry,
               h_bry, theta_s, theta_b, hc, N, vtransform,
               Nzgoodmin, depthp, depthb, angle_bry,
               LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
               LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
               LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
               LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    #
    #
    ################################################################
    # function interp_bry(obctype,ncglo,tndx_glo,ncbry,tndx_bry,\
    #              h_bry,theta_s,theta_b,hc,N,vtransform,\
    #	           Nzgoodmin,depth,angle_bry,\
    #              LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,\
    #              LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,\
    #	           LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry)
    #
    #
    #   Interpolates all the variable for one boundary
    #
    #   Input:
    #
    # obctype,ncglo,tndx_glo,ncbry,tndx_bry,
    # h_bry,theta_s,theta_b,hc,N,vtransform,
    # Nzgoodmin,depth,angle_bry,
    # LonT_bry,LatT_bry,iminT_bry,imaxT_bry,jminT_bry,jmaxT_bry,elemT_bry,coefT_bry,
    # LonU_bry,LatU_bry,iminU_bry,imaxU_bry,jminU_bry,jmaxU_bry,elemU_bry,coefU_bry,
    # LonV_bry,LatV_bry,iminV_bry,imaxV_bry,jminV_bry,jmaxV_bry,elemV_bry,coefV_bry
    #
    #   Output:
    #
    #   ncbry    Netcdf file structure
    #
    ################################################################
    #
    #

    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthp, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthp, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Nitrate
    #
    #

    print('Interpolate Nitrate...')

    no3_bry = glor.interp3d(ncbgco, 'no3', tndx_glo, Nzgoodmin, depthb, z_rho,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Orthophosphate
    #
    #

    print('Interpolate Orthophosphate...')

    po4_bry = glor.interp3d(ncbgco, 'po4', tndx_glo, Nzgoodmin, depthb, z_rho,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Silicate
    #
    #

    print('Interpolate Silicate...')

    si_bry = glor.interp3d(ncbgco, 'si', tndx_glo, Nzgoodmin, depthb, z_rho,
                           iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                           LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Dissolved Oxygen
    #
    #

    print('Interpolate Dissolved Oxygen...')

    o2_bry = glor.interp3d(ncbgco, 'o2', tndx_glo, Nzgoodmin, depthb, z_rho,
                           iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                           LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Dissolved Inorganic Carbon
    #
    #
    #
    # print('Interpolate Dissolved Inorganic Carbon...')
    #
    # dic_bry = glor.interp3d(ncbgco, 'dissic', tndx_glo, Nzgoodmin, depthb, z_rho,
    #                         iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
    #                         LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Ammonium
    #
    #
    #
    # print('Interpolate Ammonium...')
    #
    # nh4_bry = glor.interp3d(ncbgco, 'nh4', tndx_glo, Nzgoodmin, depth, z_rho,
    #                         iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
    #                         LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Iron
    #
    #

    print('Interpolate Iron...')

    fe_bry = glor.interp3d(ncbgco, 'fe', tndx_glo, Nzgoodmin, depthb, z_rho,
                           iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                           LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthp, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if obctype == 's':

        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]
        ncbry['NO3_south'][tndx_bry, :, :] = no3_bry[:, 0, :]
        ncbry['PO4_south'][tndx_bry, :, :] = po4_bry[:, 0, :]
        ncbry['Si_south'][tndx_bry, :, :] = si_bry[:, 0, :]
        ncbry['O2_south'][tndx_bry, :, :] = o2_bry[:, 0, :]
        ncbry['DIC_south'][tndx_bry, :, :] = 2150. * np.ones_like(o2_bry[:, 0, :])
        ncbry['TALK_south'][tndx_bry, :, :] = 2350. * np.ones_like(o2_bry[:, 0, :])
        ncbry['FER_south'][tndx_bry, :, :] = fe_bry[:, 0, :]

    elif obctype == 'n':

        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]
        ncbry['NO3_north'][tndx_bry, :, :] = no3_bry[:, -1, :]
        ncbry['PO4_north'][tndx_bry, :, :] = po4_bry[:, -1, :]
        ncbry['Si_north'][tndx_bry, :, :] = si_bry[:, -1, :]
        ncbry['O2_north'][tndx_bry, :, :] = o2_bry[:, -1, :]
        ncbry['DIC_north'][tndx_bry, :, :] = 2150. * np.ones_like(o2_bry[:, -1, :])
        ncbry['TALK_north'][tndx_bry, :, :] = 2350. * np.ones_like(o2_bry[:, -1, :])
        ncbry['FER_north'][tndx_bry, :, :] = fe_bry[:, -1, :]

    elif obctype == 'e':

        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]
        ncbry['NO3_east'][tndx_bry, :, :] = no3_bry[:, :, 0]
        ncbry['PO4_east'][tndx_bry, :, :] = po4_bry[:, :, 0]
        ncbry['Si_east'][tndx_bry, :, :] = si_bry[:, :, 0]
        ncbry['O2_east'][tndx_bry, :, :] = o2_bry[:, :, 0]
        ncbry['DIC_east'][tndx_bry, :, :] = 2150. * np.ones_like(o2_bry[:, :, 0])
        ncbry['TALK_east'][tndx_bry, :, :] = 2350. * np.ones_like(o2_bry[:, :, 0])
        ncbry['FER_east'][tndx_bry, :, :] = fe_bry[:, :, 0]

    elif obctype == 'w':

        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]
        ncbry['NO3_west'][tndx_bry, :, :] = no3_bry[:, :, -1]
        ncbry['PO4_west'][tndx_bry, :, :] = po4_bry[:, :, -1]
        ncbry['Si_west'][tndx_bry, :, :] = si_bry[:, :, -1]
        ncbry['O2_west'][tndx_bry, :, :] = o2_bry[:, :, -1]
        ncbry['DIC_west'][tndx_bry, :, :] = 2150. * np.ones_like(o2_bry[:, :, -1])
        ncbry['TALK_west'][tndx_bry, :, :] = 2350. * np.ones_like(o2_bry[:, :, -1])
        ncbry['FER_west'][tndx_bry, :, :] = fe_bry[:, :, -1]

    return ncbry


#
#
# #####################################################################
# #### FUNCTION INTERP_BRY ############################################
# #####################################################################
#
#


def interp_bry_PISCES_NORESM(obctype,
                             PNI, Tidxn,
                             ncglo,
                             ncpism1o, ncpiso, ncpisp1o,
                             NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
                             tndx_glo, ncbry, tndx_bry,
                             h_bry, theta_s, theta_b,
                             hc, N, vtransform, Nzgoodmin, depthg, depthp, depthn, angle_bry,
                             LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                             LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                             LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                             LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
                             LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    # interp_bry_PISCES_NORESM('s',
    #
    #                          PISCES_NORESM_interpd, Tininxp,
    #
    #                          ncglo,
    #
    #                          ncpisprio, ncpiscuro, ncpisposo,
    #                          NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
    #
    #                          tndx_glo, ncbry, tndx_bry,
    #                          h_bry, theta_s, theta_b,
    #                          hc, N, vtransform, Nzgoodmin, depth, angle_bry,
    #                          LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
    #                          LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
    #                          LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
    #                          LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
    #                          LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)

    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, 'add_offset')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, 'add_offset')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        print('Interpolate Nitrate...')

        # ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
        # ncnornio = ncnorni.variables
        #
        # no3_bry1 = glor.interp3d(ncnornio, 'no3no2', Tidxn[0], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # no3_bry2 = glor.interp3d(ncnornio, 'no3no2', Tidxn[1], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # no3_bry3 = glor.interp3d(ncnornio, 'no3no2', Tidxn[2], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        # ncnorni.close()

        no3_bry1 = glor.interp3d(ncpism1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        no3_bry2 = glor.interp3d(ncpiso, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        no3_bry3 = glor.interp3d(ncpisp1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        #
        # 3: Orthophosphate
        #
        #

        print('Interpolate Orthophosphate...')

        # ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
        # ncnorpoo = ncnorpo.variables
        #
        # po4_bry1 = glor.interp3d(ncnorpoo, 'po4', Tidxn[0], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # po4_bry2 = glor.interp3d(ncnorpoo, 'po4', Tidxn[1], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # po4_bry3 = glor.interp3d(ncnorpoo, 'po4', Tidxn[2], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        # ncnorpo.close()

        po4_bry1 = glor.interp3d(ncpism1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        po4_bry2 = glor.interp3d(ncpiso, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        po4_bry3 = glor.interp3d(ncpisp1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        #
        # 3: Silicate
        #
        #

        print('Interpolate Silicate...')

        # ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
        # ncnorsio = ncnorsi.variables
        #
        # si_bry1 = glor.interp3d(ncnorsio, 'si', Tidxn[0], Nzgoodmin, depthn, z_rho,
        #                         iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                         LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # si_bry2 = glor.interp3d(ncnorsio, 'si', Tidxn[1], Nzgoodmin, depthn, z_rho,
        #                         iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                         LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # si_bry3 = glor.interp3d(ncnorsio, 'si', Tidxn[2], Nzgoodmin, depthn, z_rho,
        #                         iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                         LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        # ncnorsi.close()

        si_bry1 = glor.interp3d(ncpism1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        si_bry2 = glor.interp3d(ncpiso, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        si_bry3 = glor.interp3d(ncpisp1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        print('Interpolate Dissolved Oxygen...')

        # ncnordo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[3] + NORESM_ending, 'r')
        # ncnordoo = ncnordo.variables
        #
        # o2_bry1 = glor.interp3d(ncnordoo, 'o2', Tidxn[0], Nzgoodmin, depthn, z_rho,
        #                         iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                         LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # o2_bry2 = glor.interp3d(ncnordoo, 'o2', Tidxn[1], Nzgoodmin, depthn, z_rho,
        #                         iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                         LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # o2_bry3 = glor.interp3d(ncnordoo, 'o2', Tidxn[2], Nzgoodmin, depthn, z_rho,
        #                         iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                         LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        # ncnordo.close()

        o2_bry1 = glor.interp3d(ncpism1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        o2_bry2 = glor.interp3d(ncpiso, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        o2_bry3 = glor.interp3d(ncpisp1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        print('Interpolate Dissolved Inorganic Carbon...')

        # ncnordic = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[4] + NORESM_ending, 'r')
        # ncnordico = ncnordic.variables
        #
        # dic_bry1 = glor.interp3d(ncnordico, 'dissic', Tidxn[0], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # dic_bry2 = glor.interp3d(ncnordico, 'dissic', Tidxn[1], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # dic_bry3 = glor.interp3d(ncnordico, 'dissic', Tidxn[2], Nzgoodmin, depthn, z_rho,
        #                          iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                          LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        # ncnordic.close()

        dic_bry1 = glor.interp3d(ncpism1o, 'dic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dic_bry2 = glor.interp3d(ncpiso, 'dic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dic_bry3 = glor.interp3d(ncpisp1o, 'dic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        #
        # 3: Total Alkalinity
        #
        #

        print('Interpolate Total Alkalinity...')

        # ncnoralkalini = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[5] + NORESM_ending, 'r')
        # ncnoralkalinio = ncnoralkalini.variables
        #
        # talk_bry1 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[0], Nzgoodmin, depthn, z_rho,
        #                           iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                           LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # talk_bry2 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[1], Nzgoodmin, depthn, z_rho,
        #                           iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                           LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        #
        # talk_bry3 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[2], Nzgoodmin, depthn, z_rho,
        #                           iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
        #                           LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        # ncnoralkalini.close()

        talk_bry1 = glor.interp3d(ncpism1o, 'alk', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        talk_bry2 = glor.interp3d(ncpiso, 'alk', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        talk_bry3 = glor.interp3d(ncpisp1o, 'alk', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        # PISCES
        #
        #
        # 3: Ammonium
        #
        #

        print('Interpolate Calcite...')

        calc_bry1 = glor.interp3d(ncpism1o, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        calc_bry2 = glor.interp3d(ncpiso, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        calc_bry3 = glor.interp3d(ncpisp1o, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5b: Particulate Organic Carbon from IBI PISCES
        #

        print('Interpolate Particulate Organic Carbon...')

        poc_bry1 = glor.interp3d(ncpism1o, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        poc_bry2 = glor.interp3d(ncpiso, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        poc_bry3 = glor.interp3d(ncpisp1o, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5c: Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Nanophytoplankton...')

        phy_bry1 = glor.interp3d(ncpism1o, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy_bry2 = glor.interp3d(ncpiso, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy_bry3 = glor.interp3d(ncpism1o, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5d: Microzooplankton from IBI PISCES
        #

        print('Interpolate Microzooplankton...')

        zoo_bry1 = glor.interp3d(ncpism1o, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo_bry2 = glor.interp3d(ncpiso, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo_bry3 = glor.interp3d(ncpisp1o, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5e: Dissolved Organic Carbon from IBI PISCES
        #

        print('Interpolate Dissolved Organic Carbon...')

        doc_bry1 = glor.interp3d(ncpism1o, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        doc_bry2 = glor.interp3d(ncpiso, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        doc_bry3 = glor.interp3d(ncpisp1o, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5f: Diatom from IBI PISCES
        #

        print('Interpolate Diatom...')

        phy2_bry1 = glor.interp3d(ncpism1o, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy2_bry2 = glor.interp3d(ncpiso, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy2_bry3 = glor.interp3d(ncpisp1o, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5g: Mesozooplankton from IBI PISCES
        #

        print('Interpolate Mesozooplankton...')

        zoo2_bry1 = glor.interp3d(ncpism1o, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo2_bry2 = glor.interp3d(ncpiso, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo2_bry3 = glor.interp3d(ncpisp1o, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5h: Biogenic Silica from IBI PISCES
        #

        print('Interpolate Biogenic Silica...')

        gsi_bry1 = glor.interp3d(ncpism1o, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        gsi_bry2 = glor.interp3d(ncpiso, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        gsi_bry3 = glor.interp3d(ncpisp1o, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5j: Big Particle Iron from IBI PISCES
        #

        print('Interpolate Big Particle Iron...')

        bfe_bry1 = glor.interp3d(ncpism1o, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        bfe_bry2 = glor.interp3d(ncpiso, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        bfe_bry3 = glor.interp3d(ncpisp1o, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5k: Big Particulate Organic Carbon from IBI PISCES
        #

        print('Interpolate Big Particulate Organic Carbon...')

        goc_bry1 = glor.interp3d(ncpism1o, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        goc_bry2 = glor.interp3d(ncpiso, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        goc_bry3 = glor.interp3d(ncpisp1o, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5l: Iron in the small particles from IBI PISCES
        #

        print('Interpolate Iron in the small particles...')

        sfe_bry1 = glor.interp3d(ncpism1o, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        sfe_bry2 = glor.interp3d(ncpiso, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        sfe_bry3 = glor.interp3d(ncpisp1o, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5m: Iron content of the diatoms from IBI PISCES
        #

        print('Interpolate Iron content of the diatoms...')

        dfe_bry1 = glor.interp3d(ncpism1o, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dfe_bry2 = glor.interp3d(ncpiso, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dfe_bry3 = glor.interp3d(ncpisp1o, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5n: Silicon content of the Diatoms from IBI PISCES
        #

        print('Interpolate Silicon content of the Diatoms...')

        dsi_bry1 = glor.interp3d(ncpism1o, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dsi_bry2 = glor.interp3d(ncpiso, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dsi_bry3 = glor.interp3d(ncpisp1o, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5o: Iron content of the Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Iron content of the Nanophytoplankton...')

        nfe_bry1 = glor.interp3d(ncpism1o, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nfe_bry2 = glor.interp3d(ncpiso, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nfe_bry3 = glor.interp3d(ncpisp1o, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5p: Nanophytoplankton Chlorophyll from IBI PISCES
        #

        print('Interpolate Nanophytoplankton Chlorophyll...')

        nchl_bry1 = glor.interp3d(ncpism1o, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nchl_bry2 = glor.interp3d(ncpiso, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nchl_bry3 = glor.interp3d(ncpisp1o, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5q: Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Diatom Chlorophyll...')

        dchl_bry1 = glor.interp3d(ncpism1o, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dchl_bry2 = glor.interp3d(ncpiso, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dchl_bry3 = glor.interp3d(ncpisp1o, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        if obctype == 's':
            # # NORESM previous month
            # ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]*1000
            # ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]*1000
            # ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]*1000
            # ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]*1000
            # ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :]*1000
            # ncbry['TALK_south'][0, :, :] = talk_bry1[:, 0, :]*1000
            # PISCES previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]
            ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :]
            ncbry['TALK_south'][0, :, :] = talk_bry1[:, 0, :]

            # PISCES previous month
            ncbry['DCHL_south'][0, :, :] = dchl_bry1[:, 0, :]
            ncbry['NCHL_south'][0, :, :] = nchl_bry1[:, 0, :]
            ncbry['NFE_south'][0, :, :] = nfe_bry1[:, 0, :]
            ncbry['DSI_south'][0, :, :] = dsi_bry1[:, 0, :]
            ncbry['DFE_south'][0, :, :] = dfe_bry1[:, 0, :]
            ncbry['SFE_south'][0, :, :] = sfe_bry1[:, 0, :]
            ncbry['GOC_south'][0, :, :] = goc_bry1[:, 0, :]
            ncbry['BFE_south'][0, :, :] = bfe_bry1[:, 0, :]
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['BSI_south'][0, :, :] = gsi_bry1[:, 0, :]
            ncbry['MESO_south'][0, :, :] = zoo2_bry1[:, 0, :]
            ncbry['DIA_south'][0, :, :] = phy2_bry1[:, 0, :]
            ncbry['DOC_south'][0, :, :] = doc_bry1[:, 0, :]
            ncbry['ZOO_south'][0, :, :] = zoo_bry1[:, 0, :]
            ncbry['NANO_south'][0, :, :] = phy_bry1[:, 0, :]
            ncbry['CACO3_south'][0, :, :] = calc_bry1[:, 0, :]
            ncbry['POC_south'][0, :, :] = poc_bry1[:, 0, :]
            ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]

            # # NORESM current month
            # ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]*1000
            # ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]*1000
            # ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]*1000
            # ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]*1000
            # ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :]*1000
            # ncbry['TALK_south'][1, :, :] = talk_bry2[:, 0, :]*1000
            # PISCES current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]
            ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :]
            ncbry['TALK_south'][1, :, :] = talk_bry2[:, 0, :]

            # PISCES current month
            ncbry['DCHL_south'][1, :, :] = dchl_bry2[:, 0, :]
            ncbry['NCHL_south'][1, :, :] = nchl_bry2[:, 0, :]
            ncbry['NFE_south'][1, :, :] = nfe_bry2[:, 0, :]
            ncbry['DSI_south'][1, :, :] = dsi_bry2[:, 0, :]
            ncbry['DFE_south'][1, :, :] = dfe_bry2[:, 0, :]
            ncbry['SFE_south'][1, :, :] = sfe_bry2[:, 0, :]
            ncbry['GOC_south'][1, :, :] = goc_bry2[:, 0, :]
            ncbry['BFE_south'][1, :, :] = bfe_bry2[:, 0, :]
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['BSI_south'][1, :, :] = gsi_bry2[:, 0, :]
            ncbry['MESO_south'][1, :, :] = zoo2_bry2[:, 0, :]
            ncbry['DIA_south'][1, :, :] = phy2_bry2[:, 0, :]
            ncbry['DOC_south'][1, :, :] = doc_bry2[:, 0, :]
            ncbry['ZOO_south'][1, :, :] = zoo_bry2[:, 0, :]
            ncbry['NANO_south'][1, :, :] = phy_bry2[:, 0, :]
            ncbry['CACO3_south'][1, :, :] = calc_bry2[:, 0, :]
            ncbry['POC_south'][1, :, :] = poc_bry2[:, 0, :]
            ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]

            # # NORESM next month
            # ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]*1000
            # ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]*1000
            # ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]*1000
            # ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]*1000
            # ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :]*1000
            # ncbry['TALK_south'][2, :, :] = talk_bry3[:, 0, :]*1000
            # PISCES next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]
            ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :]
            ncbry['TALK_south'][2, :, :] = talk_bry3[:, 0, :]

            # PISCES next month
            ncbry['DCHL_south'][2, :, :] = dchl_bry3[:, 0, :]
            ncbry['NCHL_south'][2, :, :] = nchl_bry3[:, 0, :]
            ncbry['NFE_south'][2, :, :] = nfe_bry3[:, 0, :]
            ncbry['DSI_south'][2, :, :] = dsi_bry3[:, 0, :]
            ncbry['DFE_south'][2, :, :] = dfe_bry3[:, 0, :]
            ncbry['SFE_south'][2, :, :] = sfe_bry3[:, 0, :]
            ncbry['GOC_south'][2, :, :] = goc_bry3[:, 0, :]
            ncbry['BFE_south'][2, :, :] = bfe_bry3[:, 0, :]
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['BSI_south'][2, :, :] = gsi_bry3[:, 0, :]
            ncbry['MESO_south'][2, :, :] = zoo2_bry3[:, 0, :]
            ncbry['DIA_south'][2, :, :] = phy2_bry3[:, 0, :]
            ncbry['DOC_south'][2, :, :] = doc_bry3[:, 0, :]
            ncbry['ZOO_south'][2, :, :] = zoo_bry3[:, 0, :]
            ncbry['NANO_south'][2, :, :] = phy_bry3[:, 0, :]
            ncbry['CACO3_south'][2, :, :] = calc_bry3[:, 0, :]
            ncbry['POC_south'][2, :, :] = poc_bry3[:, 0, :]
            ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]

        elif obctype == 'n':
            # # NORESM previous month
            # ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]*1000
            # ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]*1000
            # ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]*1000
            # ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]*1000
            # ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :]*1000
            # ncbry['TALK_north'][0, :, :] = talk_bry1[:, -1, :]*1000
            # PISCES previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]
            ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :]
            ncbry['TALK_north'][0, :, :] = talk_bry1[:, -1, :]

            # PISCES previous month
            ncbry['DCHL_north'][0, :, :] = dchl_bry1[:, -1, :]
            ncbry['NCHL_north'][0, :, :] = nchl_bry1[:, -1, :]
            ncbry['NFE_north'][0, :, :] = nfe_bry1[:, -1, :]
            ncbry['DSI_north'][0, :, :] = dsi_bry1[:, -1, :]
            ncbry['DFE_north'][0, :, :] = dfe_bry1[:, -1, :]
            ncbry['SFE_north'][0, :, :] = sfe_bry1[:, -1, :]
            ncbry['GOC_north'][0, :, :] = goc_bry1[:, -1, :]
            ncbry['BFE_north'][0, :, :] = bfe_bry1[:, -1, :]
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['BSI_north'][0, :, :] = gsi_bry1[:, -1, :]
            ncbry['MESO_north'][0, :, :] = zoo2_bry1[:, -1, :]
            ncbry['DIA_north'][0, :, :] = phy2_bry1[:, -1, :]
            ncbry['DOC_north'][0, :, :] = doc_bry1[:, -1, :]
            ncbry['ZOO_north'][0, :, :] = zoo_bry1[:, -1, :]
            ncbry['NANO_north'][0, :, :] = phy_bry1[:, -1, :]
            ncbry['CACO3_north'][0, :, :] = calc_bry1[:, -1, :]
            ncbry['POC_north'][0, :, :] = poc_bry1[:, -1, :]
            ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]

            # # NORESM current month
            # ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]*1000
            # ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]*1000
            # ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]*1000
            # ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]*1000
            # ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :]*1000
            # ncbry['TALK_north'][1, :, :] = talk_bry2[:, -1, :]*1000
            # PISCES current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]
            ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :]
            ncbry['TALK_north'][1, :, :] = talk_bry2[:, -1, :]

            # PISCES current month
            ncbry['DCHL_north'][1, :, :] = dchl_bry2[:, -1, :]
            ncbry['NCHL_north'][1, :, :] = nchl_bry2[:, -1, :]
            ncbry['NFE_north'][1, :, :] = nfe_bry2[:, -1, :]
            ncbry['DSI_north'][1, :, :] = dsi_bry2[:, -1, :]
            ncbry['DFE_north'][1, :, :] = dfe_bry2[:, -1, :]
            ncbry['SFE_north'][1, :, :] = sfe_bry2[:, -1, :]
            ncbry['GOC_north'][1, :, :] = goc_bry2[:, -1, :]
            ncbry['BFE_north'][1, :, :] = bfe_bry2[:, -1, :]
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['BSI_north'][1, :, :] = gsi_bry2[:, -1, :]
            ncbry['MESO_north'][1, :, :] = zoo2_bry2[:, -1, :]
            ncbry['DIA_north'][1, :, :] = phy2_bry2[:, -1, :]
            ncbry['DOC_north'][1, :, :] = doc_bry2[:, -1, :]
            ncbry['ZOO_north'][1, :, :] = zoo_bry2[:, -1, :]
            ncbry['NANO_north'][1, :, :] = phy_bry2[:, -1, :]
            ncbry['CACO3_north'][1, :, :] = calc_bry2[:, -1, :]
            ncbry['POC_north'][1, :, :] = poc_bry2[:, -1, :]
            ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]

            # # NORESM next month
            # ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]*1000
            # ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]*1000
            # ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]*1000
            # ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]*1000
            # ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :]*1000
            # ncbry['TALK_north'][2, :, :] = talk_bry3[:, -1, :]*1000
            # PISCES next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]
            ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :]
            ncbry['TALK_north'][2, :, :] = talk_bry3[:, -1, :]

            # PISCES next month
            ncbry['DCHL_north'][2, :, :] = dchl_bry3[:, -1, :]
            ncbry['NCHL_north'][2, :, :] = nchl_bry3[:, -1, :]
            ncbry['NFE_north'][2, :, :] = nfe_bry3[:, -1, :]
            ncbry['DSI_north'][2, :, :] = dsi_bry3[:, -1, :]
            ncbry['DFE_north'][2, :, :] = dfe_bry3[:, -1, :]
            ncbry['SFE_north'][2, :, :] = sfe_bry3[:, -1, :]
            ncbry['GOC_north'][2, :, :] = goc_bry3[:, -1, :]
            ncbry['BFE_north'][2, :, :] = bfe_bry3[:, -1, :]
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['BSI_north'][2, :, :] = gsi_bry3[:, -1, :]
            ncbry['MESO_north'][2, :, :] = zoo2_bry3[:, -1, :]
            ncbry['DIA_north'][2, :, :] = phy2_bry3[:, -1, :]
            ncbry['DOC_north'][2, :, :] = doc_bry3[:, -1, :]
            ncbry['ZOO_north'][2, :, :] = zoo_bry3[:, -1, :]
            ncbry['NANO_north'][2, :, :] = phy_bry3[:, -1, :]
            ncbry['CACO3_north'][2, :, :] = calc_bry3[:, -1, :]
            ncbry['POC_north'][2, :, :] = poc_bry3[:, -1, :]
            ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]

        elif obctype == 'e':
            # # NORESM previous month
            # ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]*1000
            # ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]*1000
            # ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]*1000
            # ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]*1000
            # ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0]*1000
            # ncbry['TALK_east'][0, :, :] = talk_bry1[:, :, 0]*1000
            # PISCES previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]
            ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0]
            ncbry['TALK_east'][0, :, :] = talk_bry1[:, :, 0]

            # PISCES previous month
            ncbry['DCHL_east'][0, :, :] = dchl_bry1[:, :, 0]
            ncbry['NCHL_east'][0, :, :] = nchl_bry1[:, :, 0]
            ncbry['NFE_east'][0, :, :] = nfe_bry1[:, :, 0]
            ncbry['DSI_east'][0, :, :] = dsi_bry1[:, :, 0]
            ncbry['DFE_east'][0, :, :] = dfe_bry1[:, :, 0]
            ncbry['SFE_east'][0, :, :] = sfe_bry1[:, :, 0]
            ncbry['GOC_east'][0, :, :] = goc_bry1[:, :, 0]
            ncbry['BFE_east'][0, :, :] = bfe_bry1[:, :, 0]
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['BSI_east'][0, :, :] = gsi_bry1[:, :, 0]
            ncbry['MESO_east'][0, :, :] = zoo2_bry1[:, :, 0]
            ncbry['DIA_east'][0, :, :] = phy2_bry1[:, :, 0]
            ncbry['DOC_east'][0, :, :] = doc_bry1[:, :, 0]
            ncbry['ZOO_east'][0, :, :] = zoo_bry1[:, :, 0]
            ncbry['NANO_east'][0, :, :] = phy_bry1[:, :, 0]
            ncbry['CACO3_east'][0, :, :] = calc_bry1[:, :, 0]
            ncbry['POC_east'][0, :, :] = poc_bry1[:, :, 0]
            ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]

            # # NORESM current month
            # ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]*1000
            # ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]*1000
            # ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]*1000
            # ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]*1000
            # ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0]*1000
            # ncbry['TALK_east'][1, :, :] = talk_bry2[:, :, 0]*1000
            # PISCES current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]
            ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0]
            ncbry['TALK_east'][1, :, :] = talk_bry2[:, :, 0]

            # PISCES current month
            ncbry['DCHL_east'][1, :, :] = dchl_bry2[:, :, 0]
            ncbry['NCHL_east'][1, :, :] = nchl_bry2[:, :, 0]
            ncbry['NFE_east'][1, :, :] = nfe_bry2[:, :, 0]
            ncbry['DSI_east'][1, :, :] = dsi_bry2[:, :, 0]
            ncbry['DFE_east'][1, :, :] = dfe_bry2[:, :, 0]
            ncbry['SFE_east'][1, :, :] = sfe_bry2[:, :, 0]
            ncbry['GOC_east'][1, :, :] = goc_bry2[:, :, 0]
            ncbry['BFE_east'][1, :, :] = bfe_bry2[:, :, 0]
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['BSI_east'][1, :, :] = gsi_bry2[:, :, 0]
            ncbry['MESO_east'][1, :, :] = zoo2_bry2[:, :, 0]
            ncbry['DIA_east'][1, :, :] = phy2_bry2[:, :, 0]
            ncbry['DOC_east'][1, :, :] = doc_bry2[:, :, 0]
            ncbry['ZOO_east'][1, :, :] = zoo_bry2[:, :, 0]
            ncbry['NANO_east'][1, :, :] = phy_bry2[:, :, 0]
            ncbry['CACO3_east'][1, :, :] = calc_bry2[:, :, 0]
            ncbry['POC_east'][1, :, :] = poc_bry2[:, :, 0]
            ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]

            # # NORESM next month
            # ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]*1000
            # ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]*1000
            # ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]*1000
            # ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]*1000
            # ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0]*1000
            # ncbry['TALK_east'][2, :, :] = talk_bry3[:, :, 0]*1000
            # PISCES next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]
            ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0]
            ncbry['TALK_east'][2, :, :] = talk_bry3[:, :, 0]

            # PISCES next month
            ncbry['DCHL_east'][2, :, :] = dchl_bry3[:, :, 0]
            ncbry['NCHL_east'][2, :, :] = nchl_bry3[:, :, 0]
            ncbry['NFE_east'][2, :, :] = nfe_bry3[:, :, 0]
            ncbry['DSI_east'][2, :, :] = dsi_bry3[:, :, 0]
            ncbry['DFE_east'][2, :, :] = dfe_bry3[:, :, 0]
            ncbry['SFE_east'][2, :, :] = sfe_bry3[:, :, 0]
            ncbry['GOC_east'][2, :, :] = goc_bry3[:, :, 0]
            ncbry['BFE_east'][2, :, :] = bfe_bry3[:, :, 0]
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['BSI_east'][2, :, :] = gsi_bry3[:, :, 0]
            ncbry['MESO_east'][2, :, :] = zoo2_bry3[:, :, 0]
            ncbry['DIA_east'][2, :, :] = phy2_bry3[:, :, 0]
            ncbry['DOC_east'][2, :, :] = doc_bry3[:, :, 0]
            ncbry['ZOO_east'][2, :, :] = zoo_bry3[:, :, 0]
            ncbry['NANO_east'][2, :, :] = phy_bry3[:, :, 0]
            ncbry['CACO3_east'][2, :, :] = calc_bry3[:, :, 0]
            ncbry['POC_east'][2, :, :] = poc_bry3[:, :, 0]
            ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]

        elif obctype == 'w':
            # # NORESM previous month
            # ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]*1000
            # ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]*1000
            # ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]*1000
            # ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]*1000
            # ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1]*1000
            # ncbry['TALK_west'][0, :, :] = talk_bry1[:, :, -1]*1000
            # PISCES previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]
            ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1]
            ncbry['TALK_west'][0, :, :] = talk_bry1[:, :, -1]

            # PISCES previous month
            ncbry['DCHL_west'][0, :, :] = dchl_bry1[:, :, -1]
            ncbry['NCHL_west'][0, :, :] = nchl_bry1[:, :, -1]
            ncbry['NFE_west'][0, :, :] = nfe_bry1[:, :, -1]
            ncbry['DSI_west'][0, :, :] = dsi_bry1[:, :, -1]
            ncbry['DFE_west'][0, :, :] = dfe_bry1[:, :, -1]
            ncbry['SFE_west'][0, :, :] = sfe_bry1[:, :, -1]
            ncbry['GOC_west'][0, :, :] = goc_bry1[:, :, -1]
            ncbry['BFE_west'][0, :, :] = bfe_bry1[:, :, -1]
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['BSI_west'][0, :, :] = gsi_bry1[:, :, -1]
            ncbry['MESO_west'][0, :, :] = zoo2_bry1[:, :, -1]
            ncbry['DIA_west'][0, :, :] = phy2_bry1[:, :, -1]
            ncbry['DOC_west'][0, :, :] = doc_bry1[:, :, -1]
            ncbry['ZOO_west'][0, :, :] = zoo_bry1[:, :, -1]
            ncbry['NANO_west'][0, :, :] = phy_bry1[:, :, -1]
            ncbry['CACO3_west'][0, :, :] = calc_bry1[:, :, -1]
            ncbry['POC_west'][0, :, :] = poc_bry1[:, :, -1]
            ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]

            # # NORESM current month
            # ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]*1000
            # ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]*1000
            # ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]*1000
            # ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]*1000
            # ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1]*1000
            # ncbry['TALK_west'][1, :, :] = talk_bry2[:, :, -1]*1000
            # PISCES current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]
            ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1]
            ncbry['TALK_west'][1, :, :] = talk_bry2[:, :, -1]

            # PISCES current month
            ncbry['DCHL_west'][1, :, :] = dchl_bry2[:, :, -1]
            ncbry['NCHL_west'][1, :, :] = nchl_bry2[:, :, -1]
            ncbry['NFE_west'][1, :, :] = nfe_bry2[:, :, -1]
            ncbry['DSI_west'][1, :, :] = dsi_bry2[:, :, -1]
            ncbry['DFE_west'][1, :, :] = dfe_bry2[:, :, -1]
            ncbry['SFE_west'][1, :, :] = sfe_bry2[:, :, -1]
            ncbry['GOC_west'][1, :, :] = goc_bry2[:, :, -1]
            ncbry['BFE_west'][1, :, :] = bfe_bry2[:, :, -1]
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['BSI_west'][1, :, :] = gsi_bry2[:, :, -1]
            ncbry['MESO_west'][1, :, :] = zoo2_bry2[:, :, -1]
            ncbry['DIA_west'][1, :, :] = phy2_bry2[:, :, -1]
            ncbry['DOC_west'][1, :, :] = doc_bry2[:, :, -1]
            ncbry['ZOO_west'][1, :, :] = zoo_bry2[:, :, -1]
            ncbry['NANO_west'][1, :, :] = phy_bry2[:, :, -1]
            ncbry['CACO3_west'][1, :, :] = calc_bry2[:, :, -1]
            ncbry['POC_west'][1, :, :] = poc_bry2[:, :, -1]
            ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]

            # # NORESM next month
            # ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]*1000
            # ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]*1000
            # ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]*1000
            # ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]*1000
            # ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1]*1000
            # ncbry['TALK_west'][2, :, :] = talk_bry3[:, :, -1]*1000
            # PISCES next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]
            ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1]
            ncbry['TALK_west'][2, :, :] = talk_bry3[:, :, -1]

            # PISCES next month
            ncbry['DCHL_west'][2, :, :] = dchl_bry3[:, :, -1]
            ncbry['NCHL_west'][2, :, :] = nchl_bry3[:, :, -1]
            ncbry['NFE_west'][2, :, :] = nfe_bry3[:, :, -1]
            ncbry['DSI_west'][2, :, :] = dsi_bry3[:, :, -1]
            ncbry['DFE_west'][2, :, :] = dfe_bry3[:, :, -1]
            ncbry['SFE_west'][2, :, :] = sfe_bry3[:, :, -1]
            ncbry['GOC_west'][2, :, :] = goc_bry3[:, :, -1]
            ncbry['BFE_west'][2, :, :] = bfe_bry3[:, :, -1]
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['BSI_west'][2, :, :] = gsi_bry3[:, :, -1]
            ncbry['MESO_west'][2, :, :] = zoo2_bry3[:, :, -1]
            ncbry['DIA_west'][2, :, :] = phy2_bry3[:, :, -1]
            ncbry['DOC_west'][2, :, :] = doc_bry3[:, :, -1]
            ncbry['ZOO_west'][2, :, :] = zoo_bry3[:, :, -1]
            ncbry['NANO_west'][2, :, :] = phy_bry3[:, :, -1]
            ncbry['CACO3_west'][2, :, :] = calc_bry3[:, :, -1]
            ncbry['POC_west'][2, :, :] = poc_bry3[:, :, -1]
            ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry


def interp_bry_BGC_PHY(obctype,
                       PNI, Tidxn,
                       ncglo,
                       ncpism1o, ncpiso, ncpisp1o,
                       ncphymm1o, ncphymo, ncphymp1o,
                       NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
                       tndx_glo, ncbry, tndx_bry,
                       h_bry, theta_s, theta_b,
                       hc, N, vtransform, Nzgoodmin, depthg, depthp, depthn, angle_bry,
                       LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                       LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                       LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                       LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
                       LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(ncpism1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry2 = glor.interp3d(ncpiso, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry3 = glor.interp3d(ncpisp1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Orthophosphate
        #
        #

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(ncpism1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry2 = glor.interp3d(ncpiso, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry3 = glor.interp3d(ncpisp1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Silicate
        #
        #

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(ncpism1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry2 = glor.interp3d(ncpiso, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry3 = glor.interp3d(ncpisp1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(ncpism1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry2 = glor.interp3d(ncpiso, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry3 = glor.interp3d(ncpisp1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        print('Interpolate Dissolved Inorganic Carbon...')

        dic_bry1 = glor.interp3d(ncpism1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        dic_bry2 = glor.interp3d(ncpiso, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        dic_bry3 = glor.interp3d(ncpisp1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: pH (to get TALK later via CO2SYS)
        #
        #

        print('Interpolate Total Alkalinity...')

        ph_bry1 = glor.interp3d(ncpism1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry2 = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry3 = glor.interp3d(ncpisp1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # PISCES
        #

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # Monthly Salinity and Temperature from IBI for calculation of TALK from DIC, pH, temp and salinity
        #

        salm_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        salm_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        salm_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        if obctype == 's':
            # PISCES previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]
            ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
            ncbry['PH_south'][0, :, :] = ph_bry1[:, 0, :]
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]
            ncbry['tmpm_south'][0, :, :] = tmpm_bry1[:, 0, :]
            ncbry['salm_south'][0, :, :] = salm_bry1[:, 0, :]

            # PISCES current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]
            ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
            ncbry['PH_south'][1, :, :] = ph_bry2[:, 0, :]
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]
            ncbry['tmpm_south'][1, :, :] = tmpm_bry2[:, 0, :]
            ncbry['salm_south'][1, :, :] = salm_bry2[:, 0, :]

            # PISCES next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]
            ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
            ncbry['PH_south'][2, :, :] = ph_bry3[:, 0, :]
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]
            ncbry['tmpm_south'][2, :, :] = tmpm_bry3[:, 0, :]
            ncbry['salm_south'][2, :, :] = salm_bry3[:, 0, :]

        elif obctype == 'n':

            # PISCES previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]
            ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
            ncbry['PH_north'][0, :, :] = ph_bry1[:, -1, :]
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]
            ncbry['tmpm_north'][0, :, :] = tmpm_bry1[:, -1, :]
            ncbry['salm_north'][0, :, :] = salm_bry1[:, -1, :]

            # PISCES current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]
            ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
            ncbry['PH_north'][1, :, :] = ph_bry2[:, -1, :]
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]
            ncbry['tmpm_north'][1, :, :] = tmpm_bry2[:, -1, :]
            ncbry['salm_north'][1, :, :] = salm_bry2[:, -1, :]

            # PISCES next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]
            ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
            ncbry['PH_north'][2, :, :] = ph_bry3[:, -1, :]
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]
            ncbry['tmpm_north'][2, :, :] = tmpm_bry3[:, -1, :]
            ncbry['salm_north'][2, :, :] = salm_bry3[:, -1, :]

        elif obctype == 'e':
            # PISCES previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]
            ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
            ncbry['PH_east'][0, :, :] = ph_bry1[:, :, 0]
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]
            ncbry['tmpm_east'][0, :, :] = tmpm_bry1[:, :, 0]
            ncbry['salm_east'][0, :, :] = salm_bry1[:, :, 0]

            # PISCES current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]
            ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
            ncbry['PH_east'][1, :, :] = ph_bry2[:, :, 0]
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]
            ncbry['tmpm_east'][1, :, :] = tmpm_bry2[:, :, 0]
            ncbry['salm_east'][1, :, :] = salm_bry2[:, :, 0]

            # PISCES next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]
            ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
            ncbry['PH_east'][2, :, :] = ph_bry3[:, :, 0]
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]
            ncbry['tmpm_east'][2, :, :] = tmpm_bry3[:, :, 0]
            ncbry['salm_east'][2, :, :] = salm_bry3[:, :, 0]

        elif obctype == 'w':
            # PISCES previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]
            ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
            ncbry['PH_west'][0, :, :] = ph_bry1[:, :, -1]
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]
            ncbry['tmpm_west'][0, :, :] = tmpm_bry1[:, :, -1]
            ncbry['salm_west'][0, :, :] = salm_bry1[:, :, -1]

            # PISCES current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]
            ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
            ncbry['PH_west'][1, :, :] = ph_bry2[:, :, -1]
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]
            ncbry['tmpm_west'][1, :, :] = tmpm_bry2[:, :, -1]
            ncbry['salm_west'][1, :, :] = salm_bry2[:, :, -1]

            # PISCES next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]
            ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
            ncbry['PH_west'][2, :, :] = ph_bry3[:, :, -1]
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]
            ncbry['tmpm_west'][2, :, :] = tmpm_bry3[:, :, -1]
            ncbry['salm_west'][2, :, :] = salm_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry


def interp_bry_mth_BGC_PHY(obctype,
                           Tidxn,
                           ncpism1o, ncpiso, ncpisp1o,
                           ncphymm1o, ncphymo, ncphymp1o,
                           NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
                           tndx_glo, ncbry,
                           h_bry, theta_s, theta_b,
                           hc, N, vtransform, Nzgoodmin, depthg, depthp, depthn, angle_bry,
                           LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                           LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                           LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                           LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
                           LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry1, NzGood1) = glor.interp_tracers(ncphymm1o, 'zos', tndx_glo, -1,
                                               iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                               LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    (zeta_bry2, NzGood2) = glor.interp_tracers(ncphymo, 'zos', tndx_glo, -1,
                                               iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                               LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    (zeta_bry3, NzGood3) = glor.interp_tracers(ncphymp1o, 'zos', tndx_glo, -1,
                                               iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                               LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho1 = vgrd.zlevs(h_bry, zeta_bry1, theta_s, theta_b, hc, N, 'r', vtransform)
    z_w1 = vgrd.zlevs(h_bry, zeta_bry1, theta_s, theta_b, hc, N, 'w', vtransform)

    z_rho2 = vgrd.zlevs(h_bry, zeta_bry2, theta_s, theta_b, hc, N, 'r', vtransform)
    z_w2 = vgrd.zlevs(h_bry, zeta_bry2, theta_s, theta_b, hc, N, 'w', vtransform)

    z_rho3 = vgrd.zlevs(h_bry, zeta_bry3, theta_s, theta_b, hc, N, 'r', vtransform)
    z_w3 = vgrd.zlevs(h_bry, zeta_bry3, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    temp_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    temp_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    salt_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    salt_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry1, v_bry1] = glor.interp3d_uv(ncphymm1o, tndx_glo, Nzgoodmin, depthg, z_rho1, cosa, sina,
                                        iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                        LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                        iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                        LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    [u_bry2, v_bry2] = glor.interp3d_uv(ncphymo, tndx_glo, Nzgoodmin, depthg, z_rho2, cosa, sina,
                                        iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                        LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                        iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                        LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    [u_bry3, v_bry3] = glor.interp3d_uv(ncphymp1o, tndx_glo, Nzgoodmin, depthg, z_rho3, cosa, sina,
                                        iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                        LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                        iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                        LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry1, h0_1) = vgrd.vintegr(u_bry1, rho2u_3d(z_w1), rho2u_3d(z_rho1), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry1, h0_1) = vgrd.vintegr(v_bry1, rho2v_3d(z_w1), rho2v_3d(z_rho1), np.nan, np.nan) / rho2v_2d(h_bry)

    (ubar_bry2, h0_2) = vgrd.vintegr(u_bry2, rho2u_3d(z_w2), rho2u_3d(z_rho2), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry2, h0_2) = vgrd.vintegr(v_bry2, rho2v_3d(z_w2), rho2v_3d(z_rho2), np.nan, np.nan) / rho2v_2d(h_bry)

    (ubar_bry3, h0_3) = vgrd.vintegr(u_bry3, rho2u_3d(z_w3), rho2u_3d(z_rho3), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry3, h0_3) = vgrd.vintegr(v_bry3, rho2v_3d(z_w3), rho2v_3d(z_rho3), np.nan, np.nan) / rho2v_2d(h_bry)

    # if PNI is 1:
    # # NORESM
    # #
    # #
    # # 3: Nitrate
    # #
    # #

    print('Interpolate Nitrate...')

    no3_bry1 = glor.interp3d(ncpism1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho1,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    no3_bry2 = glor.interp3d(ncpiso, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho2,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    no3_bry3 = glor.interp3d(ncpisp1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho3,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Orthophosphate
    #
    #

    print('Interpolate Orthophosphate...')

    po4_bry1 = glor.interp3d(ncpism1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho1,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    po4_bry2 = glor.interp3d(ncpiso, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho2,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    po4_bry3 = glor.interp3d(ncpisp1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho3,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Silicate
    #
    #

    print('Interpolate Silicate...')

    si_bry1 = glor.interp3d(ncpism1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho1,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    si_bry2 = glor.interp3d(ncpiso, 'si', tndx_glo, Nzgoodmin, depthp, z_rho2,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    si_bry3 = glor.interp3d(ncpisp1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho3,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Dissolved Oxygen
    #
    #

    print('Interpolate Dissolved Oxygen...')

    o2_bry1 = glor.interp3d(ncpism1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho1,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    o2_bry2 = glor.interp3d(ncpiso, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho2,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    o2_bry3 = glor.interp3d(ncpisp1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho3,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: Dissolved Inorganic Carbon
    #
    #

    print('Interpolate Dissolved Inorganic Carbon...')

    dic_bry1 = glor.interp3d(ncpism1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho1,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    dic_bry2 = glor.interp3d(ncpiso, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho2,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    dic_bry3 = glor.interp3d(ncpisp1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho3,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    #
    # 3: pH (to get TALK later via CO2SYS)
    #
    #

    print('Interpolate Total Alkalinity...')

    ph_bry1 = glor.interp3d(ncpism1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho1,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    ph_bry2 = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho2,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    ph_bry3 = glor.interp3d(ncpisp1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho3,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    # PISCES
    #

    #
    # 5i: Dissolved Iron from IBI PISCES
    #

    print('Interpolate Dissolved Iron...')

    fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho1,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho2,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho3,
                            iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                            LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    # 5r: Ammonium from IBI PISCES
    #

    print('Interpolate Ammonium...')

    nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho1,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho2,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho3,
                             iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                             LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

    #
    # Monthly Salinity and Temperature from IBI for calculation of TALK from DIC, pH, temp and salinity
    #

    salm_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    salm_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    salm_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    tmpm_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    tmpm_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    tmpm_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    if obctype == 's':
        # PISCES previous month
        ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]
        ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]
        ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]
        ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]
        ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
        ncbry['PH_south'][0, :, :] = ph_bry1[:, 0, :]
        ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
        ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]
        ncbry['tmpm_south'][0, :, :] = tmpm_bry1[:, 0, :]
        ncbry['salm_south'][0, :, :] = salm_bry1[:, 0, :]

        # PISCES current month
        ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]
        ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]
        ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]
        ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]
        ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
        ncbry['PH_south'][1, :, :] = ph_bry2[:, 0, :]
        ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
        ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]
        ncbry['tmpm_south'][1, :, :] = tmpm_bry2[:, 0, :]
        ncbry['salm_south'][1, :, :] = salm_bry2[:, 0, :]

        # PISCES next month
        ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]
        ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]
        ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]
        ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]
        ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
        ncbry['PH_south'][2, :, :] = ph_bry3[:, 0, :]
        ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
        ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]
        ncbry['tmpm_south'][2, :, :] = tmpm_bry3[:, 0, :]
        ncbry['salm_south'][2, :, :] = salm_bry3[:, 0, :]

    elif obctype == 'n':

        # PISCES previous month
        ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]
        ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]
        ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]
        ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]
        ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
        ncbry['PH_north'][0, :, :] = ph_bry1[:, -1, :]
        ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
        ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]
        ncbry['tmpm_north'][0, :, :] = tmpm_bry1[:, -1, :]
        ncbry['salm_north'][0, :, :] = salm_bry1[:, -1, :]

        # PISCES current month
        ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]
        ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]
        ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]
        ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]
        ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
        ncbry['PH_north'][1, :, :] = ph_bry2[:, -1, :]
        ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
        ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]
        ncbry['tmpm_north'][1, :, :] = tmpm_bry2[:, -1, :]
        ncbry['salm_north'][1, :, :] = salm_bry2[:, -1, :]

        # PISCES next month
        ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]
        ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]
        ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]
        ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]
        ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
        ncbry['PH_north'][2, :, :] = ph_bry3[:, -1, :]
        ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
        ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]
        ncbry['tmpm_north'][2, :, :] = tmpm_bry3[:, -1, :]
        ncbry['salm_north'][2, :, :] = salm_bry3[:, -1, :]

    elif obctype == 'e':
        # PISCES previous month
        ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]
        ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]
        ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]
        ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]
        ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
        ncbry['PH_east'][0, :, :] = ph_bry1[:, :, 0]
        ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
        ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]
        ncbry['tmpm_east'][0, :, :] = tmpm_bry1[:, :, 0]
        ncbry['salm_east'][0, :, :] = salm_bry1[:, :, 0]

        # PISCES current month
        ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]
        ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]
        ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]
        ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]
        ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
        ncbry['PH_east'][1, :, :] = ph_bry2[:, :, 0]
        ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
        ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]
        ncbry['tmpm_east'][1, :, :] = tmpm_bry2[:, :, 0]
        ncbry['salm_east'][1, :, :] = salm_bry2[:, :, 0]

        # PISCES next month
        ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]
        ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]
        ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]
        ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]
        ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
        ncbry['PH_east'][2, :, :] = ph_bry3[:, :, 0]
        ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
        ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]
        ncbry['tmpm_east'][2, :, :] = tmpm_bry3[:, :, 0]
        ncbry['salm_east'][2, :, :] = salm_bry3[:, :, 0]

    elif obctype == 'w':
        # PISCES previous month
        ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]
        ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]
        ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]
        ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]
        ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
        ncbry['PH_west'][0, :, :] = ph_bry1[:, :, -1]
        ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
        ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]
        ncbry['tmpm_west'][0, :, :] = tmpm_bry1[:, :, -1]
        ncbry['salm_west'][0, :, :] = salm_bry1[:, :, -1]

        # PISCES current month
        ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]
        ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]
        ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]
        ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]
        ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
        ncbry['PH_west'][1, :, :] = ph_bry2[:, :, -1]
        ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
        ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]
        ncbry['tmpm_west'][1, :, :] = tmpm_bry2[:, :, -1]
        ncbry['salm_west'][1, :, :] = salm_bry2[:, :, -1]

        # PISCES next month
        ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]
        ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]
        ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]
        ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]
        ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
        ncbry['PH_west'][2, :, :] = ph_bry3[:, :, -1]
        ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
        ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]
        ncbry['tmpm_west'][2, :, :] = tmpm_bry3[:, :, -1]
        ncbry['salm_west'][2, :, :] = salm_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][0, :] = zeta_bry1[0, :]
        ncbry['temp_south'][0, :, :] = temp_bry1[:, 0, :]
        ncbry['salt_south'][0, :, :] = salt_bry1[:, 0, :]
        ncbry['u_south'][0, :, :] = u_bry1[:, 0, :]
        ncbry['v_south'][0, :, :] = v_bry1[:, 0, :]
        ncbry['ubar_south'][0, :] = ubar_bry1[0, :]
        ncbry['vbar_south'][0, :] = vbar_bry1[0, :]

        ncbry['zeta_south'][1, :] = zeta_bry2[0, :]
        ncbry['temp_south'][1, :, :] = temp_bry2[:, 0, :]
        ncbry['salt_south'][1, :, :] = salt_bry2[:, 0, :]
        ncbry['u_south'][1, :, :] = u_bry2[:, 0, :]
        ncbry['v_south'][1, :, :] = v_bry2[:, 0, :]
        ncbry['ubar_south'][1, :] = ubar_bry2[0, :]
        ncbry['vbar_south'][1, :] = vbar_bry2[0, :]

        ncbry['zeta_south'][2, :] = zeta_bry3[0, :]
        ncbry['temp_south'][2, :, :] = temp_bry3[:, 0, :]
        ncbry['salt_south'][2, :, :] = salt_bry3[:, 0, :]
        ncbry['u_south'][2, :, :] = u_bry3[:, 0, :]
        ncbry['v_south'][2, :, :] = v_bry3[:, 0, :]
        ncbry['ubar_south'][2, :] = ubar_bry3[0, :]
        ncbry['vbar_south'][2, :] = vbar_bry3[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][0, :] = zeta_bry1[-1, :]
        ncbry['temp_north'][0, :, :] = temp_bry1[:, -1, :]
        ncbry['salt_north'][0, :, :] = salt_bry1[:, -1, :]
        ncbry['u_north'][0, :, :] = u_bry1[:, -1, :]
        ncbry['v_north'][0, :, :] = v_bry1[:, -1, :]
        ncbry['ubar_north'][0, :] = ubar_bry1[-1, :]
        ncbry['vbar_north'][0, :] = vbar_bry1[-1, :]

        ncbry['zeta_north'][1, :] = zeta_bry2[-1, :]
        ncbry['temp_north'][1, :, :] = temp_bry2[:, -1, :]
        ncbry['salt_north'][1, :, :] = salt_bry2[:, -1, :]
        ncbry['u_north'][1, :, :] = u_bry2[:, -1, :]
        ncbry['v_north'][1, :, :] = v_bry2[:, -1, :]
        ncbry['ubar_north'][1, :] = ubar_bry2[-1, :]
        ncbry['vbar_north'][1, :] = vbar_bry2[-1, :]

        ncbry['zeta_north'][2, :] = zeta_bry3[-1, :]
        ncbry['temp_north'][2, :, :] = temp_bry3[:, -1, :]
        ncbry['salt_north'][2, :, :] = salt_bry3[:, -1, :]
        ncbry['u_north'][2, :, :] = u_bry3[:, -1, :]
        ncbry['v_north'][2, :, :] = v_bry3[:, -1, :]
        ncbry['ubar_north'][2, :] = ubar_bry3[-1, :]
        ncbry['vbar_north'][2, :] = vbar_bry3[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][0, :] = zeta_bry1[:, 0]
        ncbry['temp_east'][0, :, :] = temp_bry1[:, :, 0]
        ncbry['salt_east'][0, :, :] = salt_bry1[:, :, 0]
        ncbry['u_east'][0, :, :] = u_bry1[:, :, 0]
        ncbry['v_east'][0, :, :] = v_bry1[:, :, 0]
        ncbry['ubar_east'][0, :] = ubar_bry1[:, 0]
        ncbry['vbar_east'][0, :] = vbar_bry1[:, 0]

        ncbry['zeta_east'][1, :] = zeta_bry2[:, 0]
        ncbry['temp_east'][1, :, :] = temp_bry2[:, :, 0]
        ncbry['salt_east'][1, :, :] = salt_bry2[:, :, 0]
        ncbry['u_east'][1, :, :] = u_bry2[:, :, 0]
        ncbry['v_east'][1, :, :] = v_bry2[:, :, 0]
        ncbry['ubar_east'][1, :] = ubar_bry2[:, 0]
        ncbry['vbar_east'][1, :] = vbar_bry2[:, 0]

        ncbry['zeta_east'][2, :] = zeta_bry3[:, 0]
        ncbry['temp_east'][2, :, :] = temp_bry3[:, :, 0]
        ncbry['salt_east'][2, :, :] = salt_bry3[:, :, 0]
        ncbry['u_east'][2, :, :] = u_bry3[:, :, 0]
        ncbry['v_east'][2, :, :] = v_bry3[:, :, 0]
        ncbry['ubar_east'][2, :] = ubar_bry3[:, 0]
        ncbry['vbar_east'][2, :] = vbar_bry3[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][0, :] = zeta_bry1[:, -1]
        ncbry['temp_west'][0, :, :] = temp_bry1[:, :, -1]
        ncbry['salt_west'][0, :, :] = salt_bry1[:, :, -1]
        ncbry['u_west'][0, :, :] = u_bry1[:, :, -1]
        ncbry['v_west'][0, :, :] = v_bry1[:, :, -1]
        ncbry['ubar_west'][0, :] = ubar_bry1[:, -1]
        ncbry['vbar_west'][0, :] = vbar_bry1[:, -1]

        ncbry['zeta_west'][1, :] = zeta_bry2[:, -1]
        ncbry['temp_west'][1, :, :] = temp_bry2[:, :, -1]
        ncbry['salt_west'][1, :, :] = salt_bry2[:, :, -1]
        ncbry['u_west'][1, :, :] = u_bry2[:, :, -1]
        ncbry['v_west'][1, :, :] = v_bry2[:, :, -1]
        ncbry['ubar_west'][1, :] = ubar_bry2[:, -1]
        ncbry['vbar_west'][1, :] = vbar_bry2[:, -1]

        ncbry['zeta_west'][2, :] = zeta_bry3[:, -1]
        ncbry['temp_west'][2, :, :] = temp_bry3[:, :, -1]
        ncbry['salt_west'][2, :, :] = salt_bry3[:, :, -1]
        ncbry['u_west'][2, :, :] = u_bry3[:, :, -1]
        ncbry['v_west'][2, :, :] = v_bry3[:, :, -1]
        ncbry['ubar_west'][2, :] = ubar_bry3[:, -1]
        ncbry['vbar_west'][2, :] = vbar_bry3[:, -1]

    return ncbry

def interp_bry_BFM(obctype,
                   PNI, Tidxn,
                   ncglo,
                   ncpism1o, ncpiso, ncpisp1o,
                   ncphymm1o, ncphymo, ncphymp1o,
                   NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
                   tndx_glo, ncbry, tndx_bry,
                   h_bry, theta_s, theta_b,
                   hc, N, vtransform, Nzgoodmin, depthg, depthp, depthn, angle_bry,
                   LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                   LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                   LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                   LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
                   LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(ncpism1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry2 = glor.interp3d(ncpiso, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry3 = glor.interp3d(ncpisp1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Orthophosphate
        #
        #

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(ncpism1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry2 = glor.interp3d(ncpiso, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry3 = glor.interp3d(ncpisp1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Silicate
        #
        #

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(ncpism1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry2 = glor.interp3d(ncpiso, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry3 = glor.interp3d(ncpisp1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(ncpism1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry2 = glor.interp3d(ncpiso, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry3 = glor.interp3d(ncpisp1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        print('Interpolate Dissolved Inorganic Carbon...')

        dic_bry1 = glor.interp3d(ncpism1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        dic_bry2 = glor.interp3d(ncpiso, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        dic_bry3 = glor.interp3d(ncpisp1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: pH (to get TALK later via CO2SYS)
        #
        #

        print('Interpolate Total Alkalinity...')

        ph_bry1 = glor.interp3d(ncpism1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry2 = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry3 = glor.interp3d(ncpisp1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # PISCES
        #

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # Monthly Salinity and Temperature from IBI for calculation of TALK from DIC, pH, temp and salinity
        #

        salm_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        salm_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        salm_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        if obctype == 's':
            # PISCES previous month
            ncbry['N3n_south'][0, :, :] = no3_bry1[:, 0, :]
            ncbry['N1p_south'][0, :, :] = po4_bry1[:, 0, :]
            ncbry['N5s_south'][0, :, :] = si_bry1[:, 0, :]
            ncbry['O2o_south'][0, :, :] = o2_bry1[:, 0, :]
            ncbry['O3c_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
            ncbry['PH_south'][0, :, :] = ph_bry1[:, 0, :]
            ncbry['N7f_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['N4n_south'][0, :, :] = nh4_bry1[:, 0, :]
            ncbry['tmpm_south'][0, :, :] = tmpm_bry1[:, 0, :]
            ncbry['salm_south'][0, :, :] = salm_bry1[:, 0, :]

            # PISCES current month
            ncbry['N3n_south'][1, :, :] = no3_bry2[:, 0, :]
            ncbry['N1p_south'][1, :, :] = po4_bry2[:, 0, :]
            ncbry['N5s_south'][1, :, :] = si_bry2[:, 0, :]
            ncbry['O2o_south'][1, :, :] = o2_bry2[:, 0, :]
            ncbry['O3c_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
            ncbry['PH_south'][1, :, :] = ph_bry2[:, 0, :]
            ncbry['N7f_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['N4n_south'][1, :, :] = nh4_bry2[:, 0, :]
            ncbry['tmpm_south'][1, :, :] = tmpm_bry2[:, 0, :]
            ncbry['salm_south'][1, :, :] = salm_bry2[:, 0, :]

            # PISCES next month
            ncbry['N3n_south'][2, :, :] = no3_bry3[:, 0, :]
            ncbry['N1p_south'][2, :, :] = po4_bry3[:, 0, :]
            ncbry['N5s_south'][2, :, :] = si_bry3[:, 0, :]
            ncbry['O2o_south'][2, :, :] = o2_bry3[:, 0, :]
            ncbry['O3c_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
            ncbry['PH_south'][2, :, :] = ph_bry3[:, 0, :]
            ncbry['N7f_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['N4n_south'][2, :, :] = nh4_bry3[:, 0, :]
            ncbry['tmpm_south'][2, :, :] = tmpm_bry3[:, 0, :]
            ncbry['salm_south'][2, :, :] = salm_bry3[:, 0, :]

        elif obctype == 'n':

            # PISCES previous month
            ncbry['N3n_north'][0, :, :] = no3_bry1[:, -1, :]
            ncbry['N1p_north'][0, :, :] = po4_bry1[:, -1, :]
            ncbry['N5s_north'][0, :, :] = si_bry1[:, -1, :]
            ncbry['O2o_north'][0, :, :] = o2_bry1[:, -1, :]
            ncbry['O3c_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
            ncbry['PH_north'][0, :, :] = ph_bry1[:, -1, :]
            ncbry['N7f_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['N4n_north'][0, :, :] = nh4_bry1[:, -1, :]
            ncbry['tmpm_north'][0, :, :] = tmpm_bry1[:, -1, :]
            ncbry['salm_north'][0, :, :] = salm_bry1[:, -1, :]

            # PISCES current month
            ncbry['N3n_north'][1, :, :] = no3_bry2[:, -1, :]
            ncbry['N1p_north'][1, :, :] = po4_bry2[:, -1, :]
            ncbry['N5s_north'][1, :, :] = si_bry2[:, -1, :]
            ncbry['O2o_north'][1, :, :] = o2_bry2[:, -1, :]
            ncbry['O3c_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
            ncbry['PH_north'][1, :, :] = ph_bry2[:, -1, :]
            ncbry['N7f_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['N4n_north'][1, :, :] = nh4_bry2[:, -1, :]
            ncbry['tmpm_north'][1, :, :] = tmpm_bry2[:, -1, :]
            ncbry['salm_north'][1, :, :] = salm_bry2[:, -1, :]

            # PISCES next month
            ncbry['N3n_north'][2, :, :] = no3_bry3[:, -1, :]
            ncbry['N1p_north'][2, :, :] = po4_bry3[:, -1, :]
            ncbry['N5s_north'][2, :, :] = si_bry3[:, -1, :]
            ncbry['O2o_north'][2, :, :] = o2_bry3[:, -1, :]
            ncbry['O3c_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
            ncbry['PH_north'][2, :, :] = ph_bry3[:, -1, :]
            ncbry['N7f_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['N4n_north'][2, :, :] = nh4_bry3[:, -1, :]
            ncbry['tmpm_north'][2, :, :] = tmpm_bry3[:, -1, :]
            ncbry['salm_north'][2, :, :] = salm_bry3[:, -1, :]

        elif obctype == 'e':
            # PISCES previous month
            ncbry['N3n_east'][0, :, :] = no3_bry1[:, :, 0]
            ncbry['N1p_east'][0, :, :] = po4_bry1[:, :, 0]
            ncbry['N5s_east'][0, :, :] = si_bry1[:, :, 0]
            ncbry['O2o_east'][0, :, :] = o2_bry1[:, :, 0]
            ncbry['O3c_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
            ncbry['PH_east'][0, :, :] = ph_bry1[:, :, 0]
            ncbry['N7f_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['N4n_east'][0, :, :] = nh4_bry1[:, :, 0]
            ncbry['tmpm_east'][0, :, :] = tmpm_bry1[:, :, 0]
            ncbry['salm_east'][0, :, :] = salm_bry1[:, :, 0]

            # PISCES current month
            ncbry['N3n_east'][1, :, :] = no3_bry2[:, :, 0]
            ncbry['N1p_east'][1, :, :] = po4_bry2[:, :, 0]
            ncbry['N5s_east'][1, :, :] = si_bry2[:, :, 0]
            ncbry['O2o_east'][1, :, :] = o2_bry2[:, :, 0]
            ncbry['O3c_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
            ncbry['PH_east'][1, :, :] = ph_bry2[:, :, 0]
            ncbry['N7f_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['N4n_east'][1, :, :] = nh4_bry2[:, :, 0]
            ncbry['tmpm_east'][1, :, :] = tmpm_bry2[:, :, 0]
            ncbry['salm_east'][1, :, :] = salm_bry2[:, :, 0]

            # PISCES next month
            ncbry['N3n_east'][2, :, :] = no3_bry3[:, :, 0]
            ncbry['N1p_east'][2, :, :] = po4_bry3[:, :, 0]
            ncbry['N5s_east'][2, :, :] = si_bry3[:, :, 0]
            ncbry['O2o_east'][2, :, :] = o2_bry3[:, :, 0]
            ncbry['O3c_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
            ncbry['PH_east'][2, :, :] = ph_bry3[:, :, 0]
            ncbry['N7f_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['N4n_east'][2, :, :] = nh4_bry3[:, :, 0]
            ncbry['tmpm_east'][2, :, :] = tmpm_bry3[:, :, 0]
            ncbry['salm_east'][2, :, :] = salm_bry3[:, :, 0]

        elif obctype == 'w':
            # PISCES previous month
            ncbry['N3n_west'][0, :, :] = no3_bry1[:, :, -1]
            ncbry['N1p_west'][0, :, :] = po4_bry1[:, :, -1]
            ncbry['N5s_west'][0, :, :] = si_bry1[:, :, -1]
            ncbry['O2o_west'][0, :, :] = o2_bry1[:, :, -1]
            ncbry['O3c_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
            ncbry['PH_west'][0, :, :] = ph_bry1[:, :, -1]
            ncbry['N7f_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['N4n_west'][0, :, :] = nh4_bry1[:, :, -1]
            ncbry['tmpm_west'][0, :, :] = tmpm_bry1[:, :, -1]
            ncbry['salm_west'][0, :, :] = salm_bry1[:, :, -1]

            # PISCES current month
            ncbry['N3n_west'][1, :, :] = no3_bry2[:, :, -1]
            ncbry['N1p_west'][1, :, :] = po4_bry2[:, :, -1]
            ncbry['N5s_west'][1, :, :] = si_bry2[:, :, -1]
            ncbry['O2o_west'][1, :, :] = o2_bry2[:, :, -1]
            ncbry['O3c_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
            ncbry['PH_west'][1, :, :] = ph_bry2[:, :, -1]
            ncbry['N7f_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['N4n_west'][1, :, :] = nh4_bry2[:, :, -1]
            ncbry['tmpm_west'][1, :, :] = tmpm_bry2[:, :, -1]
            ncbry['salm_west'][1, :, :] = salm_bry2[:, :, -1]

            # PISCES next month
            ncbry['N3n_west'][2, :, :] = no3_bry3[:, :, -1]
            ncbry['N1p_west'][2, :, :] = po4_bry3[:, :, -1]
            ncbry['N5s_west'][2, :, :] = si_bry3[:, :, -1]
            ncbry['O2o_west'][2, :, :] = o2_bry3[:, :, -1]
            ncbry['O3c_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
            ncbry['PH_west'][2, :, :] = ph_bry3[:, :, -1]
            ncbry['N7f_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['N4n_west'][2, :, :] = nh4_bry3[:, :, -1]
            ncbry['tmpm_west'][2, :, :] = tmpm_bry3[:, :, -1]
            ncbry['salm_west'][2, :, :] = salm_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry

def interp_bry_IBIclim(obctype,
                       PNI, Tidxn,
                       ncglo,
                       ncpism1o, ncpiso, ncpisp1o,
                       ncphymm1o, ncphymo, ncphymp1o,
                       tndx_glo, ncbry, tndx_bry,
                       h_bry, theta_s, theta_b,
                       hc, N, vtransform, Nzgoodmin, depthg, depthp, angle_bry,
                       LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                       LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                       LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                       LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(ncpism1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry2 = glor.interp3d(ncpiso, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry3 = glor.interp3d(ncpisp1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Orthophosphate
        #
        #

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(ncpism1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry2 = glor.interp3d(ncpiso, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry3 = glor.interp3d(ncpisp1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Silicate
        #
        #

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(ncpism1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry2 = glor.interp3d(ncpiso, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry3 = glor.interp3d(ncpisp1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(ncpism1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry2 = glor.interp3d(ncpiso, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry3 = glor.interp3d(ncpisp1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        print('Interpolate Dissolved Inorganic Carbon...')

        dic_bry1 = glor.interp3d(ncpism1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        dic_bry2 = glor.interp3d(ncpiso, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        dic_bry3 = glor.interp3d(ncpisp1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: pH (to get TALK later via CO2SYS)
        #
        #

        print('Interpolate Total Alkalinity...')

        ph_bry1 = glor.interp3d(ncpism1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry2 = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry3 = glor.interp3d(ncpisp1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # PISCES
        #

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # Monthly Salinity and Temperature from IBI for calculation of TALK from DIC, pH, temp and salinity
        #

        salm_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        salm_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        salm_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        tmpm_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        tmpm_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        tmpm_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        if obctype == 's':
            # PISCES previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]
            ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
            ncbry['PH_south'][0, :, :] = ph_bry1[:, 0, :]
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]
            ncbry['tmpm_south'][0, :, :] = tmpm_bry1[:, 0, :]
            ncbry['salm_south'][0, :, :] = salm_bry1[:, 0, :]

            # PISCES current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]
            ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
            ncbry['PH_south'][1, :, :] = ph_bry2[:, 0, :]
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]
            ncbry['tmpm_south'][1, :, :] = tmpm_bry2[:, 0, :]
            ncbry['salm_south'][1, :, :] = salm_bry2[:, 0, :]

            # PISCES next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]
            ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
            ncbry['PH_south'][2, :, :] = ph_bry3[:, 0, :]
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]
            ncbry['tmpm_south'][2, :, :] = tmpm_bry3[:, 0, :]
            ncbry['salm_south'][2, :, :] = salm_bry3[:, 0, :]

        elif obctype == 'n':

            # PISCES previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]
            ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
            ncbry['PH_north'][0, :, :] = ph_bry1[:, -1, :]
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]
            ncbry['tmpm_north'][0, :, :] = tmpm_bry1[:, -1, :]
            ncbry['salm_north'][0, :, :] = salm_bry1[:, -1, :]

            # PISCES current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]
            ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
            ncbry['PH_north'][1, :, :] = ph_bry2[:, -1, :]
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]
            ncbry['tmpm_north'][1, :, :] = tmpm_bry2[:, -1, :]
            ncbry['salm_north'][1, :, :] = salm_bry2[:, -1, :]

            # PISCES next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]
            ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
            ncbry['PH_north'][2, :, :] = ph_bry3[:, -1, :]
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]
            ncbry['tmpm_north'][2, :, :] = tmpm_bry3[:, -1, :]
            ncbry['salm_north'][2, :, :] = salm_bry3[:, -1, :]

        elif obctype == 'e':
            # PISCES previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]
            ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
            ncbry['PH_east'][0, :, :] = ph_bry1[:, :, 0]
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]
            ncbry['tmpm_east'][0, :, :] = tmpm_bry1[:, :, 0]
            ncbry['salm_east'][0, :, :] = salm_bry1[:, :, 0]

            # PISCES current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]
            ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
            ncbry['PH_east'][1, :, :] = ph_bry2[:, :, 0]
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]
            ncbry['tmpm_east'][1, :, :] = tmpm_bry2[:, :, 0]
            ncbry['salm_east'][1, :, :] = salm_bry2[:, :, 0]

            # PISCES next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]
            ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
            ncbry['PH_east'][2, :, :] = ph_bry3[:, :, 0]
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]
            ncbry['tmpm_east'][2, :, :] = tmpm_bry3[:, :, 0]
            ncbry['salm_east'][2, :, :] = salm_bry3[:, :, 0]

        elif obctype == 'w':
            # PISCES previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]
            ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
            ncbry['PH_west'][0, :, :] = ph_bry1[:, :, -1]
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]
            ncbry['tmpm_west'][0, :, :] = tmpm_bry1[:, :, -1]
            ncbry['salm_west'][0, :, :] = salm_bry1[:, :, -1]

            # PISCES current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]
            ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
            ncbry['PH_west'][1, :, :] = ph_bry2[:, :, -1]
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]
            ncbry['tmpm_west'][1, :, :] = tmpm_bry2[:, :, -1]
            ncbry['salm_west'][1, :, :] = salm_bry2[:, :, -1]

            # PISCES next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]
            ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
            ncbry['PH_west'][2, :, :] = ph_bry3[:, :, -1]
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]
            ncbry['tmpm_west'][2, :, :] = tmpm_bry3[:, :, -1]
            ncbry['salm_west'][2, :, :] = salm_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry

def interp_bry_WOU_BROU(obctype,
                        PNI, Tidxn,
                        ncglo,
                        ncpism1o, ncpiso, ncpisp1o,
                        woan_preo, woan_curo, woan_poso,
                        woap_preo, woap_curo, woap_poso,
                        woas_preo, woas_curo, woas_poso,
                        woao_preo, woao_curo, woao_poso,
                        ncbroudo, ncbrouto, brou_pre, brou_cur, brou_pos,
                        ncphymm1o, ncphymo, ncphymp1o,
                        tndx_glo, ncbry, tndx_bry,
                        h_bry, theta_s, theta_b,
                        hc, N, vtransform, Nzgoodmin,
                        depthg, depthp, depthwn, depthwo, depthbd,
                        angle_bry,
                        LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                        LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                        LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                        LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry,
                        LonW_bry, LatW_bry, iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry, elemW_bry, coefW_bry,
                        LonB_bry, LatB_bry, iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry, elemB_bry, coefB_bry):
    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(woan_preo, 'n_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                 iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                 LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        no3_bry2 = glor.interp3d(woan_curo, 'n_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                 iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                 LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        no3_bry3 = glor.interp3d(woan_poso, 'n_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                 iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                 LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        #
        #
        # 3: Orthophosphate
        #
        #

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(woap_preo, 'p_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                 iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                 LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        po4_bry2 = glor.interp3d(woap_curo, 'p_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                 iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                 LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        po4_bry3 = glor.interp3d(woap_poso, 'p_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                 iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                 LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        #
        #
        # 3: Silicate
        #
        #

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(woas_preo, 'i_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        si_bry2 = glor.interp3d(woas_curo, 'i_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        si_bry3 = glor.interp3d(woas_poso, 'i_an', tndx_glo, Nzgoodmin, depthwn, z_rho,
                                iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(woao_preo, 'o_an', tndx_glo, Nzgoodmin, depthwo, z_rho,
                                iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        o2_bry2 = glor.interp3d(woao_curo, 'o_an', tndx_glo, Nzgoodmin, depthwo, z_rho,
                                iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        o2_bry3 = glor.interp3d(woao_poso, 'o_an', tndx_glo, Nzgoodmin, depthwo, z_rho,
                                iminW_bry, imaxW_bry, jminW_bry, jmaxW_bry,
                                LonW_bry, LatW_bry, coefW_bry, elemW_bry, '_FillValue')

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        print('Interpolate Dissolved Inorganic Carbon...')

        dic_bry1 = glor.interp3d(ncbroudo, 'TCO2_NNGv2LDEO', brou_pre, Nzgoodmin, depthbd, z_rho,
                                 iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry,
                                 LonB_bry, LatB_bry, coefB_bry, elemB_bry, [])

        dic_bry2 = glor.interp3d(ncbroudo, 'TCO2_NNGv2LDEO', brou_cur, Nzgoodmin, depthbd, z_rho,
                                 iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry,
                                 LonB_bry, LatB_bry, coefB_bry, elemB_bry, [])

        dic_bry3 = glor.interp3d(ncbroudo, 'TCO2_NNGv2LDEO', brou_pos, Nzgoodmin, depthbd, z_rho,
                                 iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry,
                                 LonB_bry, LatB_bry, coefB_bry, elemB_bry, [])

        #
        #
        # 3: Total Alkalinity
        #
        #

        print('Interpolate Total Alkalinity...')

        tlk_bry1 = glor.interp3d(ncbrouto, 'AT_NNGv2', brou_pre, Nzgoodmin, depthbd, z_rho,
                                 iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry,
                                 LonB_bry, LatB_bry, coefB_bry, elemB_bry, [])

        tlk_bry2 = glor.interp3d(ncbrouto, 'AT_NNGv2', brou_cur, Nzgoodmin, depthbd, z_rho,
                                 iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry,
                                 LonB_bry, LatB_bry, coefB_bry, elemB_bry, [])

        tlk_bry3 = glor.interp3d(ncbrouto, 'AT_NNGv2', brou_pos, Nzgoodmin, depthbd, z_rho,
                                 iminB_bry, imaxB_bry, jminB_bry, jmaxB_bry,
                                 LonB_bry, LatB_bry, coefB_bry, elemB_bry, [])

        #
        #
        # 3: pH (to get TALK later via CO2SYS)
        #
        #

        print('Interpolate Total Alkalinity...')

        ph_bry1 = glor.interp3d(ncpism1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry2 = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        ph_bry3 = glor.interp3d(ncpisp1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # PISCES
        #

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # Monthly Salinity and Temperature from IBI for calculation of TALK from DIC, pH, temp and salinity
        #

        salm_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        salm_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        salm_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        tmpm_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                                  iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                  LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        if obctype == 's':
            # PISCES previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]*1.025
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]*1.025
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]*1.025
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]*1.025
            ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :]*1.025
            ncbry['TALK_south'][0, :, :] = tlk_bry1[:, 0, :]*1.025
            ncbry['PH_south'][0, :, :] = ph_bry1[:, 0, :]
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]
            ncbry['tmpm_south'][0, :, :] = tmpm_bry1[:, 0, :]
            ncbry['salm_south'][0, :, :] = salm_bry1[:, 0, :]

            # PISCES current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]*1.025
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]*1.025
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]*1.025
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]*1.025
            ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :]*1.025
            ncbry['TALK_south'][1, :, :] = tlk_bry2[:, 0, :]*1.025
            ncbry['PH_south'][1, :, :] = ph_bry2[:, 0, :]
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]
            ncbry['tmpm_south'][1, :, :] = tmpm_bry2[:, 0, :]
            ncbry['salm_south'][1, :, :] = salm_bry2[:, 0, :]

            # PISCES next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]*1.025
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]*1.025
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]*1.025
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]*1.025
            ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :]*1.025
            ncbry['TALK_south'][2, :, :] = tlk_bry3[:, 0, :]*1.025
            ncbry['PH_south'][2, :, :] = ph_bry3[:, 0, :]
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]
            ncbry['tmpm_south'][2, :, :] = tmpm_bry3[:, 0, :]
            ncbry['salm_south'][2, :, :] = salm_bry3[:, 0, :]

        elif obctype == 'n':

            # PISCES previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]*1.025
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]*1.025
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]*1.025
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]*1.025
            ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :]*1.025
            ncbry['TALK_north'][0, :, :] = tlk_bry1[:, -1, :]*1.025
            ncbry['PH_north'][0, :, :] = ph_bry1[:, -1, :]
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]
            ncbry['tmpm_north'][0, :, :] = tmpm_bry1[:, -1, :]
            ncbry['salm_north'][0, :, :] = salm_bry1[:, -1, :]

            # PISCES current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]*1.025
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]*1.025
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]*1.025
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]*1.025
            ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :]*1.025
            ncbry['TALK_north'][1, :, :] = tlk_bry2[:, -1, :]*1.025
            ncbry['PH_north'][1, :, :] = ph_bry2[:, -1, :]
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]
            ncbry['tmpm_north'][1, :, :] = tmpm_bry2[:, -1, :]
            ncbry['salm_north'][1, :, :] = salm_bry2[:, -1, :]

            # PISCES next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]*1.025
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]*1.025
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]*1.025
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]*1.025
            ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :]*1.025
            ncbry['TALK_north'][2, :, :] = tlk_bry3[:, -1, :]*1.025
            ncbry['PH_north'][2, :, :] = ph_bry3[:, -1, :]
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]
            ncbry['tmpm_north'][2, :, :] = tmpm_bry3[:, -1, :]
            ncbry['salm_north'][2, :, :] = salm_bry3[:, -1, :]

        elif obctype == 'e':
            # PISCES previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]*1.025
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]*1.025
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]*1.025
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]*1.025
            ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0]*1.025
            ncbry['TALK_east'][0, :, :] = tlk_bry1[:, :, 0]*1.025
            ncbry['PH_east'][0, :, :] = ph_bry1[:, :, 0]
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]
            ncbry['tmpm_east'][0, :, :] = tmpm_bry1[:, :, 0]
            ncbry['salm_east'][0, :, :] = salm_bry1[:, :, 0]

            # PISCES current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]*1.025
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]*1.025
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]*1.025
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]*1.025
            ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0]*1.025
            ncbry['TALK_east'][1, :, :] = tlk_bry2[:, :, 0]*1.025
            ncbry['PH_east'][1, :, :] = ph_bry2[:, :, 0]
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]
            ncbry['tmpm_east'][1, :, :] = tmpm_bry2[:, :, 0]
            ncbry['salm_east'][1, :, :] = salm_bry2[:, :, 0]

            # PISCES next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]*1.025
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]*1.025
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]*1.025
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]*1.025
            ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0]*1.025
            ncbry['TALK_east'][2, :, :] = tlk_bry3[:, :, 0]*1.025
            ncbry['PH_east'][2, :, :] = ph_bry3[:, :, 0]
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]
            ncbry['tmpm_east'][2, :, :] = tmpm_bry3[:, :, 0]
            ncbry['salm_east'][2, :, :] = salm_bry3[:, :, 0]

        elif obctype == 'w':
            # PISCES previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]*1.025
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]*1.025
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]*1.025
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]*1.025
            ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1]*1.025
            ncbry['TALK_west'][0, :, :] = tlk_bry1[:, :, -1]*1.025
            ncbry['PH_west'][0, :, :] = ph_bry1[:, :, -1]
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]
            ncbry['tmpm_west'][0, :, :] = tmpm_bry1[:, :, -1]
            ncbry['salm_west'][0, :, :] = salm_bry1[:, :, -1]

            # PISCES current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]*1.025
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]*1.025
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]*1.025
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]*1.025
            ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1]*1.025
            ncbry['TALK_west'][1, :, :] = tlk_bry2[:, :, -1]*1.025
            ncbry['PH_west'][1, :, :] = ph_bry2[:, :, -1]
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]
            ncbry['tmpm_west'][1, :, :] = tmpm_bry2[:, :, -1]
            ncbry['salm_west'][1, :, :] = salm_bry2[:, :, -1]

            # PISCES next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]*1.025
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]*1.025
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]*1.025
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]*1.025
            ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1]*1.025
            ncbry['TALK_west'][2, :, :] = tlk_bry3[:, :, -1]*1.025
            ncbry['PH_west'][2, :, :] = ph_bry3[:, :, -1]
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]
            ncbry['tmpm_west'][2, :, :] = tmpm_bry3[:, :, -1]
            ncbry['salm_west'][2, :, :] = salm_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry


def interp_bry_ESM(obctype,
                   zfiles_3m, tempfiles_3m, saltfiles_3m, ufiles_3m, vfiles_3m,
                   nitfiles_3m, po4files_3m, sifiles_3m, o2files_3m, fefiles_3m,
                   dicfiles_3m, talkfiles_3m,
                   zfiles, tempfiles, saltfiles, ufiles, vfiles,
                   nitfiles, po4files, sifiles, o2files, fefiles,
                   dicfiles, talkfiles,
                   iyear, imonth,
                   ncbry,
                   h_bry, theta_s, theta_b,
                   hc, N, vtransform, Nzgoodmin, depth, angle_bry,
                   LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry):
    #
    #
    # 1: SSH
    #
    #

    Tinin2 = datetime(iyear, imonth, 15)

    if imonth is 1:
        Tinin1 = datetime(iyear-1, 12, 15)
        Tinin3 = datetime(iyear, imonth+1, 15)
    elif imonth is 12:
        Tinin1 = datetime(iyear, imonth-1 , 15)
        Tinin3 = datetime(iyear+1, 1 , 15)
    else:
        Tinin1 = datetime(iyear, imonth-1, 15)
        Tinin3 = datetime(iyear, imonth+1, 15)

    z1 = netcdf(zfiles[zfiles_3m[0]], 'r')
    z1o = z1.variables
    zininx1 = d2i(Tinin1, z1o['time'], select='nearest')
    z2 = netcdf(zfiles[zfiles_3m[1]], 'r')
    z2o = z2.variables
    zininx2 = d2i(Tinin2, z2o['time'], select='nearest')
    z3 = netcdf(zfiles[zfiles_3m[2]], 'r')
    z3o = z3.variables
    zininx3 = d2i(Tinin3, z3o['time'], select='nearest')

    print('Interpolate SSH...')
    (zeta_bry1, NzGood) = glor.interp_tracers(z1o, 'zos', zininx1, -1,
                                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    (zeta_bry2, NzGood) = glor.interp_tracers(z2o, 'zos', zininx2, -1,
                                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    (zeta_bry3, NzGood) = glor.interp_tracers(z3o, 'zos', zininx3, -1,
                                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #
    z_rho1 = vgrd.zlevs(h_bry, zeta_bry1, theta_s, theta_b, hc, N, 'r', vtransform)
    z_w1 = vgrd.zlevs(h_bry, zeta_bry1, theta_s, theta_b, hc, N, 'w', vtransform)
    z_rho2 = vgrd.zlevs(h_bry, zeta_bry2, theta_s, theta_b, hc, N, 'r', vtransform)
    z_w2 = vgrd.zlevs(h_bry, zeta_bry2, theta_s, theta_b, hc, N, 'w', vtransform)
    z_rho3 = vgrd.zlevs(h_bry, zeta_bry3, theta_s, theta_b, hc, N, 'r', vtransform)
    z_w3 = vgrd.zlevs(h_bry, zeta_bry3, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #
    t1 = netcdf(tempfiles[tempfiles_3m[0]], 'r')
    t1o = t1.variables
    tininx1 = d2i(Tinin1, t1o['time'], select='nearest')
    t2 = netcdf(tempfiles[tempfiles_3m[1]], 'r')
    t2o = t2.variables
    tininx2 = d2i(Tinin2, t2o['time'], select='nearest')
    t3 = netcdf(tempfiles[tempfiles_3m[2]], 'r')
    t3o = t3.variables
    tininx3 = d2i(Tinin3, t3o['time'], select='nearest')

    print('Interpolate Temperature...')
    temp_bry1 = glor.interp3d(t1o, 'thetao', tininx1, Nzgoodmin, depth, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    temp_bry2 = glor.interp3d(t2o, 'thetao', tininx2, Nzgoodmin, depth, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    temp_bry3 = glor.interp3d(t3o, 'thetao', tininx3, Nzgoodmin, depth, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #
    s1 = netcdf(saltfiles[saltfiles_3m[0]], 'r')
    s1o = s1.variables
    sininx1 = d2i(Tinin1, s1o['time'], select='nearest')
    s2 = netcdf(saltfiles[saltfiles_3m[1]], 'r')
    s2o = s2.variables
    sininx2 = d2i(Tinin2, s2o['time'], select='nearest')
    s3 = netcdf(saltfiles[saltfiles_3m[2]], 'r')
    s3o = s3.variables
    sininx3 = d2i(Tinin3, s3o['time'], select='nearest')

    print('Interpolate Salinity...')

    salt_bry1 = glor.interp3d(s1o, 'so', sininx1, Nzgoodmin, depth, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    salt_bry2 = glor.interp3d(s2o, 'so', sininx2, Nzgoodmin, depth, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    salt_bry3 = glor.interp3d(s3o, 'so', sininx3, Nzgoodmin, depth, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #
    u1 = netcdf(ufiles[ufiles_3m[0]], 'r')
    u1o = u1.variables
    uininx1 = d2i(Tinin1, u1o['time'], select='nearest')
    u2 = netcdf(ufiles[ufiles_3m[1]], 'r')
    u2o = u2.variables
    uininx2 = d2i(Tinin2, u2o['time'], select='nearest')
    u3 = netcdf(ufiles[ufiles_3m[2]], 'r')
    u3o = u3.variables
    uininx3 = d2i(Tinin3, u3o['time'], select='nearest')
    v1 = netcdf(vfiles[vfiles_3m[0]], 'r')
    v1o = v1.variables
    v2 = netcdf(vfiles[vfiles_3m[1]], 'r')
    v2o = v2.variables
    v3 = netcdf(vfiles[vfiles_3m[2]], 'r')
    v3o = v3.variables

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry1, v_bry1] = glor.interp3d_uv_ESM(u1o, v1o, uininx1, Nzgoodmin, depth, z_rho1, cosa, sina,
                                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                            LonT_bry, LatT_bry, coefT_bry, elemT_bry,
                                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                            LonT_bry, LatT_bry, coefT_bry, elemT_bry)
    [u_bry2, v_bry2] = glor.interp3d_uv_ESM(u2o, v2o, uininx2, Nzgoodmin, depth, z_rho2, cosa, sina,
                                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                            LonT_bry, LatT_bry, coefT_bry, elemT_bry,
                                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                            LonT_bry, LatT_bry, coefT_bry, elemT_bry)
    [u_bry3, v_bry3] = glor.interp3d_uv_ESM(u3o, v3o, uininx3, Nzgoodmin, depth, z_rho3, cosa, sina,
                                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                            LonT_bry, LatT_bry, coefT_bry, elemT_bry,
                                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                            LonT_bry, LatT_bry, coefT_bry, elemT_bry)
    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    (ubar_bry1, h0) = vgrd.vintegr(u_bry1, rho2u_3d(z_w1), rho2u_3d(z_rho1), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry1, h0) = vgrd.vintegr(v_bry1, rho2v_3d(z_w1), rho2v_3d(z_rho1), np.nan, np.nan) / rho2v_2d(h_bry)
    (ubar_bry2, h0) = vgrd.vintegr(u_bry2, rho2u_3d(z_w2), rho2u_3d(z_rho2), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry2, h0) = vgrd.vintegr(v_bry2, rho2v_3d(z_w2), rho2v_3d(z_rho2), np.nan, np.nan) / rho2v_2d(h_bry)
    (ubar_bry3, h0) = vgrd.vintegr(u_bry3, rho2u_3d(z_w3), rho2u_3d(z_rho3), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry3, h0) = vgrd.vintegr(v_bry3, rho2v_3d(z_w3), rho2v_3d(z_rho3), np.nan, np.nan) / rho2v_2d(h_bry)

    #
    #
    # 3: Nitrate
    #
    #
    n1 = netcdf(nitfiles[nitfiles_3m[0]], 'r')
    n1o = n1.variables
    nininx1 = d2i(Tinin1, n1o['time'], select='nearest')
    n2 = netcdf(nitfiles[nitfiles_3m[1]], 'r')
    n2o = n2.variables
    nininx2 = d2i(Tinin2, n2o['time'], select='nearest')
    n3 = netcdf(nitfiles[nitfiles_3m[2]], 'r')
    n3o = n3.variables
    nininx3 = d2i(Tinin3, n3o['time'], select='nearest')

    print('Interpolate Nitrate...')

    no3_bry1 = glor.interp3d(n1o, 'no3', nininx1, Nzgoodmin, depth, z_rho1,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    no3_bry2 = glor.interp3d(n2o, 'no3', nininx2, Nzgoodmin, depth, z_rho2,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    no3_bry3 = glor.interp3d(n3o, 'no3', nininx3, Nzgoodmin, depth, z_rho3,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Orthophosphate
    #
    #
    p1 = netcdf(po4files[po4files_3m[0]], 'r')
    p1o = p1.variables
    pininx1 = d2i(Tinin1, p1o['time'], select='nearest')
    p2 = netcdf(po4files[po4files_3m[1]], 'r')
    p2o = p2.variables
    pininx2 = d2i(Tinin2, p2o['time'], select='nearest')
    p3 = netcdf(po4files[po4files_3m[2]], 'r')
    p3o = p3.variables
    pininx3 = d2i(Tinin3, p3o['time'], select='nearest')

    print('Interpolate Orthophosphate...')

    po4_bry1 = glor.interp3d(p1o, 'po4', pininx1, Nzgoodmin, depth, z_rho1,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    po4_bry2 = glor.interp3d(p2o, 'po4', pininx2, Nzgoodmin, depth, z_rho2,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    po4_bry3 = glor.interp3d(p3o, 'po4', pininx3, Nzgoodmin, depth, z_rho3,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Silicate
    #
    #
    si1 = netcdf(sifiles[sifiles_3m[0]], 'r')
    si1o = si1.variables
    siininx1 = d2i(Tinin1, si1o['time'], select='nearest')
    si2 = netcdf(sifiles[sifiles_3m[1]], 'r')
    si2o = si2.variables
    siininx2 = d2i(Tinin2, si2o['time'], select='nearest')
    si3 = netcdf(sifiles[sifiles_3m[2]], 'r')
    si3o = si3.variables
    siininx3 = d2i(Tinin3, si3o['time'], select='nearest')

    print('Interpolate Silicate...')

    si_bry1 = glor.interp3d(si1o, 'si', siininx1, Nzgoodmin, depth, z_rho1,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    si_bry2 = glor.interp3d(si2o, 'si', siininx2, Nzgoodmin, depth, z_rho2,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    si_bry3 = glor.interp3d(si3o, 'si', siininx3, Nzgoodmin, depth, z_rho3,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Dissolved Oxygen
    #
    #
    o21 = netcdf(o2files[o2files_3m[0]], 'r')
    o21o = o21.variables
    o2ininx1 = d2i(Tinin1, o21o['time'], select='nearest')
    o22 = netcdf(o2files[o2files_3m[1]], 'r')
    o22o = o22.variables
    o2ininx2 = d2i(Tinin2, o22o['time'], select='nearest')
    o23 = netcdf(o2files[o2files_3m[2]], 'r')
    o23o = o23.variables
    o2ininx3 = d2i(Tinin3, o23o['time'], select='nearest')

    print('Interpolate Dissolved Oxygen...')

    o2_bry1 = glor.interp3d(o21o, 'o2', o2ininx1, Nzgoodmin, depth, z_rho1,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    o2_bry2 = glor.interp3d(o22o, 'o2', o2ininx2, Nzgoodmin, depth, z_rho2,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    o2_bry3 = glor.interp3d(o23o, 'o2', o2ininx3, Nzgoodmin, depth, z_rho3,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Dissolved Inorganic Carbon
    #
    #
    dic1 = netcdf(dicfiles[dicfiles_3m[0]], 'r')
    dic1o = dic1.variables
    dicininx1 = d2i(Tinin1, dic1o['time'], select='nearest')
    dic2 = netcdf(dicfiles[dicfiles_3m[1]], 'r')
    dic2o = dic2.variables
    dicininx2 = d2i(Tinin2, dic2o['time'], select='nearest')
    dic3 = netcdf(dicfiles[dicfiles_3m[2]], 'r')
    dic3o = dic3.variables
    dicininx3 = d2i(Tinin3, dic3o['time'], select='nearest')

    print('Interpolate Dissolved Inorganic Carbon...')

    dic_bry1 = glor.interp3d(dic1o, 'dissic', dicininx1, Nzgoodmin, depth, z_rho1,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    dic_bry2 = glor.interp3d(dic2o, 'dissic', dicininx2, Nzgoodmin, depth, z_rho2,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    dic_bry3 = glor.interp3d(dic3o, 'dissic', dicininx3, Nzgoodmin, depth, z_rho3,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: pH (to get TALK later via CO2SYS)
    #
    #
    talk1 = netcdf(talkfiles[talkfiles_3m[0]], 'r')
    talk1o = talk1.variables
    talkininx1 = d2i(Tinin1, talk1o['time'], select='nearest')
    talk2 = netcdf(talkfiles[talkfiles_3m[1]], 'r')
    talk2o = talk2.variables
    talkininx2 = d2i(Tinin2, talk2o['time'], select='nearest')
    talk3 = netcdf(talkfiles[talkfiles_3m[2]], 'r')
    talk3o = talk3.variables
    talkininx3 = d2i(Tinin3, talk3o['time'], select='nearest')

    print('Interpolate Total Alkalinity...')

    talk_bry1 = glor.interp3d(talk1o, 'talk', talkininx1, Nzgoodmin, depth, z_rho1,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    talk_bry2 = glor.interp3d(talk2o, 'talk', talkininx2, Nzgoodmin, depth, z_rho2,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    talk_bry3 = glor.interp3d(talk3o, 'talk', talkininx3, Nzgoodmin, depth, z_rho3,
                              iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                              LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # PISCES
    #

    #
    # 5i: Dissolved Iron from IBI PISCES
    #
    fe1 = netcdf(fefiles[fefiles_3m[0]], 'r')
    fe1o = fe1.variables
    feininx1 = d2i(Tinin1, fe1o['time'], select='nearest')
    fe2 = netcdf(fefiles[fefiles_3m[1]], 'r')
    fe2o = fe2.variables
    feininx2 = d2i(Tinin2, fe2o['time'], select='nearest')
    fe3 = netcdf(fefiles[fefiles_3m[2]], 'r')
    fe3o = fe3.variables
    feininx3 = d2i(Tinin3, fe3o['time'], select='nearest')

    print('Interpolate Dissolved Iron...')

    fe_bry1 = glor.interp3d(fe1o, 'dfe', feininx1, Nzgoodmin, depth, z_rho1,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    fe_bry2 = glor.interp3d(fe2o, 'dfe', feininx2, Nzgoodmin, depth, z_rho2,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    fe_bry3 = glor.interp3d(fe3o, 'dfe', feininx3, Nzgoodmin, depth, z_rho3,
                            iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                            LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][0, :] = zeta_bry1[0, :]
        ncbry['temp_south'][0, :, :] = temp_bry1[:, 0, :]
        ncbry['salt_south'][0, :, :] = salt_bry1[:, 0, :]
        ncbry['u_south'][0, :, :] = u_bry1[:, 0, :]
        ncbry['v_south'][0, :, :] = v_bry1[:, 0, :]
        ncbry['ubar_south'][0, :] = ubar_bry1[0, :]
        ncbry['vbar_south'][0, :] = vbar_bry1[0, :]

        ncbry['zeta_south'][1, :] = zeta_bry2[0, :]
        ncbry['temp_south'][1, :, :] = temp_bry2[:, 0, :]
        ncbry['salt_south'][1, :, :] = salt_bry2[:, 0, :]
        ncbry['u_south'][1, :, :] = u_bry2[:, 0, :]
        ncbry['v_south'][1, :, :] = v_bry2[:, 0, :]
        ncbry['ubar_south'][1, :] = ubar_bry2[0, :]
        ncbry['vbar_south'][1, :] = vbar_bry2[0, :]

        ncbry['zeta_south'][2, :] = zeta_bry3[0, :]
        ncbry['temp_south'][2, :, :] = temp_bry3[:, 0, :]
        ncbry['salt_south'][2, :, :] = salt_bry3[:, 0, :]
        ncbry['u_south'][2, :, :] = u_bry3[:, 0, :]
        ncbry['v_south'][2, :, :] = v_bry3[:, 0, :]
        ncbry['ubar_south'][2, :] = ubar_bry3[0, :]
        ncbry['vbar_south'][2, :] = vbar_bry3[0, :]

        # PISCES previous month
        ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :] * 1000
        ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :] * 1000
        ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :] * 1000
        ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :] * 1000
        ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
        ncbry['TALK_south'][0, :, :] = talk_bry1[:, 0, :] * 1000
        ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :] * 1000

        # PISCES current month
        ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :] * 1000
        ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :] * 1000
        ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :] * 1000
        ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :] * 1000
        ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
        ncbry['TALK_south'][1, :, :] = talk_bry2[:, 0, :] * 1000
        ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :] * 1000

        # PISCES next month
        ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :] * 1000
        ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :] * 1000
        ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :] * 1000
        ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :] * 1000
        ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
        ncbry['TALK_south'][2, :, :] = talk_bry3[:, 0, :] * 1000
        ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :] * 1000

    elif obctype == 'n':

        # PHYSICS
        ncbry['zeta_north'][0, :] = zeta_bry1[-1, :]
        ncbry['temp_north'][0, :, :] = temp_bry1[:, -1, :]
        ncbry['salt_north'][0, :, :] = salt_bry1[:, -1, :]
        ncbry['u_north'][0, :, :] = u_bry1[:, -1, :]
        ncbry['v_north'][0, :, :] = v_bry1[:, -1, :]
        ncbry['ubar_north'][0, :] = ubar_bry1[-1, :]
        ncbry['vbar_north'][0, :] = vbar_bry1[-1, :]

        ncbry['zeta_north'][1, :] = zeta_bry2[-1, :]
        ncbry['temp_north'][1, :, :] = temp_bry2[:, -1, :]
        ncbry['salt_north'][1, :, :] = salt_bry2[:, -1, :]
        ncbry['u_north'][1, :, :] = u_bry2[:, -1, :]
        ncbry['v_north'][1, :, :] = v_bry2[:, -1, :]
        ncbry['ubar_north'][1, :] = ubar_bry2[-1, :]
        ncbry['vbar_north'][1, :] = vbar_bry2[-1, :]

        ncbry['zeta_north'][2, :] = zeta_bry3[-1, :]
        ncbry['temp_north'][2, :, :] = temp_bry3[:, -1, :]
        ncbry['salt_north'][2, :, :] = salt_bry3[:, -1, :]
        ncbry['u_north'][2, :, :] = u_bry3[:, -1, :]
        ncbry['v_north'][2, :, :] = v_bry3[:, -1, :]
        ncbry['ubar_north'][2, :] = ubar_bry3[-1, :]
        ncbry['vbar_north'][2, :] = vbar_bry3[-1, :]

        # PISCES previous month
        ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :] * 1000
        ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :] * 1000
        ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :] * 1000
        ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :] * 1000
        ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
        ncbry['TALK_north'][0, :, :] = talk_bry1[:, -1, :] * 1000
        ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :] * 1000

        # PISCES current month
        ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :] * 1000
        ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :] * 1000
        ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :] * 1000
        ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :] * 1000
        ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
        ncbry['TALK_north'][1, :, :] = talk_bry2[:, -1, :] * 1000
        ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :] * 1000

        # PISCES next month
        ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :] * 1000
        ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :] * 1000
        ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :] * 1000
        ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :] * 1000
        ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
        ncbry['TALK_north'][2, :, :] = talk_bry3[:, -1, :] * 1000
        ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :] * 1000

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][0, :] = zeta_bry1[:, 0]
        ncbry['temp_east'][0, :, :] = temp_bry1[:, :, 0]
        ncbry['salt_east'][0, :, :] = salt_bry1[:, :, 0]
        ncbry['u_east'][0, :, :] = u_bry1[:, :, 0]
        ncbry['v_east'][0, :, :] = v_bry1[:, :, 0]
        ncbry['ubar_east'][0, :] = ubar_bry1[:, 0]
        ncbry['vbar_east'][0, :] = vbar_bry1[:, 0]

        ncbry['zeta_east'][1, :] = zeta_bry2[:, 0]
        ncbry['temp_east'][1, :, :] = temp_bry2[:, :, 0]
        ncbry['salt_east'][1, :, :] = salt_bry2[:, :, 0]
        ncbry['u_east'][1, :, :] = u_bry2[:, :, 0]
        ncbry['v_east'][1, :, :] = v_bry2[:, :, 0]
        ncbry['ubar_east'][1, :] = ubar_bry2[:, 0]
        ncbry['vbar_east'][1, :] = vbar_bry2[:, 0]

        ncbry['zeta_east'][2, :] = zeta_bry3[:, 0]
        ncbry['temp_east'][2, :, :] = temp_bry3[:, :, 0]
        ncbry['salt_east'][2, :, :] = salt_bry3[:, :, 0]
        ncbry['u_east'][2, :, :] = u_bry3[:, :, 0]
        ncbry['v_east'][2, :, :] = v_bry3[:, :, 0]
        ncbry['ubar_east'][2, :] = ubar_bry3[:, 0]
        ncbry['vbar_east'][2, :] = vbar_bry3[:, 0]

        # PISCES previous month
        ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0] * 1000
        ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0] * 1000
        ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0] * 1000
        ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0] * 1000
        ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
        ncbry['TALK_east'][0, :, :] = talk_bry1[:, :, 0] * 1000
        ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0] * 1000

        # PISCES current month
        ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0] * 1000
        ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0] * 1000
        ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0] * 1000
        ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0] * 1000
        ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
        ncbry['TALK_east'][1, :, :] = talk_bry2[:, :, 0] * 1000
        ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0] * 1000

        # PISCES next month
        ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0] * 1000
        ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0] * 1000
        ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0] * 1000
        ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0] * 1000
        ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
        ncbry['TALK_east'][2, :, :] = talk_bry3[:, :, 0] * 1000
        ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0] * 1000

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][0, :] = zeta_bry1[:, -1]
        ncbry['temp_west'][0, :, :] = temp_bry1[:, :, -1]
        ncbry['salt_west'][0, :, :] = salt_bry1[:, :, -1]
        ncbry['u_west'][0, :, :] = u_bry1[:, :, -1]
        ncbry['v_west'][0, :, :] = v_bry1[:, :, -1]
        ncbry['ubar_west'][0, :] = ubar_bry1[:, -1]
        ncbry['vbar_west'][0, :] = vbar_bry1[:, -1]

        ncbry['zeta_west'][1, :] = zeta_bry2[:, -1]
        ncbry['temp_west'][1, :, :] = temp_bry2[:, :, -1]
        ncbry['salt_west'][1, :, :] = salt_bry2[:, :, -1]
        ncbry['u_west'][1, :, :] = u_bry2[:, :, -1]
        ncbry['v_west'][1, :, :] = v_bry2[:, :, -1]
        ncbry['ubar_west'][1, :] = ubar_bry2[:, -1]
        ncbry['vbar_west'][1, :] = vbar_bry2[:, -1]

        ncbry['zeta_west'][2, :] = zeta_bry3[:, -1]
        ncbry['temp_west'][2, :, :] = temp_bry3[:, :, -1]
        ncbry['salt_west'][2, :, :] = salt_bry3[:, :, -1]
        ncbry['u_west'][2, :, :] = u_bry3[:, :, -1]
        ncbry['v_west'][2, :, :] = v_bry3[:, :, -1]
        ncbry['ubar_west'][2, :] = ubar_bry3[:, -1]
        ncbry['vbar_west'][2, :] = vbar_bry3[:, -1]

        # PISCES previous month
        ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1] * 1000
        ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1] * 1000
        ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1] * 1000
        ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1] * 1000
        ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
        ncbry['TALK_west'][0, :, :] = talk_bry1[:, :, -1] * 1000
        ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1] * 1000

        # PISCES current month
        ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1] * 1000
        ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1] * 1000
        ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1] * 1000
        ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1] * 1000
        ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
        ncbry['TALK_west'][1, :, :] = talk_bry2[:, :, -1] * 1000
        ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1] * 1000

        # PISCES next month
        ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1] * 1000
        ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1] * 1000
        ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1] * 1000
        ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1] * 1000
        ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
        ncbry['TALK_west'][2, :, :] = talk_bry3[:, :, -1] * 1000
        ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1] * 1000

    return ncbry


def interp_bry_GLO_BGC_PHY(obctype,
                           PNI,
                           ncglo,
                           ncpism1o, ncpiso, ncpisp1o,
                           tndx_glo, ncbry, tndx_bry,
                           h_bry, theta_s, theta_b,
                           hc, N, vtransform, Nzgoodmin, depthg, depthp, angle_bry,
                           LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                           LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                           LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                           LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(ncpism1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry2 = glor.interp3d(ncpiso, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        no3_bry3 = glor.interp3d(ncpisp1o, 'no3', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Orthophosphate
        #
        #

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(ncpism1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry2 = glor.interp3d(ncpiso, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        po4_bry3 = glor.interp3d(ncpisp1o, 'po4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Silicate
        #
        #

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(ncpism1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry2 = glor.interp3d(ncpiso, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        si_bry3 = glor.interp3d(ncpisp1o, 'si', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(ncpism1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry2 = glor.interp3d(ncpiso, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        o2_bry3 = glor.interp3d(ncpisp1o, 'o2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #
        #
        # print('Interpolate Dissolved Inorganic Carbon...')
        #
        # dic_bry1 = glor.interp3d(ncpism1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                          iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                          LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
        #
        # dic_bry2 = glor.interp3d(ncpiso, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                          iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                          LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
        #
        # dic_bry3 = glor.interp3d(ncpisp1o, 'dissic', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                          iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                          LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        #
        # 3: pH (to get TALK later via CO2SYS)
        #
        #
        #
        # print('Interpolate Total Alkalinity...')
        #
        # ph_bry1 = glor.interp3d(ncpism1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                         iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                         LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
        #
        # ph_bry2 = glor.interp3d(ncpiso, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                         iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                         LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
        #
        # ph_bry3 = glor.interp3d(ncpisp1o, 'ph', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                         iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                         LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # PISCES
        #

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # 5r: Ammonium from IBI PISCES
        #
        #
        # print('Interpolate Ammonium...')
        #
        # nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                          iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                          LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
        #
        # nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                          iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                          LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')
        #
        # nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
        #                          iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
        #                          LonP_bry, LatP_bry, coefP_bry, elemP_bry, '_FillValue')

        #
        # Monthly Salinity and Temperature from IBI for calculation of TALK from DIC, pH, temp and salinity
        #
        #
        # salm_bry1 = glor.interp3d(ncphymm1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
        #                           iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
        #                           LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
        #
        # salm_bry2 = glor.interp3d(ncphymo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
        #                           iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
        #                           LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
        #
        # salm_bry3 = glor.interp3d(ncphymp1o, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
        #                           iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
        #                           LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
        #
        # tmpm_bry1 = glor.interp3d(ncphymm1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
        #                           iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
        #                           LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
        #
        # tmpm_bry2 = glor.interp3d(ncphymo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
        #                           iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
        #                           LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')
        #
        # tmpm_bry3 = glor.interp3d(ncphymp1o, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
        #                           iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
        #                           LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

        if obctype == 's':
            # PISCES previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :]
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :]
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :]
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :]
            ncbry['DIC_south'][0, :, :] = 2150. * np.ones_like(o2_bry1[:, 0, :])
            ncbry['TALK_south'][0, :, :] = 2350. * np.ones_like(o2_bry1[:, 0, :])
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]

            # PISCES current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :]
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :]
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :]
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :]
            ncbry['DIC_south'][1, :, :] = 2150. * np.ones_like(o2_bry2[:, 0, :])
            ncbry['TALK_south'][1, :, :] = 2350. * np.ones_like(o2_bry2[:, 0, :])
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]

            # PISCES next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :]
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :]
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :]
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :]
            ncbry['DIC_south'][2, :, :] = 2150. * np.ones_like(o2_bry3[:, 0, :])
            ncbry['TALK_south'][2, :, :] = 2350. * np.ones_like(o2_bry3[:, 0, :])
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]

        elif obctype == 'n':

            # PISCES previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :]
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :]
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :]
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :]
            ncbry['DIC_north'][0, :, :] = 2150. * np.ones_like(o2_bry1[:, -1, :])
            ncbry['TALK_north'][0, :, :] = 2350. * np.ones_like(o2_bry1[:, -1, :])
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]

            # PISCES current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :]
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :]
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :]
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :]
            ncbry['DIC_north'][1, :, :] = 2150. * np.ones_like(o2_bry2[:, -1, :])
            ncbry['TALK_north'][1, :, :] = 2350. * np.ones_like(o2_bry2[:, -1, :])
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]

            # PISCES next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :]
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :]
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :]
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :]
            ncbry['DIC_north'][2, :, :] = 2150. * np.ones_like(o2_bry3[:, -1, :])
            ncbry['TALK_north'][2, :, :] = 2350. * np.ones_like(o2_bry3[:, -1, :])
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]

        elif obctype == 'e':
            # PISCES previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0]
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0]
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0]
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0]
            ncbry['DIC_east'][0, :, :] = 2150. * np.ones_like(o2_bry1[:, :, 0])
            ncbry['TALK_east'][0, :, :] = 2350. * np.ones_like(o2_bry1[:, :, 0])
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]

            # PISCES current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0]
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0]
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0]
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0]
            ncbry['DIC_east'][1, :, :] = 2150. * np.ones_like(o2_bry2[:, :, 0])
            ncbry['TALK_east'][1, :, :] = 2350. * np.ones_like(o2_bry2[:, :, 0])
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]

            # PISCES next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0]
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0]
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0]
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0]
            ncbry['DIC_east'][2, :, :] = 2150. * np.ones_like(o2_bry3[:, :, 0])
            ncbry['TALK_east'][2, :, :] = 2350. * np.ones_like(o2_bry3[:, :, 0])
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]

        elif obctype == 'w':
            # PISCES previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1]
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1]
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1]
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1]
            ncbry['DIC_west'][0, :, :] = 2150. * np.ones_like(o2_bry1[:, :, -1])
            ncbry['TALK_west'][0, :, :] = 2350. * np.ones_like(o2_bry1[:, :, -1])
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]

            # PISCES current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1]
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1]
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1]
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1]
            ncbry['DIC_west'][1, :, :] = 2150. * np.ones_like(o2_bry2[:, :, -1])
            ncbry['TALK_west'][1, :, :] = 2350. * np.ones_like(o2_bry2[:, :, -1])
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]

            # PISCES next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1]
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1]
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1]
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1]
            ncbry['DIC_west'][2, :, :] = 2150. * np.ones_like(o2_bry3[:, :, -1])
            ncbry['TALK_west'][2, :, :] = 2350. * np.ones_like(o2_bry3[:, :, -1])
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry


def interp_bry_bc_PISCES_NORESM(obctype,
                                PNI, Tidxn,
                                ncglo,
                                ncpism1o, ncpiso, ncpisp1o,
                                NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
                                tndx_glo, ncbry, tndx_bry,
                                h_bry, theta_s, theta_b,
                                hc, N, vtransform, Nzgoodmin, depthg, depthp, depthn, angle_bry,
                                LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                                LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                                LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                                LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
                                LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    # interp_bry_PISCES_NORESM('s',
    #
    #                          PISCES_NORESM_interpd, Tininxp,
    #
    #                          ncglo,
    #
    #                          ncpisprio, ncpiscuro, ncpisposo,
    #                          NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
    #
    #                          tndx_glo, ncbry, tndx_bry,
    #                          h_bry, theta_s, theta_b,
    #                          hc, N, vtransform, Nzgoodmin, depth, angle_bry,
    #                          LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
    #                          LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
    #                          LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
    #                          LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
    #                          LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)

    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao_bc', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, 'add_offset')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so_bc', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, 'add_offset')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
        ncnornio = ncnorni.variables

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(ncnornio, 'no3no2', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        no3_bry2 = glor.interp3d(ncnornio, 'no3no2', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        no3_bry3 = glor.interp3d(ncnornio, 'no3no2', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnorni.close()

        #
        #
        # 3: Orthophosphate
        #
        #

        ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
        ncnorpoo = ncnorpo.variables

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(ncnorpoo, 'po4', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        po4_bry2 = glor.interp3d(ncnorpoo, 'po4', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        po4_bry3 = glor.interp3d(ncnorpoo, 'po4', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnorpo.close()

        #
        #
        # 3: Silicate
        #
        #

        ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
        ncnorsio = ncnorsi.variables

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(ncnorsio, 'si', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        si_bry2 = glor.interp3d(ncnorsio, 'si', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        si_bry3 = glor.interp3d(ncnorsio, 'si', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnorsi.close()

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        ncnordo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[3] + NORESM_ending, 'r')
        ncnordoo = ncnordo.variables

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(ncnordoo, 'o2', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        o2_bry2 = glor.interp3d(ncnordoo, 'o2', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        o2_bry3 = glor.interp3d(ncnordoo, 'o2', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnordo.close()

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        ncnordic = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[4] + NORESM_ending, 'r')
        ncnordico = ncnordic.variables

        print('Interpolate Dissolved Inorganic Carbon...')

        dic_bry1 = glor.interp3d(ncnordico, 'dissic', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        dic_bry2 = glor.interp3d(ncnordico, 'dissic', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        dic_bry3 = glor.interp3d(ncnordico, 'dissic', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnordic.close()

        #
        #
        # 3: Total Alkalinity
        #
        #

        ncnoralkalini = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[5] + NORESM_ending, 'r')
        ncnoralkalinio = ncnoralkalini.variables

        print('Interpolate Total Alkalinity...')

        talk_bry1 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                  iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                  LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        talk_bry2 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                  iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                  LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        talk_bry3 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                  iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                  LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnoralkalini.close()

        # PISCES
        #
        #
        # 3: Ammonium
        #
        #

        print('Interpolate Calcite...')

        calc_bry1 = glor.interp3d(ncpism1o, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        calc_bry2 = glor.interp3d(ncpiso, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        calc_bry3 = glor.interp3d(ncpisp1o, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5b: Particulate Organic Carbon from IBI PISCES
        #

        print('Interpolate Particulate Organic Carbon...')

        poc_bry1 = glor.interp3d(ncpism1o, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        poc_bry2 = glor.interp3d(ncpiso, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        poc_bry3 = glor.interp3d(ncpisp1o, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5c: Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Nanophytoplankton...')

        phy_bry1 = glor.interp3d(ncpism1o, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy_bry2 = glor.interp3d(ncpiso, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy_bry3 = glor.interp3d(ncpism1o, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5d: Microzooplankton from IBI PISCES
        #

        print('Interpolate Microzooplankton...')

        zoo_bry1 = glor.interp3d(ncpism1o, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo_bry2 = glor.interp3d(ncpiso, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo_bry3 = glor.interp3d(ncpisp1o, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5e: Dissolved Organic Carbon from IBI PISCES
        #

        print('Interpolate Dissolved Organic Carbon...')

        doc_bry1 = glor.interp3d(ncpism1o, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        doc_bry2 = glor.interp3d(ncpiso, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        doc_bry3 = glor.interp3d(ncpisp1o, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5f: Diatom from IBI PISCES
        #

        print('Interpolate Diatom...')

        phy2_bry1 = glor.interp3d(ncpism1o, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy2_bry2 = glor.interp3d(ncpiso, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy2_bry3 = glor.interp3d(ncpisp1o, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5g: Mesozooplankton from IBI PISCES
        #

        print('Interpolate Mesozooplankton...')

        zoo2_bry1 = glor.interp3d(ncpism1o, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo2_bry2 = glor.interp3d(ncpiso, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo2_bry3 = glor.interp3d(ncpisp1o, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5h: Biogenic Silica from IBI PISCES
        #

        print('Interpolate Biogenic Silica...')

        gsi_bry1 = glor.interp3d(ncpism1o, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        gsi_bry2 = glor.interp3d(ncpiso, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        gsi_bry3 = glor.interp3d(ncpisp1o, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5j: Big Particle Iron from IBI PISCES
        #

        print('Interpolate Big Particle Iron...')

        bfe_bry1 = glor.interp3d(ncpism1o, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        bfe_bry2 = glor.interp3d(ncpiso, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        bfe_bry3 = glor.interp3d(ncpisp1o, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5k: Big Particulate Organic Carbon from IBI PISCES
        #

        print('Interpolate Big Particulate Organic Carbon...')

        goc_bry1 = glor.interp3d(ncpism1o, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        goc_bry2 = glor.interp3d(ncpiso, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        goc_bry3 = glor.interp3d(ncpisp1o, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5l: Iron in the small particles from IBI PISCES
        #

        print('Interpolate Iron in the small particles...')

        sfe_bry1 = glor.interp3d(ncpism1o, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        sfe_bry2 = glor.interp3d(ncpiso, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        sfe_bry3 = glor.interp3d(ncpisp1o, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5m: Iron content of the diatoms from IBI PISCES
        #

        print('Interpolate Iron content of the diatoms...')

        dfe_bry1 = glor.interp3d(ncpism1o, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dfe_bry2 = glor.interp3d(ncpiso, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dfe_bry3 = glor.interp3d(ncpisp1o, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5n: Silicon content of the Diatoms from IBI PISCES
        #

        print('Interpolate Silicon content of the Diatoms...')

        dsi_bry1 = glor.interp3d(ncpism1o, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dsi_bry2 = glor.interp3d(ncpiso, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dsi_bry3 = glor.interp3d(ncpisp1o, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5o: Iron content of the Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Iron content of the Nanophytoplankton...')

        nfe_bry1 = glor.interp3d(ncpism1o, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nfe_bry2 = glor.interp3d(ncpiso, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nfe_bry3 = glor.interp3d(ncpisp1o, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5p: Nanophytoplankton Chlorophyll from IBI PISCES
        #

        print('Interpolate Nanophytoplankton Chlorophyll...')

        nchl_bry1 = glor.interp3d(ncpism1o, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nchl_bry2 = glor.interp3d(ncpiso, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nchl_bry3 = glor.interp3d(ncpisp1o, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5q: Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Diatom Chlorophyll...')

        dchl_bry1 = glor.interp3d(ncpism1o, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dchl_bry2 = glor.interp3d(ncpiso, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dchl_bry3 = glor.interp3d(ncpisp1o, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        if obctype == 's':
            # NORESM previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :] * 1000
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :] * 1000
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :] * 1000
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :] * 1000
            ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
            ncbry['TALK_south'][0, :, :] = talk_bry1[:, 0, :] * 1000

            # PISCES previous month
            ncbry['DCHL_south'][0, :, :] = dchl_bry1[:, 0, :]
            ncbry['NCHL_south'][0, :, :] = nchl_bry1[:, 0, :]
            ncbry['NFE_south'][0, :, :] = nfe_bry1[:, 0, :]
            ncbry['DSI_south'][0, :, :] = dsi_bry1[:, 0, :]
            ncbry['DFE_south'][0, :, :] = dfe_bry1[:, 0, :]
            ncbry['SFE_south'][0, :, :] = sfe_bry1[:, 0, :]
            ncbry['GOC_south'][0, :, :] = goc_bry1[:, 0, :]
            ncbry['BFE_south'][0, :, :] = bfe_bry1[:, 0, :]
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['BSI_south'][0, :, :] = gsi_bry1[:, 0, :]
            ncbry['MESO_south'][0, :, :] = zoo2_bry1[:, 0, :]
            ncbry['DIA_south'][0, :, :] = phy2_bry1[:, 0, :]
            ncbry['DOC_south'][0, :, :] = doc_bry1[:, 0, :]
            ncbry['ZOO_south'][0, :, :] = zoo_bry1[:, 0, :]
            ncbry['NANO_south'][0, :, :] = phy_bry1[:, 0, :]
            ncbry['CACO3_south'][0, :, :] = calc_bry1[:, 0, :]
            ncbry['POC_south'][0, :, :] = poc_bry1[:, 0, :]
            ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]

            # NORESM current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :] * 1000
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :] * 1000
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :] * 1000
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :] * 1000
            ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
            ncbry['TALK_south'][1, :, :] = talk_bry2[:, 0, :] * 1000

            # PISCES current month
            ncbry['DCHL_south'][1, :, :] = dchl_bry2[:, 0, :]
            ncbry['NCHL_south'][1, :, :] = nchl_bry2[:, 0, :]
            ncbry['NFE_south'][1, :, :] = nfe_bry2[:, 0, :]
            ncbry['DSI_south'][1, :, :] = dsi_bry2[:, 0, :]
            ncbry['DFE_south'][1, :, :] = dfe_bry2[:, 0, :]
            ncbry['SFE_south'][1, :, :] = sfe_bry2[:, 0, :]
            ncbry['GOC_south'][1, :, :] = goc_bry2[:, 0, :]
            ncbry['BFE_south'][1, :, :] = bfe_bry2[:, 0, :]
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['BSI_south'][1, :, :] = gsi_bry2[:, 0, :]
            ncbry['MESO_south'][1, :, :] = zoo2_bry2[:, 0, :]
            ncbry['DIA_south'][1, :, :] = phy2_bry2[:, 0, :]
            ncbry['DOC_south'][1, :, :] = doc_bry2[:, 0, :]
            ncbry['ZOO_south'][1, :, :] = zoo_bry2[:, 0, :]
            ncbry['NANO_south'][1, :, :] = phy_bry2[:, 0, :]
            ncbry['CACO3_south'][1, :, :] = calc_bry2[:, 0, :]
            ncbry['POC_south'][1, :, :] = poc_bry2[:, 0, :]
            ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]

            # NORESM next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :] * 1000
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :] * 1000
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :] * 1000
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :] * 1000
            ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
            ncbry['TALK_south'][2, :, :] = talk_bry3[:, 0, :] * 1000

            # PISCES next month
            ncbry['DCHL_south'][2, :, :] = dchl_bry3[:, 0, :]
            ncbry['NCHL_south'][2, :, :] = nchl_bry3[:, 0, :]
            ncbry['NFE_south'][2, :, :] = nfe_bry3[:, 0, :]
            ncbry['DSI_south'][2, :, :] = dsi_bry3[:, 0, :]
            ncbry['DFE_south'][2, :, :] = dfe_bry3[:, 0, :]
            ncbry['SFE_south'][2, :, :] = sfe_bry3[:, 0, :]
            ncbry['GOC_south'][2, :, :] = goc_bry3[:, 0, :]
            ncbry['BFE_south'][2, :, :] = bfe_bry3[:, 0, :]
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['BSI_south'][2, :, :] = gsi_bry3[:, 0, :]
            ncbry['MESO_south'][2, :, :] = zoo2_bry3[:, 0, :]
            ncbry['DIA_south'][2, :, :] = phy2_bry3[:, 0, :]
            ncbry['DOC_south'][2, :, :] = doc_bry3[:, 0, :]
            ncbry['ZOO_south'][2, :, :] = zoo_bry3[:, 0, :]
            ncbry['NANO_south'][2, :, :] = phy_bry3[:, 0, :]
            ncbry['CACO3_south'][2, :, :] = calc_bry3[:, 0, :]
            ncbry['POC_south'][2, :, :] = poc_bry3[:, 0, :]
            ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]

        elif obctype == 'n':
            # NORESM previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :] * 1000
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :] * 1000
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :] * 1000
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :] * 1000
            ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
            ncbry['TALK_north'][0, :, :] = talk_bry1[:, -1, :] * 1000

            # PISCES previous month
            ncbry['DCHL_north'][0, :, :] = dchl_bry1[:, -1, :]
            ncbry['NCHL_north'][0, :, :] = nchl_bry1[:, -1, :]
            ncbry['NFE_north'][0, :, :] = nfe_bry1[:, -1, :]
            ncbry['DSI_north'][0, :, :] = dsi_bry1[:, -1, :]
            ncbry['DFE_north'][0, :, :] = dfe_bry1[:, -1, :]
            ncbry['SFE_north'][0, :, :] = sfe_bry1[:, -1, :]
            ncbry['GOC_north'][0, :, :] = goc_bry1[:, -1, :]
            ncbry['BFE_north'][0, :, :] = bfe_bry1[:, -1, :]
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['BSI_north'][0, :, :] = gsi_bry1[:, -1, :]
            ncbry['MESO_north'][0, :, :] = zoo2_bry1[:, -1, :]
            ncbry['DIA_north'][0, :, :] = phy2_bry1[:, -1, :]
            ncbry['DOC_north'][0, :, :] = doc_bry1[:, -1, :]
            ncbry['ZOO_north'][0, :, :] = zoo_bry1[:, -1, :]
            ncbry['NANO_north'][0, :, :] = phy_bry1[:, -1, :]
            ncbry['CACO3_north'][0, :, :] = calc_bry1[:, -1, :]
            ncbry['POC_north'][0, :, :] = poc_bry1[:, -1, :]
            ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]

            # NORESM current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :] * 1000
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :] * 1000
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :] * 1000
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :] * 1000
            ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
            ncbry['TALK_north'][1, :, :] = talk_bry2[:, -1, :] * 1000

            # PISCES current month
            ncbry['DCHL_north'][1, :, :] = dchl_bry2[:, -1, :]
            ncbry['NCHL_north'][1, :, :] = nchl_bry2[:, -1, :]
            ncbry['NFE_north'][1, :, :] = nfe_bry2[:, -1, :]
            ncbry['DSI_north'][1, :, :] = dsi_bry2[:, -1, :]
            ncbry['DFE_north'][1, :, :] = dfe_bry2[:, -1, :]
            ncbry['SFE_north'][1, :, :] = sfe_bry2[:, -1, :]
            ncbry['GOC_north'][1, :, :] = goc_bry2[:, -1, :]
            ncbry['BFE_north'][1, :, :] = bfe_bry2[:, -1, :]
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['BSI_north'][1, :, :] = gsi_bry2[:, -1, :]
            ncbry['MESO_north'][1, :, :] = zoo2_bry2[:, -1, :]
            ncbry['DIA_north'][1, :, :] = phy2_bry2[:, -1, :]
            ncbry['DOC_north'][1, :, :] = doc_bry2[:, -1, :]
            ncbry['ZOO_north'][1, :, :] = zoo_bry2[:, -1, :]
            ncbry['NANO_north'][1, :, :] = phy_bry2[:, -1, :]
            ncbry['CACO3_north'][1, :, :] = calc_bry2[:, -1, :]
            ncbry['POC_north'][1, :, :] = poc_bry2[:, -1, :]
            ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]

            # NORESM next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :] * 1000
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :] * 1000
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :] * 1000
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :] * 1000
            ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
            ncbry['TALK_north'][2, :, :] = talk_bry3[:, -1, :] * 1000

            # PISCES next month
            ncbry['DCHL_north'][2, :, :] = dchl_bry3[:, -1, :]
            ncbry['NCHL_north'][2, :, :] = nchl_bry3[:, -1, :]
            ncbry['NFE_north'][2, :, :] = nfe_bry3[:, -1, :]
            ncbry['DSI_north'][2, :, :] = dsi_bry3[:, -1, :]
            ncbry['DFE_north'][2, :, :] = dfe_bry3[:, -1, :]
            ncbry['SFE_north'][2, :, :] = sfe_bry3[:, -1, :]
            ncbry['GOC_north'][2, :, :] = goc_bry3[:, -1, :]
            ncbry['BFE_north'][2, :, :] = bfe_bry3[:, -1, :]
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['BSI_north'][2, :, :] = gsi_bry3[:, -1, :]
            ncbry['MESO_north'][2, :, :] = zoo2_bry3[:, -1, :]
            ncbry['DIA_north'][2, :, :] = phy2_bry3[:, -1, :]
            ncbry['DOC_north'][2, :, :] = doc_bry3[:, -1, :]
            ncbry['ZOO_north'][2, :, :] = zoo_bry3[:, -1, :]
            ncbry['NANO_north'][2, :, :] = phy_bry3[:, -1, :]
            ncbry['CACO3_north'][2, :, :] = calc_bry3[:, -1, :]
            ncbry['POC_north'][2, :, :] = poc_bry3[:, -1, :]
            ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]

        elif obctype == 'e':
            # NORESM previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0] * 1000
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0] * 1000
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0] * 1000
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0] * 1000
            ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
            ncbry['TALK_east'][0, :, :] = talk_bry1[:, :, 0] * 1000

            # PISCES previous month
            ncbry['DCHL_east'][0, :, :] = dchl_bry1[:, :, 0]
            ncbry['NCHL_east'][0, :, :] = nchl_bry1[:, :, 0]
            ncbry['NFE_east'][0, :, :] = nfe_bry1[:, :, 0]
            ncbry['DSI_east'][0, :, :] = dsi_bry1[:, :, 0]
            ncbry['DFE_east'][0, :, :] = dfe_bry1[:, :, 0]
            ncbry['SFE_east'][0, :, :] = sfe_bry1[:, :, 0]
            ncbry['GOC_east'][0, :, :] = goc_bry1[:, :, 0]
            ncbry['BFE_east'][0, :, :] = bfe_bry1[:, :, 0]
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['BSI_east'][0, :, :] = gsi_bry1[:, :, 0]
            ncbry['MESO_east'][0, :, :] = zoo2_bry1[:, :, 0]
            ncbry['DIA_east'][0, :, :] = phy2_bry1[:, :, 0]
            ncbry['DOC_east'][0, :, :] = doc_bry1[:, :, 0]
            ncbry['ZOO_east'][0, :, :] = zoo_bry1[:, :, 0]
            ncbry['NANO_east'][0, :, :] = phy_bry1[:, :, 0]
            ncbry['CACO3_east'][0, :, :] = calc_bry1[:, :, 0]
            ncbry['POC_east'][0, :, :] = poc_bry1[:, :, 0]
            ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]

            # NORESM current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0] * 1000
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0] * 1000
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0] * 1000
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0] * 1000
            ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
            ncbry['TALK_east'][1, :, :] = talk_bry2[:, :, 0] * 1000

            # PISCES current month
            ncbry['DCHL_east'][1, :, :] = dchl_bry2[:, :, 0]
            ncbry['NCHL_east'][1, :, :] = nchl_bry2[:, :, 0]
            ncbry['NFE_east'][1, :, :] = nfe_bry2[:, :, 0]
            ncbry['DSI_east'][1, :, :] = dsi_bry2[:, :, 0]
            ncbry['DFE_east'][1, :, :] = dfe_bry2[:, :, 0]
            ncbry['SFE_east'][1, :, :] = sfe_bry2[:, :, 0]
            ncbry['GOC_east'][1, :, :] = goc_bry2[:, :, 0]
            ncbry['BFE_east'][1, :, :] = bfe_bry2[:, :, 0]
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['BSI_east'][1, :, :] = gsi_bry2[:, :, 0]
            ncbry['MESO_east'][1, :, :] = zoo2_bry2[:, :, 0]
            ncbry['DIA_east'][1, :, :] = phy2_bry2[:, :, 0]
            ncbry['DOC_east'][1, :, :] = doc_bry2[:, :, 0]
            ncbry['ZOO_east'][1, :, :] = zoo_bry2[:, :, 0]
            ncbry['NANO_east'][1, :, :] = phy_bry2[:, :, 0]
            ncbry['CACO3_east'][1, :, :] = calc_bry2[:, :, 0]
            ncbry['POC_east'][1, :, :] = poc_bry2[:, :, 0]
            ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]

            # NORESM next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0] * 1000
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0] * 1000
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0] * 1000
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0] * 1000
            ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
            ncbry['TALK_east'][2, :, :] = talk_bry3[:, :, 0] * 1000

            # PISCES next month
            ncbry['DCHL_east'][2, :, :] = dchl_bry3[:, :, 0]
            ncbry['NCHL_east'][2, :, :] = nchl_bry3[:, :, 0]
            ncbry['NFE_east'][2, :, :] = nfe_bry3[:, :, 0]
            ncbry['DSI_east'][2, :, :] = dsi_bry3[:, :, 0]
            ncbry['DFE_east'][2, :, :] = dfe_bry3[:, :, 0]
            ncbry['SFE_east'][2, :, :] = sfe_bry3[:, :, 0]
            ncbry['GOC_east'][2, :, :] = goc_bry3[:, :, 0]
            ncbry['BFE_east'][2, :, :] = bfe_bry3[:, :, 0]
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['BSI_east'][2, :, :] = gsi_bry3[:, :, 0]
            ncbry['MESO_east'][2, :, :] = zoo2_bry3[:, :, 0]
            ncbry['DIA_east'][2, :, :] = phy2_bry3[:, :, 0]
            ncbry['DOC_east'][2, :, :] = doc_bry3[:, :, 0]
            ncbry['ZOO_east'][2, :, :] = zoo_bry3[:, :, 0]
            ncbry['NANO_east'][2, :, :] = phy_bry3[:, :, 0]
            ncbry['CACO3_east'][2, :, :] = calc_bry3[:, :, 0]
            ncbry['POC_east'][2, :, :] = poc_bry3[:, :, 0]
            ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]

        elif obctype == 'w':
            # NORESM previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1] * 1000
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1] * 1000
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1] * 1000
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1] * 1000
            ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
            ncbry['TALK_west'][0, :, :] = talk_bry1[:, :, -1] * 1000

            # PISCES previous month
            ncbry['DCHL_west'][0, :, :] = dchl_bry1[:, :, -1]
            ncbry['NCHL_west'][0, :, :] = nchl_bry1[:, :, -1]
            ncbry['NFE_west'][0, :, :] = nfe_bry1[:, :, -1]
            ncbry['DSI_west'][0, :, :] = dsi_bry1[:, :, -1]
            ncbry['DFE_west'][0, :, :] = dfe_bry1[:, :, -1]
            ncbry['SFE_west'][0, :, :] = sfe_bry1[:, :, -1]
            ncbry['GOC_west'][0, :, :] = goc_bry1[:, :, -1]
            ncbry['BFE_west'][0, :, :] = bfe_bry1[:, :, -1]
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['BSI_west'][0, :, :] = gsi_bry1[:, :, -1]
            ncbry['MESO_west'][0, :, :] = zoo2_bry1[:, :, -1]
            ncbry['DIA_west'][0, :, :] = phy2_bry1[:, :, -1]
            ncbry['DOC_west'][0, :, :] = doc_bry1[:, :, -1]
            ncbry['ZOO_west'][0, :, :] = zoo_bry1[:, :, -1]
            ncbry['NANO_west'][0, :, :] = phy_bry1[:, :, -1]
            ncbry['CACO3_west'][0, :, :] = calc_bry1[:, :, -1]
            ncbry['POC_west'][0, :, :] = poc_bry1[:, :, -1]
            ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]

            # NORESM current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1] * 1000
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1] * 1000
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1] * 1000
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1] * 1000
            ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
            ncbry['TALK_west'][1, :, :] = talk_bry2[:, :, -1] * 1000

            # PISCES current month
            ncbry['DCHL_west'][1, :, :] = dchl_bry2[:, :, -1]
            ncbry['NCHL_west'][1, :, :] = nchl_bry2[:, :, -1]
            ncbry['NFE_west'][1, :, :] = nfe_bry2[:, :, -1]
            ncbry['DSI_west'][1, :, :] = dsi_bry2[:, :, -1]
            ncbry['DFE_west'][1, :, :] = dfe_bry2[:, :, -1]
            ncbry['SFE_west'][1, :, :] = sfe_bry2[:, :, -1]
            ncbry['GOC_west'][1, :, :] = goc_bry2[:, :, -1]
            ncbry['BFE_west'][1, :, :] = bfe_bry2[:, :, -1]
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['BSI_west'][1, :, :] = gsi_bry2[:, :, -1]
            ncbry['MESO_west'][1, :, :] = zoo2_bry2[:, :, -1]
            ncbry['DIA_west'][1, :, :] = phy2_bry2[:, :, -1]
            ncbry['DOC_west'][1, :, :] = doc_bry2[:, :, -1]
            ncbry['ZOO_west'][1, :, :] = zoo_bry2[:, :, -1]
            ncbry['NANO_west'][1, :, :] = phy_bry2[:, :, -1]
            ncbry['CACO3_west'][1, :, :] = calc_bry2[:, :, -1]
            ncbry['POC_west'][1, :, :] = poc_bry2[:, :, -1]
            ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]

            # NORESM next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1] * 1000
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1] * 1000
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1] * 1000
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1] * 1000
            ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
            ncbry['TALK_west'][2, :, :] = talk_bry3[:, :, -1] * 1000

            # PISCES next month
            ncbry['DCHL_west'][2, :, :] = dchl_bry3[:, :, -1]
            ncbry['NCHL_west'][2, :, :] = nchl_bry3[:, :, -1]
            ncbry['NFE_west'][2, :, :] = nfe_bry3[:, :, -1]
            ncbry['DSI_west'][2, :, :] = dsi_bry3[:, :, -1]
            ncbry['DFE_west'][2, :, :] = dfe_bry3[:, :, -1]
            ncbry['SFE_west'][2, :, :] = sfe_bry3[:, :, -1]
            ncbry['GOC_west'][2, :, :] = goc_bry3[:, :, -1]
            ncbry['BFE_west'][2, :, :] = bfe_bry3[:, :, -1]
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['BSI_west'][2, :, :] = gsi_bry3[:, :, -1]
            ncbry['MESO_west'][2, :, :] = zoo2_bry3[:, :, -1]
            ncbry['DIA_west'][2, :, :] = phy2_bry3[:, :, -1]
            ncbry['DOC_west'][2, :, :] = doc_bry3[:, :, -1]
            ncbry['ZOO_west'][2, :, :] = zoo_bry3[:, :, -1]
            ncbry['NANO_west'][2, :, :] = phy_bry3[:, :, -1]
            ncbry['CACO3_west'][2, :, :] = calc_bry3[:, :, -1]
            ncbry['POC_west'][2, :, :] = poc_bry3[:, :, -1]
            ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry


def interp_bry_GLO_PISCES_NORESM(obctype,
                                 PNI, Tidxn,
                                 ncglo,
                                 ncpism1o, ncpiso, ncpisp1o,
                                 NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
                                 tndx_glo, ncbry, tndx_bry,
                                 h_bry, theta_s, theta_b,
                                 hc, N, vtransform, Nzgoodmin, depthg, depthp, depthn, angle_bry,
                                 LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
                                 LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
                                 LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
                                 LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
                                 LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry):
    # interp_bry_PISCES_NORESM('s',
    #
    #                          PISCES_NORESM_interpd, Tininxp,
    #
    #                          ncglo,
    #
    #                          ncpisprio, ncpiscuro, ncpisposo,
    #                          NORESMfiles_dir, NORESM_prefix, NIVAvars, NORESM_ending,
    #
    #                          tndx_glo, ncbry, tndx_bry,
    #                          h_bry, theta_s, theta_b,
    #                          hc, N, vtransform, Nzgoodmin, depth, angle_bry,
    #                          LonT_bry, LatT_bry, iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry, elemT_bry, coefT_bry,
    #                          LonU_bry, LatU_bry, iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry, elemU_bry, coefU_bry,
    #                          LonV_bry, LatV_bry, iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry, elemV_bry, coefV_bry,
    #                          LonN_bry, LatN_bry, iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry, elemN_bry, coefN_bry,
    #                          LonP_bry, LatP_bry, iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry, elemP_bry, coefP_bry)

    #
    #
    # 1: SSH
    #
    #

    print('Interpolate SSH...')

    (zeta_bry, NzGood) = glor.interp_tracers(ncglo, 'zos', tndx_glo, -1,
                                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    # Get CROCO sigma coordinate at rho and w points (using zeta)
    #

    z_rho = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'r', vtransform)

    z_w = vgrd.zlevs(h_bry, zeta_bry, theta_s, theta_b, hc, N, 'w', vtransform)

    #
    #
    # 2: Temperature
    #
    #

    print('Interpolate Temperature...')

    temp_bry = glor.interp3d(ncglo, 'thetao', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 3: Salinity
    #
    #

    print('Interpolate Salinity...')

    salt_bry = glor.interp3d(ncglo, 'so', tndx_glo, Nzgoodmin, depthg, z_rho,
                             iminT_bry, imaxT_bry, jminT_bry, jmaxT_bry,
                             LonT_bry, LatT_bry, coefT_bry, elemT_bry, '_FillValue')

    #
    #
    # 4: U and V
    #
    # (interpolate on z levels at rho points - rotate to align with the grid -
    #  put to u and v points - vertical interpolation to sigma grid)
    #
    #

    cosa = np.cos(angle_bry)
    sina = np.sin(angle_bry)

    [u_bry, v_bry] = glor.interp3d_uv(ncglo, tndx_glo, Nzgoodmin, depthg, z_rho, cosa, sina,
                                      iminU_bry, imaxU_bry, jminU_bry, jmaxU_bry,
                                      LonU_bry, LatU_bry, coefU_bry, elemU_bry,
                                      iminV_bry, imaxV_bry, jminV_bry, jmaxV_bry,
                                      LonV_bry, LatV_bry, coefV_bry, elemV_bry)

    #
    #
    # 5: UBAR and VBAR
    #
    # Here it could be nice to get the barotropic transport from GLORYS, to put it on CROCO grid,
    # and to correct u and v accordingly...
    # But let's start simple and just integrate u and v on CROCO grid.
    #
    #

    (ubar_bry, h0) = vgrd.vintegr(u_bry, rho2u_3d(z_w), rho2u_3d(z_rho), np.nan, np.nan) / rho2u_2d(h_bry)
    (vbar_bry, h0) = vgrd.vintegr(v_bry, rho2v_3d(z_w), rho2v_3d(z_rho), np.nan, np.nan) / rho2v_2d(h_bry)

    if PNI is 1:
        # NORESM
        #
        #
        # 3: Nitrate
        #
        #

        ncnorni = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[0] + NORESM_ending, 'r')
        ncnornio = ncnorni.variables

        print('Interpolate Nitrate...')

        no3_bry1 = glor.interp3d(ncnornio, 'no3no2', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        no3_bry2 = glor.interp3d(ncnornio, 'no3no2', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        no3_bry3 = glor.interp3d(ncnornio, 'no3no2', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnorni.close()

        #
        #
        # 3: Orthophosphate
        #
        #

        ncnorpo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[1] + NORESM_ending, 'r')
        ncnorpoo = ncnorpo.variables

        print('Interpolate Orthophosphate...')

        po4_bry1 = glor.interp3d(ncnorpoo, 'po4', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        po4_bry2 = glor.interp3d(ncnorpoo, 'po4', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        po4_bry3 = glor.interp3d(ncnorpoo, 'po4', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnorpo.close()

        #
        #
        # 3: Silicate
        #
        #

        ncnorsi = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[2] + NORESM_ending, 'r')
        ncnorsio = ncnorsi.variables

        print('Interpolate Silicate...')

        si_bry1 = glor.interp3d(ncnorsio, 'si', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        si_bry2 = glor.interp3d(ncnorsio, 'si', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        si_bry3 = glor.interp3d(ncnorsio, 'si', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnorsi.close()

        #
        #
        # 3: Dissolved Oxygen
        #
        #

        ncnordo = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[3] + NORESM_ending, 'r')
        ncnordoo = ncnordo.variables

        print('Interpolate Dissolved Oxygen...')

        o2_bry1 = glor.interp3d(ncnordoo, 'o2', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        o2_bry2 = glor.interp3d(ncnordoo, 'o2', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        o2_bry3 = glor.interp3d(ncnordoo, 'o2', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnordo.close()

        #
        #
        # 3: Dissolved Inorganic Carbon
        #
        #

        ncnordic = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[4] + NORESM_ending, 'r')
        ncnordico = ncnordic.variables

        print('Interpolate Dissolved Inorganic Carbon...')

        dic_bry1 = glor.interp3d(ncnordico, 'dissic', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        dic_bry2 = glor.interp3d(ncnordico, 'dissic', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        dic_bry3 = glor.interp3d(ncnordico, 'dissic', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                 iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                 LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnordic.close()

        #
        #
        # 3: Total Alkalinity
        #
        #

        ncnoralkalini = netcdf(NORESMfiles_dir + NORESM_prefix + NIVAvars[5] + NORESM_ending, 'r')
        ncnoralkalinio = ncnoralkalini.variables

        print('Interpolate Total Alkalinity...')

        talk_bry1 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[0], Nzgoodmin, depthn, z_rho,
                                  iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                  LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        talk_bry2 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[1], Nzgoodmin, depthn, z_rho,
                                  iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                  LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')

        talk_bry3 = glor.interp3d(ncnoralkalinio, 'talk', Tidxn[2], Nzgoodmin, depthn, z_rho,
                                  iminN_bry, imaxN_bry, jminN_bry, jmaxN_bry,
                                  LonN_bry, LatN_bry, coefN_bry, elemN_bry, 'add_offset')
        ncnoralkalini.close()

        # PISCES
        #
        #
        # 3: Ammonium
        #
        #

        print('Interpolate Calcite...')

        calc_bry1 = glor.interp3d(ncpism1o, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        calc_bry2 = glor.interp3d(ncpiso, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        calc_bry3 = glor.interp3d(ncpisp1o, 'calc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5b: Particulate Organic Carbon from IBI PISCES
        #

        print('Interpolate Particulate Organic Carbon...')

        poc_bry1 = glor.interp3d(ncpism1o, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        poc_bry2 = glor.interp3d(ncpiso, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        poc_bry3 = glor.interp3d(ncpisp1o, 'poc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5c: Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Nanophytoplankton...')

        phy_bry1 = glor.interp3d(ncpism1o, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy_bry2 = glor.interp3d(ncpiso, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy_bry3 = glor.interp3d(ncpism1o, 'phy', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5d: Microzooplankton from IBI PISCES
        #

        print('Interpolate Microzooplankton...')

        zoo_bry1 = glor.interp3d(ncpism1o, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo_bry2 = glor.interp3d(ncpiso, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo_bry3 = glor.interp3d(ncpisp1o, 'zoo', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5e: Dissolved Organic Carbon from IBI PISCES
        #

        print('Interpolate Dissolved Organic Carbon...')

        doc_bry1 = glor.interp3d(ncpism1o, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        doc_bry2 = glor.interp3d(ncpiso, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        doc_bry3 = glor.interp3d(ncpisp1o, 'doc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5f: Diatom from IBI PISCES
        #

        print('Interpolate Diatom...')

        phy2_bry1 = glor.interp3d(ncpism1o, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy2_bry2 = glor.interp3d(ncpiso, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        phy2_bry3 = glor.interp3d(ncpisp1o, 'phy2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5g: Mesozooplankton from IBI PISCES
        #

        print('Interpolate Mesozooplankton...')

        zoo2_bry1 = glor.interp3d(ncpism1o, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo2_bry2 = glor.interp3d(ncpiso, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        zoo2_bry3 = glor.interp3d(ncpisp1o, 'zoo2', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5h: Biogenic Silica from IBI PISCES
        #

        print('Interpolate Biogenic Silica...')

        gsi_bry1 = glor.interp3d(ncpism1o, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        gsi_bry2 = glor.interp3d(ncpiso, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        gsi_bry3 = glor.interp3d(ncpisp1o, 'gsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5i: Dissolved Iron from IBI PISCES
        #

        print('Interpolate Dissolved Iron...')

        fe_bry1 = glor.interp3d(ncpism1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        fe_bry2 = glor.interp3d(ncpiso, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        fe_bry3 = glor.interp3d(ncpisp1o, 'fe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5j: Big Particle Iron from IBI PISCES
        #

        print('Interpolate Big Particle Iron...')

        bfe_bry1 = glor.interp3d(ncpism1o, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        bfe_bry2 = glor.interp3d(ncpiso, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        bfe_bry3 = glor.interp3d(ncpisp1o, 'bfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5k: Big Particulate Organic Carbon from IBI PISCES
        #

        print('Interpolate Big Particulate Organic Carbon...')

        goc_bry1 = glor.interp3d(ncpism1o, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        goc_bry2 = glor.interp3d(ncpiso, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        goc_bry3 = glor.interp3d(ncpisp1o, 'goc', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5l: Iron in the small particles from IBI PISCES
        #

        print('Interpolate Iron in the small particles...')

        sfe_bry1 = glor.interp3d(ncpism1o, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        sfe_bry2 = glor.interp3d(ncpiso, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        sfe_bry3 = glor.interp3d(ncpisp1o, 'sfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5m: Iron content of the diatoms from IBI PISCES
        #

        print('Interpolate Iron content of the diatoms...')

        dfe_bry1 = glor.interp3d(ncpism1o, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dfe_bry2 = glor.interp3d(ncpiso, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dfe_bry3 = glor.interp3d(ncpisp1o, 'dfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5n: Silicon content of the Diatoms from IBI PISCES
        #

        print('Interpolate Silicon content of the Diatoms...')

        dsi_bry1 = glor.interp3d(ncpism1o, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dsi_bry2 = glor.interp3d(ncpiso, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dsi_bry3 = glor.interp3d(ncpisp1o, 'dsi', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5o: Iron content of the Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Iron content of the Nanophytoplankton...')

        nfe_bry1 = glor.interp3d(ncpism1o, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nfe_bry2 = glor.interp3d(ncpiso, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nfe_bry3 = glor.interp3d(ncpisp1o, 'nfe', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5p: Nanophytoplankton Chlorophyll from IBI PISCES
        #

        print('Interpolate Nanophytoplankton Chlorophyll...')

        nchl_bry1 = glor.interp3d(ncpism1o, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nchl_bry2 = glor.interp3d(ncpiso, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nchl_bry3 = glor.interp3d(ncpisp1o, 'nchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5q: Nanophytoplankton from IBI PISCES
        #

        print('Interpolate Diatom Chlorophyll...')

        dchl_bry1 = glor.interp3d(ncpism1o, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dchl_bry2 = glor.interp3d(ncpiso, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        dchl_bry3 = glor.interp3d(ncpisp1o, 'dchl', tndx_glo, Nzgoodmin, depthp, z_rho,
                                  iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                  LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        #
        # 5r: Ammonium from IBI PISCES
        #

        print('Interpolate Ammonium...')

        nh4_bry1 = glor.interp3d(ncpism1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nh4_bry2 = glor.interp3d(ncpiso, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        nh4_bry3 = glor.interp3d(ncpisp1o, 'nh4', tndx_glo, Nzgoodmin, depthp, z_rho,
                                 iminP_bry, imaxP_bry, jminP_bry, jmaxP_bry,
                                 LonP_bry, LatP_bry, coefP_bry, elemP_bry, 'add_offset')

        if obctype == 's':
            # NORESM previous month
            ncbry['NO3_south'][0, :, :] = no3_bry1[:, 0, :] * 1000
            ncbry['PO4_south'][0, :, :] = po4_bry1[:, 0, :] * 1000
            ncbry['Si_south'][0, :, :] = si_bry1[:, 0, :] * 1000
            ncbry['O2_south'][0, :, :] = o2_bry1[:, 0, :] * 1000
            ncbry['DIC_south'][0, :, :] = dic_bry1[:, 0, :] * 1000
            ncbry['TALK_south'][0, :, :] = talk_bry1[:, 0, :] * 1000

            # PISCES previous month
            ncbry['DCHL_south'][0, :, :] = dchl_bry1[:, 0, :]
            ncbry['NCHL_south'][0, :, :] = nchl_bry1[:, 0, :]
            ncbry['NFE_south'][0, :, :] = nfe_bry1[:, 0, :]
            ncbry['DSI_south'][0, :, :] = dsi_bry1[:, 0, :]
            ncbry['DFE_south'][0, :, :] = dfe_bry1[:, 0, :]
            ncbry['SFE_south'][0, :, :] = sfe_bry1[:, 0, :]
            ncbry['GOC_south'][0, :, :] = goc_bry1[:, 0, :]
            ncbry['BFE_south'][0, :, :] = bfe_bry1[:, 0, :]
            ncbry['FER_south'][0, :, :] = fe_bry1[:, 0, :]
            ncbry['BSI_south'][0, :, :] = gsi_bry1[:, 0, :]
            ncbry['MESO_south'][0, :, :] = zoo2_bry1[:, 0, :]
            ncbry['DIA_south'][0, :, :] = phy2_bry1[:, 0, :]
            ncbry['DOC_south'][0, :, :] = doc_bry1[:, 0, :]
            ncbry['ZOO_south'][0, :, :] = zoo_bry1[:, 0, :]
            ncbry['NANO_south'][0, :, :] = phy_bry1[:, 0, :]
            ncbry['CACO3_south'][0, :, :] = calc_bry1[:, 0, :]
            ncbry['POC_south'][0, :, :] = poc_bry1[:, 0, :]
            ncbry['NH4_south'][0, :, :] = nh4_bry1[:, 0, :]

            # NORESM current month
            ncbry['NO3_south'][1, :, :] = no3_bry2[:, 0, :] * 1000
            ncbry['PO4_south'][1, :, :] = po4_bry2[:, 0, :] * 1000
            ncbry['Si_south'][1, :, :] = si_bry2[:, 0, :] * 1000
            ncbry['O2_south'][1, :, :] = o2_bry2[:, 0, :] * 1000
            ncbry['DIC_south'][1, :, :] = dic_bry2[:, 0, :] * 1000
            ncbry['TALK_south'][1, :, :] = talk_bry2[:, 0, :] * 1000

            # PISCES current month
            ncbry['DCHL_south'][1, :, :] = dchl_bry2[:, 0, :]
            ncbry['NCHL_south'][1, :, :] = nchl_bry2[:, 0, :]
            ncbry['NFE_south'][1, :, :] = nfe_bry2[:, 0, :]
            ncbry['DSI_south'][1, :, :] = dsi_bry2[:, 0, :]
            ncbry['DFE_south'][1, :, :] = dfe_bry2[:, 0, :]
            ncbry['SFE_south'][1, :, :] = sfe_bry2[:, 0, :]
            ncbry['GOC_south'][1, :, :] = goc_bry2[:, 0, :]
            ncbry['BFE_south'][1, :, :] = bfe_bry2[:, 0, :]
            ncbry['FER_south'][1, :, :] = fe_bry2[:, 0, :]
            ncbry['BSI_south'][1, :, :] = gsi_bry2[:, 0, :]
            ncbry['MESO_south'][1, :, :] = zoo2_bry2[:, 0, :]
            ncbry['DIA_south'][1, :, :] = phy2_bry2[:, 0, :]
            ncbry['DOC_south'][1, :, :] = doc_bry2[:, 0, :]
            ncbry['ZOO_south'][1, :, :] = zoo_bry2[:, 0, :]
            ncbry['NANO_south'][1, :, :] = phy_bry2[:, 0, :]
            ncbry['CACO3_south'][1, :, :] = calc_bry2[:, 0, :]
            ncbry['POC_south'][1, :, :] = poc_bry2[:, 0, :]
            ncbry['NH4_south'][1, :, :] = nh4_bry2[:, 0, :]

            # NORESM next month
            ncbry['NO3_south'][2, :, :] = no3_bry3[:, 0, :] * 1000
            ncbry['PO4_south'][2, :, :] = po4_bry3[:, 0, :] * 1000
            ncbry['Si_south'][2, :, :] = si_bry3[:, 0, :] * 1000
            ncbry['O2_south'][2, :, :] = o2_bry3[:, 0, :] * 1000
            ncbry['DIC_south'][2, :, :] = dic_bry3[:, 0, :] * 1000
            ncbry['TALK_south'][2, :, :] = talk_bry3[:, 0, :] * 1000

            # PISCES next month
            ncbry['DCHL_south'][2, :, :] = dchl_bry3[:, 0, :]
            ncbry['NCHL_south'][2, :, :] = nchl_bry3[:, 0, :]
            ncbry['NFE_south'][2, :, :] = nfe_bry3[:, 0, :]
            ncbry['DSI_south'][2, :, :] = dsi_bry3[:, 0, :]
            ncbry['DFE_south'][2, :, :] = dfe_bry3[:, 0, :]
            ncbry['SFE_south'][2, :, :] = sfe_bry3[:, 0, :]
            ncbry['GOC_south'][2, :, :] = goc_bry3[:, 0, :]
            ncbry['BFE_south'][2, :, :] = bfe_bry3[:, 0, :]
            ncbry['FER_south'][2, :, :] = fe_bry3[:, 0, :]
            ncbry['BSI_south'][2, :, :] = gsi_bry3[:, 0, :]
            ncbry['MESO_south'][2, :, :] = zoo2_bry3[:, 0, :]
            ncbry['DIA_south'][2, :, :] = phy2_bry3[:, 0, :]
            ncbry['DOC_south'][2, :, :] = doc_bry3[:, 0, :]
            ncbry['ZOO_south'][2, :, :] = zoo_bry3[:, 0, :]
            ncbry['NANO_south'][2, :, :] = phy_bry3[:, 0, :]
            ncbry['CACO3_south'][2, :, :] = calc_bry3[:, 0, :]
            ncbry['POC_south'][2, :, :] = poc_bry3[:, 0, :]
            ncbry['NH4_south'][2, :, :] = nh4_bry3[:, 0, :]

        elif obctype == 'n':
            # NORESM previous month
            ncbry['NO3_north'][0, :, :] = no3_bry1[:, -1, :] * 1000
            ncbry['PO4_north'][0, :, :] = po4_bry1[:, -1, :] * 1000
            ncbry['Si_north'][0, :, :] = si_bry1[:, -1, :] * 1000
            ncbry['O2_north'][0, :, :] = o2_bry1[:, -1, :] * 1000
            ncbry['DIC_north'][0, :, :] = dic_bry1[:, -1, :] * 1000
            ncbry['TALK_north'][0, :, :] = talk_bry1[:, -1, :] * 1000

            # PISCES previous month
            ncbry['DCHL_north'][0, :, :] = dchl_bry1[:, -1, :]
            ncbry['NCHL_north'][0, :, :] = nchl_bry1[:, -1, :]
            ncbry['NFE_north'][0, :, :] = nfe_bry1[:, -1, :]
            ncbry['DSI_north'][0, :, :] = dsi_bry1[:, -1, :]
            ncbry['DFE_north'][0, :, :] = dfe_bry1[:, -1, :]
            ncbry['SFE_north'][0, :, :] = sfe_bry1[:, -1, :]
            ncbry['GOC_north'][0, :, :] = goc_bry1[:, -1, :]
            ncbry['BFE_north'][0, :, :] = bfe_bry1[:, -1, :]
            ncbry['FER_north'][0, :, :] = fe_bry1[:, -1, :]
            ncbry['BSI_north'][0, :, :] = gsi_bry1[:, -1, :]
            ncbry['MESO_north'][0, :, :] = zoo2_bry1[:, -1, :]
            ncbry['DIA_north'][0, :, :] = phy2_bry1[:, -1, :]
            ncbry['DOC_north'][0, :, :] = doc_bry1[:, -1, :]
            ncbry['ZOO_north'][0, :, :] = zoo_bry1[:, -1, :]
            ncbry['NANO_north'][0, :, :] = phy_bry1[:, -1, :]
            ncbry['CACO3_north'][0, :, :] = calc_bry1[:, -1, :]
            ncbry['POC_north'][0, :, :] = poc_bry1[:, -1, :]
            ncbry['NH4_north'][0, :, :] = nh4_bry1[:, -1, :]

            # NORESM current month
            ncbry['NO3_north'][1, :, :] = no3_bry2[:, -1, :] * 1000
            ncbry['PO4_north'][1, :, :] = po4_bry2[:, -1, :] * 1000
            ncbry['Si_north'][1, :, :] = si_bry2[:, -1, :] * 1000
            ncbry['O2_north'][1, :, :] = o2_bry2[:, -1, :] * 1000
            ncbry['DIC_north'][1, :, :] = dic_bry2[:, -1, :] * 1000
            ncbry['TALK_north'][1, :, :] = talk_bry2[:, -1, :] * 1000

            # PISCES current month
            ncbry['DCHL_north'][1, :, :] = dchl_bry2[:, -1, :]
            ncbry['NCHL_north'][1, :, :] = nchl_bry2[:, -1, :]
            ncbry['NFE_north'][1, :, :] = nfe_bry2[:, -1, :]
            ncbry['DSI_north'][1, :, :] = dsi_bry2[:, -1, :]
            ncbry['DFE_north'][1, :, :] = dfe_bry2[:, -1, :]
            ncbry['SFE_north'][1, :, :] = sfe_bry2[:, -1, :]
            ncbry['GOC_north'][1, :, :] = goc_bry2[:, -1, :]
            ncbry['BFE_north'][1, :, :] = bfe_bry2[:, -1, :]
            ncbry['FER_north'][1, :, :] = fe_bry2[:, -1, :]
            ncbry['BSI_north'][1, :, :] = gsi_bry2[:, -1, :]
            ncbry['MESO_north'][1, :, :] = zoo2_bry2[:, -1, :]
            ncbry['DIA_north'][1, :, :] = phy2_bry2[:, -1, :]
            ncbry['DOC_north'][1, :, :] = doc_bry2[:, -1, :]
            ncbry['ZOO_north'][1, :, :] = zoo_bry2[:, -1, :]
            ncbry['NANO_north'][1, :, :] = phy_bry2[:, -1, :]
            ncbry['CACO3_north'][1, :, :] = calc_bry2[:, -1, :]
            ncbry['POC_north'][1, :, :] = poc_bry2[:, -1, :]
            ncbry['NH4_north'][1, :, :] = nh4_bry2[:, -1, :]

            # NORESM next month
            ncbry['NO3_north'][2, :, :] = no3_bry3[:, -1, :] * 1000
            ncbry['PO4_north'][2, :, :] = po4_bry3[:, -1, :] * 1000
            ncbry['Si_north'][2, :, :] = si_bry3[:, -1, :] * 1000
            ncbry['O2_north'][2, :, :] = o2_bry3[:, -1, :] * 1000
            ncbry['DIC_north'][2, :, :] = dic_bry3[:, -1, :] * 1000
            ncbry['TALK_north'][2, :, :] = talk_bry3[:, -1, :] * 1000

            # PISCES next month
            ncbry['DCHL_north'][2, :, :] = dchl_bry3[:, -1, :]
            ncbry['NCHL_north'][2, :, :] = nchl_bry3[:, -1, :]
            ncbry['NFE_north'][2, :, :] = nfe_bry3[:, -1, :]
            ncbry['DSI_north'][2, :, :] = dsi_bry3[:, -1, :]
            ncbry['DFE_north'][2, :, :] = dfe_bry3[:, -1, :]
            ncbry['SFE_north'][2, :, :] = sfe_bry3[:, -1, :]
            ncbry['GOC_north'][2, :, :] = goc_bry3[:, -1, :]
            ncbry['BFE_north'][2, :, :] = bfe_bry3[:, -1, :]
            ncbry['FER_north'][2, :, :] = fe_bry3[:, -1, :]
            ncbry['BSI_north'][2, :, :] = gsi_bry3[:, -1, :]
            ncbry['MESO_north'][2, :, :] = zoo2_bry3[:, -1, :]
            ncbry['DIA_north'][2, :, :] = phy2_bry3[:, -1, :]
            ncbry['DOC_north'][2, :, :] = doc_bry3[:, -1, :]
            ncbry['ZOO_north'][2, :, :] = zoo_bry3[:, -1, :]
            ncbry['NANO_north'][2, :, :] = phy_bry3[:, -1, :]
            ncbry['CACO3_north'][2, :, :] = calc_bry3[:, -1, :]
            ncbry['POC_north'][2, :, :] = poc_bry3[:, -1, :]
            ncbry['NH4_north'][2, :, :] = nh4_bry3[:, -1, :]

        elif obctype == 'e':
            # NORESM previous month
            ncbry['NO3_east'][0, :, :] = no3_bry1[:, :, 0] * 1000
            ncbry['PO4_east'][0, :, :] = po4_bry1[:, :, 0] * 1000
            ncbry['Si_east'][0, :, :] = si_bry1[:, :, 0] * 1000
            ncbry['O2_east'][0, :, :] = o2_bry1[:, :, 0] * 1000
            ncbry['DIC_east'][0, :, :] = dic_bry1[:, :, 0] * 1000
            ncbry['TALK_east'][0, :, :] = talk_bry1[:, :, 0] * 1000

            # PISCES previous month
            ncbry['DCHL_east'][0, :, :] = dchl_bry1[:, :, 0]
            ncbry['NCHL_east'][0, :, :] = nchl_bry1[:, :, 0]
            ncbry['NFE_east'][0, :, :] = nfe_bry1[:, :, 0]
            ncbry['DSI_east'][0, :, :] = dsi_bry1[:, :, 0]
            ncbry['DFE_east'][0, :, :] = dfe_bry1[:, :, 0]
            ncbry['SFE_east'][0, :, :] = sfe_bry1[:, :, 0]
            ncbry['GOC_east'][0, :, :] = goc_bry1[:, :, 0]
            ncbry['BFE_east'][0, :, :] = bfe_bry1[:, :, 0]
            ncbry['FER_east'][0, :, :] = fe_bry1[:, :, 0]
            ncbry['BSI_east'][0, :, :] = gsi_bry1[:, :, 0]
            ncbry['MESO_east'][0, :, :] = zoo2_bry1[:, :, 0]
            ncbry['DIA_east'][0, :, :] = phy2_bry1[:, :, 0]
            ncbry['DOC_east'][0, :, :] = doc_bry1[:, :, 0]
            ncbry['ZOO_east'][0, :, :] = zoo_bry1[:, :, 0]
            ncbry['NANO_east'][0, :, :] = phy_bry1[:, :, 0]
            ncbry['CACO3_east'][0, :, :] = calc_bry1[:, :, 0]
            ncbry['POC_east'][0, :, :] = poc_bry1[:, :, 0]
            ncbry['NH4_east'][0, :, :] = nh4_bry1[:, :, 0]

            # NORESM current month
            ncbry['NO3_east'][1, :, :] = no3_bry2[:, :, 0] * 1000
            ncbry['PO4_east'][1, :, :] = po4_bry2[:, :, 0] * 1000
            ncbry['Si_east'][1, :, :] = si_bry2[:, :, 0] * 1000
            ncbry['O2_east'][1, :, :] = o2_bry2[:, :, 0] * 1000
            ncbry['DIC_east'][1, :, :] = dic_bry2[:, :, 0] * 1000
            ncbry['TALK_east'][1, :, :] = talk_bry2[:, :, 0] * 1000

            # PISCES current month
            ncbry['DCHL_east'][1, :, :] = dchl_bry2[:, :, 0]
            ncbry['NCHL_east'][1, :, :] = nchl_bry2[:, :, 0]
            ncbry['NFE_east'][1, :, :] = nfe_bry2[:, :, 0]
            ncbry['DSI_east'][1, :, :] = dsi_bry2[:, :, 0]
            ncbry['DFE_east'][1, :, :] = dfe_bry2[:, :, 0]
            ncbry['SFE_east'][1, :, :] = sfe_bry2[:, :, 0]
            ncbry['GOC_east'][1, :, :] = goc_bry2[:, :, 0]
            ncbry['BFE_east'][1, :, :] = bfe_bry2[:, :, 0]
            ncbry['FER_east'][1, :, :] = fe_bry2[:, :, 0]
            ncbry['BSI_east'][1, :, :] = gsi_bry2[:, :, 0]
            ncbry['MESO_east'][1, :, :] = zoo2_bry2[:, :, 0]
            ncbry['DIA_east'][1, :, :] = phy2_bry2[:, :, 0]
            ncbry['DOC_east'][1, :, :] = doc_bry2[:, :, 0]
            ncbry['ZOO_east'][1, :, :] = zoo_bry2[:, :, 0]
            ncbry['NANO_east'][1, :, :] = phy_bry2[:, :, 0]
            ncbry['CACO3_east'][1, :, :] = calc_bry2[:, :, 0]
            ncbry['POC_east'][1, :, :] = poc_bry2[:, :, 0]
            ncbry['NH4_east'][1, :, :] = nh4_bry2[:, :, 0]

            # NORESM next month
            ncbry['NO3_east'][2, :, :] = no3_bry3[:, :, 0] * 1000
            ncbry['PO4_east'][2, :, :] = po4_bry3[:, :, 0] * 1000
            ncbry['Si_east'][2, :, :] = si_bry3[:, :, 0] * 1000
            ncbry['O2_east'][2, :, :] = o2_bry3[:, :, 0] * 1000
            ncbry['DIC_east'][2, :, :] = dic_bry3[:, :, 0] * 1000
            ncbry['TALK_east'][2, :, :] = talk_bry3[:, :, 0] * 1000

            # PISCES next month
            ncbry['DCHL_east'][2, :, :] = dchl_bry3[:, :, 0]
            ncbry['NCHL_east'][2, :, :] = nchl_bry3[:, :, 0]
            ncbry['NFE_east'][2, :, :] = nfe_bry3[:, :, 0]
            ncbry['DSI_east'][2, :, :] = dsi_bry3[:, :, 0]
            ncbry['DFE_east'][2, :, :] = dfe_bry3[:, :, 0]
            ncbry['SFE_east'][2, :, :] = sfe_bry3[:, :, 0]
            ncbry['GOC_east'][2, :, :] = goc_bry3[:, :, 0]
            ncbry['BFE_east'][2, :, :] = bfe_bry3[:, :, 0]
            ncbry['FER_east'][2, :, :] = fe_bry3[:, :, 0]
            ncbry['BSI_east'][2, :, :] = gsi_bry3[:, :, 0]
            ncbry['MESO_east'][2, :, :] = zoo2_bry3[:, :, 0]
            ncbry['DIA_east'][2, :, :] = phy2_bry3[:, :, 0]
            ncbry['DOC_east'][2, :, :] = doc_bry3[:, :, 0]
            ncbry['ZOO_east'][2, :, :] = zoo_bry3[:, :, 0]
            ncbry['NANO_east'][2, :, :] = phy_bry3[:, :, 0]
            ncbry['CACO3_east'][2, :, :] = calc_bry3[:, :, 0]
            ncbry['POC_east'][2, :, :] = poc_bry3[:, :, 0]
            ncbry['NH4_east'][2, :, :] = nh4_bry3[:, :, 0]

        elif obctype == 'w':
            # NORESM previous month
            ncbry['NO3_west'][0, :, :] = no3_bry1[:, :, -1] * 1000
            ncbry['PO4_west'][0, :, :] = po4_bry1[:, :, -1] * 1000
            ncbry['Si_west'][0, :, :] = si_bry1[:, :, -1] * 1000
            ncbry['O2_west'][0, :, :] = o2_bry1[:, :, -1] * 1000
            ncbry['DIC_west'][0, :, :] = dic_bry1[:, :, -1] * 1000
            ncbry['TALK_west'][0, :, :] = talk_bry1[:, :, -1] * 1000

            # PISCES previous month
            ncbry['DCHL_west'][0, :, :] = dchl_bry1[:, :, -1]
            ncbry['NCHL_west'][0, :, :] = nchl_bry1[:, :, -1]
            ncbry['NFE_west'][0, :, :] = nfe_bry1[:, :, -1]
            ncbry['DSI_west'][0, :, :] = dsi_bry1[:, :, -1]
            ncbry['DFE_west'][0, :, :] = dfe_bry1[:, :, -1]
            ncbry['SFE_west'][0, :, :] = sfe_bry1[:, :, -1]
            ncbry['GOC_west'][0, :, :] = goc_bry1[:, :, -1]
            ncbry['BFE_west'][0, :, :] = bfe_bry1[:, :, -1]
            ncbry['FER_west'][0, :, :] = fe_bry1[:, :, -1]
            ncbry['BSI_west'][0, :, :] = gsi_bry1[:, :, -1]
            ncbry['MESO_west'][0, :, :] = zoo2_bry1[:, :, -1]
            ncbry['DIA_west'][0, :, :] = phy2_bry1[:, :, -1]
            ncbry['DOC_west'][0, :, :] = doc_bry1[:, :, -1]
            ncbry['ZOO_west'][0, :, :] = zoo_bry1[:, :, -1]
            ncbry['NANO_west'][0, :, :] = phy_bry1[:, :, -1]
            ncbry['CACO3_west'][0, :, :] = calc_bry1[:, :, -1]
            ncbry['POC_west'][0, :, :] = poc_bry1[:, :, -1]
            ncbry['NH4_west'][0, :, :] = nh4_bry1[:, :, -1]

            # NORESM current month
            ncbry['NO3_west'][1, :, :] = no3_bry2[:, :, -1] * 1000
            ncbry['PO4_west'][1, :, :] = po4_bry2[:, :, -1] * 1000
            ncbry['Si_west'][1, :, :] = si_bry2[:, :, -1] * 1000
            ncbry['O2_west'][1, :, :] = o2_bry2[:, :, -1] * 1000
            ncbry['DIC_west'][1, :, :] = dic_bry2[:, :, -1] * 1000
            ncbry['TALK_west'][1, :, :] = talk_bry2[:, :, -1] * 1000

            # PISCES current month
            ncbry['DCHL_west'][1, :, :] = dchl_bry2[:, :, -1]
            ncbry['NCHL_west'][1, :, :] = nchl_bry2[:, :, -1]
            ncbry['NFE_west'][1, :, :] = nfe_bry2[:, :, -1]
            ncbry['DSI_west'][1, :, :] = dsi_bry2[:, :, -1]
            ncbry['DFE_west'][1, :, :] = dfe_bry2[:, :, -1]
            ncbry['SFE_west'][1, :, :] = sfe_bry2[:, :, -1]
            ncbry['GOC_west'][1, :, :] = goc_bry2[:, :, -1]
            ncbry['BFE_west'][1, :, :] = bfe_bry2[:, :, -1]
            ncbry['FER_west'][1, :, :] = fe_bry2[:, :, -1]
            ncbry['BSI_west'][1, :, :] = gsi_bry2[:, :, -1]
            ncbry['MESO_west'][1, :, :] = zoo2_bry2[:, :, -1]
            ncbry['DIA_west'][1, :, :] = phy2_bry2[:, :, -1]
            ncbry['DOC_west'][1, :, :] = doc_bry2[:, :, -1]
            ncbry['ZOO_west'][1, :, :] = zoo_bry2[:, :, -1]
            ncbry['NANO_west'][1, :, :] = phy_bry2[:, :, -1]
            ncbry['CACO3_west'][1, :, :] = calc_bry2[:, :, -1]
            ncbry['POC_west'][1, :, :] = poc_bry2[:, :, -1]
            ncbry['NH4_west'][1, :, :] = nh4_bry2[:, :, -1]

            # NORESM next month
            ncbry['NO3_west'][2, :, :] = no3_bry3[:, :, -1] * 1000
            ncbry['PO4_west'][2, :, :] = po4_bry3[:, :, -1] * 1000
            ncbry['Si_west'][2, :, :] = si_bry3[:, :, -1] * 1000
            ncbry['O2_west'][2, :, :] = o2_bry3[:, :, -1] * 1000
            ncbry['DIC_west'][2, :, :] = dic_bry3[:, :, -1] * 1000
            ncbry['TALK_west'][2, :, :] = talk_bry3[:, :, -1] * 1000

            # PISCES next month
            ncbry['DCHL_west'][2, :, :] = dchl_bry3[:, :, -1]
            ncbry['NCHL_west'][2, :, :] = nchl_bry3[:, :, -1]
            ncbry['NFE_west'][2, :, :] = nfe_bry3[:, :, -1]
            ncbry['DSI_west'][2, :, :] = dsi_bry3[:, :, -1]
            ncbry['DFE_west'][2, :, :] = dfe_bry3[:, :, -1]
            ncbry['SFE_west'][2, :, :] = sfe_bry3[:, :, -1]
            ncbry['GOC_west'][2, :, :] = goc_bry3[:, :, -1]
            ncbry['BFE_west'][2, :, :] = bfe_bry3[:, :, -1]
            ncbry['FER_west'][2, :, :] = fe_bry3[:, :, -1]
            ncbry['BSI_west'][2, :, :] = gsi_bry3[:, :, -1]
            ncbry['MESO_west'][2, :, :] = zoo2_bry3[:, :, -1]
            ncbry['DIA_west'][2, :, :] = phy2_bry3[:, :, -1]
            ncbry['DOC_west'][2, :, :] = doc_bry3[:, :, -1]
            ncbry['ZOO_west'][2, :, :] = zoo_bry3[:, :, -1]
            ncbry['NANO_west'][2, :, :] = phy_bry3[:, :, -1]
            ncbry['CACO3_west'][2, :, :] = calc_bry3[:, :, -1]
            ncbry['POC_west'][2, :, :] = poc_bry3[:, :, -1]
            ncbry['NH4_west'][2, :, :] = nh4_bry3[:, :, -1]

    if obctype == 's':
        # PHYSICS
        ncbry['zeta_south'][tndx_bry, :] = zeta_bry[0, :]
        ncbry['temp_south'][tndx_bry, :, :] = temp_bry[:, 0, :]
        ncbry['salt_south'][tndx_bry, :, :] = salt_bry[:, 0, :]
        ncbry['u_south'][tndx_bry, :, :] = u_bry[:, 0, :]
        ncbry['v_south'][tndx_bry, :, :] = v_bry[:, 0, :]
        ncbry['ubar_south'][tndx_bry, :] = ubar_bry[0, :]
        ncbry['vbar_south'][tndx_bry, :] = vbar_bry[0, :]

    elif obctype == 'n':
        # PHYSICS
        ncbry['zeta_north'][tndx_bry, :] = zeta_bry[-1, :]
        ncbry['temp_north'][tndx_bry, :, :] = temp_bry[:, -1, :]
        ncbry['salt_north'][tndx_bry, :, :] = salt_bry[:, -1, :]
        ncbry['u_north'][tndx_bry, :, :] = u_bry[:, -1, :]
        ncbry['v_north'][tndx_bry, :, :] = v_bry[:, -1, :]
        ncbry['ubar_north'][tndx_bry, :] = ubar_bry[-1, :]
        ncbry['vbar_north'][tndx_bry, :] = vbar_bry[-1, :]

    elif obctype == 'e':
        # PHYSICS
        ncbry['zeta_east'][tndx_bry, :] = zeta_bry[:, 0]
        ncbry['temp_east'][tndx_bry, :, :] = temp_bry[:, :, 0]
        ncbry['salt_east'][tndx_bry, :, :] = salt_bry[:, :, 0]
        ncbry['u_east'][tndx_bry, :, :] = u_bry[:, :, 0]
        ncbry['v_east'][tndx_bry, :, :] = v_bry[:, :, 0]
        ncbry['ubar_east'][tndx_bry, :] = ubar_bry[:, 0]
        ncbry['vbar_east'][tndx_bry, :] = vbar_bry[:, 0]

    elif obctype == 'w':
        # PHYSICS
        ncbry['zeta_west'][tndx_bry, :] = zeta_bry[:, -1]
        ncbry['temp_west'][tndx_bry, :, :] = temp_bry[:, :, -1]
        ncbry['salt_west'][tndx_bry, :, :] = salt_bry[:, :, -1]
        ncbry['u_west'][tndx_bry, :, :] = u_bry[:, :, -1]
        ncbry['v_west'][tndx_bry, :, :] = v_bry[:, :, -1]
        ncbry['ubar_west'][tndx_bry, :] = ubar_bry[:, -1]
        ncbry['vbar_west'][tndx_bry, :] = vbar_bry[:, -1]

    return ncbry
