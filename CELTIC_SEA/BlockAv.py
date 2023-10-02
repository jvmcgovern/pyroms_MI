import pandas as pd
import os
import glob
import numpy as np
import numpy.matlib
import platform
from netCDF4 import Dataset, num2date
import time
import h5py
import gzip
import shutil
import tarfile
import pickle
import pathlib

# Purpose: average finer dataset onto the grid of coarser dataset.
# Input: Fine gridded dataset, coarse gridded dataset.

# Setting up dictionary of netCDF items (subject to revision)
VarDict = \
{'chl': 'chlorophyll', 'CHL': 'chlorophyll', 'Chlorophyll': 'chlorophyll',
 'npp': 'NPP', 'nppv': 'NPP', 'PRP': 'NPP', 'PPmasked': 'NPP',
 'zeu': 'eup', 'ZEU_mean': 'eup', 'EUP': 'eup',
 'n_an': 'no3', 'no3': 'no3', 'NO3': 'no3',
 'p_an': 'po4', 'po4': 'po4', 'PO4': 'po4',
 'o_an': 'o2', 'o2': 'o2', 'Oxygen': 'o2', 'OXI': 'o2',
 'i_an': 'si', 'si': 'si', 'SIL': 'si',
 'ph': 'ph',
 'fco2_ave_weighted': 'spco2', 'spco2_raw': 'spco2', 'spco2': 'spco2',
 'dic': 'dic',
 'latitude': 'latitude', 'nav_lat': 'latitude', 'lat': 'latitude', 'ylat': 'latitude', 'LATITUDE': 'latitude',
 'longitude': 'longitude', 'nav_lon': 'longitude', 'lon': 'longitude', 'xlon': 'longitude', 'LONGITUDE': 'longitude',
 'depth': 'depth', 'deptht': 'depth', 'z': 'depth'}


TimeDict = {'time_counter': 'time', 'time': 'time', 'tmnth': 'time', 'date': 'time',  # Climatologies
            'date_time': 'time', 'TIME': 'time'}                                      # ERDDAP, gliders etc.

platchk = platform.system()
if platchk == 'Linux':
    slash = '/'
else:  # platchk == 'Windows'
    slash = '\\'


# Finding the list of file names, paths, and the start and end date of the file sequence
def freview(fpath, fwild, fdatype, tvars):
    fllist = glob.glob(os.path.join(fpath, fwild))
    flpath = pd.DataFrame(fllist, columns=['FullPath'])
    flpath[['foldpath', 'filNCnam']] = flpath['FullPath'].str.rsplit(pat=slash, n=1, expand=True)
    # from file name, determine dates
    if len(flpath) == 1:  # single file, contains either multiple times or single time, with time variable
        onefile = True
        fil_nc = Dataset(os.path.join(flpath['foldpath'].values[0], flpath['filNCnam'].values[0]), 'r')
        if 'clm' not in fdatype:  # dealing with data, IBI model outputs never presented as single file
            tlist, _ = keysearch(fil_nc.variables, tvars)
            tQ = fil_nc.variables[tlist[0]]
            t_unit = fil_nc.variables[tlist[0]].units
            try:
                t_cal = fil_nc.variables[tlist[0]].calendar
            except:
                t_cal = "standard"
            flstrd = num2date(tQ[0], units=t_unit, calendar=t_cal, only_use_cftime_datetimes=False)  # 1st time
            flendd = num2date(tQ[-1], units=t_unit, calendar=t_cal, only_use_cftime_datetimes=False) # Last time
            dtcalc = abs(flendd-flstrd).days
            dtdelta = dtcalc/len(tQ)
        elif 'clm' in fdatype:  # dealing with a climatology, monthly, annual, or (rarely,) daily climatology
            flstrd = []
            flendd = []
            if fdatype == 'dclm':
                dtdelta = 1
            elif fdatype == 'mclm':
                dtdelta = 30
            elif fdatype == 'yclm':
                dtdelta = 365
    elif len(flpath) == 12 or (len(flpath) % 12 == 0):  # probably monthly files, set date as 15th of each month
        onefile = False
        if 'clm' not in fdatype:  # dealing with data
            flpath[['date6']] = flpath['filNCnam'].str.extract(r'(\d{6})')
            flpath[['Day']] = '15'
            flpath['DateStr'] = pd.to_datetime(flpath['date6'] + flpath['Day'], format='%Y%m%d')
            flpath.sort_values(by=['DateStr'], inplace=True, ascending=True)
            flstrd = flpath['DateStr'].iloc[0]
            flendd = flpath['DateStr'].iloc[-1]
            dtcalc = abs(flendd-flstrd).days
            dtdelta = dtcalc/len(flpath)
        elif 'clm' in fdatype:  # dealing with a climatology, monthly, annual, or (rarely,) daily climatology
            flstrd = []
            flendd = []
            if fdatype == 'dclm':
                dtdelta = 1
            elif fdatype == 'mclm':
                dtdelta = 30
            elif fdatype == 'yclm':
                dtdelta = 365
            flpath[['digits']] = flpath['filNCnam'].str.findall(r'\d+')
            flpath.sort_values(by=['digits'], inplace=True, ascending=True)
            flpath.reset_index(drop=True, inplace=True)
    else:  # probably daily files with date in filename
        onefile = False
        flpath[['date8']] = flpath['filNCnam'].str.extract(r'(\d{8})')
        flpath['DateStr'] = pd.to_datetime(flpath['date8'], format='%Y%m%d')
        flpath.sort_values(by=['DateStr'], inplace=True, ascending=True)
        flstrd = flpath['DateStr'].iloc[0]
        flendd = flpath['DateStr'].iloc[-1]
        dtcalc = abs(flendd-flstrd).days
        dtdelta = dtcalc/len(flpath)
    if 'clm' in fdatype:
        if fdatype == 'yclm':
            dt = 'zy'
        elif fdatype == 'mclm':
            dt = 'm'
        elif fdatype == 'dclm':
            dt = 'd'
        else:
            if 0.9 < dtdelta < 1.1:
                dt = 'd'
            elif 25 < dtdelta < 35:
                dt = 'm'
            else:
                dt = 'zy'
    else:
        if 0.9 < dtdelta < 1.1:
            dt = 'd'
        elif 25 < dtdelta < 35:
            dt = 'm'
        elif flstrd == flendd or flstrd and flendd == []:
            dt = 'zy'  # deltadate is zero (z), or, annual climatology (y)
        else:
            dt = 'zy'

    return fllist, flpath, flstrd, flendd, dt, onefile


def freview_OSU(fpath, fwild, fdatype):
    fllist_hg = glob.glob(os.path.join(fpath, fwild + '.hdf.gz'))
    fllist_h = glob.glob(os.path.join(fpath, fwild + '.hdf'))
    if (len(fllist_hg) > 0) & (len(fllist_hg) > len(fllist_h)):
        flpath_hg = pd.DataFrame(fllist_hg, columns=['FullPath'])
        flpath_hg[['foldpath', 'filnam']] = flpath_hg['FullPath'].str.rsplit(pat=slash, n=1, expand=True)
        for hgi in range(len(flpath_hg)):
            with gzip.open(os.path.join(flpath_hg['foldpath'].values[hgi],
                                        flpath_hg['filnam'].values[hgi]), 'rb') as f_in:
                with open(os.path.join(flpath_hg['foldpath'].values[hgi],
                                       flpath_hg['filnam'].values[hgi][:-3]), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
    elif len(fllist_hg) == 0:
        fllist_t = glob.glob(os.path.join(fpath, '*.tar'))
        flpath_t = pd.DataFrame(fllist_t, columns=['FullPath'])
        flpath_t[['foldpath', 'filnam']] = flpath_t['FullPath'].str.rsplit(pat=slash, n=1, expand=True)
        for t in range(len(flpath_t)):
            tar = tarfile.open(os.path.join(flpath_t['foldpath'].values[t], flpath_t['filnam'].values[t]))
            tar.extractall()
            tar.close()
        fllist_hg = glob.glob(os.path.join(fpath, fwild + '.hdf.gz'))
        flpath_hg = pd.DataFrame(fllist_hg, columns=['FullPath'])
        flpath_hg[['foldpath', 'filnam']] = flpath_hg['FullPath'].str.rsplit(pat=slash, n=1, expand=True)
        for hgi in range(len(flpath_hg)):
            with gzip.open(os.path.join(flpath_hg['foldpath'].values[hgi],
                                        flpath_hg['filnam'].values[hgi]), 'rb') as f_in:
                with open(os.path.join(flpath_hg['foldpath'].values[hgi],
                                       flpath_hg['filnam'].values[hgi][:-3]), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

    fllist_h = glob.glob(os.path.join(fpath, fwild + '.hdf'))
    flpath_h = pd.DataFrame(fllist_h, columns=['FullPath'])
    flpath_h[['foldpath', 'filnam']] = flpath_h['FullPath'].str.rsplit(pat=slash, n=1, expand=True)

    # from file name, determine dates
    if len(flpath_h) == 12 or (len(flpath_h) % 12 == 0):  # probably monthly files, set date as 15th of each month
        onefile = False
        if 'clm' not in fdatype:  # dealing with data
            flpath_h[['date7']] = flpath_h['filnam'].str.extract(r'(\d{7})')
            flpath_h[['Day']] = '15'
            flpath_h['DateStr'] = pd.to_datetime(flpath_h['date7'], format='%Y%j')
            flpath_h.sort_values(by=['DateStr'], inplace=True, ascending=True)
            flstrd = flpath_h['DateStr'].iloc[0]
            flendd = flpath_h['DateStr'].iloc[-1]
            dtcalc = abs(flendd-flstrd).days
            dtdelta = dtcalc/len(flpath_h)
        elif 'clm' in fdatype:  # dealing with a climatology, monthly, annual, or (rarely,) daily climatology
            flstrd = []
            flendd = []
            if fdatype == 'dclm':
                dtdelta = 1
            elif fdatype == 'mclm':
                dtdelta = 30
            elif fdatype == 'yclm':
                dtdelta = 365
            flpath_h[['digits']] = flpath_h['filnam'].str.findall(r'\d+')
            flpath_h.sort_values(by=['digits'], inplace=True, ascending=True)
            flpath_h.reset_index(drop=True, inplace=True)
    else:  # probably daily files with date in filename
        onefile = False
        flpath_h[['date7']] = flpath_h['filnam'].str.extract(r'(\d{7})')
        flpath_h['DateStr'] = pd.to_datetime(flpath_h['date7'], format='%Y%j')
        flpath_h.sort_values(by=['DateStr'], inplace=True, ascending=True)
        flstrd = flpath_h['DateStr'].iloc[0]
        flendd = flpath_h['DateStr'].iloc[-1]
        dtcalc = abs(flendd-flstrd).days
        dtdelta = dtcalc/len(flpath_h)
    if 'clm' in fdatype:
        if fdatype == 'yclm':
            dt = 'zy'
        elif fdatype == 'mclm':
            dt = 'm'
        elif fdatype == 'dclm':
            dt = 'd'
        else:
            if 0.9 < dtdelta < 1.1:
                dt = 'd'
            elif 25 < dtdelta < 35:
                dt = 'm'
            else:
                dt = 'zy'
    else:
        if 0.9 < dtdelta < 1.1:
            dt = 'd'
        elif 25 < dtdelta < 35:
            dt = 'm'
        elif flstrd == flendd or flstrd and flendd == []:
            dt = 'zy'  # deltadate is zero (z), or, annual climatology (y)
        else:
            dt = 'zy'

    return fllist_h, flpath_h, flstrd, flendd, dt, onefile


# Finding which keys in the netcdf are common to the dictionary of terms
# List is outputted as latitude, longitude, variable, and depth if possible
def keysearch(n, f):
    syek = [] # keys
    seulav = [] # values
    for item in f.keys():
        if item in n.keys():
            syek.append(item)
            seulav.append(f.get(item))
    return syek, seulav


# Outputting the common variables in order, plus whether full depth or surface/2D
def parcommon(eno, owt, foc):
    eno_nc = Dataset(os.path.join(eno['foldpath'].values[0], eno['filNCnam'].values[0]), 'r')
    owt_nc = Dataset(os.path.join(owt['foldpath'].values[0], owt['filNCnam'].values[0]), 'r')
    inter_eno, vals_eno = keysearch(eno_nc.variables, foc)
    inter_owt, vals_owt = keysearch(owt_nc.variables, foc)
    commonvars = [value for value in vals_eno if value in vals_owt]
    try:
        depthidx_eno = [vals_eno.index('depth')]
        depth_eno = [inter_eno[i] for i in depthidx_eno]
    except:
        pass
    try:
        depthidx_owt = [vals_owt.index('depth')]
        depth_owt = [inter_owt[i] for i in depthidx_owt]
    except:
        pass
    indices_eno = [vals_eno.index(x) for x in commonvars]
    indices_owt = [vals_owt.index(x) for x in commonvars]
    inter_eno = [inter_eno[i] for i in indices_eno]
    inter_owt = [inter_owt[i] for i in indices_owt]
    if 'depth' not in vals_eno and 'depth' in vals_owt:
        dfl_eno = 'surf'      # file type 1 is at surface
        dfl_owt = 'depth'     # file type 2 is full depth
        inter_owt.append(depth_owt[0])
    elif 'depth' in vals_eno and 'depth' not in vals_owt:
        dfl_owt = 'surf'      # file type 2 is at surface
        dfl_eno = 'depth'     # file type 1 is full depth
        inter_eno.append(depth_eno[0])
    elif 'depth' in vals_eno and 'depth' in vals_owt:
        dfl_eno = 'depth'     # file type 1 is full depth
        dfl_owt = 'depth'     # file type 2 is full depth
    else:  # both are surface or 2d
        dfl_eno = 'surf'     # file type 1 is at surface
        dfl_owt = 'surf'     # file type 2 is at surface
    return inter_eno, inter_owt, dfl_eno, dfl_owt


# Outputting the common variables in order, plus whether full depth or surface/2D
def parcommon_OSU(eno, owt, foc):
    eno_nc = Dataset(os.path.join(eno['foldpath'].values[0], eno['filNCnam'].values[0]), 'r')
    owt_h5 = h5py.File(os.path.join(owt['foldpath'].values[0], owt['filnam'].values[0]), 'r')
    inter_eno, vals_eno = keysearch(eno_nc.variables, foc)
    inter_owt, vals_owt = keysearch(owt_h5, foc)
    commonvars = [value for value in vals_eno if value in vals_owt]
    indices_eno = [vals_eno.index(x) for x in commonvars]
    indices_owt = [vals_owt.index(x) for x in commonvars]
    inter_eno = [inter_eno[i] for i in indices_eno]
    inter_owt = [inter_owt[i] for i in indices_owt]
    if 'depth' not in vals_eno and 'depth' in vals_owt:
        dfl_eno = 'surf'      # file type 1 is at surface
        dfl_owt = 'depth'     # file type 2 is full depth
    elif 'depth' in vals_eno and 'depth' not in vals_owt:
        dfl_owt = 'surf'      # file type 2 is at surface
        dfl_eno = 'depth'     # file type 1 is full depth
    elif 'depth' in vals_eno and 'depth' in vals_owt:
        dfl_eno = 'depth'     # file type 1 is full depth
        dfl_owt = 'depth'     # file type 2 is full depth
    else:  # both are surface or 2d
        dfl_eno = 'surf'     # file type 1 is at surface
        dfl_owt = 'surf'     # file type 2 is at surface
    return inter_eno, inter_owt, dfl_eno, dfl_owt


# If using for establishing grid indices, produces lat, lon, sample data etc.
# If using in loop, just extract data
def varhandling(qpath, vlist, depq, dpos, dep_other, isloop):
    fdat = Dataset(os.path.join(qpath['foldpath'].values[0], qpath['filNCnam'].values[0]), 'r')
    if isloop == True:
        if depq == dep_other:  # if no == 1 and depQ == 's1' or no == 2 and depQ == 's2' or depQ == 'same':
            dat = fdat.variables[vlist[0]]
        else:
            dat = surfext(fdat.variables[vlist[0]], dpos)
        return dat
    else:
        if depq == dep_other:  # if no == 1 and depQ == 's1' or no == 2 and depQ == 's2' or depQ == 'same':
            dat = fdat.variables[vlist[0]]
        else:
            dat = surfext(fdat.variables[vlist[0]], dpos)
        datlat = fdat.variables[vlist[1]][:]
        datlon = fdat.variables[vlist[2]][:]
        if len(vlist) == 4:
            dep = fdat.variables[vlist[3]]
        else:
            dep = numpy.empty((0, 0), dtype=numpy.int8)

        if datlat.shape == datlon.shape:
            pass
        else:
            latpos = dat.shape.index(datlat.shape[0])  # Find which dimension of data is latitude
            lonpos = dat.shape.index(datlon.shape[0])  # Find which dimension of data is latitude
            datlat = np.matlib.repmat(datlat, dat.shape[lonpos], 1)
            if datlat.shape[0] != dat.shape[latpos]:
                datlat = np.transpose(datlat)
            datlon = np.matlib.repmat(datlon, dat.shape[latpos], 1)
            if datlat.shape != datlon.shape:
                datlon = np.transpose(datlon)

        return datlat, datlon, dat, dep


# If using for establishing grid indices, produces lat, lon, sample data etc.
# If using in loop, just extract data
def varhandling_h5(qpath, vlist, depq, dpos, dep_other, isloop):
    fdat = h5py.File(os.path.join(qpath['foldpath'].values[0], qpath['filnam'].values[0]), 'r')
    if isloop == True:
        if depq == dep_other:  # if no == 1 and depQ == 's1' or no == 2 and depQ == 's2' or depQ == 'same':
            dat = fdat[vlist[0]]
        else:
            dat = surfext(fdat[vlist[0]], dpos)
        return dat
    else:
        if depq == dep_other:  # if no == 1 and depQ == 's1' or no == 2 and depQ == 's2' or depQ == 'same':
            dat = fdat[vlist[0]]
        else:
            dat = surfext(fdat[vlist[0]], dpos)
        latsize = min(dat.shape)
        datlat = np.arange((90-(1/(latsize/90))), (-90+(1/(latsize/90))), ((latsize*2)/90))
        lonsize = max(dat.shape)
        datlon = np.arange((90-(1/(lonsize/180))), (-90+(1/(lonsize/180))), ((lonsize*2)/180))
        if len(vlist) == 4:
            dep = fdat[vlist[3]]
        else:
            dep = numpy.empty((0, 0), dtype=numpy.int8)
        if datlat.shape == datlon.shape:
            pass
        else:
            datlat = np.matlib.repmat(datlat, lonsize, 1)
            if datlat.shape[0] != dat.shape[0]:
                datlat = np.transpose(datlat)
            datlon = np.matlib.repmat(datlon, latsize, 1)
            if datlat.shape != datlon.shape:
                datlon = np.transpose(datlon)

        return datlat, datlon, dat, dep


# get grid resolution
def gres(lati, longi):
    latr = max(abs(lati[1, 1] - lati[2, 1]), abs(lati[1, 1] - lati[1, 2]))
    lonr = max(abs(longi[1, 1] - longi[2, 1]), abs(longi[1, 1] - longi[1, 2]))
    return latr, lonr


# establish the grid indices for averaging fine data to coarse resolution
def gridx(lt1, ln1, lt2, ln2, dt2, dloca, tloca, ltr2, lnr2, dattyp2, picklef):
    pathstr = str(pathlib.Path(__file__).parent.absolute())
    if os.path.exists(pathstr + slash + picklef):
        with open(pathstr + slash + picklef, 'rb') as f:
            crs_fin_grd = pickle.load(f)
            lt2 = pickle.load(f)
            ltd2 = pickle.load(f)
            ltid2 = pickle.load(f)
            ln2 = pickle.load(f)
            lnd2 = pickle.load(f)
            lnid2 = pickle.load(f)
            trim = pickle.load(f)
        return crs_fin_grd, lt2, ltd2, ltid2, ln2, lnd2, lnid2, trim

    else:
        ltdelta = 0.5 * ltr2
        ltmin = lt2 - ltdelta
        ltmax = lt2 + ltdelta

        lndelta = 0.5 * lnr2
        lnmin = ln2 - lndelta
        lnmax = ln2 + lndelta

        lt2, ltd2, ltid2, ln2, lnd2, lnid2, trim = IBItrim(lt2, ln2, dattyp2)

        crs_fin_grd = np.empty((lt2.shape[0], ln2.shape[1]), dtype=object)

        if len(dt2.shape) > 2:
            if not dloca:
                if dt2.mask:
                    dt2 = np.ma.masked_where(np.abs(dt2) > 1e10, dt2)
                    dt2 = dt2.sum(axis=tloca)
            else:
                dt2 = surfext(dt2, dloca)

        if ltd2 == 'r' and lnd2 == 'c':
            dt2 = dt2[ltid2, :][:, lnid2]
            ltmin = ltmin[ltid2, :][:, lnid2]
            ltmax = ltmax[ltid2, :][:, lnid2]
            lnmin = lnmin[ltid2, :][:, lnid2]
            lnmax = lnmax[ltid2, :][:, lnid2]
        else:
            dt2 = dt2[lnid2, :][:, ltid2]
            ltmin = ltmin[lnid2, :][:, ltid2]
            ltmax = ltmax[lnid2, :][:, ltid2]
            lnmin = lnmin[lnid2, :][:, ltid2]
            lnmax = lnmax[lnid2, :][:, ltid2]

        for lts in range(0, lt2.shape[0]):
            sta = time.time()
            for lns in range(0, ln2.shape[1]):
                if not dt2.mask[lts, lns]:
                    crs_fin_grd[lts, lns] = np.argwhere(((lt1 >= ltmin[lts, lns]) &
                                                         (lt1 <= ltmax[lts, lns])) &
                                                        ((ln1 >= lnmin[lts, lns]) &
                                                         (ln1 <= lnmax[lts, lns])))
            dne = time.time()
            print(dne-sta)

        with open(pathstr + slash + picklef, 'wb') as f:
            pickle.dump(crs_fin_grd, f)
            pickle.dump(lt2, f)
            pickle.dump(ltd2, f)
            pickle.dump(ltid2, f)
            pickle.dump(ln2, f)
            pickle.dump(lnd2, f)
            pickle.dump(lnid2, f)
            pickle.dump(trim, f)

        return crs_fin_grd, lt2, ltd2, ltid2, ln2, lnd2, lnid2, trim


# Establish from filepath whether the main parameter has a time or depth dimension, and whether either could be
# squeezed out
# Also establish which dimension is latitude and longitude, for later indexing
def ncquery(bath, ncvars, list_t, typ_dat):  # path (rhymes with bath), list of variables, list of time variables
    dat = Dataset(os.path.join(bath['foldpath'].values[0], bath['filNCnam'].values[0]), 'r')
    v_nar = ncvars
    t_nar, _ = keysearch(dat.variables, list_t)
    mvar = np.empty(len(v_nar), dtype=object)
    dsq = False
    tsq = False
    for v in range(0, len(v_nar)):
        mvar[v] = dat.variables[v_nar[v]].shape  # Find shape of each variable
    if len(v_nar) == 3:  # file has no depth
        deppos = []
    else:
        deppos = mvar[0].index(mvar[3][0])  # Find which dimension of data is depth

    if mvar[1] == mvar[2]:
        lat = dat.variables[v_nar[1]][:]
        lon = dat.variables[v_nar[2]][:]
        _, latdim, _, _, londim, _, _ = IBItrim(lat, lon, typ_dat)
        if latdim == 'r' and londim == 'c':
            latpos = mvar[0].index(mvar[1][0])
            lonpos = latpos + 1
        else:
            lonpos = mvar[0].index(mvar[1][0])
            latpos = lonpos + 1
    else:
        latpos = mvar[0].index(mvar[1][0])  # Find which dimension of data is latitude
        lonpos = mvar[0].index(mvar[2][0])  # Find which dimension of data is latitude
        if len(v_nar) == 3:  # file has no depth
            if len(mvar[0]) == 3:
                deppos = list(range(len(mvar[0])))
                deppos.remove(latpos)
                deppos.remove(lonpos)
                deppos = deppos[0]
    if not t_nar:
        timpos = []
    else:
        tvar = dat.variables[t_nar[0]].shape
        timpos = mvar[0].index(tvar[0]) # find position of time in data shape
        if len(v_nar) == 4 and tvar[0] == mvar[3]:
            # file has 1 depth, 1 time and both dimensions can be squeezed out
            dsq = True
            tsq = True
        elif tvar[0] == 1:
            tsq = True
         # shape of the main variable
              # location of depth dimension
                      # squeeze depth (Y/N)
                           # location of time dimension
                                   # squeeze time (Y/N)
    if timpos == deppos:
        deppos = []
    return mvar[0], deppos, dsq, timpos, tsq, latpos, lonpos


# Establish from filepath whether the main parameter has a time or depth dimension, and whether either could be
# squeezed out
# Also establish which dimension is latitude and longitude, for later indexing
def h5query(bath, ncvars, list_t, typ_dat): # path (rhymes with bath), list of variables, list of time variables
    dat = h5py.File(os.path.join(bath['foldpath'].values[0], bath['filnam'].values[0]), 'r')
    v_nar = ncvars
    t_nar = keysearch(dat, list_t)
    mvar = np.empty(len(v_nar), dtype=object)
    dsq = False
    tsq = False
    for v in range(0, len(v_nar)):
        mvar[v] = dat[v_nar[v]].shape  # Find shape of each variable
    if len(v_nar) == 3:  # file has no depth
        deppos = []
    else:
        deppos = mvar[0].index(mvar[3][0])  # Find which dimension of data is depth

    if mvar[1] == mvar[2]:
        lat = dat[v_nar[1]][:]
        lon = dat[v_nar[2]][:]
        _, latdim, _, _, londim, _, _ = IBItrim(lat, lon, typ_dat)
        if latdim == 'r' and londim == 'c':
            latpos = mvar[0].index(mvar[1][0])
            lonpos = latpos + 1
        else:
            lonpos = mvar[0].index(mvar[1][0])
            latpos = lonpos + 1
    else:
        latpos = mvar[0].index(mvar[1][0])  # Find which dimension of data is latitude
        lonpos = mvar[0].index(mvar[2][0])  # Find which dimension of data is latitude
        if len(v_nar) == 3:  # file has no depth
            if len(mvar[0]) == 3:
                deppos = list(range(len(mvar[0])))
                deppos.remove(latpos)
                deppos.remove(lonpos)
                deppos = deppos[0]
    if not t_nar:
        timpos = []
    else:
        tvar = dat[t_nar[0]].shape
        timpos = mvar[0].index(tvar[0]) # find position of time in data shape
        if len(v_nar) == 4 and tvar[0] == mvar[3]:
            # file has 1 depth, 1 time and both dimensions can be squeezed out
            dsq = True
            tsq = True
        elif tvar[0] == 1:
            tsq = True
         # shape of the main variable
              # location of depth dimension
                      # squeeze depth (Y/N)
                           # location of time dimension
                                   # squeeze time (Y/N)
    if timpos == deppos:
        deppos = []
    return mvar[0], deppos, dsq, timpos, tsq, latpos, lonpos


# Restrict single file to dates for which
def shorten(htap, firstd, lastd, timelist, varlist, dim_t):
    fine = Dataset(os.path.join(htap['foldpath'].values[0], htap['filNCnam'].values[0]), 'r')
    tname, _ = keysearch(fine.variables, timelist)
    vnames = varlist
    tQ = fine.variables[tname[0]]
    t_unit = fine.variables[tname[0]].units
    try:
        t_cal = fine.variables[tname[0]].calendar
    except:
        t_cal = "gregorian"
    fdates = num2date(tQ, units=t_unit, calendar=t_cal)  # Dates
    firstd2 = fdates[0]
    firstd2 = firstd2.replace(day=firstd.day)
    firstd2 = firstd2.replace(month=firstd.month)
    firstd2 = firstd2.replace(year=firstd.year)
    lastd2 = fdates[-1]
    lastd2 = lastd2.replace(day=lastd.day)
    lastd2 = lastd2.replace(month=lastd.month)
    lastd2 = lastd2.replace(year=lastd.year)
    firstd = firstd2
    lastd = lastd2
    ID = np.logical_and(firstd <= fdates, fdates <= lastd)
    Data = tslice(fine.variables[vnames[0]], ID, dim_t)
    return Data


# Restrict lat/lon dimensions to the extent of the IBI service domain
def IBItrim(lat, lon, datyp):
    latlim = [26, 56]
    latchkr = abs(lat[0, 0] - lat[1, 0])
    latchkc = abs(lat[0, 1] - lat[0, 0])
    if latchkr > latchkc:
        latdim = 'r'
        if datyp != 'IBI':
            latID = np.logical_and(latlim[0] <= lat[:, 0], lat[:, 0] <= latlim[1])
            latrim = sum(latID) < np.shape(lat)[0]
        else:
            latID = np.ones_like(lat[:, 0], dtype=bool)
            latrim = True
    else:
        latdim = 'c'
        if datyp != 'IBI':
            latID = np.logical_and(latlim[0] <= lat[0, :], lat[0, :] <= latlim[1])
            latrim = sum(latID) < np.shape(lat)[1]
        else:
            latID = np.ones_like(lat[0, :], dtype=bool)
            latrim = True

    lonlim = [-19, 5]
    lonchkr = abs(lon[0, 0] - lon[1, 0])
    lonchkc = abs(lon[0, 1] - lon[0, 0])
    if lonchkr > lonchkc:
        londim = 'r'
        if datyp != 'IBI':
            lonID = np.logical_and(lon[:, 0] >= lonlim[0], lon[:, 0] <= lonlim[1])
            lotrim = sum(lonID) < np.shape(lon)[0]
        else:
            lonID = np.ones_like(lon[:, 0], dtype=bool)
            lotrim = True
    else:
        londim = 'c'
        if datyp != 'IBI':
            lonID = np.logical_and(lon[0, :] >= lonlim[0], lon[0, :] <= lonlim[1])
            lotrim = sum(lonID) < np.shape(lon)[1]
        else:
            lonID = np.ones_like(lon[0, :], dtype=bool)
            lotrim = True
    if latdim == 'r' and londim == 'c':
        lat = lat[latID, :][:, lonID]
        lon = lon[latID, :][:, lonID]
    else:
        lat = lat[lonID, :][:, latID]
        lon = lon[lonID, :][:, latID]
    if lotrim == False and latrim == False:
        trim = False
    else:
        trim = True
    return lat, latdim, latID, lon, londim, lonID, trim


# Use to limit single file to times that are common to the coarse and fine datasets (see shorten)
def tslice(atad, DI, tdim):
    s_array = atad.shape
    if len(s_array) == 4 and tdim == 0:
        atad_o = atad[DI, :, :, :]
    elif len(s_array) == 3 and tdim == 0:
        atad_o = atad[DI, :, :]
    elif len(s_array) == 4 and tdim == 1:
        atad_o = atad[:, DI, :, :]
    else:
        atad_o = atad
    return atad_o


def surfext(atad, d_in_shape):
    s_array = atad.shape
    if len(s_array) == 4 and d_in_shape == 0:
        atad_o = np.squeeze(atad[0, :, :, :])
    elif len(s_array) == 3 and d_in_shape == 0:
        atad_o = np.squeeze(atad[0, :, :])
    elif len(s_array) == 4 and d_in_shape == 1:
        atad_o = np.squeeze(atad[:, 0, :, :])
    else:
        atad_o = atad
    return atad_o


# If data has dimensions to squeeze, squeeze them out, turn off flags for requirement to squeeze (*_squ), and empty the
# position of dimensions that were squeeze (*pos)
def dimsquash(dartray, dpos, d_squ, tpos, t_squ):
    dartray_o = np.squeeze(dartray)
    if d_squ == True:
        d_squ = False
        dpos = []
    if t_squ == True:
        t_squ = False
        tpos = []
    return dartray_o, dpos, d_squ, tpos, t_squ


def tbuild(tflag1, tflag2, date1, date2):
    # no time, as assembling climatology (summing the long dataset to the year and averaging)
    if tflag1 == 'zy' or tflag2 == 'zy':
        tarray = []
        tflag_both = 'zy'
    # both at monthly resolution
    elif tflag1 == 'm' or tflag2 == 'm':
        # times are at monthly resolution from date1 to date2
        tarray = pd.date_range(date1, date2, freq='M')
        tarray = pd.date_range(date1, periods=len(tarray.month.unique())*len(tarray.year.unique()), freq='M')
        tarray = tarray - numpy.timedelta64(15, 'D')
        tflag_both = 'm'
    # both at daily resolution
    else:
        # times are at daily resolution from date1 to date2
        tarray = pd.date_range(date1, date2, freq='D')
        tflag_both = 'd'
    # datetime string array for full
    return tarray, tflag_both


# Output the smaller data array, with latitude and longitude narrowed to IBI extents
def cutter(atad, DI, dim):
    s_array = atad.shape
    if len(s_array) == 4:
        if dim == 0:
            atad_o = atad[DI, :, :, :]
        elif dim == 1:
            atad_o = atad[:, DI, :, :]
        elif dim == 2:
            atad_o = atad[:, :, DI, :]
        elif dim == 3:
            atad_o = atad[:, :, :, DI]
        else:
            atad_o = atad
    elif len(s_array) == 3:
        if dim == 0:
            atad_o = atad[DI, :, :]
        elif dim == 1:
            atad_o = atad[:, DI, :]
        elif dim == 2:
            atad_o = atad[:, :, DI]
        else:
            atad_o = atad
    elif len(s_array) == 2:
        if dim == 0:
            atad_o = atad[DI, :]
        elif dim == 1:
            atad_o = atad[:, DI]
        else:
            atad_o = atad
    else:
        atad_o = atad
    return atad_o


def vertlevels(datain, data_tofill, depin, depout, latin, lonin, dep_limit):
    data_ismask = (data_tofill == data_tofill)
    if isinstance(latin, int):
        for lons in range(0, lonin.shape[1]):
            if dep_limit[0, lons] > 0:
                data_tofill[0, lons] = np.interp(depout, depin[:int(dep_limit[0, lons])],
                                                    datain[:int(dep_limit[0, lons]), 0, lons])
                data_ismask[0, lons] = False
    elif isinstance(lonin, int):
        for lats in range(0, latin.shape[0]):
            if dep_limit[lats, 0] > 0:
                data_tofill[lats, 0] = np.interp(depout, depin[:int(dep_limit[lats, 0])],
                                                    datain[:int(dep_limit[lats, 0]), lats, 0])
                data_ismask[lats, 0] = False
    else:
        for lats in range(0, latin.shape[0]):
            for lons in range(0, lonin.shape[1]):
                if dep_limit[lats, lons] > 0:
                    doutlim = np.argwhere(depout[:] >= depin[dep_limit[lats, lons]])
                    if len(doutlim) > 1:
                        doutlim = int(doutlim[0] - 1)
                    elif len(doutlim) == 1:
                        doutlim = int(doutlim - 1)
                    else:
                        doutlim = int(len(depout))
                    data_tofill[:doutlim, lats, lons] = np.interp(depout[:doutlim], depin[:int(dep_limit[lats, lons])],
                                                                  datain[:int(dep_limit[lats, lons]), lats, lons])
                    data_ismask[:doutlim, lats, lons] = False
    return data_tofill, data_ismask


def aggregate(path, t, tdelta, dchkboth, dpos, param, trim, latpos, latIdx, lonpos, lonIdx):
    if tdelta == 'zy':
        nupath = path[(path['DateStr'].dt.year == t.year)]
    elif tdelta == 'm':
        nupath = path[(path['DateStr'].dt.year == t.year) & (path['DateStr'].dt.month == t.month)]
    else:  # must be daily resolution, therefore aggregating subdaily data
        nupath = path
    for f in range(len(nupath)):
        if f == 0:
            data  = data_read(nupath, f, dchkboth, dpos, param, trim, latpos, latIdx, lonpos, lonIdx)
        else:
            data1 = data_read(nupath, f, dchkboth, dpos, param, trim, latpos, latIdx, lonpos, lonIdx)
            data = np.add(data, data1)
    data = np.divide(data, len(nupath))
    data = np.ma.masked_where(data > 10**4, data)
    return data


def aggregate_h5(path, t, tdelta, dchkboth, dpos, param, trim, latpos, latIdx, lonpos, lonIdx):
    if tdelta == 'zy':
        nupath = path[(path['DateStr'].dt.year == t.year)]
    elif tdelta == 'm':
        nupath = path[(path['DateStr'].dt.year == t.year) & (path['DateStr'].dt.month == t.month)]
    else:  # must be daily resolution, therefore aggregating subdaily data
        nupath = path
    for f in range(len(nupath)):
        if f == 0:
            data  = data_read_h5(nupath, f, dchkboth, dpos, param, trim, latpos, latIdx, lonpos, lonIdx)
        else:
            data1 = data_read_h5(nupath, f, dchkboth, dpos, param, trim, latpos, latIdx, lonpos, lonIdx)
            data = np.add(data, data1)
    data = np.divide(data, len(nupath))
    data = np.ma.masked_where(data > 10**4, data)
    return data


def data_read(htap, ind, dflag, d_loc, parname, trim, latpos, latIdx, lonpos, lonIdx):
    finfil = Dataset(os.path.join(htap['foldpath'].values[ind], htap['filNCnam'].values[ind]), 'r')
    if dflag == 'depth':
        rdata = finfil.variables[parname]
    elif dflag == 'surf' and d_loc == []:
        rdata = finfil.variables[parname]
    else: # dflag == 'surf' and d_loc is:
        rdata = surfext(finfil.variables[parname], d_loc)
    if trim:
        rdata = cutter(rdata, latIdx, latpos)
        rdata = cutter(rdata, lonIdx, lonpos)
    return rdata


def data_read_h5(htap, ind, dflag, d_loc, parname, trim, latpos, latIdx, lonpos, lonIdx):
    finfil = h5py.File(os.path.join(htap['foldpath'].values[ind], htap['filnam'].values[ind]), 'r')
    if dflag == 'depth':
        rdata = finfil[parname]
    elif dflag == 'surf' and d_loc == []:
        rdata = finfil[parname]
    else:  # dflag == 'surf' and d_loc is:
        rdata = surfext(finfil[parname], d_loc)
    if trim:
        rdata = cutter(rdata, latIdx, latpos)
        rdata = cutter(rdata, lonIdx, lonpos)
    return rdata


def blkav(data, data_tofill, lat, lon, datsamp, index, dflag):
    data_ismask = (data_tofill == data_tofill)
    if dflag == 'surf':
        for lats in range(0, lat.shape[0]):
            for lons in range(0, lon.shape[1]):
                try:
                    if not datsamp.mask[lats, lons]:
                        try:
                            if index[lats, lons] is None:
                                data_ismask[0, lats, lons] = True
                            elif index[lats, lons].any():
                                msk_counting = np.ma.count_masked(data[0, index[lats, lons][:, 0],
                                                                       index[lats, lons][:, 1]])
                                if msk_counting == 0:
                                    mndat = np.ma.mean(data[0, index[lats, lons][:, 0], index[lats, lons][:, 1]])
                                    if mndat < 10**4:
                                        data_tofill[0, lats, lons] = mndat
                                        data_ismask[0, lats, lons] = False
                                    else:
                                        data_tofill[0, lats, lons] = 0
                                        data_ismask[0, lats, lons] = True
                                else:
                                    data_ismask[0, lats, lons] = True
                            else:
                                data_ismask[0, lats, lons] = True
                        except:
                            if index[lats, lons] is None:
                                data_ismask[0, lats, lons] = True
                            elif index[lats, lons].any():
                                msk_counting = np.ma.count_masked(data[index[lats, lons][:, 0],
                                                                       index[lats, lons][:, 1]])
                                if msk_counting == 0:
                                    mndat = np.ma.mean(data[index[lats, lons][:, 0], index[lats, lons][:, 1]])
                                    if mndat < 10**4:
                                        data_tofill[0, lats, lons] = mndat
                                        data_ismask[0, lats, lons] = False
                                    else:
                                        data_tofill[0, lats, lons] = 0
                                        data_ismask[0, lats, lons] = True
                                else:
                                    data_ismask[0, lats, lons] = True
                            else:
                                data_ismask[0, lats, lons] = True
                except:
                    try:
                        if index[lats, lons] is None:
                            data_ismask[0, lats, lons] = True
                        elif index[lats, lons].any():
                            msk_counting = np.ma.count_masked(data[0, index[lats, lons][:, 0],
                                                                   index[lats, lons][:, 1]])
                            if msk_counting == 0:
                                mndat = np.ma.mean(data[0, index[lats, lons][:, 0], index[lats, lons][:, 1]])
                                if mndat < 10**4:
                                    data_tofill[0, lats, lons] = mndat
                                    data_ismask[0, lats, lons] = False
                                else:
                                    data_tofill[0, lats, lons] = 0
                                    data_ismask[0, lats, lons] = True
                            else:
                                data_ismask[0, lats, lons] = True
                        else:
                            data_ismask[0, lats, lons] = True
                    except:
                        if index[lats, lons] is None:
                            data_ismask[0, lats, lons] = True
                        elif index[lats, lons].any():
                            msk_counting = np.ma.count_masked(data[index[lats, lons][:, 0], index[lats, lons][:, 1]])
                            if msk_counting == 0:
                                mndat = np.ma.mean(data[index[lats, lons][:, 0], index[lats, lons][:, 1]])
                                if mndat < 10**4:
                                    data_tofill[0, lats, lons] = mndat
                                    data_ismask[0, lats, lons] = False
                                else:
                                    data_tofill[0, lats, lons] = 0
                                    data_ismask[0, lats, lons] = True
                            else:
                                data_ismask[0, lats, lons] = True
                        else:
                            data_ismask[0, lats, lons] = True
    else:  # dflag == 'depth'
        for d in range(0, data.shape[0]):
            for lats in range(0, lat.shape[0]):
                for lons in range(0, lon.shape[1]):
                    try:
                        if not datsamp.mask[lats, lons]:
                            if index[lats, lons] is None:
                                data_ismask[0, lats, lons] = True
                            elif index[lats, lons].any():
                                if np.ma.count_masked(data[d, index[lats, lons][:, 0],
                                                           index[lats, lons][:, 1]]) == 0:
                                    mndat = np.ma.mean(data[d, index[lats, lons][:, 0], index[lats, lons][:, 1]])
                                    if mndat < 10**4:
                                        data_tofill[d, lats, lons] = mndat
                                        data_ismask[d, lats, lons] = False
                                    else:
                                        data_tofill[d, lats, lons] = 0
                                        data_ismask[d, lats, lons] = True
                                else:
                                    data_ismask[d, lats, lons] = True
                            else:
                                data_ismask[d, lats, lons] = True
                    except:
                        if index[lats, lons] is None:
                            data_ismask[0, lats, lons] = True
                        elif index[lats, lons].any():
                            if np.ma.count_masked(data[d, index[lats, lons][:, 0], index[lats, lons][:, 1]]) == 0:
                                mndat = np.ma.mean(data[d, index[lats, lons][:, 0], index[lats, lons][:, 1]])
                                if mndat < 10**4:
                                    data_tofill[d, lats, lons] = mndat
                                    data_ismask[d, lats, lons] = False
                                else:
                                    data_tofill[d, lats, lons] = 0
                                    data_ismask[d, lats, lons] = True
                            else:
                                data_ismask[d, lats, lons] = True
                        else:
                            data_ismask[d, lats, lons] = True
    return data_tofill, data_ismask


def BLAV(path1, wild1, datyp1, path2, wild2, datyp2, unit_coeff, Picklename):
    # datyp1, datyp2 can be 'IBI', 'data', 'dclm', 'mclm' or 'yclm'
    if datyp1 == 'IBI':
        pathibi = path1
        wildibi = wild1
        pathdat = path2
        wilddat = wild2
        typedat = datyp2
    else:  # datyp2 == 'IBI':
        pathibi = path2
        wildibi = wild2
        pathdat = path1
        wilddat = wild1
        typedat = datyp1
    # list all files in folder using wildcards
    listibi, pathibi, strdibi, enddibi, dtibi, ofibi = freview(pathibi, wildibi, 'IBI',   TimeDict)
    listdat, pathdat, strddat, endddat, dtdat, ofdat = freview(pathdat, wilddat, typedat, TimeDict)
    if not strdibi:
        strdibi = strddat
    if not strddat:
        strddat = strdibi
    if not enddibi:
        enddibi = endddat
    if not endddat:
        endddat = enddibi
    strddat2 = strdibi
    strddat2 = strddat2.replace(day=strddat.day)
    strddat2 = strddat2.replace(month=strddat.month)
    strddat2 = strddat2.replace(year=strddat.year)
    endddat2 = enddibi
    endddat2 = endddat2.replace(day=endddat.day)
    endddat2 = endddat2.replace(month=endddat.month)
    endddat2 = endddat2.replace(year=endddat.year)
    strddat = strddat2
    endddat = endddat2

    # having list of files and corresponding dates, determine max of the start dates and minima of the end dates (across
    # both datasets)
    strd = max(strdibi, strddat)
    endd = min(enddibi, endddat)

    # limit processing to files between the two dates
    # From file strings (get listing using wildcards), and ensure that processing is only
    # applied to overlapping period i.e. max of start times for 2 datasets, min of end times for two datasets.
    if ofibi == False:
        pathibi = pathibi[(pathibi['DateStr'] >= strd) & (pathibi['DateStr'] <= endd)]
    if ofdat == False and typedat == 'data':
        pathdat = pathdat[(pathdat['DateStr'] >= strd) & (pathdat['DateStr'] <= endd)]

    # get common variables from both params
    interibi, interdat, dchkibi, dchkdat = parcommon(pathibi, pathdat, VarDict)

    # establish shape of data, which dimension is depth, whether it can be squeezed, which dim is time,
    # squeezed, and which dimension is latitude and longitude
    vardims_ibi, dpos_ibi, dsq_ibi, tpos_ibi, tsq_ibi, lapos_ibi, lopos_ibi = \
        ncquery(pathibi, interibi, TimeDict, 'IBI')
    vardims_dat, dpos_dat, dsq_dat, tpos_dat, tsq_dat, lapos_dat, lopos_dat = \
        ncquery(pathdat, interdat, TimeDict, typedat)


    # read latitude and longitude of the first IBI file in list
    latibi, lonibi, datibi, cdepibi = varhandling(pathibi, interibi, dchkibi, dpos_ibi, dchkdat, False)

    # read latitude and longitude of the first data file in each list
    latdat, londat, datdat, cdepdat = varhandling(pathdat, interdat, dchkdat, dpos_dat, dchkibi, False)

    # find fine grid and coarse grid resolution
    latibir, lonibir = gres(latibi, lonibi)
    latdatr, londatr = gres(latdat, londat)
    # get area of each grid cell (decimal degree squared)
    gibi_2 = latibir * lonibir
    gdat_2 = latdatr * londatr

    # establish grid indices for averaging finer spatial data to coarser grid structure
    if gibi_2 < gdat_2:  # IBI grid is finer than data
        finidx, latdat, latdatdim, latdatID, londat, londatdim, londatID, trimflag = gridx(latibi, lonibi, latdat,
                                                                                           londat, datdat, dpos_dat,
                                                                                           tpos_dat, latdatr, londatr,
                                                                                           typedat, Picklename)
        # Latitude and longitude taken from coarser grid, and restricted to IBI domain
        latlen = sum(latdatID)
        lonlen = sum(londatID)
        if trimflag == True:
            trim_ID = 'data'
            lapos = lapos_dat
            lopos = lopos_dat
            latID  = latdatID
            lonID  = londatID
        else:
            trim_ID = 'none'
            lapos = []
            lopos = []
            latID  = []
            lonID  = []
        coarse = 'data'
    else:  # data/climatology is finer than IBI
        finidx, latibi, latibidim, latibiID, lonibi, lonibidim, lonibiID, trimflag = gridx(latdat, londat, latibi,
                                                                                           lonibi, datibi, dpos_ibi,
                                                                                           tpos_ibi, latibir, lonibir,
                                                                                           'IBI', Picklename)
        # Latitude and longitude taken from coarser grid, and restricted to IBI domain
        latlen = sum(latibiID)
        lonlen = sum(lonibiID)
        if trimflag == True:
            trim_ID = 'IBI'
            lapos = lapos_ibi
            lopos = lopos_ibi
            latID  = latibiID
            lonID  = lonibiID
        else:
            trim_ID = 'none'
            lapos = []
            lopos = []
            latID  = []
            lonID  = []
        coarse = 'IBI'

    # Time dimension - lowest of two frequencies if different
    times, dtboth = tbuild(dtibi, dtdat, strd, endd)
    tlen = len(times)
    # Depth dimension - either surface or full depth, depending on data characteristics or manual input
    if dchkibi and dchkdat == 'depth':
        dlen_ibi = vardims_ibi[dpos_ibi]
        dlen_dat = vardims_dat[dpos_dat]
        dchkboth = 'depth'
    else:  # either one or both data are surface only, which means only surface can be compared
        dlen_ibi = 1  # just taking surface, can be "squashed out" afterwards
        dlen_dat = 1
        dchkboth = 'surf'

    if 'clm' in typedat:
        pathdat = pd.concat([pathdat] * int((tlen/len(pathdat))), ignore_index=True)

    # Set up empty array for coarse and fine grid, such that both will be equivalent in size, shape, etc.
    # What size will the common array be? lat, lon, depth, time?
    datV = np.zeros((tlen, dlen_dat, latlen, lonlen))
    datV_msk = (datV == datV)
    datM = np.zeros((tlen, dlen_dat, latlen, lonlen))
    datM_msk = (datM == datM)
    tpos = 0
    dpos = 1
    latpos = 2
    lonpos = 3
    # if the coarser grid is non-IBI (data or climatology), restrict latitudes and longitudes (i) before indexing,
    # and (ii) during the construction of the two final arrays.

    if trim_ID == 'IBI':
        trimIBI = True
    else:
        trimIBI = False

    if ofibi == True:
        datM = shorten(pathibi, strd, endd, TimeDict, VarDict, tpos_ibi)
        datM = cutter(datM, latibiID, lapos_ibi)
        datM = cutter(datM, lonibiID, lopos_ibi)
    else:  # multiple IBI files, need processing
        if coarse == 'IBI':
            if dtibi != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate(pathibi, times[ff], dtboth, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos,
                                       latID, lopos, lonID)
                    datagg = np.ma.masked_array.squeeze(datagg)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(datagg.mask)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm, maskchk = vertlevels(datagg, datfrm, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm = np.ma.masked_where(maskchk == True, datfrm)
                        datM[ff, :, :, :] = datfrm
                        datM_msk[ff, :, :, :] = maskchk
                    else:
                        datM[ff, :, :, :] = datagg
                        datM_msk[ff, :, :, :] = datagg.mask
            else:  # if not, just read IBI data, stack into array
                for ff in range(tlen):
                    data = data_read(pathibi, ff, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos, latID, lopos, lonID)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(data.mask)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm, maskchk = vertlevels(data, datfrm, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm = np.ma.masked_where(maskchk == True, datfrm)
                        datM[ff, :, :, :] = datfrm
                        datM_msk[ff, :, :, :] = maskchk
                    else:
                        datM[ff, :, :, :] = data
                        datM_msk[ff, :, :, :] = data.mask
        else:  # data is coarser than IBI and so, IBI needs block averaging
            if dtibi != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate(pathibi, times[ff], dtboth, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos,
                                       latID, lopos, lonID)
                    datagg = np.ma.masked_array.squeeze(datagg)
                    datfrm1 = np.zeros((dlen_ibi, latlen, lonlen))
                    if ofdat == True:
                        datdatII = shorten(pathdat, strd, endd, TimeDict, interdat, tpos_dat)
                        datdatII = cutter(datdatII, latdatID, lapos_dat)
                        datdatII = cutter(datdatII, londatID, lopos_dat)
                        datdatII = cutter(datdatII, ff, tpos_dat)
                    elif dpos_dat != []:
                        datdatII = surfext(datdat, dpos_dat)
                        if latdatdim == 'r' and londatdim == 'c':
                            datdatII = datdatII[latdatID, :][:, londatID]
                        else:
                            datdatII = datdatII[londatID, :][:, latdatID]
                    else: # dpos_dat == []:
                        if latdatdim == 'r' and londatdim == 'c':
                            datdatII = datdat[latdatID, :][:, londatID]
                        else:
                            datdatII = datdat[londatID, :][:, latdatID]
                    datfrm1, maskchk = blkav(datagg, datfrm1, latdat, londat, datdatII, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(maskchk)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm2 = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm2, maskchkII = vertlevels(datfrm1, datfrm2, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm2 = np.ma.masked_where(maskchkII == True, datfrm2)
                        datM[ff, :, :, :] = datfrm2
                        datM_msk[ff, :, :, :] = maskchkII
                    else:
                        datM[ff, :, :, :] = datfrm1
                        datM_msk[ff, :, :, :] = maskchk
            else:  # if not, just read IBI data, stack into array
                for ff in range(0, tlen):
                    data = data_read(pathibi, ff, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos, latID, lopos, lonID)
                    try:
                        data = np.ma.masked_array.squeeze(data)
                    except:
                        data = np.squeeze(data)
                    datfrm1 = np.zeros((dlen_ibi, latlen, lonlen))
                    datdatII = surfext(datdat, dpos_dat)
                    if latdatdim == 'r' and londatdim == 'c':
                        datdatII = datdatII[latdatID, :][:, londatID]
                    else:
                        datdatII = datdatII[londatID, :][:, latdatID]
                    datfrm1, maskchk = blkav(data, datfrm1, latdat, londat, datdatII, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(maskchk)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm2 = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm2, maskchkII = vertlevels(datfrm1, datfrm2, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm2 = np.ma.masked_where(maskchkII == True, datfrm2)
                        datM[ff, :, :, :] = datfrm2
                        datM_msk[ff, :, :, :] = maskchkII
                    else:
                        datM[ff, :, :, :] = datfrm1
                        datM_msk[ff, :, :, :] = maskchk

    if trim_ID == 'data':
        trimdat = True
    else:
        trimdat = False

    if ofdat == True:
        datV = shorten(pathdat, strd, endd, TimeDict, interdat, tpos_dat)
        datV = cutter(datV, latdatID, lapos_dat)
        datV = cutter(datV, londatID, lopos_dat)
    else:  # multiple data/climatology files, need processing
        if coarse != 'IBI':  # data is coarse, IBI needs block averaging
            if dtdat != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate(pathdat, times[ff], dtboth, dchkboth, dpos_dat, interdat[0], trimdat, lapos,
                                       latID, lopos, lonID)
                    datagg = np.ma.masked_array.squeeze(datagg)
                    datV[ff, :, :, :] = datagg
                    datV_msk[ff, :, :, :] = datagg.mask
            else:  # if not, just read data, stack into array
                for ff in range(tlen):
                    data = data_read(pathdat, ff, dchkboth, dpos_dat, interdat[0], trimdat, lapos, latID, lopos, lonID)
                    datV[ff, :, :, :] = data
                    datV_msk[ff, :, :, :] = data.mask
        else:  # IBI is coarser than data/climatology and so, data/climatology needs block averaging
            if dtdat != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate(pathdat, times[ff], dtboth, dchkboth, dpos_dat, interdat[0], trimdat, lapos,
                                       latID, lopos, lonID)
                    datfrm1 = np.zeros((dlen_dat, latlen, lonlen))
                    datfrm1, maskchk = blkav(datagg, datfrm1, latdat, londat, datdat, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    datV[ff, :, :, :] = datfrm1
                    datV_msk[ff, :, :, :] = maskchk
            else:  # if not, just read data, stack into array
                for ff in range(tlen):
                    sta = time.time()
                    data = data_read(pathdat, ff, dchkboth, dpos_dat, interdat[0], trimdat, lapos, latID, lopos, lonID)
                    datfrm1 = np.zeros((dlen_dat, latlen, lonlen))
                    datfrm1, maskchk = blkav(data, datfrm1, latibi, lonibi,  datibi, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    datV[ff, :, :, :] = datfrm1
                    datV_msk[ff, :, :, :] = maskchk
                    end = time.time()
                    print(end-sta)
    datV = datV * unit_coeff
    datV_datM_msk = datM_msk | datV_msk
    datM = np.ma.masked_where(datV_datM_msk == True, datM)
    datV = np.ma.masked_where(datV_datM_msk == True, datV)
    return datV, datM, times, tpos, cdepdat, dpos, dchkboth, latdat, latpos, londat, lonpos



def BLAV_OSU(path1, wild1, datyp1, path2, wild2, datyp2, Picklename):
    # datyp1, datyp2 can be 'IBI', 'data', 'dclm', 'mclm' or 'yclm'
    if datyp1 == 'IBI':
        pathibi = path1
        wildibi = wild1
        pathdat = path2
        wilddat = wild2
        typedat = datyp2
    else:  # datyp2 == 'IBI':
        pathibi = path2
        wildibi = wild2
        pathdat = path1
        wilddat = wild1
        typedat = datyp1
    # list all files in folder using wildcards
    listibi, pathibi, strdibi, enddibi, dtibi, ofibi = freview(pathibi, wildibi, 'IBI',   TimeDict)
    listdat, pathdat, strddat, endddat, dtdat, ofdat = freview_OSU(pathdat, wilddat, typedat)
    if not strdibi:
        strdibi = strddat
    if not strddat:
        strddat = strdibi
    if not enddibi:
        enddibi = endddat
    if not endddat:
        endddat = enddibi

    # having list of files and corresponding dates, determine max of the start dates and minima of the end dates (across
    # both datasets)
    strd = max(strdibi, strddat)
    endd = min(enddibi, endddat)

    # limit processing to files between the two dates
    # From file strings (get listing using wildcards), and ensure that processing is only
    # applied to overlapping period i.e. max of start times for 2 datasets, min of end times for two datasets.
    if ofibi == False:
        pathibi = pathibi[(pathibi['DateStr'] >= strd) & (pathibi['DateStr'] <= endd)]
    if ofdat == False and typedat == 'data':
        pathdat = pathdat[(pathdat['DateStr'] >= strd) & (pathdat['DateStr'] <= endd)]

    # get common variables from both params
    interibi, interdat, dchkibi, dchkdat = parcommon_OSU(pathibi, pathdat, VarDict)

    # establish shape of data, which dimension is depth, whether it can be squeezed, which dim is time, whether it can be
    # squeezed, and which dimension is latitude and longitude
    vardims_ibi, dpos_ibi, dsq_ibi, tpos_ibi, tsq_ibi, lapos_ibi, lopos_ibi = ncquery(pathibi, VarDict, TimeDict)
    vardims_dat, dpos_dat, dsq_dat, tpos_dat, tsq_dat, lapos_dat, lopos_dat = h5query(pathdat, VarDict, TimeDict)


    # read latitude and longitude of the first IBI file in list
    latibi, lonibi, datibi, cdepibi = varhandling(pathibi, interibi, dchkibi, dpos_ibi, dchkdat, False)

    # read latitude and longitude of the first data file in each list
    latdat, londat, datdat, cdepdat = varhandling_h5(pathdat, interdat, dchkdat, dpos_dat, dchkibi, False)

    # find fine grid and coarse grid resolution
    latibir, lonibir = gres(latibi, lonibi)
    latdatr, londatr = gres(latdat, londat)
    # get area of each grid cell (decimal degree squared)
    gibi_2 = latibir * lonibir
    gdat_2 = latdatr * londatr

    # establish grid indices for averaging finer spatial data to coarser grid structure
    if gibi_2 < gdat_2:  # IBI grid is finer than data
        finidx, latdat, latdatdim, latdatID, londat, londatdim, londatID, trimflag = gridx(latibi, lonibi, latdat,
                                                                                           londat, datdat, dpos_dat,
                                                                                           tpos_dat, latdatr, londatr,
                                                                                           Picklename)
        # Latitude and longitude taken from coarser grid, and restricted to IBI domain
        latlen = sum(latdatID)
        lonlen = sum(londatID)
        if trimflag == True:
            trim_ID = 'data'
            lapos = lapos_dat
            lopos = lopos_dat
            latID  = latdatID
            lonID  = londatID
        else:
            trim_ID = 'none'
            lapos = []
            lopos = []
            latID  = []
            lonID  = []
        coarse = 'data'
    else:  # data/climatology is finer than IBI
        finidx, latibi, latibidim, latibiID, lonibi, lonibidim, lonibiID, trimflag = gridx(latdat, londat, latibi,
                                                                                           lonibi, datibi, dpos_ibi,
                                                                                           tpos_ibi, latibir, lonibir,
                                                                                           Picklename)
        # Latitude and longitude taken from coarser grid, and restricted to IBI domain
        latlen = sum(latibiID)
        lonlen = sum(lonibiID)
        if trimflag == True:
            trim_ID = 'IBI'
            lapos = lapos_ibi
            lopos = lopos_ibi
            latID  = latibiID
            lonID  = lonibiID
        else:
            trim_ID = 'none'
            lapos = []
            lopos = []
            latID  = []
            lonID  = []
        coarse = 'IBI'

    # Time dimension - lowest of two frequencies if different
    times, dtboth = tbuild(dtibi, dtdat, strd, endd)
    tlen = len(times)
    # Depth dimension - either surface or full depth, depending on data characteristics or manual input
    if dchkibi and dchkdat == 'depth':
        dlen_ibi = vardims_ibi[dpos_ibi]
        dlen_dat = vardims_dat[dpos_dat]
        dchkboth = 'depth'
    else:  # either one or both data are surface only, which means only surface can be compared
        dlen_ibi = 1  # just taking surface, can be "squashed out" afterwards
        dlen_dat = 1
        dchkboth = 'surf'

    if 'clm' in typedat:
        pathdat = pd.concat([pathdat] * int((tlen/len(pathdat))), ignore_index=True)

    # Set up empty array for coarse and fine grid, such that both will be equivalent in size, shape, etc.
    # What size will the common array be? lat, lon, depth, time?
    datV = np.zeros((tlen, dlen_dat, latlen, lonlen))
    datV_msk = (datV == datV)
    datM = np.zeros((tlen, dlen_dat, latlen, lonlen))
    datM_msk = (datM == datM)
    tpos = 0
    dpos = 1
    latpos = 2
    lonpos = 3
    # if the coarser grid is non-IBI (data or climatology), restrict latitudes and longitudes (i) before indexing,
    # and (ii) during the construction of the two final arrays.

    if trim_ID == 'IBI':
        trimIBI = True
    else:
        trimIBI = False

    if ofibi == True:
        datM = shorten(pathibi, strd, endd, TimeDict, interibi, tpos_ibi)
        datM = cutter(datM, latibiID, lapos_ibi)
        datM = cutter(datM, lonibiID, lopos_ibi)
    else:  # multiple IBI files, need processing
        if coarse == 'IBI':
            if dtibi != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate(pathibi, times[ff], dtboth, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos,
                                       latID, lopos, lonID)
                    datagg = np.ma.masked_array.squeeze(datagg)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(datagg.mask)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm, maskchk = vertlevels(datagg, datfrm, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm = np.ma.masked_where(maskchk == True, datfrm)
                        datM[ff, :, :, :] = datfrm
                        datM_msk[ff, :, :, :] = maskchk
                    else:
                        datM[ff, :, :, :] = datagg
                        datM_msk[ff, :, :, :] = datagg.mask
            else:  # if not, just read IBI data, stack into array
                for ff in range(tlen):
                    data = data_read(pathibi, ff, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos, latID, lopos, lonID)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(data.mask)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm, maskchk = vertlevels(data, datfrm, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm = np.ma.masked_where(maskchk == True, datfrm)
                        datM[ff, :, :, :] = datfrm
                        datM_msk[ff, :, :, :] = maskchk
                    else:
                        datM[ff, :, :, :] = data
                        datM_msk[ff, :, :, :] = maskchk
        else:  # data is coarser than IBI and so, IBI needs block averaging
            if dtibi != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate(pathibi, times[ff], dtboth, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos,
                                       latID, lopos, lonID)
                    datagg = np.ma.masked_array.squeeze(datagg)
                    datfrm1 = np.zeros((dlen_ibi, latlen, lonlen))
                    datdatII = surfext(datdat, dpos_dat)
                    if latdatdim == 'r' and londatdim == 'c':
                        datdatII = datdatII[latdatID, :][:, londatID]
                    else:
                        datdatII = datdatII[londatID, :][:, latdatID]
                    datfrm1, maskchk = blkav(datagg, datfrm1, latdat, londat, datdatII, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(maskchk)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm2 = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm2, maskchkII = vertlevels(datfrm1, datfrm2, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm2 = np.ma.masked_where(maskchkII == True, datfrm2)
                        datM[ff, :, :, :] = datfrm2
                        datM_msk[ff, :, :, :] = maskchkII
                    else:
                        datM[ff, :, :, :] = datfrm1
                        datM_msk[ff, :, :, :] = maskchk
            else:  # if not, just read IBI data, stack into array
                for ff in range(0, tlen):
                    data = data_read(pathibi, ff, dchkboth, dpos_ibi, interibi[0], trimIBI, lapos, latID, lopos, lonID)
                    data = np.ma.masked_array.squeeze(data)
                    datfrm1 = np.zeros((dlen_ibi, latlen, lonlen))
                    datdatII = surfext(datdat, dpos_dat)
                    if latdatdim == 'r' and londatdim == 'c':
                        datdatII = datdatII[latdatID, :][:, londatID]
                    else:
                        datdatII = datdatII[londatID, :][:, latdatID]
                    datfrm1, maskchk = blkav(data, datfrm1, latdat, londat, datdatII, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    if dchkboth == 'depth' and dlen_ibi != dlen_dat:
                        d_count = sum(maskchk)
                        d_limit = (dlen_ibi * np.ones_like(d_count))-d_count
                        datfrm2 = np.zeros((dlen_dat, latlen, lonlen))
                        datfrm2, maskchkII = vertlevels(datfrm1, datfrm2, cdepibi, cdepdat, latdat, londat, d_limit)
                        datfrm2 = np.ma.masked_where(maskchkII == True, datfrm2)
                        datM[ff, :, :, :] = datfrm2
                        datM_msk[ff, :, :, :] = maskchkII
                    else:
                        datM[ff, :, :, :] = datfrm1
                        datM_msk[ff, :, :, :] = maskchk

    if trim_ID == 'data':
        trimdat = True
    else:
        trimdat = False

    if ofdat == True:
        datV = shorten(pathdat, strd, endd, TimeDict, interdat, tpos_dat)
        datV = cutter(datV, latdatID, lapos_dat)
        datV = cutter(datV, londatID, lopos_dat)
    else:  # multiple data/climatology files, need processing
        if coarse != 'IBI':  # data is coarse, IBI needs block averaging
            if dtdat != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate_h5(pathdat, times[ff], dtboth, dchkboth, dpos_dat, interdat[0], trimdat, lapos,
                                       latID, lopos, lonID)
                    datagg = np.ma.masked_array.squeeze(datagg)
                    datV[ff, :, :, :] = datagg
                    datV_msk[ff, :, :, :] = datagg.mask
            else:  # if not, just read data, stack into array
                for ff in range(tlen):
                    data = data_read_h5(pathdat, ff, dchkboth, dpos_dat, interdat[0], trimdat, lapos, latID, lopos, lonID)
                    datV[ff, :, :, :] = data
                    datV_msk[ff, :, :, :] = data.mask
        else:  # IBI is coarser than data/climatology and so, data/climatology needs block averaging
            if dtdat != dtboth:  # check if aggregation needed
                # if so, do aggregation
                for ff in range(tlen):
                    datagg = aggregate_h5(pathdat, times[ff], dtboth, dchkboth, dpos_dat, interdat[0], trimdat, lapos,
                                       latID, lopos, lonID)
                    datfrm1 = np.zeros((dlen_dat, latlen, lonlen))
                    datfrm1, maskchk = blkav(datagg, datfrm1, latdat, londat, datdat, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    datV[ff, :, :, :] = datfrm1
                    datV_msk[ff, :, :, :] = maskchk
            else:  # if not, just read data, stack into array
                for ff in range(tlen):
                    data = data_read_h5(pathdat, ff, dchkboth, dpos_dat, interdat[0], trimdat, lapos, latID, lopos, lonID)
                    datfrm1 = np.zeros((dlen_dat, latlen, lonlen))
                    datfrm1, maskchk = blkav(data, datfrm1, latdat, londat,  datdat, finidx, dchkboth)
                    datfrm1 = np.ma.masked_where(maskchk == True, datfrm1)
                    datV[ff, :, :, :] = datfrm1
                    datV_msk[ff, :, :, :] = maskchk
    datV_datM_msk = datM_msk | datV_msk
    datM = np.ma.masked_where(datV_datM_msk == True, datM)
    datV = np.ma.masked_where(datV_datM_msk == True, datV)
    return datV, datM, times, tpos, cdepdat, dpos, dchkboth, latdat, latpos, londat, lonpos


def BLAV_ncout(ofile, varname, valdat, moddat, times, tpos, depths, dpos, latdat, latpos, londat, lonpos):
    ncout1 = Dataset(ofile, 'w', format='NETCDF4')
    ndep = valdat.shape[dpos]
    nlat = valdat.shape[latpos]
    nlon = valdat.shape[lonpos]
    ntim = valdat.shape[tpos]

    par_shp = []
    par_shp.insert(tpos, 'depth')
    par_shp.insert(tpos, 'time')
    # par_shp.insert(latpos, 'lat')
    # par_shp.insert(lonpos, 'lon')
    par_shp.insert(latpos, 'y')
    par_shp.insert(lonpos, 'x')
    shp_tup = tuple(par_shp)

    # define axis size
    ncout1.createDimension('time', ntim)
    ncout1.createDimension('depth', ndep)
    ncout1.createDimension('lat', nlat)
    ncout1.createDimension('lon', nlon)
    ncout1.createDimension('y', nlat)
    ncout1.createDimension('x', nlon)

    # create time axis
    time = ncout1.createVariable('time', np.dtype('double').char, ('time', ))
    time.long_name = 'time'
    time.units = 'hours since 1950-01-01 00:00:00'  # input units: days, output units: hours
    time.calendar = 'standard'
    # time.axis = 'T'

    # create depth axis
    depth = ncout1.createVariable('dep', np.dtype('double').char, ('depth',))
    depth.long_name = 'Depth'
    depth.units = 'm'
    # depth.axis = 'Z'

    # create latitude axis
    lat = ncout1.createVariable('lat', np.dtype('double').char, ('y', 'x'))
    lat.standard_name = 'latitude'
    lat.long_name = 'latitude'
    lat.units = 'degrees_north'
    # lat.axis = 'Y'

    # create longitude axis
    lon = ncout1.createVariable('lon', np.dtype('double').char, ('y', 'x'))
    lon.standard_name = 'longitude'
    lon.long_name = 'longitude'
    lon.units = 'degrees_east'
    # lon.axis = 'X'

    # create model array
    modelled = ncout1.createVariable(varname + '_MOD', np.dtype('double').char, shp_tup)
    modelled.standard_name = varname + '_IBI'
    modelled.long_name = varname + '_IBI'
    modelled.units = '-'
    # modelled.axis = 'Z'

    # create modelled array
    validate = ncout1.createVariable(varname + '_VAL', np.dtype('double').char, shp_tup)
    validate.standard_name = varname + '_val'
    validate.long_name = varname + '_val'
    validate.units = '-'
    # validate.axis = 'Z'

    # copy axis from original dataset
    time[:] = times[:]
    depth[:] = depths[:]
    lat[:] = latdat[:]
    lon[:] = londat[:]
    modelled[:] = moddat[:]
    validate[:] = valdat[:]

    # close files
    ncout1.close()
    return


def IBIsection(datval, datmod, timeposit, depths, latorlon, sect, lat, ltpos, lon, lnpos):
        s_array = datval.shape
        if latorlon == 'lat':
            dm1 = ltpos
            latlon = lat
            remainingdim = lon
            dm2 = lnpos
        elif latorlon == 'lon':
            dm1 = lnpos
            latlon = lon
            remainingdim = lat
            dm2 = ltpos
        if len(latlon.shape) == 2:
            if abs(latlon[1, 0] - latlon[0, 0]) > abs(latlon[0, 1] - latlon[0, 0]):
                cri = np.where(min(abs(latlon - sect)[:, 0]) == abs(latlon - sect)[:, 0])
                remainingdim = np.mean(remainingdim[cri, :], axis=1)
            else:
                cri = np.where(min(abs(latlon - sect)[0, :]) == abs(latlon - sect)[0, :])
                remainingdim = np.mean(remainingdim[:, cri], axis=0)
        elif len(latlon.shape) == 1:
            cri = np.where(min(abs(latlon - sect)) == abs(latlon - sect))
        if len(s_array) == 4:
            if dm1 == 0:
                datmod_Xsn = np.mean(datmod[cri, :, :, :], axis=dm2)
                datval_Xsn = np.mean(datval[cri, :, :, :], axis=dm2)
            elif dm1 == 1:
                datmod_Xsn = np.mean(datmod[:, cri, :, :], axis=dm2)
                datval_Xsn = np.mean(datval[:, cri, :, :], axis=dm2)
            elif dm1 == 2:
                datmod_Xsn = np.mean(datmod[:, :, cri, :], axis=dm2)
                datval_Xsn = np.mean(datval[:, :, cri, :], axis=dm2)
            elif dm1 == 3:
                datmod_Xsn = np.mean(datmod[:, :, :, cri], axis=dm2)
                datval_Xsn = np.mean(datval[:, :, :, cri], axis=dm2)
            else:
                datmod_Xsn = []
                datval_Xsn = []
        elif len(s_array) == 3:
            if dm1 == 0:
                datmod_Xsn = np.mean(datmod[cri, :, :], axis=dm2)
                datval_Xsn = np.mean(datval[cri, :, :], axis=dm2)
            elif dm1 == 1:
                datmod_Xsn = np.mean(datmod[:, cri, :], axis=dm2)
                datval_Xsn = np.mean(datval[:, cri, :], axis=dm2)
            elif dm1 == 2:
                datmod_Xsn = np.mean(datmod[:, :, cri], axis=dm2)
                datval_Xsn = np.mean(datval[:, :, cri], axis=dm2)
            else:
                datmod_Xsn = []
                datval_Xsn = []
        elif len(s_array) == 2:
            if dm1 == 0:
                datmod_Xsn = np.mean(datmod[cri, :], axis=dm2)
                datval_Xsn = np.mean(datval[cri, :], axis=dm2)
            elif dm1 == 1:
                datmod_Xsn = np.mean(datmod[:, cri], axis=dm2)
                datval_Xsn = np.mean(datval[:, cri], axis=dm2)
            else:
                datmod_Xsn = []
                datval_Xsn = []
        else:
                datmod_Xsn = []
                datval_Xsn = []
        datmod_Xsn = np.ma.mean(datmod_Xsn, axis=timeposit)
        datmod_Xsn = np.ma.squeeze(datmod_Xsn)
        datval_Xsn = np.ma.mean(datval_Xsn, axis=timeposit)
        datval_Xsn = np.ma.squeeze(datval_Xsn)

        return datval_Xsn, datmod_Xsn, remainingdim


def Hovsection(datval, datmod, timeposit, depths, depos, latorlon, sect, lat, ltpos, lon, lnpos, Hov_depth):
        s_array = datval.shape
        if latorlon == 'lat':
            dm1 = ltpos
            latlon = lat
            remainingdim = lon
            dm2 = lnpos
        elif latorlon == 'lon':
            dm1 = lnpos
            latlon = lon
            remainingdim = lat
            dm2 = ltpos
        if len(latlon.shape) == 2:
            if abs(latlon[1, 0] - latlon[0, 0]) > abs(latlon[0, 1] - latlon[0, 0]):
                cri = np.where(min(abs(latlon - sect)[:, 0]) == abs(latlon - sect)[:, 0])
                remainingdim = np.mean(remainingdim[cri, :], axis=1)
            else:
                cri = np.where(min(abs(latlon - sect)[0, :]) == abs(latlon - sect)[0, :])
                remainingdim = np.mean(remainingdim[:, cri], axis=0)
        elif len(latlon.shape) == 1:
            cri = np.where(min(abs(latlon - sect)) == abs(latlon - sect))
        if len(s_array) == 4:
            if dm1 == 0:
                datmod_Xsn = np.mean(datmod[cri, :, :, :], axis=dm2)
                datval_Xsn = np.mean(datval[cri, :, :, :], axis=dm2)
            elif dm1 == 1:
                datmod_Xsn = np.mean(datmod[:, cri, :, :], axis=dm2)
                datval_Xsn = np.mean(datval[:, cri, :, :], axis=dm2)
            elif dm1 == 2:
                datmod_Xsn = np.mean(datmod[:, :, cri, :], axis=dm2)
                datval_Xsn = np.mean(datval[:, :, cri, :], axis=dm2)
            elif dm1 == 3:
                datmod_Xsn = np.mean(datmod[:, :, :, cri], axis=dm2)
                datval_Xsn = np.mean(datval[:, :, :, cri], axis=dm2)
            else:
                datmod_Xsn = []
                datval_Xsn = []
        elif len(s_array) == 3:
            if dm1 == 0:
                datmod_Xsn = np.mean(datmod[cri, :, :], axis=dm2)
                datval_Xsn = np.mean(datval[cri, :, :], axis=dm2)
            elif dm1 == 1:
                datmod_Xsn = np.mean(datmod[:, cri, :], axis=dm2)
                datval_Xsn = np.mean(datval[:, cri, :], axis=dm2)
            elif dm1 == 2:
                datmod_Xsn = np.mean(datmod[:, :, cri], axis=dm2)
                datval_Xsn = np.mean(datval[:, :, cri], axis=dm2)
            else:
                datmod_Xsn = []
                datval_Xsn = []
        elif len(s_array) == 2:
            if dm1 == 0:
                datmod_Xsn = np.mean(datmod[cri, :], axis=dm2)
                datval_Xsn = np.mean(datval[cri, :], axis=dm2)
            elif dm1 == 1:
                datmod_Xsn = np.mean(datmod[:, cri], axis=dm2)
                datval_Xsn = np.mean(datval[:, cri], axis=dm2)
            else:
                datmod_Xsn = []
                datval_Xsn = []
        else:
                datmod_Xsn = []
                datval_Xsn = []

        if Hov_depth == 0:
            datmod_Xsn = surfext(datmod_Xsn, depos)
            datval_Xsn = surfext(datval_Xsn, depos)
        else:

            modstencil = np.zeros((datmod_Xsn.shape[timeposit], 0, datmod_Xsn.shape[ltpos], datmod_Xsn.shape[lnpos]))
            modstencil_msk = np.zeros((datmod_Xsn.shape[timeposit], 0, datmod_Xsn.shape[ltpos],
                                       datmod_Xsn.shape[lnpos]))
            valstencil = np.zeros((datval_Xsn.shape[timeposit], 0, datval_Xsn.shape[ltpos], datval_Xsn.shape[lnpos]))
            valstencil_msk = np.zeros((datval_Xsn.shape[timeposit], 0, datval_Xsn.shape[ltpos],
                                       datval_Xsn.shape[lnpos]))
            for ts in range(datmod.shape[timeposit]):
                mod_subset = datmod_Xsn[ts, :, :, :]
                d_count = sum(mod_subset.mask)
                d_limit = (datmod_Xsn.shape[depos] * np.ones_like(d_count))-d_count

                modstencil_ts = np.zeros((datmod_Xsn.shape[ltpos], datmod_Xsn.shape[lnpos]))
                if latorlon == 'lat':
                    datmod_Xsn_ts, datmod_msk = vertlevels(mod_subset, modstencil_ts, depths,
                                                           Hov_depth, sect, remainingdim, d_limit)
                else:
                    datmod_Xsn_ts, datmod_msk = vertlevels(mod_subset, modstencil_ts, depths,
                                                           Hov_depth, remainingdim, sect, d_limit)
                datmod_Xsn_ts = np.ma.masked_where(datmod_msk == True, datmod_Xsn_ts)
                modstencil[ts, :, :, :] = datmod_Xsn_ts
                modstencil_msk[ts, :, :, :] = datmod_msk

                val_subset = datval_Xsn[ts, :, :, :]
                d_count = sum(val_subset.mask)
                d_limit = (datval_Xsn.shape[depos] * np.ones_like(d_count))-d_count
                valstencil_ts = np.zeros((datval_Xsn.shape[ltpos], datval_Xsn.shape[lnpos]))
                if latorlon == 'lat':
                    datval_Xsn_ts, datval_msk = vertlevels(val_subset, valstencil_ts, depths,
                                                           Hov_depth, sect, remainingdim, d_limit)
                else:
                    datval_Xsn_ts, datval_msk = vertlevels(val_subset, valstencil_ts, depths,
                                                           Hov_depth, remainingdim, sect, d_limit)
                datval_Xsn_ts = np.ma.masked_where(datval_msk == True, datval_Xsn_ts)
                valstencil[ts, :, :, :] = datval_Xsn_ts
                valstencil_msk[ts, :, :, :] = datval_msk
            modstencil = np.ma.masked_where(modstencil_msk == True, modstencil)
            valstencil = np.ma.masked_where(valstencil_msk == True, valstencil)
            datmod_Xsn = np.ma.squeeze(modstencil)
            datval_Xsn = np.ma.squeeze(valstencil)

        return datval_Xsn, datmod_Xsn, remainingdim

