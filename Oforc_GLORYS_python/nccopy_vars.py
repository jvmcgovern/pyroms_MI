import numpy as np


def nccopy_vars(varnamev, oldnc, newnc, opt):
    # function nccopy_vars(varnamev,oldnc,newnc,opt)
    # Copy NetCDF variables varnamev [1*nvar], including all attributes,
    # from file oldnc to newnc.
    # Variable data is written for variables with opt.write_data = 1 [1*nvar] (def = 1)
    # Phil Wallhead 23/03/2019
    # Joe McGovern 26/09/2023 -- convert to python

    if not iscell(varnamev):
        varnamev = np.array([varnamev])

    nvarv = len(varnamev)
    if len(varargin) < 4:
        opt = []

    if hasattr(opt, 'verbose') == 1:
        verbose = opt.verbose
    else:
        verbose = 0

    if hasattr(opt, 'write_data') == 1:
        write_data = opt.write_data
    else:
        write_data = np.ones((1, nvarv))

    if hasattr(opt, 'write_atts') == 1:
        write_atts = opt.write_atts
    else:
        write_atts = np.ones((1, nvarv))

    if hasattr(opt, 'newvarnamev') == 1:
        newvarnamev = opt.newvarnamev
    else:
        newvarnamev = []

    if hasattr(opt, 'newdims') == 1:
        newdims = opt.newdims
    else:
        newdims = []

    if hasattr(opt, 'dimrename') == 1:
        dimrename = opt.dimrename
    else:
        dimrename = []

    if hasattr(opt, 'newdimname') == 1:
        newdimname = opt.newdimname
    else:
        newdimname = []

    if hasattr(opt, 'dimsubselname') == 1:
        dimsubselname = opt.dimsubselname
    else:
        dimsubselname = []

    if hasattr(opt, 'subsel') == 1:
        subsel = opt.subsel
    else:
        subsel = []

    if hasattr(opt, 'newfillvalue') == 1:
        newfillvalue = opt.newfillvalue
    else:
        newfillvalue = np.nan * np.ones((1, nvarv))

    if hasattr(opt, 'newmissing_value') == 1:
        newmissing_value = opt.newmissing_value
    else:
        newmissing_value = np.nan * np.ones((1, nvarv))

    if hasattr(opt, 'newdatatype') == 1:
        newdatatype = opt.newdatatype
    else:
        newdatatype = []

    if hasattr(opt, 'newformat') == 1:
        newformat = opt.newformat
    else:
        newformat = []

    if hasattr(opt, 'newunits') == 1:
        newunits = opt.newunits
    else:
        newunits = []

    if hasattr(opt, 'newcycle_length') == 1:
        newcycle_length = opt.newcycle_length
    else:
        newcycle_length = np.nan * np.ones((1, nvarv))

    if hasattr(opt, 'newlong_name') == 1:
        newlong_name = opt.newlong_name
    else:
        newlong_name = []

    if hasattr(opt, 'newfield') == 1:
        newfield = opt.newfield
    else:
        newfield = []

    if hasattr(opt, 'newfile') == 1:
        newfile = opt.newfile
    else:
        newfile = 1

    if len(newvarnamev) == 0 == 1:
        newvarnamev = varnamev

    if len(write_data) == 1:
        write_data = write_data * np.ones((1, nvarv))

    if len(write_atts) == 1:
        write_atts = write_atts * np.ones((1, nvarv))

    if len(newfillvalue) == 1:
        newfillvalue = newfillvalue * np.ones((1, nvarv))

    if len(newmissing_value) == 1:
        newmissing_value = newmissing_value * np.ones((1, nvarv))

    if len(newcycle_length) == 1:
        newcycle_length = newcycle_length * np.ones((1, nvarv))

    if len(newfield) == 1:
        newfield = newfield * np.ones((1, nvarv))

    samev = cell(1, nvarv)
    for i in np.arange(1, nvarv + 1).reshape(-1):
        samev[i] = 'same'

    if len(newdatatype) == 0 == 1:
        newdatatype = samev

    if len(newunits) == 0 == 1:
        newunits = samev

    if len(newlong_name) == 0 == 1:
        newlong_name = samev

    if len(newfield) == 0 == 1:
        newfield = samev

    if len(subsel) == 0 == 0:
        nsubsel = len(subsel)

    if newfile == 0:
        finfonew = ncinfo(newnc)
        varNamesnew = np.array([finfonew.Variables.Name])
    else:
        varNamesnew = []

    for j in np.arange(1, nvarv + 1).reshape(-1):
        vinfo0 = ncinfo(oldnc, varnamev[j])
        # Gather attributes and fillvalue from old file
        if len(vinfo0.Attributes) == 0 == 0:
            attNames0 = np.array([vinfo0.Attributes.Name])
            attValues0 = np.array([vinfo0.Attributes.Value])
            if np.isnan(newfillvalue(j)):
                fillvalue1 = []
                for i in np.arange(1, len(attNames0) + 1).reshape(-1):
                    if strcmpi(attNames0[i], '_FillValue') == 1:
                        fillvalue1 = attValues0[i]
            else:
                fillvalue1 = newfillvalue(j)
        # Gather dimensions from old file (or rename dimensions using newdims)
        if len(newdims) == 0 == 1:
            dims0 = vinfo0.Dimensions
            ndims0 = len(dims0)
            dims1 = cell(1, 2 * ndims0)
            for i in np.arange(1, ndims0 + 1).reshape(-1):
                dims11 = dims0(i).Name
                if str(dims11) == str(dimrename) == 1:
                    dims11 = newdimname
                dims1[[i - 1] * 2 + 1] = dims11
                if vinfo0.Dimensions(i).Unlimited == 1:
                    dims1[i * 2] = np.inf
                else:
                    if str(dims11) == str(dimsubselname) == 1:
                        dims1[i * 2] = nsubsel
                        dimsubsel = i
                    else:
                        dims1[i * 2] = dims0(i).Length
        else:
            dims1 = newdims[j]
        # Modify datatype of format here if req'd
        if strcmpi(newdatatype[j], 'same') == 0:
            newdatatype1 = newdatatype[j]
        else:
            newdatatype1 = vinfo0.Datatype
        if len(newformat) == 0 == 1:
            newformat1 = vinfo0.Format
        else:
            newformat1 = newformat
        # Create the NetCDF variable
        if not np.any(str(varNamesnew) == str(newvarnamev[j])):
            if not len(fillvalue1) == 0 and str(newformat1) == str('64bit') == 0:
                nccreate(newnc, newvarnamev[j], 'Dimensions', dims1, 'DataType', newdatatype1, 'FillValue', fillvalue1,
                         'DeflateLevel', vinfo0.DeflateLevel, 'Shuffle', vinfo0.Shuffle, 'Format', newformat1)
            else:
                nccreate(newnc, newvarnamev[j], 'Dimensions', dims1, 'DataType', newdatatype1, 'DeflateLevel',
                         vinfo0.DeflateLevel, 'Shuffle', vinfo0.Shuffle, 'Format', newformat1)
            if verbose == 1:
                print(np.array(['Created netcdf variable ', newvarnamev[j]]))
        # Write the data if req'd
        if write_data(j) == 1:
            X1 = ncread(oldnc, varnamev[j])
            if len(dimsubselname) == 0 == 0:
                indstr = '('
                for i in np.arange(1, dimsubsel - 1 + 1).reshape(-1):
                    indstr = np.array([indstr, ':,'])
                indstr = np.array([indstr, 'subsel'])
                for i in np.arange(1, ndims0 - dimsubsel + 1).reshape(-1):
                    indstr = np.array([indstr, ',:'])
                indstr = np.array([indstr, ')'])
                eval(np.array(['X1c = X1', indstr, ';']))
                ncwrite(newnc, newvarnamev[j], X1c)
            else:
                ncwrite(newnc, newvarnamev[j], X1)
            if verbose == 1:
                print(np.array(['Written data for netcdf variable ', newvarnamev[j]]))
        # Write the attributes if req'd
        if len(vinfo0.Attributes) == 0 == 0 and write_atts(j) == 1:
            attNames1 = attNames0
            attValues1 = attValues0
            for i in np.arange(1, len(attNames1) + 1).reshape(-1):
                if strcmpi(attNames1[i], '_FillValue') == 0:
                    ncwriteatt(newnc, newvarnamev[j], attNames1[i], attValues1[i])
            if verbose == 1:
                print(np.array(['Written attributes for netcdf variable ', newvarnamev[j]]))
        # Modify attributes as req'd
        if strcmpi(newunits[j], 'same') == 0:
            fileattrib(newnc, '+w')
            ncwriteatt(newnc, newvarnamev[j], 'units', newunits[j])
            if verbose == 1:
                print(np.array(['Modified units attribute for netcdf variable ', newvarnamev[j]]))
        if not np.isnan(newcycle_length(j)):
            if newcycle_length(j) == np.inf:
                ncid = netcdf.open(newnc, 'NC_WRITE')
                varid = netcdf.inqVarID(ncid, newvarnamev[j])
                netcdf.reDef(ncid)
                netcdf.delAtt(ncid, varid, 'cycle_length')
                netcdf.close(ncid)
            else:
                fileattrib(newnc, '+w')
                ncwriteatt(newnc, newvarnamev[j], 'cycle_length', newcycle_length(j))
            if verbose == 1:
                print(np.array(['Modified cycle_length attribute for netcdf variable ', newvarnamev[j]]))
        if strcmpi(newlong_name[j], 'same') == 0:
            fileattrib(newnc, '+w')
            ncwriteatt(newnc, newvarnamev[j], 'long_name', newlong_name[j])
            if verbose == 1:
                print(np.array(['Modified long_name attribute for netcdf variable ', newvarnamev[j]]))
        if not np.isnan(newmissing_value(j)):
            fileattrib(newnc, '+w')
            ncwriteatt(newnc, newvarnamev[j], 'missing_value', newmissing_value(j))
            if verbose == 1:
                print(np.array(['Modified missing_value attribute for netcdf variable ', newvarnamev[j]]))
        if strcmpi(newfield[j], 'same') == 0:
            fileattrib(newnc, '+w')
            ncwriteatt(newnc, newvarnamev[j], 'field', newfield[j])
            if verbose == 1:
                print(np.array(['Modified field attribute for netcdf variable ', newvarnamev[j]]))
