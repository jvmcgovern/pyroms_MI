# NOTE: This is a simplified script to demonstrate/test calculation of bias-corrected projections using
#       the matlab function fn_biascorrect_projections.m, applied to a test dataset of NORESM hindcast
#       and projection data covering the North Sea.
#       For an example of a real script used to make bias-corrected datasets and various analysis plots
#       see example_scripts/make_model_biascorrected_forcings.m.

# v1   -- matlab  -- Phil Wallhead, NIVA Bergen
# v1.1 -- python3 -- Joe McGovern, Marine Institute, Galway -- 26th September 2023


import numpy as np
from netCDF4 import Dataset as netcdf
import time
from datetime import date
import matplotlib.pyplot as plt
from fn_biascorrect_projections import fn_biascorrect_projections as fnbp
from fndatenum_to_decyr import fndatenum_to_decyr
from fn_spherical_distm import fn_spherical_distm


class opt:
    def __init__(self, seasonal, use_month, use_subdelta, fractional,
                 correct_iavar, dyear_iavar, XpcL, XpcH, legacy,
                 minX, tclimc, tol_dtclim, tol_dtclimp, use_XpcH_SWI,
                 latp, lonp, use_SWI_hav, yearpminv, yearpmaxv, Xhf, tdhf,
                 remove_deltascale_var, correct_substd, correct_subquantiles,
                 match_subdelta_deltamean, match_subdelta_hourly, fractional_subdelta,
                 minXsub, XpcfL, XpcfH, use_XpcfH_SWI, frcrit_capSWI, recalculate_Xpc, optSWI):
        self.seasonal = seasonal
        self.use_month = use_month
        self.use_subdelta = use_subdelta
        self.fractional = fractional
        self.correct_iavar = correct_iavar
        self.dyear_iavar = dyear_iavar
        self.XpcL = XpcL
        self.XpcH = XpcH
        self.legacy = legacy
        self.minX = minX
        self.tclimc = tclimc
        self.tol_dtclim = tol_dtclim
        self.tol_dtclimp = tol_dtclimp
        self.use_XpcH_SWI = use_XpcH_SWI
        self.latp = latp
        self.lonp = lonp
        self.use_SWI_hav = use_SWI_hav
        self.yearpminv = yearpminv
        self.yearpmaxv = yearpmaxv
        self.Xhf = Xhf
        self.tdhf = tdhf
        self.remove_deltascale_var = remove_deltascale_var
        self.correct_substd = correct_substd
        self.correct_subquantiles = correct_subquantiles
        self.match_subdelta_deltamean = match_subdelta_deltamean
        self.match_subdelta_hourly = match_subdelta_hourly
        self.fractional_subdelta = fractional_subdelta
        self.minXsub = minXsub
        self.XpcfL = XpcfL
        self.XpcfH = XpcfH
        self.use_XpcfH_SWI = use_XpcfH_SWI
        self.frcrit_capSWI = frcrit_capSWI
        self.recalculate_Xpc = recalculate_Xpc
        self.optSWI = optSWI


def create_new_netcdf(source_ncfile_path, new_ncfile_path, new_dimensions, history):
    # Example usage:
    # create_new_netcdf('source.nc', 'new.nc', {'time': None})

    # Step 1: Read the Source NetCDF File
    with netcdf(source_ncfile_path, 'r') as src_ncfile:
        # Reading global attributes
        global_attrs = {attr: getattr(src_ncfile, attr) for attr in src_ncfile.ncattrs()}

        # Reading dimensions
        dimensions = {dim: len(obj) for dim, obj in src_ncfile.dimensions.items()}

        # Reading variables and their attributes
        variables = {}
        for var, obj in src_ncfile.variables.items():
            variables[var] = {
                'datatype': obj.dtype,
                'dimensions': obj.dimensions,
                'attrs': {attr: getattr(obj, attr) for attr in obj.ncattrs()}
            }

    # Step 2: Modify the Dimensions
    dimensions.update(new_dimensions)

    # Step 3: Create a New NetCDF File
    with netcdf(new_ncfile_path, 'w') as new_ncfile:
        # Writing global attributes
        for attr, value in global_attrs.items():
            setattr(new_ncfile, attr, value)

        # Adding history as a global attribute
        setattr(new_ncfile, 'history', history)

        # Writing dimensions

        # Writing dimensions
        for dim, size in dimensions.items():
            new_ncfile.createDimension(dim, size)

        # Writing variables and their attributes
        for var, props in variables.items():
            var_obj = new_ncfile.createVariable(var, props['datatype'], props['dimensions'])
            for attr, value in props['attrs'].items():
                setattr(var_obj, attr, value)


varstrp = 'no3'

varstrh = varstrp
if varstrp == 'no3':
    varstrh = 'no3no2'

# Read in the raw projection data from netcdf.
fnamep = np.array(['../testdata/NorthSea/', varstrp,
                   '_Omon_NorESM2-MM_historical_ssp585_r1i1p1f1_gr_198001-210012_NorthSea.nc'])
ncp = netcdf(fnamep, 'r+')
ncpo = ncp.variables
Xp = np.array(ncpo[varstrp][:])
ncp.close()

print(np.array(['Loaded ', varstrp, ' from: ', fnamep]))
yearp1 = np.transpose((np.arange(1980, 2099 + 1)))

# Projection data times are in a strange unit -- easier to just redefine the times rather than read from file:
yearp = np.kron(yearp1, np.ones((12, 1)))
monthp = np.kron(np.ones((len(yearp1), 1)), np.transpose((np.arange(1, 13))))
dayp = 15 * np.ones((len(yearp1) * 12, 1))
tdp = np.zeros_like(yearp)
for tp in np.arange(0, len(yearp)):
    tdp[tp] = date.toordinal(date(yearp[tp], monthp[tp], dayp[tp]))
tyrp = fndatenum_to_decyr(tdp)
ntp = len(tdp)

# Read in the hindcast data from netcdf.
fnameh = np.array(['../testdata/NorthSea/NORESMOCv1p2_40N_1980_2020_', varstrh, '_reanal_NorthSea.nc'])
nch = netcdf(fnameh, 'r+')
ncho = nch.variables
Xh = np.array(ncho[varstrh][:])

print(np.array(['Loaded ', varstrh, ' from: ', fnameh]))
# No need to reread the spatial coordinate variables --- these are the same as in the projection data
tdh = date.toordinal(date(1970, 1, 1)) + np.array(ncho['time'][:])

yearh = np.zeros_like(tdh)
monthh = np.zeros_like(tdh)
dayh = np.zeros_like(tdh)
for ti in np.arange(0, len(tdh)):
    dfo = date.fromordinal(tdh[ti])
    yearh[ti] = dfo.year
    monthh[ti] = dfo.month
    dayh[ti] = dfo.day

tyrh = fndatenum_to_decyr(tdh)
nth = len(tdh)
# We take the spatial coordinates from the hindcast file --- the coordinates in the projection file
# are the same except that the range 0 to 360 is used for longitude rather than -180 to 180 as in the hindcast file.
lat1 = np.array(ncho['lat'][:])
lon1 = np.array(ncho['lon'][:])
z1 = np.array(ncho['depth'][:])
nch.close()

nz1 = len(z1)
# Reorganize 4D arrays (Xp,Xh) into 2D arrays [ntp,nth x ns] for use in fn_biascorrect_projections.m.
nx1 = Xp.shape[1 - 1]
ny1 = Xp.shape[2 - 1]
ns0 = nz1 * nx1 * ny1

Xhc = Xh
Xpc = Xp
Xh = np.nan * np.ones((nth, ns0))
Xp = np.nan * np.ones((ntp, ns0))
lat = np.nan * np.ones((1, ns0))
lon = lat
z = lat
q = 0
for i in np.arange(0, nz1).reshape(-1):
    for j in np.arange(0, ny1).reshape(-1):
        for k in np.arange(0, nx1).reshape(-1):
            q = q + 1
            Xh[np.arange(0, nth), q] = np.squeeze(Xhc[k, j, i, np.arange(1, nth + 1)])
            Xp[np.arange(0, ntp), q] = np.squeeze(Xpc[k, j, i, np.arange(1, ntp + 1)])
            lat[q] = lat1[k, j]
            lon[q] = lon1[k, j]
            z[q] = z1[i]

del Xhc, Xpc
#  Parse out the land mask, recording wet grid selection indices (needed later for writing to netcdf).
selan = np.argwhere(np.logical_and(not np.isnan(Xh[0, :]), not np.isnan(Xp[0, :])))[:, 0]

# Note: In this example both hindcast and projection data have the exact same land mask
# (check: sum(isnan(Xh(1,:))~=isnan(Xp(1,:)))=0).
#       If this is not the case it may be desirable to fill with nearest grid cells for data missing
#       in projection and/or hindcast datasets (filling land masks where one of (projection,hindcast) is available).
Xh = Xh[:, selan]
Xp = Xp[:, selan]
lat = lat[selan]
lon = lon[selan]
z = z[selan]
ns = len(selan)

# Set options for fn_biascorrect_projections.m.
# Note: these comprise a few basic 'essential' option parameters (Xp, tdp, Xh, tdh, yearminref, yearmaxref, method)
#       plus various extra option parameters passed via the structure variable 'opt'
yearminref = 2000
yearmaxref = 2019

method = 0

# 1 for Empirical Quantile Mapping (EQM), e.g. Cannon et al. (2015)
# 2 (def) for Repeated Hindcast Correction (RHC), e.g. Willems and Vrac (2011)
# Can also input as string ('DC','EQM,'RHC')
# See guidance text within fn_biascorrect_projections.m for a descriptive and discussion of the different methods.
opt.seasonal = 1
opt.use_month = 1
opt.fractional = 1
opt.tol_dtclim = 0
opt.qsel = np.arange(0.05, 0.95 + 0.1, 0.1)
opt.yearpminv = np.arange(1980, 2080 + 20, 20)
opt.yearpmaxv = np.arange(1999, 2099 + 20, 20)

# Note: The reference parameter choices were based on some limited testing of different value combinations
#       with the aim of optimizing the skill in reproducing monthly quantiles during the validation period (1980-1999),
#       considering *the whole pan-Arctic domain* i.e. at all latitudes >40N. These values might not be optimal
#       for the North Sea, and might not necessarily produce most-plausible results for many decades in the future.
#       A general finding seems to be, however, that the RHC method (2) is a lot less sensitive to these choices
#       of parameter values, and might therefore be considered a more robust method.
if method == 0:
    mstr1 = 'DC'

if method == 1:
    mstr1 = 'EQM'

if method == 2:
    mstr1 = 'RHC'

# Compute quantile-corrected projections (Xpc) covering 120 years (ntp=1440) and ns=7605 wet grid cells.

Xpc, out = fnbp(Xp, tdp, Xh, tdh, yearminref, yearmaxref, method, opt)
# DC method takes ~0.3 seconds
# EQM method takes ~4.2 seconds (~6.3 seconds with opt.use_quantile = 0)
# RHC method takes ~17 seconds (~23 seconds with opt.use_quantile = 0)
# Times recorded on Dell Precision 5530 Intel(R) Core(TM) i7-8850H CPU @ 2.60GHz, 32GB RAM, running Matlab R2019b

# Compare with reference solution
fnamepcr = np.array(['../testdata/NorthSea/reference_solutions/', varstrh,
                     '_NorESM2-NORESMreanalv1_1980-2099_NorthSea_', mstr1, '_v1_ref.nc'])
ncpcr = netcdf(fnamepcr, 'r+')
ncpcro = ncpcr.variables
Xpcrc = np.array(ncpcro[varstrh][:])
ncpcr.close()

Xpcr = np.nan * np.ones((ntp, ns0))
q = 0
for i in np.arange(1, nz1 + 1).reshape(-1):
    for j in np.arange(1, ny1 + 1).reshape(-1):
        for k in np.arange(1, nx1 + 1).reshape(-1):
            q = q + 1
            Xpcr[:, q] = np.squeeze(Xpcrc[k, j, i, :])

del Xpcrc
Xpcr = Xpcr[:, selan]
print(np.array(['Maximum discrepancy with reference solution = ',
                str(np.amax(np.amax(np.abs(Xpcr - Xpc)))), ' mol/m3']))

# Plot the solution!
# It is recommended to check/analyse the solution with several different plots, of which
# the quantile-quantile plots using the validation hindcast dataset (1980-1999) are especially
# recommended and may form the basis of a choise between different bias correction methods.
# See example_scripts/make_model_biascorrected_forcings.m for examples.
# Here we limit ourselves to perhaps the most basic analysis which is to compare the hindcast,
# projection, and bias-corrected projection time series at a selected point.


doplot_timeseries = 1

if doplot_timeseries == 1:
    tfontsize = 11
    lfontsize = 13
    lfontweight = 'bold'
    interpreter1 = 'tex'
    showtitle = 1
    varstr1 = varstrh
    tolstr = np.array([', tol=', str(opt.tol_dtclim), ' month(s)'])
    if method == 0:
        qstr1 = []
    else:
        qstr1 = np.array([', nq=', str(len(opt.qsel))])
    if method == 0:
        cstr1 = 'climatology'
    else:
        cstr1 = 'quantile'
    # latp1 = 59.24; lonp1 = 2.92; #N. North Sea: both EQM and RHC reference solutions look reasonable here
    latp1 = 54.88
    lonp1 = 2.85
    dist1 = fn_spherical_distm(lat, lon, latp1, lonp1)
    sels = np.argwhere(dist1 == np.amin(dist1))[:, 0]
    sels = sels(z(sels) == 5)
    if ns == 1:
        sels = 1
    plt.figure(40)
    plt.clf()
    show_Xp = 1
    show_Xpc = 1
    Xh_on_top = 1
    tyrminshow = 1980

    selh = np.argwhere(yearh >= tyrminshow)[:, 0]
    selp = np.argwhere(yearp >= tyrminshow)[:, 0]

    zstr1 = np.array([str(np.round(z[sels])), 'm'])
    tstr1 = np.array(
        [varstr1, ' monthly averages at (', str(lat[sels], '%3.1f'), 'N,', str(lon[sels], '%3.1f'), 'E,', zstr1,
         ') from NORESM reanalysis v1 (black)'])
    if show_Xp == 1:
        plt.plot(tyrp[selp], Xp[selp, sels], 'r-')
        tstr1 = np.array([tstr1, '; NORESM2-SSP585 (red)'])
    plt.show()
    plt.plot(tyrh[selh], Xh[selh, sels], 'k-')
    if show_Xpc == 1:
        plt.show()
        plt.plot(tyrp[selp], Xpc[selp, sels], 'c-')
        tstr1 = np.array(
            [tstr1, ';\n NORESM2-SSP585 with monthly ', cstr1, ' correction to NORESM reanalysis v1 (cyan)'])
    if Xh_on_top == 1:
        plt.show()
        plt.plot(tyrh(selh), Xh[selh, sels], 'k-')
    tstr1 = np.array([tstr1, ' (', mstr1, tolstr, qstr1, ')'])
    plt.axis('tight')
    if showtitle == 1:
        plt.title(tstr1, 'FontSize', tfontsize, 'FontWeight', 'bold', 'Interpreter', interpreter1)
    plt.xlabel('Year', 'FontSize', lfontsize, 'FontWeight', lfontweight, 'Interpreter', interpreter1)
    plt.ylabel(np.array([varstr1, ' [mol/m3]']), 'FontSize', lfontsize, 'FontWeight', lfontweight, 'Interpreter',
               interpreter1)
    plt.box(on=None)
    # Useful post-processing to remove x10 label in y ticks:
# ax = gca; ax.YAxis.Exponent = 0;

# Write to NetCDF if req'd
write_nc = 1
if write_nc == 1:
    verstr = np.array(['_', mstr1, '_v1'])
    yearwmin = 1980
    yearwmax = 2099
    modelstr1 = 'NorESM2-NORESMreanalv1'
    fnamepc = np.array([varstrh, '_', modelstr1, '_', str(yearwmin), '-', str(yearwmax), '_NorthSea', verstr, '.nc'])
    # Subset of time steps for writing
    sel1 = np.argwhere(np.logical_and(yearp >= yearwmin, yearp <= yearwmax))[:, 0]
    np1 = len(sel1)

    # Many different ways of writing to netcdf. Here we use the format from the hindcast file
    # modify the number of time steps to fit the bias-corrected projections.
    varstrt = 'time'

    # Add history variable as global attribute. This is important to document the options used for the bias correction.
    sourcehstr = 'NORESM reanalysis v1'
    history1 = np.array([date,
                         ': Bias-corrected projections calculated using make_testdata_biascorrected_forcings.m '
                         'with method=', str(method), ' (', mstr1, ')', ', yearminref=', str(yearminref), ', ',
                         'yearmaxref=', str(yearmaxref), ': ', '\n', '             ', sourcehstr, ' file(s): ',
                         fnameh, ';', '\n', '             ESM file(s): ', fnamep, '.'])

    create_new_netcdf(fnameh, fnamepc, {varstrt: np1}, history=history1)

    timep = tdp[sel1] - date.toordinal(date(1970, 1, 1))

    fpc = netcdf(fnamepc, 'r+')
    fpco = fpc.variables

    fh = netcdf(fnameh, 'r+')
    fho = fh.variables

    fpco['depth'][:] = fho['depth'][:]
    fpco['lat'][:] = fho['lat'][:]
    fpco['lon'][:] = fho['lon'][:]
    fpco[varstrt][:] = timep
    fpco['depth'][:] = z1

    setattr(fnamepc, 'long_name', varstrh)
    print(np.array(['Writing bias-corrected projections to: ', fnamepc]))
    sta = time.time()
    for j in np.arange(0, np1).reshape(-1):
        X1c = np.nan * np.ones((1, ns0))
        X1c[:, selan] = Xpc[sel1[j], :]
        X1c = np.reshape(X1c, nx1, ny1, nz1)
        fpco[varstrh][sel1[j], :, :, :] = X1c
        if np.mod(j, 100) == 0:
            print(np.array(['Done writing ', str[j], ' of out ', str[np1], ' time steps to netcdf']))
    del X1c
    print(np.array(['Done writing bias-corrected projections to ', fnamepc]))
    dne = time.time()
    print(dne - sta)

    # NOTE: Despite the fact that finfo.Variables.DeflateLevel = 5 for all variables, it seem ncwrite.m does not write
    #       compressed netcdf data. I therefore compress the netcdf afterwards using NCO
    #       (-> typically ~10x smaller file size):
    #       ncks -L 5 <fnamepc> tmp.nc
    #       mv tmp.nc <fnamepc>
