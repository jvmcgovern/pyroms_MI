#
# This file handles:
# (i)  Reanalysis hindcast input oceanic boundary condition files and bulk atmospheric forcing files for CROCO
#      Ocean bry: CMEMS IBI MY PHY, CMEMS IBI MY BGC Non-assimilative hindcast
#      Atm:       ECMWF ERA5, hourly, "online-handling" bulk data
# (ii) Climate projection input oceanic boundary condition files and atmospheric for CROCO
#      Interpolated from CMIP6 model of choice to file structure analogous to data in (i)
#
# The file prepares both of the aforementioned datasets to implement any of the three projection bias-correction
# methodologies which have been provided coded in matlab format by Phil Wallhead of NIVA in Bergen, and recoded to
# python 3 by Joe McGovern, Marine Institute, Rinville, Oranmore, Co. Galway.
#
# The methodologies are:
# (a) Delta-change
# (b) Empirical Quantile Mapping
# (c) Repeated Hindcast Correction
#
# V1.0 Joe McGovern, Marine Institute, Rinville, Oranmore, Co. Galway.      25th August 2023

import glob
import numpy as np
from netCDF4 import Dataset as netcdf
from netCDF4 import date2index as d2i
from calendar import monthrange, isleap
from datetime import date, datetime
import croco_glorys as glor
from fn_biascorrect_projections import fn_biascorrect_projections as fnbp
from interp_Cgrid import *


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


class out:
    def __init__(self,
                 qXp, minX, nseltclim, Xpclim, Xpstdclim,
                 qXpref, nseltclimh, Xhclim, qXhref, Fhref,
                 Xhstdclim, Rstdclim, Rqref, dqXp, Xpco, XpcH_SWI,
                 dXhfstdsub, dXhfminsub, dXhfmaxsub, dXhfrangesub,
                 IsortXhfsub, bs, bsres, Jmins, yearpminv, yearpmaxv,
                 correction_factor, Xpcfo, XpcfH_SWI, Xpco2, Xpcf,
                 tdpf, selhsub, selhsubhour):
        self.qXp = qXp
        self.minX = minX
        self.nseltclim = nseltclim
        self.Xpclim = Xpclim
        self.Xpstdclim = Xpstdclim
        self.qXpref = qXpref
        self.nseltclimh = nseltclimh
        self.Xhclim = Xhclim
        self.qXhref = qXhref
        self.Fhref = Fhref
        self.Xhstdclim = Xhstdclim
        self.Rstdclim = Rstdclim
        self.Rqref = Rqref
        self.dqXp = dqXp
        self.Xpco = Xpco
        self.XpcH_SWI = XpcH_SWI
        self.dXhfstdsub = dXhfstdsub
        self.dXhfminsub = dXhfminsub
        self.dXhfmaxsub = dXhfmaxsub
        self.dXhfrangesub = dXhfrangesub
        self.IsortXhfsub = IsortXhfsub
        self.bs = bs
        self.bsres = bsres
        self.Jmins = Jmins
        self.yearpminv = yearpminv
        self.yearpmaxv = yearpmaxv
        self.correction_factor = correction_factor
        self.Xpcfo = Xpcfo
        self.XpcfH_SWI = XpcfH_SWI
        self.Xpco2 = Xpco2
        self.Xpcf = Xpcf
        self.tdpf = tdpf
        self.selhsub = selhsub
        self.selhsubhour = selhsubhour


ocean_ana = 1
atm_ana = 1

# Default options
opt.use_daily = 1
interp_to_gregorian = 1
if opt.use_daily == 0:
    interp_to_gregorian = 0
test_method = 0

opt.seasonal = 1
opt.use_month = 0
if opt.use_daily == 1:
    opt.use_month = 0
opt.correct_iavar = 0

# Hindcast year range for file names
yearmin0 = 2000
yearmax0 = 2019

# Hindcast year range for file names
yearminp = 2000
yearmaxp = 2099

opt.tclimc = np.arange(0.5, 365.5, 1)

# 0 for Delta-Change (DC),
# 1 for Empirical Quantile Mapping (EQM),
# 2 for Repeated Hindcast Correction (RHC)
method = 2

opt.legacy = 0

opt.use_subdelta = 1
opt.remove_deltascale_var = 1
opt.correct_substd = 0
opt.correct_subquantiles = 1
if method == 2:
    opt.correct_subquantiles = 0
opt.match_subdelta_deltamean = 0
opt.match_subdelta_hourly = 0
opt.use_calc_surfaceSWI = 1

opt.yearpminv = np.arange(2000, 2090, 10)
opt.yearpmaxv = np.arange(2009, 2099, 10)
if method == 2:
    opt.yearpminv = np.arange(2000, 2090, 10)
    opt.yearpmaxv = np.arange(2009, 2099, 10)

# # Grid  point subselection indices(used only for atmospheric forcings)
# Whole ERA5 domain for ROHO800
latsel = np.arange(0, 24, 1)
lonsel = np.arange(0, 36, 1)

nlatsel = len(latsel)
nlonsel = len(lonsel)
ns = nlatsel*nlonsel

bcstrv = ['_west', '_north', '_east', '_south']

Mstart = 1
Mend = 12

if ocean_ana == 1:
    nbcstr = len(bcstrv)
else:
    nbcstr = 1

for ibcstr in range(0, nbcstr):
    if ocean_ana == 1:
        bcstr = bcstrv[ibcstr]

tstrv = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
         'November', 'December']

# Reference period to define climatologies / quantile correction factors
# Default values if not specified later
yearminref = 2010
yearmaxref = 2019

HC_bry_dir = '/media/dskfour/CS1KM_19932021_BCd/'
# HC_bry_gnrc = 'croco_bry_my_bc_mthMERCATOR_hc50m_MSL_ng'
HC_bry_gnrc = 'croco_bry_my_bc_MERCATOR_hc50m_MSL_ng'

PR_bry_dir_in = '/media/dskthree/UKESM1-0-LL/croco_bry_585/'
PR_bry_gnrc_in = 'croco_bry_UKESM1-0-LL_ssp585_hc50m_MSL_ng'

PR_bry_dir_out = '/media/dskthree/UKESM1-0-LL/croco_bry_bc_585/'
PR_bry_gnrc_out = 'croco_bry_bc_UKESM1-0-LL_ssp585_hc50m_MSL_ng'

HC_atm_dir = '/media/dskone/CELTIC/ECMWF/ERA5/ERA5_CELTIC_PISCES_II/'
PR_atm_dir_in = '/media/dsksix/UKESM1-0-LL/NATIVE/'
PR_atm_dir_out = '/media/dskone/CELTIC/ECMWF/ERA5/UKESM_atm_CELTIC_PISCES_II/'

# Yr_init_ana = 2010
# Yr_end_ana = 2019
Yr_init_ana = yearminref
Yr_end_ana = yearmaxref
Ana_len = Yr_end_ana - Yr_init_ana + 1

Yr_init_proj1 = 2015
Yr_end_proj1 = Yr_init_proj1 + Ana_len - 1

Yr_init_proj2 = Yr_end_proj1
Yr_end_proj2 = Yr_init_proj2 + Ana_len - 1

Yr_init_proj3 = Yr_end_proj2
Yr_end_proj3 = Yr_init_proj3 + Ana_len - 1

if ocean_ana == 1:
    # Reference period to define climatologies / quantile correction factors
    yearminref = 2010
    yearmaxref = 2019

    XpcL = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    XpcH = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    XpcfL = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    XpcfH = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    minX = [0, 0, 0, 0, 0, 0, 0, 0]
    minXsub = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    ocnvars = ['temp', 'salt', 'zeta', 'ubar', 'vbar', 'u', 'v',
               'DIC', 'TALK', 'NO3', 'PO4', 'Si', 'O2', 'NH4', 'FER']
    tol_dtclim = [50, 50, 50, 50, 50, 50, 50, 1, 1, 1, 1, 1, 1, 1, 1]
    tol_dtclimp = [50, 50, 50, 50, 50, 50, 50, 1, 1, 1, 1, 1, 1, 1, 1]
    fractional = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    fractional_subdelta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    seasonal = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    use_daily = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    use_month = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
    ocntstr = ['temp_time', 'salt_time', 'zeta_time', 'v2d_time', 'v2d_time', 'v3d_time', 'v3d_time',
               'dic_time', 'talk_time', 'no3_time', 'po4_time', 'si_time', 'o2_time', 'nh4_time', 'fer_time']
    source0str = 'IBI PHY reanalysis and climatologically bias-corrected IBI BGC reanalysis'

    for direction in range(0, len(bcstrv)):
        for ocv in range(0, len(ocnvars)):
            opt.XpcL = XpcL[ocv]
            opt.XpcH = XpcH[ocv]
            opt.XpcfL = XpcfL[ocv]
            opt.XpcfH = XpcfH[ocv]
            opt.minX = minX[ocv]
            opt.minXsub = minXsub[ocv]
            opt.tol_dtclim = tol_dtclim[ocv]
            opt.tol_dtclimp = tol_dtclimp[ocv]
            opt.fractional = fractional[ocv]
            opt.fractional_subdelta = fractional_subdelta[ocv]
            opt.seasonal = seasonal[ocv]
            opt.use_daily = use_daily[ocv]
            opt.use_month = use_month[ocv]

            for iyear in range(yearmin0, yearmax0 + 1):
                for imonth in range(Mstart, Mend + 1):
                    brynameHC = HC_bry_dir + HC_bry_gnrc + '_Y' + str(iyear) + 'M' \
                                + str(imonth).zfill(2) + '.nc'
                    ncbrh = netcdf(brynameHC, 'r+')
                    ncbrhy = ncbrh.variables
                    if (iyear == yearmin0) & (imonth == Mstart):
                        X0 = np.array(ncbrhy[ocnvars[ocv] + bcstrv[direction]][:])
                        td0 = np.array(ncbrhy[ocntstr[ocv]][:])
                    else:
                        X0 = np.concatenate((X0, np.array(ncbrhy[ocnvars[ocv] + bcstrv[direction]][:])))
                        td0 = np.concatenate((td0, np.array(ncbrhy[ocntstr[ocv]][:])))
                    ncbrh.close()

            td0_unq, td0_idx = np.unique(td0, return_index=True)
            td0 = td0_unq
            if len(X0.shape) == 4:
                X0 = X0[td0_idx, :, :, :]
            else:
                X0 = X0[td0_idx, :, :]

            for iyear in range(yearminp, yearmaxp + 1):
                for imonth in range(Mstart, Mend + 1):
                    brynamePRJ = PR_bry_dir_in + PR_bry_gnrc_in + '_Y' + str(iyear) + 'M' \
                                + str(imonth).zfill(2) + '.nc'
                    ncbrp = netcdf(brynamePRJ, 'r+')
                    ncbrpy = ncbrp.variables
                    if (iyear == yearmin0) & (imonth == Mstart):
                        XP = np.array(ncbrpy[ocnvars[ocv] + bcstrv[direction]][:])
                        tdp = np.array(ncbrpy[ocntstr[ocv]][:])
                    else:
                        XP = np.concatenate((XP, np.array(ncbrpy[ocnvars[ocv] + bcstrv[direction]][:])))
                        tdp = np.concatenate((tdp, np.array(ncbrpy[ocntstr[ocv]][:])))
                    ncbrp.close()

            tdp_unq, tdp_idx = np.unique(tdp, return_index=True)
            tdp = tdp_unq
            if len(XP.shape) == 4:
                XP = XP[tdp_idx, :, :, :]
            else:
                XP = XP[tdp_idx, :, :]
            # Convert arrays to matrices for processing purposes
            X0m = np.reshape(X0, (X0.shape[0], X0.shape[1] * X0.shape[2]), 'F')
            XPm = np.reshape(XP, (XP.shape[0], XP.shape[1] * XP.shape[2]), 'F')
            # Following expression results in zero
            # XP[:, 19, 1] - XPm[:, 39]
            Xpc, out = fnbp(XPm, tdp, X0m, td0, yearminref, yearmaxref, method, opt)

if atm_ana == 1:
    # Reference period to define climatologies / quantile correction factors
    yearminref = 2010
    yearmaxref = 2019

    XpcL = [np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, 0]
    XpcH = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    XpcfL = [np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, 0]
    XpcfH = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    minX = [0, 0, 0, 0, 0, 0, 0, 0]
    minXsub = [np.nan, np.nan, np.nan, np.nan, np.nan, 0, np.nan, -1e-4]
    atmvars = ['T2M', 'pmer', 'Q', 'U10M', 'V10M', 'SSR', 'STRD', 'TP']
    tol_dtclim = [50, 50, 50, 5, 5, 100, 100, 100]
    fractional = [0, 0, 0, 0, 0, 0, 0, 0]
    fractional_subdelta = [0, 0, 0, 0, 0, 1, 0, 1]
    seasonal = [1, 1, 1, 0, 0, 1, 1, 1]
    use_daily = [1, 1, 1, 1, 1, 1, 1, 1]
    use_month = [0, 0, 0, 0, 0, 0, 0, 0]
    atmtstr = ['time', 'time', 'time', 'time', 'time', 'time', 'time', 'time']
    frcrit_capSWI = 0.0
    source0str = 'ERA5 hindcast'

    for acv in range(0, len(atmvars)):
        opt.XpcL = XpcL[acv]
        opt.XpcH = XpcH[acv]
        opt.XpcfL = XpcfL[acv]
        opt.XpcfH = XpcfH[acv]
        opt.minX = minX[acv]
        opt.minXsub = minXsub[acv]
        opt.tol_dtclim = tol_dtclim[acv]
        opt.fractional = fractional[acv]
        opt.fractional_subdelta = fractional_subdelta[acv]
        opt.seasonal = seasonal[acv]
        opt.use_daily = use_daily[acv]
        opt.use_month = use_month[acv]

        for iyear in range(yearmin0, yearmax0 + 1):
            for imonth in range(Mstart, Mend + 1):
                atmnameHC = HC_atm_dir + atmvars[acv] + '_Y' + str(iyear) + 'M' \
                            + str(imonth).zfill(2) + '.nc'
                ncatmh = netcdf(atmnameHC, 'r+')
                ncatmho = ncatmh.variables
                if (iyear == yearmin0) & (imonth == Mstart):
                    X0 = np.array(ncatmho[atmvars[acv]][:])
                    td0 = np.array(ncatmho[atmtstr[acv]][:])
                else:
                    X0 = np.concatenate((X0, np.array(ncatmho[atmvars[acv]][:])))
                    td0 = np.concatenate((td0, np.array(ncatmho[atmtstr[acv]][:])))
                ncatmh.close()

        td0_unq, td0_idx = np.unique(td0, return_index=True)
        td0 = td0_unq
        X0 = X0[td0_idx, :, :]

        for iyear in range(yearminp, yearmaxp + 1):
            for imonth in range(Mstart, Mend + 1):
                atmnamePRJ = PR_atm_dir_in + atmvars[acv] + '_Y' + str(iyear) + 'M' \
                             + str(imonth).zfill(2) + '.nc'
                ncatmp = netcdf(atmnamePRJ, 'r+')
                ncatmpy = ncatmp.variables
                if (iyear == yearmin0) & (imonth == Mstart):
                    XP = np.array(ncatmpy[atmvars[acv]][:])
                    tdp = np.array(ncatmpy[atmtstr[acv]][:])
                else:
                    XP = np.concatenate((XP, np.array(ncatmpy[atmvars[acv]][:])))
                    tdp = np.concatenate((tdp, np.array(ncatmpy[atmtstr[acv]][:])))
                ncatmp.close()

        tdp_unq, tdp_idx = np.unique(tdp, return_index=True)
        tdp = tdp_unq
        XP = XP[tdp_idx, :, :]

        Xpc, out = fnbp(XP, tdp, X0, td0, yearminref, yearmaxref, method, opt)

#   read netcdfs on a loop
#   either over all ocean boundary files
#   setting up a new output boundary file types (empty) only once per year/month combo
#   then processing one variable at a time, and when done, opening the netcdf, inserting the bias-corrected
#   data file by file and moving on to the next ocean variable
#
#   or
#   setting up a new output atmospheric forcing file types (empty) only once per year/month combo
#   process one atmospheric variable at a time and copying bias-corrected data file by file

crocofiles_dir = '/home/pete/PycharmProjects/pyroms_MI/CELTIC_SEA/'
bryfiles_dir = '/media/dskthree/UKESM1-0-LL/'
grdname = crocofiles_dir + 'CELTIC_grd_v5_MSL_h10.nc'  # LAT to MSL included, INFOMAR block averaged hmin = 10m

N = 20
theta_s = 7.
theta_b = 0.
hc = 50.
vtransform = 2

obc = [1, 1, 1, 1]  # open boundaries (1=open , [S E N W])

time_bry = 0.
cycle_bry = 0.

Yorig = 1990  # year origin of time : days since Yorig-01-01

Ystart = 2015
Mstart = 1
Dstart = 1

Yend = 2099
Mend = 12
Dend = 31

ERA5sfold = '/media/dskone/CELTIC/ECMWF/ERA5/ERA5_CELTIC_PISCES_II/'
ERA5sfile = 'U10M_Y1999M06.nc'

CMIPfold = '/media/dsksix/UKESM1-0-LL/NATIVE/'

PHYfold = 'PHY/'
BGCfold = 'BGC/'
ATMfold = 'ATM/'

CMIP6_ATMstr  = ['rsds_rsus', 'rlds', 'tas', 'psl',  'pr',  'huss', 'uas',  'vas']
CMIP6_ATMlab  = ['rsds', 'rlds', 'tas', 'psl', 'pr', 'huss', 'uas', 'vas']
CMIP6_ATMtres = ['3hr',       '3hr',  'day', 'day',  '3hr', 'day',  'day',  'day']
ERA5_ATMstr   = ['SSR',       'STRD', 'T2M', 'pmer', 'TP',  'Q',    'U10M', 'V10M']
atm_ln        = ['surface_net_solar_radiation', 'surface_thermal_radiation_downwards', '2m_temperature',
                 'mean_sea_level_pressure', 'total_precipitation', 'specific_humidity', '10m_u_component_of_wind',
                 '10m_v_component_of_wind']
atm_units     = ['W m-2', 'W m-2', 'K', 'Pa', 'kg m-2 s-1', 'kg kg-1', 'm s-1', 'm s-1']
ATMpost       = '_ERA5g.nc'

phytracvars   = ['thetao', 'so']
bgctracvars   = ['no3', 'po4', 'si', 'o2', 'dfe', 'dissic', 'talk']
ele           = ['zos']
vec_vars      = ['uo', 'vo']

model = 'UKESM1-0-LL'

hist = 'historical'
# proj = 'ssp245'
proj = 'ssp585'
exp = 'r1i1p1f2_gn'
freq = 'Omon'
#
# get the sample extent
#

nces = netcdf(ERA5sfold + ERA5sfile, 'r')
nceso = nces.variables
lon = np.array(nces['lon'][:])
lat = np.array(nces['lat'][:])
lonlen = np.shape(lon)[0]
latlen = np.shape(lat)[0]
nces.close()
# Need to list all netcdfs for each of the variables to interpolate.
# For each file, it will be possible to know what years are covered by it by identifying the year in file name
# For each year and month, it will be necessary to identify the index of the pre-, present- and post- month in their
# respective files and pass the files to the interpolation process.
# Need to firstly sample one of each files to explore that lat/lon data and do the grid weighting

for v in range(0, len(ERA5_ATMstr)):
    if CMIP6_ATMtres[v] == '3hr':
        tpad = '_????????????-????????????'
        ysta1 = -34
        yend1 = -30
        ysta2 = -21
        yend2 = -17
        hrinit = 1
        hrend = 22
        mint = 30
        tinit = 0.0625
        tendi = 1
        tincr = 0.125

    elif CMIP6_ATMtres[v] == 'day':
        tpad = '_????????-????????'
        ysta1 = -26
        yend1 = -22
        ysta2 = -17
        yend2 = -13
        hrinit = 12
        hrend = 12
        mint = 00
        tinit = 0.5
        tendi = 1
        tincr = 1

    # CMIP6 ATM files
    CMIP6_hist = sorted(glob.glob(CMIPfold + ATMfold + CMIP6_ATMstr[v] + '_' + CMIP6_ATMtres[v] + '_' +
                                  model + '_' + hist +
                                  '_' + exp + tpad + ATMpost))
    CMIP6_proj = sorted(glob.glob(CMIPfold + ATMfold + CMIP6_ATMstr[v] + '_' + CMIP6_ATMtres[v] + '_' +
                                  model + '_' + proj +
                                  '_' + exp + tpad + ATMpost))
    CMIP6files = CMIP6_hist + CMIP6_proj
    files_yrs = np.ones((len(CMIP6files), 2), dtype=int)
    for file in range(0, len(CMIP6files)):
        files_yrs[file, 0] = int(CMIP6files[file][ysta1:yend1])
        files_yrs[file, 1] = int(CMIP6files[file][ysta2:yend2])

    for iyear in range(Ystart, Yend + 1):
        for imonth in range(Mstart, Mend + 1):

            atmname = CMIPfold + ERA5_ATMstr[v] + '_Y' + str(iyear) + 'M' + str(imonth).zfill(2) + '.nc'

            print(' ')
            print(' Making atm file: ' + atmname)
            print(' ')

            #
            # Create the CROCO boundary file
            #
            atime = 0

            glor.create_ERA5style_ESMatm(atmname, ERA5_ATMstr[v], atm_ln[v], atm_units[v], atime, latlen, lonlen)

            if imonth is 1:
                files_3m = [np.argwhere((files_yrs[:, 0] <= iyear-1) & (files_yrs[:, 1] >= iyear-1))[0][0],
                            np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                            np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0]]

            elif imonth == 12:
                files_3m = [np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                            np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                            np.argwhere((files_yrs[:, 0] <= iyear+1) & (files_yrs[:, 1] >= iyear+1))[0][0]]

            else:
                files_3m = [np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                            np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0],
                            np.argwhere((files_yrs[:, 0] <= iyear) & (files_yrs[:, 1] >= iyear))[0][0]]

            files_3m = np.unique(files_3m)

            if imonth == 1:
                Tinin_pre = datetime(iyear - 1, 12, monthrange(iyear - 1, 12)[1] - 2,
                                     hrinit, mint, 0)
                Tinin_post = datetime(iyear, imonth + 1, 3, hrend, mint, 0)
            elif imonth == 12:
                Tinin_pre = datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                     hrinit, mint, 0)
                Tinin_post = datetime(iyear + 1, 1, 3, hrend, mint, 0)
            else:
                Tinin_pre = datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                     hrinit, mint, 0)
                Tinin_post = datetime(iyear, imonth + 1, 3, hrend, mint, 0)

            for fn in range(0, len(files_3m)):
                # Get time array per file
                f1 = netcdf(CMIP6files[files_3m[fn]], 'r')
                f1o = f1.variables
                try:
                    Tininx_pre = d2i(Tinin_pre, f1o['time'], select='exact', calendar=f1o['time'].calendar)
                except:
                    Tininx_pre = 0
                try:
                    Tininx_post = d2i(Tinin_post, f1o['time'], select='exact', calendar=f1o['time'].calendar)
                    # if Tininx_pre == 0:
                    Tininx_post = Tininx_post + 1
                except:
                    Tininx_post = f1o['time'].shape[0]
                if fn == 0:
                    var_time = np.array(f1o['time'][Tininx_pre:Tininx_post])
                    var_ext  = np.array(f1o[CMIP6_ATMlab[v]][Tininx_pre:Tininx_post, :, :])
                else:
                    var_time = np.concatenate((var_time, np.array(f1o['time'][Tininx_pre:Tininx_post])))
                    var_ext  = np.concatenate((var_ext,
                                               np.array(f1o[CMIP6_ATMlab[v]][Tininx_pre:Tininx_post, :, :])))
                f1.close()

            if imonth == 1:
                Tnew_pre = date.toordinal(datetime(iyear - 1, 12, monthrange(iyear - 1, 12)[1] - 2,
                                                   hrinit, mint, 0)) - \
                           date.toordinal(datetime(Yorig, 1, 1))

                Tnew_post = date.toordinal(datetime(iyear, imonth + 1, 3, hrend, mint, 0)) - \
                            date.toordinal(datetime(Yorig, 1, 1))
            elif imonth == 12:
                Tnew_pre = date.toordinal(datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                                   hrinit, mint, 0)) - \
                           date.toordinal(datetime(Yorig, 1, 1))

                Tnew_post = date.toordinal(datetime(iyear + 1, 1, 3, hrend, mint, 0)) - \
                            date.toordinal(datetime(Yorig, 1, 1))
            else:
                Tnew_pre = date.toordinal(datetime(iyear, imonth - 1, monthrange(iyear, imonth - 1)[1] - 2,
                                                   hrinit, mint, 0)) - \
                           date.toordinal(datetime(Yorig, 1, 1))

                Tnew_post = date.toordinal(datetime(iyear, imonth + 1, 3, hrend, mint, 0)) - \
                            date.toordinal(datetime(Yorig, 1, 1))

            if CMIP6_ATMtres[v] == '3hr':
                if (imonth == 1) or (imonth == 3) or (imonth == 8):
                    var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :],
                                              var_ext[16:256, :, :], var_ext[-32:, :, :]))
                if (imonth == 4) or (imonth == 6):
                    var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :], var_ext[16:280, :, :]))
                if (imonth == 5) or (imonth == 7) or (imonth == 12) or (imonth == 10):
                    var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :], var_ext[16:288, :, :]))
                if (imonth == 9) or (imonth == 11):
                    var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :], var_ext[16:280, :, :]))
                if (imonth == 2):
                    if isleap(iyear) == 0:
                        var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :],
                                                  var_ext[16:240, :, :], var_ext[-24:, :, :]))
                    elif isleap(iyear) == 1:
                        var_ext = np.concatenate((var_ext[0:16, :, :], var_ext[8:16, :, :],
                                                  var_ext[16:248, :, :], var_ext[-24:, :, :]))

            if CMIP6_ATMtres[v] == 'day':
                if (imonth == 1) or (imonth == 3) or (imonth == 8):
                    var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0),
                                              var_ext[2:34, :, :], var_ext[-3:, :, :]))
                if (imonth == 4) or (imonth == 6):
                    var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0), var_ext[1:, :, :]))
                if (imonth == 5) or (imonth == 7) or (imonth == 12):
                    var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0), var_ext[1:, :, :]))
                if (imonth == 9) or (imonth == 10) or (imonth == 11):
                    var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0), var_ext[1:, :, :]))
                if (imonth == 2):
                    if isleap(iyear) == 0:
                        var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0),
                                                  var_ext[2:31, :, :], var_ext[-3:, :, :]))
                    elif isleap(iyear) == 1:
                        var_ext = np.concatenate((var_ext[0:1, :, :], np.expand_dims(var_ext[1, :, :], 0),
                                                  var_ext[2:32, :, :], var_ext[-3:, :, :]))

            Tnew = np.arange(Tnew_pre + tinit, Tnew_post + tendi, tincr)

            if len(Tnew) != var_ext.shape[0]:
                print('something')

            #
            # Open the atmospheric file for writing
            #

            ncatm = netcdf(atmname, 'a')
            ncatmo = ncatm.variables
            ncatmo[ERA5_ATMstr[v]][:] = var_ext
            ncatmo['time'][:] = Tnew
            ncatmo['lat'][:] = lat
            ncatmo['lon'][:] = lon

            #
            #  End loop on time
            #
            ncatm.close()
