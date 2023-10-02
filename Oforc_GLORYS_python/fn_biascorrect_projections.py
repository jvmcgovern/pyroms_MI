import numpy as np
import matplotlib.pyplot as plt
from datetime import date, datetime
import time
import ks_regress
import Xplin
import calc_surfaceSWI
import fnperiodic_distm
from numba import jit, prange


@jit(nopython=True, parallel=True)
def CDFgen(ssnl, d1, tol_dtclimp, nthref, seltclimh0, seltclimh, ns, Fhref, Xh, c):
    if ssnl == 1:
        selthref1 = np.argwhere(d1 <= tol_dtclimp)[:, 0]
    else:
        selthref1 = np.arange(0, nthref)
    selth1 = seltclimh0[selthref1]
    nh1 = len(selthref1)
    for j in prange(nh1):
        for k in prange(ns):
            Fhref[selthref1[j], k] = np.divide((np.sum(Xh[seltclimh, k] < Xh[selth1[j], k]) + 0.5), c)
    return selthref1, selth1, nh1, Fhref


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


class optJ:
    def __init__(self, Y, Xfun, nonnegb):
        self.Y = Y
        self.Xfun = Xfun
        self.nonnegb = nonnegb


class opto:
    def __init__(self, optJ):
        self.optJ = optJ


class optk:
    def __init__(self, kfunc, lr_order):
        self.kfunc = kfunc
        self.lr_order = lr_order


class optkres:
    def __init__(self, kfunc, lr_order):
        self.kfunc = kfunc
        self.lr_order = lr_order


class opt_subq:
    def __init__(self, model, logX, logY, nparb, bymonth, resmodel, doplot):
        self.model = model
        self.logX = logX
        self.logY = logY
        self.nparb = nparb
        self.bymonth = bymonth
        self.resmodel = resmodel
        self.doplot = doplot


class optSWI:
    def __init__(self, f_model, decl_model, use_eqtime, Q_model, cloud_model, ndays_year, use_tday_lag):
        self.f_model = f_model
        self.decl_model = decl_model
        self.use_eqtime = use_eqtime
        self.Q_model = Q_model
        self.cloud_model = cloud_model
        self.ndays_year = ndays_year
        self.use_tday_lag = use_tday_lag


def fn_biascorrect_projections(Xp, tdp, Xh, tdh, yearminref, yearmaxref, method, opt):
    # function [Xpc,out] = fn_biascorrect_projections(Xp,tdp,Xh,tdh,yearminref,yearmaxref,method,opt)

    # Computes bias-corrected projections (Xpc [ntimes*nseries]) given:
    #   Xp  = input projection time series [ntimes*nseries]
    #   tdp = input projection date numbers [ntimes,1]
    #   Xh  = hindcast time series [ntimesh*nseries]
    #   tdh = hindcast date numbers [ntimes,1]
    #   yearminref,yearmaxref = reference period min/max year
    #   method = bias correction method:
    #            0 = delta-change (e.g. Kay and Butenschon, 2018)
    #            1 = empirical quantile mapping (e.g. 'QM' method in Cannon et al., 2015)
    #            2 = repeated hindcast correction (RHC, e.g. Willems and Vrac, 2011)

    # Methods (0,1) start with the projection time series data (Xp) and aim to correct the statistical
    # variability level for consistency with the hindcast time series, either by correcting climatological
    # (e.g. monthly) np.mean values (for method=1, 'delta-change'), or by correcting the quantiles of variability
    # (for method=1, 'empirical quantile mapping'). Method 2 takes a different approach, starting with the
    # hindcast data, repeating it in blocks to fill the projection time points, and then perturbating these
    # repeated data using climatic trends in quantiles (e.g. decadal-scale changes) that are calculated from
    # the projection data. See text within function below for further advice about choice of method.

    # By default, reference climatologies and quantiles are defined on a monthly basis (opt.use_month = 1).
    # For other resolutions, set opt.use_month = 0 and provide a set of central climatological times opt.tclimc
    # as yeardays within range 0-365 (e.g. tclimc=0.5:365.5 for daily delta-change). We can also specify a
    # tolerance of +/- tol_dtclim (default = 0) about the central times tclimc, to enable more
    # robust calculations of  climatologies and quantile corrections. For example if we want to include
    # (December,January,February) data when calculating means/quantiles for January data we set opt.use_month = 1
    # and opt.tol_dtclim = 1

    # Input scalar opt.fractional determines whether a fractional (1) or absolute (0, default) correction is
    # applied to each series.

    # For method=2 (Repeated Hindcast Correction) we also need to input the target blocks of years
    # for which projections are required e.g. opt.yearpminv = 2000:10:2090; opt.yearpmaxv = 2009:10:2099.
    # Hindcast data from the reference period (yearminref to yearmaxref inclusive) will be inserted to fill
    # these blocks, and then corrected using calculated trends in quantiles (differences between the
    # projection blocks calculated from the projection model data). The last day of the reference period
    # will be duplicated to fill gaps resulting from different numbers of leap years per decade.

    # The function will also add 'subdelta variability' (variability on a timescale shorter than that
    # of the delta changes, i.e. the temporal resolution of the projection data), for example to achieve
    # hourly projections given daily averages from an ESM. For this we must set opt.use_subdelta = 1 (def=0)
    # and also supply fine-resolution hindcast time series data, with corresponding datenumbers
    # (opt.Xhf [ntimesf*nseries], opt.tdhf [ntimesf*1]).

    # The subdelta variability is imported from the hindcast as a repeated multiannual chunk specified by (
    # opt.yearminsub,opt.yearmaxsub), and is used to fill data for prediction chunks specified by (opt.yearpminv,
    # opt.yearpminv), e.g. opt.yearpminv = 2000:10:2090; opt.yearpmaxv = 2009:10:2099. The subdelta variability is
    # imposed either as additive corrections (opt.fractional_subdelta=0, default) or as multiplicative corrections (
    # opt.fractional_subdelta=1) to the deltascale averages, assuming the deltascale averages are either monthly or
    # daily averages. This variability may be further corrected to improve consistency with the projected deltascale
    # average values, either by inflating/deflating the subdelta variance (opt.correct_substd=1) or by correcting the
    # subdelta quantiles (opt.correct_subquantiles=1). By default, the subdelta variability is defined wrt a linear
    # interpolation from the deltascale averages to the finescale time points (opt.remove_deltascale_var=1, default),
    # but in some cases it may make more sense to define it wrt the deltascale averages as piecewise constant values
    # (e.g. for imposing hourly variability on daily average shortwave irradiance as piecewise constant values we set
    # opt.remove_deltascale_var=0).

    # If method=1 (Empirical Quantile Mapping), it is likely that the subdelta variability will need to be
    # adjusted for consistency with the projected deltascale averages (opt.correct_substd=1 or
    # opt.correct_subquantiles=1).
    # In such cases, an options substructure opt.opt_subq is used to specify modelling options for correcting
    # this subdelta variability. By default, a simple two-parameter linear regression model is a applied to
    # model the standard deviation or deviation quantiles as a function of deltascale mean values, but several
    # other options are available (consult code for details). If correcting deviation quantiles
    # (opt.correct_subquantiles=1)
    # it is also possible to model the residual variance of the deviation quantiles wrt the fitted model values,
    # using an additional 'residual model' (see code for details).

    # As an alternative to correcting the subdelta variability, the variability can be selected from the
    # hindcast deltascale time period that best matches each of the projected deltascale means
    # (opt.match_subdelta_deltamean=1).
    # However, this approach can run into problems because i) sufficiently-similar deltascale periods might not be found
    # in the hindcast dataset, when considering all spatial locations together (to preserve realistic
    # spatial variability), or  ii) the subselection procedure may generate large discontinuities between deltascale
    # periods at the subdelta time scale
    # (e.g. between 2300 and 0000 the next day).

    # If method=2 (Repeated Hindcast Correction), it is more likely that the subdelta variability can be imposed
    # without further correction, because the the projected changes between decades for each quantile should be
    # relatively small (this is one major advantage of the RHC method).

    # Corrected projections are capped at lower and/or upper limits if given as input (opt.XpcL,XpcH,XpcfL,XpcfH).
    # Scalars can be used here as well as arrays of size(Xp) (e.g. input opt.XpcL = 0 for nonnegative variables).
    # In addition, for surface shortwave irradiance, an upper limit can be calculated (and imposed) internally
    # using astronomical formulae for the clear-sky surface shortwave irradiance (see code for details).

    # By default, deltascale averages are recalculated from the finescale (subdelta) time series
    # (opt.recalculate_Xpc=1, default).
    # This can take some time so might be switched off for testing purposes.
    # Uses:
    #   fnperiodic_distm.m (calculate separations for periodic variables such as month/day-of-year)
    #   Xplin.m (simple linear interpolation with optional constant extrapolation)
    #   calc_surfaceSWI.m (calculates maximum daily/hourly surface shortwave irradiance,
    #                      for use as upper limits on projected values)
    #   ks_regress.m (used for fitting kernel-smoothing regression models if correcting subdelta variance/quantiles)
    #   J_fit_partially_linear_model.m (used for fitting partially-linear regression models
    #                                   if correcting subdelta variance/quantiles)
    #   optimize_theta.m (used for fitting partially-linear regression models if correcting subdelta variance/quantiles)
    #
    #
    # Phil Wallhead 10/05/2023
    # Joe McGovern, conversion to Python 26/09/2023
    #
    #
    # Some General Advice About Choice Of Method
    # ------------------------------------------
    #
    # In general we would recommend either:
    # method=1 (Empirical Quantile Mapping) or
    # method=2 (Repeated Hindcast Correction)
    # as these will more thoroughly correct the level of variability in the projections.
    # The simple delta-change approach
    # (method=0) may clearly under-correct the projections in cases where the main temporal variability signal is not
    # seasonal (e.g. with surface wind speeds), but might yet be optimal in cases where:
    # i) seasonal variability dominates, and
    # ii) the hindcast data are too few to allow robust estimation of quantiles
    # (but does allow reasonable mean estimates).

    # Between methods (1,2), method 1 (Empirical Quantile Mapping) might be preferred in cases where:
    # i) Subdelta variability is not a major concern;
    # ii) It is preferable to have projections that are closer to the original projections in character;
    # iii) Artifacts due to interpolation of original projection data (leading to spatial oversmoothing),
    # potential loss of multivariate correlation structure due to the univariate implementation, and
    # potential distortion of climatic trends due to quantile correction, are not such of a concern.
    # Note that more elaborate implementations of quantile mapping have been developed to address the issues in
    # iii) separately, though to our knowledge not yet to address all together.
    # iv) It is desired to preserve the contribution of internal variability to the uncertainty in climatic change
    # between decades.

    # Alternatively, method 2 (Repeated Hindcast Correction) might be preferred in cases where:
    # i) Subdelta variability is significant and strongly dependent on deltascale mean values. RHC is preferable
    # over EQM is such cases because the perturbations to deltascale mean values (e.g. monthly means, or daily means)
    # that are imposed by RHC due to interdecadal/climatic changes are often smaller than those imposed
    # by EQM to correct
    # the statistical variability level of the original projection data. Therefore it is more likely that the subdelta
    # variability can be imported from the hindcast dataset with no corrections (or more minor corrections) to account
    # for changes in the deltascale mean values.
    # ii) It is preferable to have projections that are closer to the hindcast data in character. The RHC method
    # puts less 'faith' in the model used to produce the original projection data, in the sense that only certain
    # aspects of the original model data --- i.e. the climatic/decadal-scale variability --- is used in the final
    # corrected projections.
    # iii) Artifacts due to spatial interpolation, inaccuracy in the multivariate (cross)-correlation structure of the
    # original projection data, and loss of cross-correlations or distortion of climatic trends due to the simple EQM
    # approach provided herein, are a major cause for concern. The RHC method preserves the full multivariate
    # spatio-temporal variability and (cross) correlation structure of the hindcast data, at the cost of introducing
    # some moderate temporal discontinuity between repetition blocks (in our experience not very noticeable),
    # and some decadal-scale
    # periodicity in the corrected projections that may be considered artifactual/spurious.
    # iv) Loss of uncertainty due to varying internal variability between (multi)decadal repetition block is not a major
    # concern. This may be the case when the corrected projections are to be used for comparing two (multi)decadal
    # periods
    # (present-day vs. future). In such cases, the fact that the internal year-to-year variation is identical
    # in both periods may actually be an advantage, allowing a more precise estimate of climatic change
    # (i.e. closer to what would be obtained by averaging over many realizations of the internal variability).

    ntp, ns = Xp.shape
    out = []

    # if len(varargin) < 8:
    #     opt = []

    if hasattr(opt, 'verbose'):
        verbose = opt.verbose
    else:
        verbose = 1

    if hasattr(opt, 'seasonal'):
        seasonal = opt.seasonal
    else:
        seasonal = 1

    if hasattr(opt, 'use_month'):
        use_month = opt.use_month
    else:
        use_month = 1

    # 1  to use month-of-year (1-12) rather than yrday (day-of-year 0-365) to compute climatological statistics
    # and deltas.

    if hasattr(opt, 'tclimc'):
        tclimc = opt.tclimc
    else:
        tclimc = np.arange(0, 12)

    # Climatological time centres used to calculate climatological (means,deltas) and quantiles
    # (def = 1:12 for monthly delta change with use_month=1).
    if hasattr(opt, 'ndays_year'):
        ndays_year = opt.ndays_year
    else:
        ndays_year = 365.25

    # Best approach here seems to be to use ndays_year = 365.25; if we use ndays_year = 365 then we
    # will combine data from Dec 31st and Jan 1st when the preceding year is a leap year.
    if hasattr(opt, 'tol_dtclim'):
        tol_dtclim = opt.tol_dtclim
    else:
        tol_dtclim = 0

    # Tolerance level for climatological time used in computing climatological statistics
    # (def = 0 for simple monthly delta change).
    # Using tol_dtclim>0 may allow a more robust definition of climatology,
    # e.g. if the no. of reference years (yearmaxref-yearminref+1) is limited.
    # E.g. if aiming for a daily climatology (tclimc=(0.5:365.5)),
    # it is probably a good idea to allow e.g. tol_dtclim = 5 (i.e. +/- 5 days tolerance)
    # in order to allow robust climatological statistics and consequent deltas.
    if hasattr(opt, 'tol_dtclimp'):
        tol_dtclimp = opt.tol_dtclimp
    else:
        tol_dtclimp = 0

    # Tolerance level for climatological time used in correcting projections (def = 0 for simple monthly delta change).
    # This may be different to tol_dtclim because tol_dtclim may be expanded to allow more robust
    # climatological statistics.
    if hasattr(opt, 'fractional'):
        fractional = opt.fractional
    else:
        fractional = 0

    if hasattr(opt, 'minX'):
        minX = opt.minX
    else:
        minX = []

    # Reference minimum value to avoid potential division by ~0 in fractional methods.
    if hasattr(opt, 'allow_block_lag'):
        allow_block_lag = opt.allow_block_lag
    else:
        allow_block_lag = 1

    if Xh.shape[1] > 0 and ns == 1:
        # In case of multiple hindcast series and one projection series,
        # apply the single projection series to all hindcast climatologies
        ns = Xh.shape[1]
        Xp = Xp * np.ones((1, ns))

    if method == 0:
        if hasattr(opt, 'correct_iavar'):
            correct_iavar = opt.correct_iavar
        else:
            correct_iavar = 0
        if hasattr(opt, 'dyear_iavar'):
            dyear_iavar = opt.dyear_iavar
        else:
            dyear_iavar = 10
        if hasattr(opt, 'use_coef_iavar'):
            use_coef_iavar = opt.use_coef_iavar
        else:
            use_coef_iavar = 1
        # use_coef_iavar=1: Correct the interannual coefficient of variation rather than standard deviation,
        # only for fractional deltas (def)

    # use_coef_iavar=2: Correct the interannual coefficient of variation rather than standard deviation
    # for all deltas

    if method == 1:
        if hasattr(opt, 'qsel'):
            qsel = opt.qsel
        else:
            qsel = []
        qsel = qsel

        # Optional set of fixed quantiles to interpolate to (e.g. 0.005:0.01:0.995) ---
        # otherwise the sorted data are used.
        # Experience to date suggests that it is usually not helpful to interpolate to fixed quantiles
        # --- better to leave qsel empty.

    if method == 2:
        if hasattr(opt, 'qsel'):
            qsel = opt.qsel
        else:
            qsel = np.arange(0.05, 0.95 + 0.1, 0.1)
        qsel = qsel

    if hasattr(opt, 'XpcL'):
        XpcL = opt.XpcL
    else:
        XpcL = []

    if hasattr(opt, 'XpcH'):
        XpcH = opt.XpcH
    else:
        XpcH = []

    if hasattr(opt, 'use_XpcH_SWI'):
        use_XpcH_SWI = opt.use_XpcH_SWI
    else:
        use_XpcH_SWI = 0

    # If opt.use_XpcH_SWI=1, maximum daily-average values of shortwave irradiance are calculated using calc_surfaceSWI.m
    # and these are used to impose an upper limit on Xpc
    if use_XpcH_SWI == 1:
        if hasattr(opt, 'latp'):
            latp = opt.latp
        else:
            latp = []
        if hasattr(opt, 'lonp'):
            lonp = opt.lonp
        else:
            lonp = []
        if hasattr(opt, 'optSWI'):
            optSWI = opt.optSWI
        else:
            optSWI = []
        if len(optSWI) == 0:
            optSWI.f_model = 1
            optSWI.decl_model = 1
            optSWI.use_eqtime = 1
            optSWI.Q_model = 0
            optSWI.cloud_model = 0
            optSWI.ndays_year = 365.2422
            optSWI.use_tday_lag = 1

            # Default options for calc_surfaceSWI.m: These have been found to give best performance (minimum overshoot)
            # in tests using hourly ERA5 data for ROHO800 model (western Norway).
        if hasattr(opt, 'use_SWI_hav'):
            use_SWI_hav = opt.use_SWI_hav
        else:
            use_SWI_hav = 1
        # 1 (def) to use hourly averages rather than instantaneous values in calculation of hourly data
        # (consistent with ERA5).

    if hasattr(opt, 'yearpminv'):
        yearpminv = opt.yearpminv
    else:
        yearpminv = []

    if hasattr(opt, 'yearpmaxv'):
        yearpmaxv = opt.yearpmaxv
    else:
        yearpmaxv = []

    # Optional input of set of prediction year periods over which the subdelta variability will be resampled
    # and/or deltascale variability will be repeated (if method=2).

    if hasattr(opt, 'legacy'):
        legacy = opt.legacy
    else:
        legacy = 0

    if legacy == 1:
        allow_block_lag = 0

    if hasattr(opt, 'use_subdelta'):
        use_subdelta = opt.use_subdelta
    else:
        use_subdelta = 0

    # 1 to use subdelta variability from input finescale hindcast (opt.Xhf,opt,tdhf) to correct the projections
    if use_subdelta == 1:
        if hasattr(opt, 'tdpfmin'):
            tdpfmin = opt.tdpfmin
        else:
            tdpfmin = int(np.floor(np.amin(tdp)))
        if hasattr(opt, 'tdpfmax'):
            tdpfmax = opt.tdpfmax
        else:
            tdpfmax = np.ceil(np.amax(tdp))
        if hasattr(opt, 'Xhf'):
            Xhf = opt.Xhf
        else:
            Xhf = []
        if hasattr(opt, 'tdhf'):
            tdhf = opt.tdhf
        else:
            tdhf = []
        if hasattr(opt, 'yearminsub'):
            yearminsub = opt.yearminsub
        else:
            yearminsub = yearminref
        if hasattr(opt, 'yearmaxsub'):
            yearmaxsub = opt.yearmaxsub
        else:
            yearmaxsub = yearmaxref

        # Optional input of minimum/maximum years (yearminsub,yearmaxsub) to sample subdelta variability from
        # -- by default set equal to (yearminsub,yearmaxsub).
        # fractional_subdelta = 1 to treat the subdelta variability as a fractional perturbation
        # (better for non-negative variables)
        if hasattr(opt, 'fractional_subdelta'):
            fractional_subdelta = opt.fractional_subdelta
        else:
            fractional_subdelta = 0

        if hasattr(opt, 'minXsub'):
            minXsub = opt.minXsub
        else:
            minXsub = 0

        # Reference minimum value to avoid potential division by ~0 in fractional subdelta methods.
        if hasattr(opt, 'loop_over_projection_times'):
            loop_over_projection_times = opt.loop_over_projection_times
        else:
            loop_over_projection_times = 0

        # Using matlab's np.kron.m below gives identical results to the loop calculation, but it ~800x faster
        if hasattr(opt, 'correct_substd'):
            correct_substd = opt.correct_substd
        else:
            correct_substd = 0

        # opt.correct_subquantiles = 1 to model the subdelta-scale quantiles
        # (e.g. hourly) as functions of delta-scale means (e.g. daily means)
        if hasattr(opt, 'correct_subquantiles'):
            correct_subquantiles = opt.correct_subquantiles
        else:
            correct_subquantiles = 0

        if hasattr(opt, 'remove_deltascale_var'):
            remove_deltascale_var = opt.remove_deltascale_var
        else:
            remove_deltascale_var = 1

        if (correct_substd == 1) or (correct_subquantiles == 1):
            if hasattr(opt, 'opt_subq'):
                opt_subq = opt.opt_subq
            else:
                opt_subq = []
            # Options structure to specify the models for subdelta quantiles
            # Defaults:
            if opt_subq == []:
                if verbose > 0:
                    print('Using default options for models to correct subdelta quantiles')
                    opt_subq.model = 'polynomial'
                    opt_subq.logX = 0
                    opt_subq.logY = 0
                    opt_subq.nparb = 2
                    opt_subq.bymonth = 1
                    opt_subq.resmodel = []
                    opt_subq.doplot = 2

            if not hasattr(opt_subq, 'logX'):
                opt_subq.logX = 0

            if not hasattr(opt_subq, 'logY'):
                opt_subq.logY = 0

            if not hasattr(opt_subq, 'nonnegb'):
                opt_subq.nonnegb = 0

            if not hasattr(opt_subq, 'nr'):
                opt_subq.nr = 1

            if not hasattr(opt_subq, 'logXres'):
                opt_subq.logXres = 0

            if not hasattr(opt_subq, 'logYres'):
                opt_subq.logYres = 0

            if not hasattr(opt_subq, 'nonnegbres'):
                opt_subq.nonnegbres = 0

            if not hasattr(opt_subq, 'nrres'):
                opt_subq.nrres = 1

            if correct_subquantiles == 1:
                if not str(opt_subq.model) == str('polynomial'):
                    opt_subq.logX = 0
                    opt_subq.logY = 0
                if not str(opt_subq.resmodel) == str('polynomial'):
                    opt_subq.logXres = 0
                    opt_subq.logYres = 0
                # Log transforms only applicable to polynomial models
                if str(opt_subq.resmodel) == str('gaussian'):
                    opt_subq.nparbres = 2

        if hasattr(opt, 'match_subdelta_deltamean'):
            match_subdelta_deltamean = opt.match_subdelta_deltamean
        else:
            match_subdelta_deltamean = 0
        # 1 to select the subdelta variability from the hindcast day/month that best matches the
        # projected deltascale mean

        if (match_subdelta_deltamean == 1) or (correct_subquantiles == 1):
            loop_over_projection_times = 1
        if hasattr(opt, 'match_subdelta_hourly'):
            match_subdelta_hourly = opt.match_subdelta_hourly
        else:
            match_subdelta_hourly = 1

        if match_subdelta_deltamean == 0:
            match_subdelta_hourly = 0

        if hasattr(opt, 'XpcfL'):
            XpcfL = opt.XpcfL
        else:
            XpcfL = []

        if hasattr(opt, 'XpcfH'):
            XpcfH = opt.XpcfH
        else:
            XpcfH = []

        if hasattr(opt, 'use_XpcfH_SWI'):
            use_XpcfH_SWI = opt.use_XpcfH_SWI
        else:
            use_XpcfH_SWI = 0

        # If opt.use_XpcfH_SWI>0, maximum hourly values of shortwave irradiance are calculated using calc_surfaceSWI.m
        # and these are used to impose an upper limit on Xpcf. Then if:
        # opt.use_XpcfH_SWI = 1: All hourly projections are limited by the maximum (clear-sky) values.
        # opt.use_XpcfH_SWI = 2: Only hourly projections >=opt.frcrit_capSWI*daily maximum are limited by the
        #                        clear-sky values.
        if use_XpcfH_SWI == 2:
            if hasattr(opt, 'frcrit_capSWI'):
                frcrit_capSWI = opt.frcrit_capSWI
            else:
                frcrit_capSWI = 0.25

        if hasattr(opt, 'recalculate_Xpc'):
            recalculate_Xpc = opt.recalculate_Xpc
        else:
            recalculate_Xpc = 1

        # opt.recalculate_Xpc = 1: Recalculate the daily average projections using the
        #                         hourly projections after limits imposed.
        if hasattr(opt, 'subdelta_by_blocks'):
            subdelta_by_blocks = opt.subdelta_by_blocks
        else:
            subdelta_by_blocks = 0

        # 1 to divide subdelta variability calculation into latitude-longitude blocks, to obtain better matches to
        #  hindcast daily means (match_subdelta_deltamean=1).
        #  This is probably not a good idea unless the target variable has very limited spatial correlation between
        #  grid cells, since the blockwise calculation will not preserve spatial correlations across block boundaries.

        if subdelta_by_blocks == 1:
            if hasattr(opt, 'latp'):
                latp = opt.latp
            else:
                latp = []

            if hasattr(opt, 'lonp'):
                lonp = opt.lonp
            else:
                lonp = []

            if hasattr(opt, 'latL_blocks'):
                latL_blocks = opt.latL_blocks
            else:
                latL_blocks = []

            if hasattr(opt, 'latH_blocks'):
                latH_blocks = opt.latH_blocks
            else:
                latH_blocks = []

            if hasattr(opt, 'lonL_blocks'):
                lonL_blocks = opt.lonL_blocks
            else:
                lonL_blocks = []
            if hasattr(opt, 'lonH_blocks'):
                lonH_blocks = opt.lonH_blocks
            else:
                lonH_blocks = []

            latL_blocks1 = np.kron(np.transpose(latL_blocks), np.ones((1, len(lonL_blocks))))
            latH_blocks1 = np.kron(np.transpose(latH_blocks), np.ones((1, len(lonH_blocks))))
            lonL_blocks1 = np.kron(np.ones((1, len(latL_blocks))), np.transpose(lonL_blocks))
            lonH_blocks1 = np.kron(np.ones((1, len(latH_blocks))), np.transpose(lonH_blocks))
            nblocks = len(latL_blocks1)
        else:
            nblocks = 1
    # get yr, mth, day from tdp
    yearp = np.zeros_like(tdp, dtype=int)
    monthp = np.zeros_like(tdp, dtype=int)
    dayp = np.zeros_like(tdp, dtype=int)
    yrdayp = np.zeros_like(tdp, dtype=int)
    for tp in np.arange(0, len(tdp)):
        tpid = datetime.strptime(str(np.datetime64('1990-01-01') + tdp[tp].astype('timedelta64[D]')),
                                 '%Y-%m-%d')
        yearp[tp] = tpid.year
        monthp[tp] = tpid.month
        dayp[tp] = tpid.day
        # get day of the year
        yrdayp[tp] = date.toordinal(date(yearp[tp], monthp[tp],
                                         dayp[tp])) - date.toordinal(date(yearp[tp] - 1, 12, 31))

    if seasonal == 1:
        ntclimc = len(tclimc)
    else:
        ntclimc = 1

    if minX == []:
        minX = (np.amin(np.concatenate((Xp[:], Xh[:]))) -
                1 * (np.amax(np.concatenate((Xp[:], Xh[:]))) - np.amin(np.concatenate((Xp[:], Xh[:])))))
        # Set reference minimum to minimum minus the range --- this will avoid potential division by ~0,
        # and also avoids inaccuracy due to small variations very far from zero causing correction factors
        # very close to 1 (e.g. for pressure).
        # NOTE: minX is only used for correction when fractional = 1.
        # NOTE: Since the ratios are defined by e.g. (qXh-minX)/(qXp-minX), the above minX will restrict
        # the possible range of ratio values between 0.5 and 2.

        # out.minX = minX
    if method == 0:
        methodstr = 'delta-change'

    if method == 1:
        methodstr = 'empirical quantile mapping'

    if method == 2:
        methodstr = 'repeated hindcast correction'

    # if verbose > 0:
    #     print(np.array(['Calculating bias-corrected projections using ',methodstr,' method']))
    #     tic

    if (method == 0) or (method == 1):
        # Calculate projection climatologies / reference quantiles
        nseltclim = np.nan * np.ones((ntclimc, 1))
        if method == 0:
            Xpclim = np.nan * np.ones((ntclimc, ns))
            if correct_iavar == 1:
                Xpstdclim = Xpclim
        if method == 1:
            qXpref = np.zeros(ntclimc, 1)
        for i in np.arange(0, ntclimc + 1).reshape(-1):
            if seasonal == 1:
                # Calculate temporal separations in a periodic sense using fnperiodic_distm.m
                if use_month == 1:
                    dist1 = fnperiodic_distm.fnperiodic_distm(monthp, tclimc[i], 12)
                else:
                    dist1 = fnperiodic_distm.fnperiodic_distm(yrdayp, tclimc[i], ndays_year)
                # Choose all data within tolerance tol_dtclim of climatological time centre tclimc[i]
                selt = np.argwhere(dist1 <= tol_dtclim)
            else:
                selt = np.arange(0, ntp)
            seltclim = selt[yearp[selt] >= yearminref & yearp[selt] <= yearmaxref]
            nseltclim[i] = len(seltclim)
            if method == 0:
                Xpclim[i, :] = np.mean(Xp[seltclim, :])
                if correct_iavar == 1:
                    X1 = np.array([np.ones((nseltclim[i], 1)), yearp[seltclim]])
                    b1 = np.linalg.solve(X1, Xp[seltclim, :])
                    Xphat1 = X1 * b1
                    r1 = Xp[seltclim, :] - Xphat1
                    Xpstdclim[i, :] = np.std(r1)
                    if use_coef_iavar == 1 and fractional == 1:
                        Xpstdclim[i, :] = np.std(r1 / Xphat1)
                    if use_coef_iavar == 2:
                        Xpstdclim[i, :] = np.std(r1 / Xphat1)
            else:
                if method == 1:
                    if not len(qsel) == 0:
                        qXpref[i] = np.quantile(Xp[seltclim, :], qsel)
                    else:
                        qXpref[i] = np.sort(Xp[seltclim, :])
                    # Interpolate to fixed set of quantiles if qsel is specified, otherwise use the sorted values.
        out.nseltclim = nseltclim
        if method == 0:
            out.Xpclim = Xpclim
            if correct_iavar == 1:
                out.Xpstdclim = Xpstdclim
        if method == 1:
            out.qXpref = qXpref

    if method == 2:
        # Calculate projection model trends in quantiles between repeated hindcast periods
        nperiods = len(yearpminv)
        nq = len(qsel)
        qXp = np.nan * np.ones((ntclimc, nperiods, nq, ns))
        for i in np.arange(0, ntclimc).reshape(-1):
            if seasonal == 1:
                # Calculate temporal separations in a periodic sense using fnperiodic_distm.m
                if use_month == 1:
                    dist1 = fnperiodic_distm.fnperiodic_distm(monthp, tclimc[i], 12)
                else:
                    dist1 = fnperiodic_distm.fnperiodic_distm(yrdayp, tclimc[i], ndays_year)
                # Choose all data within tolerance tol_dtclim of climatological time centre tclimc[i]
                selt0 = np.argwhere(dist1 <= tol_dtclim)[:, 0]
            else:
                selt0 = np.arange(0, ntp)
            for j in np.arange(0, nperiods).reshape(-1):
                selt = selt0[np.argwhere(np.logical_and((yearp[selt0] >= yearpminv[j]),
                                                        (yearp[selt0] <= yearpmaxv[j])))[:, 0]]
                qXp[i, j, 0:nq, 0:ns] = np.quantile(Xp[selt, :], qsel, axis=0)
        # out.qXp = qXp

    # Apply the delta changes / quantile corrections to generate bias-corrected projections.
    # yearh,monthh,__ = datevec(tdh)
    yearh = np.zeros_like(tdh, dtype=int)
    monthh = np.zeros_like(tdh, dtype=int)
    dayh = np.zeros_like(tdh, dtype=int)
    yrdayh = np.zeros_like(tdh, dtype=int)
    for th in np.arange(0, len(tdh)):
        thid = datetime.strptime(str(np.datetime64('1990-01-01') + tdh[th].astype('timedelta64[D]')), '%Y-%m-%d')
        yearh[th] = thid.year
        monthh[th] = thid.month
        dayh[th] = thid.day
        # get day of the year
        yrdayh[th] = date.toordinal(date(yearh[th], monthh[th],
                                         dayh[th])) - date.toordinal(date(yearh[th] - 1, 12, 31))

    seltclimh0 = np.argwhere(np.logical_and(yearh >= yearminref, yearh <= yearmaxref))[:, 0]

    nthref = len(seltclimh0)
    if np.all(np.diff(yearh) == 1):
        delta_timescale = 'yearly'
    else:
        if np.logical_or(np.all(np.diff(monthh) == 1), np.all(np.diff(monthh) == - 11)):
            delta_timescale = 'monthly'
        else:
            if np.all(np.diff(tdh) == 1):
                delta_timescale = 'daily'
            else:
                delta_timescale = []

    # Preallocation for hindcast statistic arrays
    nseltclimh = np.nan * np.ones((ntclimc, 1))[:, 0]
    Xpc = np.nan * np.ones((ntp, ns))
    if method == 0:
        Xhclim = np.nan * np.ones((ntclimc, ns))
        if correct_iavar == 1:
            Xhstdclim = Xhclim
    elif method == 1:
        qXhref = np.zeros(ntclimc, 1)
    elif method == 2:
        Fhref = np.nan * np.ones((nthref, ns))

    # First we calculate the climatology (or reference quantiles) for the hindcast model.
    # For simple delta change, the corrected projections (Xpc) are also calculated within the same loop.
    for i in np.arange(0, ntclimc).reshape(-1):
        if seasonal == 1:
            if use_month == 1:
                dist1 = fnperiodic_distm.fnperiodic_distm(monthh[seltclimh0], tclimc[i], 12)
                dist1p = fnperiodic_distm.fnperiodic_distm(monthp, tclimc[i], 12)
            else:
                dist1 = fnperiodic_distm.fnperiodic_distm(yrdayh[seltclimh0], tclimc[i], ndays_year)
                dist1p = fnperiodic_distm.fnperiodic_distm(yrdayp, tclimc[i], ndays_year)

            # Choose all data within tolerance tol_dtclim of climatological time centre tclimc[i]
            seltclimh = seltclimh0[np.argwhere(dist1 <= tol_dtclim)[:, 0]]
            # Choose all data within tolerance tol_dtclimp of climatological time centre tclimc[i]
            seltp = np.argwhere(dist1p <= tol_dtclimp)[:, 0]
        else:
            seltclimh = seltclimh0
            seltp = np.arange(0, ntp)
        nseltclimh[i] = len(seltclimh)
        if method == 0:
            Xhclim[i, :] = np.mean(Xh[seltclimh, :])
            nseltp = len(seltp)
            # Here we calculate simple delta-change projections with no correction of interannual variance or quantiles.
            if fractional == 1:
                Xpc[seltp, :] = minX + np.multiply((Xp[seltp, :] - minX),
                                                   (np.ones((nseltp, 1)) * Xhclim[i, :] - minX)) / (
                                        np.ones((nseltp, 1)) * Xpclim[i, :] - minX)
            else:
                Xpc[seltp, :] = np.ones((nseltp, 1)) * Xhclim[i, :] + Xp[seltp, :] - np.ones((nseltp, 1)) * Xpclim[i, :]
            if correct_iavar == 1:
                Xh1 = np.array([np.ones((nseltclimh[i], 1)), yearh[seltclimh]])
                bh1 = np.linalg.solve(Xh1, Xh[seltclimh, :])
                Xhhat1 = Xh1 * bh1
                rh1 = Xh[seltclimh, :] - Xhhat1
                Xhstdclim[i, :] = np.std(rh1)
                if use_coef_iavar == 1 and fractional == 1:
                    Xhstdclim[i, :] = np.std(rh1 / Xhhat1)
                if use_coef_iavar == 2:
                    Xhstdclim[i, :] = np.std(rh1 / Xhhat1)
        elif method == 1:
            if not len(qsel) == 0:
                qXhref[i] = np.quantile(Xh[seltclimh, :], qsel)
            else:
                if nseltclimh[i] < nseltclim[i]:
                    seltclimhc = np.concatenate((seltclimh,
                                                 seltclimh[np.arange(0,
                                                                     (nseltclim[i] - nseltclimh[i]))]),axis=0)
                else:
                    if nseltclimh[i] > nseltclim[i]:
                        seltclimhc = seltclimh[np.arange(0, nseltclim[i])]
                    else:
                        seltclimhc = seltclimh
                qXhref[i] = np.sort(Xh[seltclimhc, :])
        elif method == 2:
            # Here we need to record the cdf value F for each time point in the repeated hindcast series
            sta = time.time()
            selthref1, selth1, nh1, Fhref = CDFgen(seasonal, dist1, tol_dtclimp, nthref, seltclimh0, seltclimh, ns,
                                                   Fhref, Xh, nseltclimh[i])
            dne = time.time()
            print(dne - sta)

            # NOTE: This is consistent with how quantiles are calculated within quantile.m.
    #      There is no need for interpolation because the values Xh(seltclimh1(j),k) will
    #      always correspond exactly to one of the values in Xh(seltclimh,k).
    # NOTE: For leap years, the final year-day (yrday=365.5) will be compared to a subset
    #      seltclimh that is equally large as for the other year days, as long as the
    #      tolerance tol_dtclim>0.

    # out.nseltclimh = nseltclimh
    # if method == 0:
    #     out.Xhclim = Xhclim
    #
    # if method == 1:
    #     out.qXhref = qXhref
    #
    # if method == 2:
    #     out.Fhref = Fhref

    # Next we apply the correction factors to the projections (if not applied already for simple delta change)
    if (method == 0) and (correct_iavar == 1):
        # Here we inflate the residuals from linear regression of ndivdec decadal/multiannual periods,
        # using the ratio of standard deviations of linear trend residuals during the baseline period.
        yearpmin = np.amin(yearp)
        yearpmax = np.amax(yearp)
        yearpspan = yearpmax - yearpmin
        ndivdec = np.round(yearpspan / dyear_iavar)
        Rstdclim = Xhstdclim / Xpstdclim
        Xpco = Xpc
        for i in np.arange(0, ntclimc).reshape(-1):
            if seasonal == 1:
                if use_month == 1:
                    dist1p = fnperiodic_distm.fnperiodic_distm(monthp, tclimc[i], 12)
                else:
                    dist1p = fnperiodic_distm.fnperiodic_distm(yrdayp, tclimc[i], ndays_year)
                # Choose all data within tolerance tol_dtclimp of climatological time centre tclimc[i]
                seltp = np.argwhere(dist1p <= tol_dtclimp)[:, 0]
            else:
                seltp = np.arange(0, ntp)
            for j in np.arange(0, ndivdec).reshape(-1):
                if j == ndivdec:
                    seltp2 = seltp[np.argwhere(yearp[seltp] >= yearpmin + (j - 1) * dyear_iavar)[:, 0]]
                    pass
                else:
                    seltp2 = seltp[np.argwhere(np.logical_and(yearp[seltp] >= yearpmin + (j - 1) * dyear_iavar,
                                                              yearp[seltp] < yearpmin + j * dyear_iavar))[:, 0]]
                    pass
                nseltp2 = len(seltp2)
                # Xp1 = np.array([np.ones((nseltp2, 1)), tdp[seltp2]])
                Xp1 = np.concatenate((np.ones((nseltp2, 1)), tdp[seltp2]), axis=1)
                bp1 = np.linalg.solve(Xp1, Xp[seltp2, :])
                Xphat1 = Xp1 * bp1
                rp1 = Xpco[seltp2, :] - Xphat1
                Xpc[seltp2, :] = Xphat1 + np.multiply((np.ones((nseltp2, 1)) * Rstdclim[i, :]), rp1)
                if (use_coef_iavar == 1) and (fractional == 1):
                    Xpc[seltp2, :] = np.multiply(Xphat1, (1 + np.multiply(Rstdclim[i, :], rp1) / Xphat1))
                if use_coef_iavar == 2:
                    Xpc[seltp2, :] = np.multiply(Xphat1, (1 + np.multiply(Rstdclim[i, :], rp1) / Xphat1))
        out.Xhstdclim = Xhstdclim
        out.Rstdclim = Rstdclim
    elif method == 1:
        # Here we adjust the projections using factors derived by interpolating the ratios/differences
        # of quantiles during the reference period. This is an implementation of simple empirical
        # quantile mapping (QM, see Cannon et al., 2015) which does not preserve climatic trends,
        # and which potentially lead to unrealistic distortion of climate change signals (e.g. Cannon et al., 2015).
        # However, such distortion effects should be minimal as long as the projected climatic trends are weak
        # relative to the full standard deviation of ESM variability during the reference period.
        Rqref = np.zeros(ntclimc, 1)
        for i in np.arange(0, ntclimc).reshape(-1):
            if seasonal == 1:
                if use_month == 1:
                    dist1p = fnperiodic_distm.fnperiodic_distm(monthp, tclimc[i], 12)
                else:
                    dist1p = fnperiodic_distm.fnperiodic_distm(yrdayp, tclimc[i], ndays_year)
                # Choose all data within tolerance tol_dtclimp of climatological time centre tclimc[i]
                seltp = np.argwhere(dist1p <= tol_dtclimp)[:, 0]
            else:
                seltp = np.arange(0, ntp)
            qXpref1 = qXpref[i]
            qXhref1 = qXhref[i]
            Rqref1 = (qXhref1 - minX) / (qXpref1 - minX)
            Rqref[i] = Rqref1
            for j in np.arange(0, ns).reshape(-1):
                qXpref11 = qXpref1[:, j]
                if fractional == 1:
                    Rqref11 = Rqref1[:, j]
                    if np.amin(np.diff(qXpref11)) == 0:
                        qXpref11o = qXpref11
                        Rqref11o = Rqref11
                        qXpref11 = np.unique(qXpref11)
                        Rqref11 = np.nan * qXpref11
                        for k in np.arange(0, len(qXpref11)).reshape(-1):
                            Rqref11[k] = np.mean(Rqref11o[np.argwhere(qXpref11o == qXpref11[k])[:, 0]])
                    Xpc[seltp, j] = minX + np.multiply((Xp[seltp, j] - minX),
                                                       (Xplin.Xplin(Xp[seltp, j], qXpref11, 1) * Rqref11))
                    # Apply quantile correction, interpolating between correction factors.
                    # NOTE: We extrapolate a constant (1 passed to Xplin.m), meaning that the correction factors
                    #       for values beyond the range of variability in the reference period are assigned the
                    #       correction factor for the lowest or highest value during the reference period.
                    # NOTE: A stability constant (minX) is used to stabilize the calculation.
                    #       The calculation is not very sensitive to the value of this constant, as long as it is
                    #       large enough to avoid ~0 values in the denominator of Rq011.
                else:
                    dqref11 = qXhref1[:, j] - qXpref11
                    if np.amin(np.diff(qXpref11)) == 0:
                        qXpref11o = qXpref11
                        dqref11o = dqref11
                        qXpref11 = np.unique(qXpref11)
                        dqref11 = np.nan * qXpref11
                        for k in np.arange(0, len(qXpref11)).reshape(-1):
                            dqref11[k] = np.mean(dqref11o[np.argwhere(qXpref11o == qXpref11[k])[:, 0]])
                    Xpc[seltp, j] = Xp[seltp, j] + Xplin.Xplin(Xp[seltp, j], qXpref11, 1) * dqref11
                    # Apply quantile correction, interpolating between additive corrections
                    # (with constant extrapolation)
                    # NOTE: Neither of these quantile corrections are strictly trend-preserving, since the
                    #      corrections depend on the original values of the projections; hence
                    #      e.g. the correction to the median will vary (slightly) over time as the median value
                    #      increases/decreases.

        out.qXhref = qXhref
        out.Rqref = Rqref
    elif method == 2:
        iref = np.argwhere(np.logical_and(yearpminv == yearminref, yearpmaxv == yearmaxref))[:, 0][0]
        dtdhref = np.zeros_like(seltclimh0, dtype=int)
        for dth in np.arange(0, len(seltclimh0)):
            dtdhref[dth] = (date.toordinal(date(yearh[seltclimh0[dth]],
                                                monthh[seltclimh0[dth]], dayh[seltclimh0[dth]])) -
                            date.toordinal(date(yearminref, 1, 1)))
        dqXp = np.nan * qXp
        for i in np.arange(0, nperiods).reshape(-1):
            sel1 = np.argwhere(np.logical_and(yearp >= yearpminv[i], yearp <= yearpmaxv[i]))[:, 0]
            nt1 = len(sel1)
            if delta_timescale == 'daily' and allow_block_lag == 1:
                dtdp1 = (date.toordinal(date(yearp[sel1[0]], monthp[sel1[0]], dayp[sel1[0]])) -
                         date.toordinal(date(yearpminv[i], 1, 1)))
                # Starting index within repeated hindcast series
                ihref1 = np.argwhere(abs(dtdhref - dtdp1) <= tol_dtclimp)[:, 0][0]
                # For example if the projection time series only runs from Jan 2nd of the projection block,
                # we should start the repeated hindcast block also from Jan 2nd.
            else:
                ihref1 = 1
            selth1 = seltclimh0[np.arange(ihref1, np.amin([nt1 + ihref1, nthref]))]
            Fhref1 = Fhref[np.arange(ihref1, np.amin([nt1 + ihref1, nthref])), :]
            if i == iref:
                Xpc[sel1, :] = Xh[seltclimh0, :]
            else:
                for j in np.arange(0, ntclimc).reshape(-1):
                    if seasonal == 1:
                        if use_month == 1:
                            dist1 = fnperiodic_distm.fnperiodic_distm(monthh[selth1], tclimc[j], 12)
                        else:
                            dist1 = fnperiodic_distm.fnperiodic_distm(yrdayh[selth1], tclimc[j], ndays_year)
                        if legacy == 1:
                            selthref1 = np.argwhere(dist1 <= tol_dtclim)[:, 0]
                            pass
                        else:
                            # Choose all data within tolerance tol_dtclimp of climatological time centre tclimc[i]
                            selthref1 = np.argwhere(dist1 <= tol_dtclimp)[:, 0]
                            # Note: This should pick up only one time index: tol_dtclimp is only to allow for small
                            #       discrepancies e.g. daily average timestamps not located exactly at midday.
                            pass
                    else:
                        selthref1 = np.arange(0, len(selth1))
                    seltp1 = (sel1[0] - 1) + selthref1
                    for l in np.arange(0, ns).reshape(-1):
                        if fractional == 1:
                            dqXp1 = (np.squeeze(qXp[j, i, :, l]) - minX) / (
                                    np.squeeze(qXp[j, iref, :, l]) - minX)
                            Xpc[seltp1, l] = minX + np.multiply((Xh[selth1[selthref1], l] - minX),
                                                                Xplin.Xplin(Fhref1[selthref1, l], qsel,
                                                                            1)) * dqXp1
                        else:
                            dqXp1 = np.squeeze(qXp[j, i, :, l]) - np.squeeze(qXp[j, iref, :, l])
                            Xpc[seltp1, l] = Xh[selth1[selthref1], l] + Xplin.Xplin(
                                Fhref1[selthref1, l], qsel, 1) * dqXp1
                        # Interpolate correction linearly from the fixed quantiles (e.g. 0.005:0.01:0.995)
                        # to the quantile/CDF value represented by the reference hindcast datum,
                        # with constant extrapolation.

                        dqXp[j, i, :, l] = dqXp1
                if delta_timescale == 'daily' and np.isnan(Xpc[sel1[-1], 0]):
                    if verbose > 0:
                        print(np.array(['Repeating last day of hindcast reference period to apply to years ',
                                        str(yearpminv[i]), ' to ', str(yearpmaxv[i])]))
                    Xpc[sel1[-1], :] = Xpc[sel1[-2], :]
                else:
                    if nt1 > (nthref + 1):
                        raise Exception(np.array(
                            ['Not yet coded re-use of hindcast data for nt1 = ', str(nt1), ' and nthref = ',
                             str(nthref)]))
            # if verbose > 0:
            #     print(np.array(['Done repeated hindcast corrected projections for years ',str(yearpminv[i]),
            #     ' through ',str(yearpmaxv(i))]))
            #     tic
        out.dqXp = dqXp

    # Impose lower/upper limits if provided
    if (XpcL != []) or (XpcH != []):
        Xpco = Xpc
        if XpcL != []:
            Xpc = np.amax(XpcL, Xpc)
        if XpcH != []:
            Xpc = np.amin(XpcH, Xpc)

        if Xpco != Xpc:
            out.Xpco = Xpco
        del Xpco

    if use_XpcH_SWI == 1:
        if verbose > 0:
            print('Calculating maximum daily-average values of shortwave irradiance')
        optSWI.calc_dav = 1
        optSWI.year = yearp
        __, outSWI = calc_surfaceSWI.calc_surfaceSWI(latp, lonp, yrdayp, optSWI)
        if verbose > 0:
            print('Done maximum daily-average values of shortwave irradiance')
        out.XpcH_SWI = outSWI.Q0dav
        Xpc = np.amin(out.XpcH_SWI, Xpc)

    # if verbose > 0:
    #     print(np.array(['Done bias-corrected projections using ',methodstr,' method']))
    #     toc

    # If required, add 'subdelta' variability at the extra temporal resolution of the hindcast (e.g. hourly)
    if (use_subdelta == 1) and (Xhf != []) and (tdhf != []):
        # if verbose > 0:
        #     print('Calculating finescale projections with subdelta variability')
        #     tic
        if correct_subquantiles == 1:
            if delta_timescale == 'daily':
                nqsub = 24
            else:
                if delta_timescale == 'monthly':
                    nqsub = 28
                else:
                    raise Exception('Subdelta variability quantiles not yet defined if '
                                    'hindcast data opt.Xh neither daily nor monthly')

        seltsub = np.argwhere(np.logical_and(yearh >= yearminsub, yearh <= yearmaxsub))[:, 0]
        Xhsub = Xh[seltsub, :]
        tdhsub = tdh[seltsub]
        ntsub = len(tdhsub)

        yearhsub = np.zeros_like(tdhsub, dtype=int)
        monthhsub = np.zeros_like(tdhsub, dtype=int)
        dayhsub = np.zeros_like(tdhsub, dtype=int)
        for tp in np.arange(0, len(tdhsub)):
            tpid = datetime.strptime(str(np.datetime64('1990-01-01') + tdhsub[tp].astype('timedelta64[D]')),
                                     '%Y-%m-%d')
            yearhsub[tp] = tpid.year
            monthhsub[tp] = tpid.month
            dayhsub[tp] = tpid.day

        yearhf = np.zeros_like(tdhf, dtype=int)
        monthhf = np.zeros_like(tdhf, dtype=int)
        dayhf = np.zeros_like(tdhf, dtype=int)
        for tp in np.arange(0, len(tdhf)):
            tpid = datetime.strptime(str(np.datetime64('1990-01-01') + tdhf[tp].astype('timedelta64[D]')),
                                     '%Y-%m-%d')
            yearhf[tp] = tpid.year
            monthhf[tp] = tpid.month
            dayhf[tp] = tpid.day

        seltfsub = np.argwhere(np.logical_and(yearhf >= yearminsub, yearhf <= yearmaxsub))[:, 0]
        Xhfsub = Xhf[seltfsub, :]
        tdhfsub = tdhf[seltfsub]
        ntfsub = len(tdhfsub)
        yearhfsub = np.zeros_like(tdhfsub, dtype=int)
        monthhfsub = np.zeros_like(tdhfsub, dtype=int)
        for tp in np.arange(0, len(tdhfsub)):
            tpid = datetime.strptime(str(np.datetime64('1990-01-01') + tdhfsub[tp].astype('timedelta64[D]')),
                                     '%Y-%m-%d')
            yearhfsub[tp] = tpid.year
            monthhfsub[tp] = tpid.month

        if remove_deltascale_var == 1:
            Xhi = np.interp(tdh, Xh, tdhf, 'linear', 'extrap')
            dXhf = Xhf - Xhi
            del Xhi
            # Note: Here we allow linear extrapolation, consistent with projection baseline
            # variability (Xpci) calculated below.
        if correct_subquantiles == 0:
            dXhfsub = np.nan * Xhfsub
        if correct_substd == 1:
            dXhfstdsub = np.nan * Xhsub
            dXhfminsub = np.nan * Xhsub
            dXhfmaxsub = np.nan * Xhsub
        if correct_subquantiles == 1:
            qXhfsub = np.nan * np.ones((ntsub, nqsub, ns))
            IsortXhfsub = qXhfsub
        m = 0
        for i in np.arange(0, ntsub).reshape(-1):
            if delta_timescale == 'daily':
                sel1 = np.argwhere(
                    np.logical_and(np.logical_and(yearhf == yearhsub[i], monthhf == monthhsub[i]),
                                   dayhf == dayhsub[i]))[:, 0]
            else:
                if delta_timescale == 'monthly':
                    sel1 = np.argwhere(np.logical_and(yearhf == yearhsub[i], monthhf == monthhsub[i]))[:, 0]
                else:
                    if delta_timescale == 'yearly':
                        sel1 = np.argwhere(yearhf == yearhsub[i])[:, 0]
            sel2 = np.arange(m + 1, m + len(sel1) + 1)
            if remove_deltascale_var == 1:
                if correct_subquantiles == 1:
                    Xhfc1 = dXhf[sel1, :]
                    # These are the finescale anomalies about the deltascale variation (day-to-day, or month-to-month).
                else:
                    Xhfc1 = dXhf[sel1, :] + Xhsub[i, :]
                    # This is the finescale variation with the signal
                    # from deltascale variation (day-to-day, or month-to-month) removed.
            else:
                Xhfc1 = Xhf[sel1, :]
            if correct_subquantiles == 1:
                qXhfsub[i, :, :], IsortXhfsub[i, :, :] = np.sort(Xhfc1)
            else:
                if fractional_subdelta == 1:
                    dXhfsub[sel2, :] = (Xhfc1 - minXsub) / (np.ones((len(sel1), 1)) * Xhsub[i, :] - minXsub)
                    # Check: mean(dXhfsub(sel2,:)) #This should be = 1 for all time series IFF
                    # remove_deltascale_var=0 and the hindcast daily/monthly averages were calculated by
                    # simple averages of all hourly/daily data within each day/month.
                else:
                    dXhfsub[sel2, :] = Xhfc1 - Xhsub[i, :]
                    # Check: mean(dXhfsub(sel2,:)) #This should be = 0 for all time series IFF
                    # remove_deltascale_var=0 and the hindcast daily/monthly averages were calculated by
                    # simple averages of all hourly/daily data within each day/month.
            # Record subdelta statistics
            if correct_substd == 1:
                dXhfstdsub[i, :] = np.std(dXhfsub[sel2, :])
                dXhfminsub[i, :] = np.amin(dXhfsub[sel2, :])
                dXhfmaxsub[i, :] = np.amax(dXhfsub[sel2, :])
            m = m + len(sel1)
        if correct_substd == 1:
            out.dXhfstdsub = dXhfstdsub
            out.dXhfminsub = dXhfminsub
            out.dXhfmaxsub = dXhfmaxsub
            out.dXhfrangesub = dXhfmaxsub - dXhfminsub
        if correct_subquantiles == 1:
            out.IsortXhfsub = IsortXhfsub
        if correct_substd == 1:
            if opt_subq.bymonth == 1:
                ndivt_subq = 12
                str1 = ' by month-of-year'
            else:
                ndivt_subq = 1
                str1 = []
            if verbose > 0:
                print(np.array(['Fitting ', opt_subq.model, ' models', str1,
                                ' to correct subdelta standard deviation as function of delta-scale mean']))
            if opt_subq.model == 'ksmooth':
                yphatm = np.nan * np.ones((ndivt_subq, len(opt_subq.xp), ns))
                optk.kfunc = opt_subq.kfunc
                optk.lr_order = opt_subq.lr_order
                if hasattr(opt_subq, 'ndclosest'):
                    optk.ndclosest = opt_subq.ndclosest
            else:
                bs = np.nan * np.ones((ndivt_subq, opt_subq.nparb, ns))
            ifig0 = 130
            for k in np.arange(0, ndivt_subq).reshape(-1):
                if opt_subq.bymonth == 1:
                    reck = np.argwhere(monthhsub == k)[:, 0]
                else:
                    reck = np.arange(0, ntsub)
                for i in np.arange(0, ns).reshape(-1):
                    rec = reck[np.argwhere(np.logical_and(not np.isnan(Xhsub[reck, i]),
                                                          not np.isnan(dXhfstdsub[reck, i])))[:, 0]]
                    if opt_subq.model == 'polynomial':
                        if opt_subq.logX == 1:
                            rec = rec[Xhsub[rec, i] > 0]
                        if opt_subq.logY == 1:
                            rec = rec[dXhfstdsub[rec, i] > 0]
                        nrec = len(rec)
                        if opt_subq.logX == 1:
                            x1 = np.log(Xhsub[rec, i])
                        else:
                            x1 = Xhsub[rec, i]
                        if opt_subq.logY == 1:
                            y1 = np.log(dXhfstdsub[rec, i])
                        else:
                            y1 = dXhfstdsub[rec, i]
                        if opt_subq.nparb == 3:
                            Xreg1 = np.concatenate((np.ones((nrec, 1)), x1, x1 ** 2), axis=1)
                        else:
                            Xreg1 = np.concatenate((np.ones((nrec, 1)), x1), axis=1)
                        # Fit model for subdelta variability std as function of daily mean values
                        # [bs,bsse,Cbs,bsint,statss] = fit_m(y1,Xreg1,struct('doplotset',1,'selxshow',1)); stop;
                        bs1 = np.linalg.solve(Xreg1, y1)
                    else:
                        if opt_subq.model == 'ksmooth':
                            x1 = Xhsub[rec, i]
                            y1 = dXhfstdsub[rec, i]
                            yphatm[k, :, i] = np.amax(0, ks_regress.ks_regress(y1, x1,
                                                                               opt_subq.xp, opt_subq.bdw, optk))
                            # NOTE: We kernel-smooth to a fixed set of prediction nodes, then interpolate
                            # these later with constant extrapolation.
                            # Re-use of the raw input data is generally too slow.
                    if opt_subq.model != 'ksmooth':
                        bs[k, :, i] = bs1
                    if (opt_subq.doplot == 1) and (i == 1):
                        if opt_subq.model == 'polynomial':
                            if opt_subq.bymonth == 1:
                                Xpcmin1 = np.amin(Xpc[np.argwhere(monthp == k)[:, 0], i])
                                Xpcmax1 = np.amax(Xpc[np.argwhere(monthp == k)[:, 0], i])
                            else:
                                Xpcmin1 = np.amin(Xpc[:, i])
                                Xpcmax1 = np.amax(Xpc[:, i])
                            xp1 = np.transpose(np.linspace(Xpcmin1, Xpcmax1, 200))
                            if opt_subq.logX == 1:
                                xp1t = np.log(xp1)
                            else:
                                xp1t = xp1
                            if opt_subq.nparb == 4:
                                yphat1t = bs1[0] + bs1[1] * xp1t + bs1[2] * xp1t ** 2 + bs1[3] * xp1t ** 3
                            else:
                                if opt_subq.nparb == 3:
                                    yphat1t = bs1[0] + bs1[1] * xp1t + bs1[2] * xp1t ** 2
                                else:
                                    yphat1t = bs1[0] + bs1[1] * xp1t
                            if opt_subq.logY == 1:
                                yphat1 = np.exp(yphat1t)
                            else:
                                yphat1 = yphat1t
                        else:
                            if opt_subq.model == 'ksmooth':
                                xp1 = opt_subq.xp
                                yphat1 = np.transpose(np.squeeze(yphatm[k, :, i]))
                        plt.figure(ifig0)
                        plt.subplot(3, 4, k)
                        plt.plot(Xhsub[rec, i], dXhfstdsub[rec, i], 'k.', xp1, yphat1, 'r-', 'LineWidth', 2)
                        plt.xlabel(opt_subq.xstr, 'FontSize', 11)
                        plt.ylabel(opt_subq.ystr, 'FontSize', 11)
                        plt.title(opt_subq.tstr[k], 'FontWeight', 'bold', 'FontSize', 11)
                        if hasattr(opt_subq, 'xlim'):
                            # set(gca,'XLim', opt_subq.xlim)
                            plt.set_xlim(opt_subq.xlim)
                if opt_subq.bymonth == 1 and verbose > 0:
                    print(np.array(['Done ', str(k), ' of 12 months']))
            if opt_subq.model != 'ksmooth':
                out.bs = bs
            del x1, y1, Xreg1
            if verbose > 0:
                print('Done fitting regression models for subdelta standard deviation')
        else:
            if correct_subquantiles == 1:
                if opt_subq.bymonth == 1:
                    ndivt_subq = 12
                    str1 = ' by month-of-year'
                else:
                    ndivt_subq = 1
                    str1 = []
                if verbose > 0:
                    print(np.array(['Fitting ', opt_subq.model, ' models', str1,
                                    ' to correct subdelta quantiles as function of delta-scale mean']))
                if opt_subq.model == 'ksmooth':
                    yphatm = np.nan * np.ones(ndivt_subq, nqsub, len(opt_subq.xp), ns)
                    optk.kfunc = opt_subq.kfunc
                    optk.lr_order = opt_subq.lr_order
                    if hasattr(opt_subq, 'ndclosest'):
                        optk.ndclosest = opt_subq.ndclosest
                else:
                    bs = np.nan * np.ones(ndivt_subq, nqsub, opt_subq.nparb, ns)
                    Jmins = np.nan * np.ones((ndivt_subq, nqsub, ns))
                if opt_subq.resmodel == 'ksmooth':
                    ypreshatm = np.nan * np.ones(ndivt_subq, nqsub, len(opt_subq.xpres), ns)
                    if hasattr(opt_subq, 'ndclosestres'):
                        optkres.ndclosest = opt_subq.ndclosestres
                else:
                    if opt_subq.resmodel:
                        bsres = np.nan * np.ones(ndivt_subq, nqsub, opt_subq.nparbres, ns)
                        Jminsres = np.nan * np.ones((ndivt_subq, nqsub, ns))
                ifig0 = 150
                ifig1 = 170
                if ns > 1:
                    opt_subq.doplot = 0
                for k in np.arange(0, ndivt_subq).reshape(-1):
                    if opt_subq.bymonth == 1:
                        reck = np.argwhere(monthhsub == k)[:, 0]
                    else:
                        reck = np.arange(0, ntsub)
                    if k > 1:
                        opt_subq.doplot = 0

                    if opt_subq.doplot == 1:
                        plt.figure(ifig0 + k)

                        if opt_subq.resmodel:
                            plt.figure(ifig1 + k)

                    for j in np.arange(0, nqsub).reshape(-1):
                        for i in np.arange(0, ns).reshape(-1):
                            rec = reck
                            nrec = len(rec)
                            if opt_subq.model == 'polynomial':
                                if opt_subq.logX == 1:
                                    rec = rec[Xhsub[rec, i] > 0]
                                if opt_subq.logY == 1:
                                    rec = rec[np.argwhere(np.squeeze(qXhfsub[rec, j, i] > 0))[:, 0]]
                                nrec = len(rec)
                                if opt_subq.logX == 1:
                                    x1 = np.log(Xhsub[rec, i])
                                else:
                                    x1 = Xhsub[rec, i]
                                if opt_subq.logY == 1:
                                    y1 = np.squeeze(np.log(qXhfsub[rec, j, i]))
                                else:
                                    y1 = np.squeeze(qXhfsub[rec, j, i])
                                if opt_subq.nparb == 4:
                                    # Xreg1 = np.array([np.ones((nrec, 1)), x1, x1 ** 2, x1 ** 3])
                                    Xreg1 = np.concatenate((np.ones((nrec, 1)), x1, x1 ** 2, x1 ** 3), axis=1)
                                else:
                                    if opt_subq.nparb == 3:
                                        # Xreg1 = np.array([np.ones((nrec, 1)), x1, x1 ** 2])
                                        Xreg1 = np.concatenate((np.ones((nrec, 1)), x1, x1 ** 2), axis=1)
                                    else:
                                        # Xreg1 = np.array([np.ones((nrec, 1)), x1])
                                        Xreg1 = np.concatenate((np.ones((nrec, 1)), x1), axis=1)

                                # Fit model for subdelta quantiles as function of daily mean values
                                # [bs1,bs1se,Cbs1,bs1int,statss1] = fit_m(y1,Xreg1,struct('doplotset',0,'selxshow',1));
                                bs1 = np.linalg.solve(Xreg1, y1)
                                if opt_subq.resmodel:
                                    yhat1 = Xreg1 * bs1
                            else:
                                if opt_subq.model == 'tanh':
                                    x1 = Xhsub[rec, i]
                                    y1 = np.squeeze(qXhfsub[rec, j, i])
                                    bs100 = 1 / np.amax(x1)
                                    if opt_subq.nparb == 3:
                                        # Xfun = lambda p=None: np.array([np.ones((nrec, 1)), np.tanh(p[1] * x1)])
                                        Xfun = lambda p=None: np.concatenate((np.ones((nrec, 1)), np.tanh(p[0] * x1)),
                                                                             axis=1)
                                    else:
                                        Xfun = lambda p=None: np.tanh(p[0] * x1)
                                    # optJ = struct('Y', y1, 'Xfun', Xfun, 'nonnegb', opt_subq.nonnegb)
                                    optJ.Y = y1
                                    optJ.Xfun = Xfun
                                    optJ.nonnegb = opt_subq.nonnegb
                                    # opto = struct('optJ', optJ)
                                    opto.optJ = optJ
                                    # phat = optimize_theta.optimize_theta(bs100, opt_subq.pL, opt_subq.pH,
                                    #                                      J_fit_partially_linear_model, opt_subq.nr,
                                    #                                      opto)
                                    # Jmins1, outJ = J_fit_partially_linear_model.J_fit_partially_linear_model(phat, optJ)
                                    phat = []
                                    Jmins1, outJ = [[], []]
                                    # bs1 = np.array([[outJ.b], [phat]])
                                    bs1 = np.concatenate(([outJ.b], [phat]), axis=1)
                                    if not len(opt_subq.resmodel) == 0:
                                        yhat1 = outJ.Ym
                                    Jmins[k, j, i] = Jmins1
                                else:
                                    if opt_subq.model == 'powerlaw':
                                        x1 = Xhsub[rec, i]
                                        y1 = np.squeeze(qXhfsub[rec, j, i])
                                        if opt_subq.nparb == 2:
                                            bs100 = np.linalg.solve(x1, y1)
                                            bs10 = [np.log(bs100), 1]
                                            Xfun = lambda p=None: x1 ** p[0]
                                            # optJ = struct('Y', y1, 'Xfun', Xfun, 'nonnegb', opt_subq.nonnegb)
                                            # opto = struct('optJ', optJ)
                                            optJ.Y = y1
                                            optJ.Xfun = Xfun
                                            optJ.nonnegb = opt_subq.nonnegb
                                            opto.optJ = optJ
                                            # phat = optimize_theta(bs10(2), opt_subq.pL, opt_subq.pH,
                                            #                       J_fit_partially_linear_model, opt_subq.nr, opto)
                                            # Jmins1, outJ = J_fit_partially_linear_model(phat, optJ)
                                            phat = []
                                            Jmins1, outJ = [[], []]
                                            bs1 = np.concatenate((outJ.b, phat), axis=0)
                                            if not len(opt_subq.resmodel) == 0:
                                                yhat1 = outJ.Ym
                                            Jmins[k, j, i] = Jmins1
                                        else:
                                            if opt_subq.nparb == 3:
                                                bs100 = np.linalg.solve(x1, y1)
                                                bs10 = np.concatenate(([np.log(bs100)], [1], [0]), axis=0)
                                                Jfun1 = lambda p=None: np.amax(0, p[0] * x1 ** p[1] + p[2]) - y1
                                                # Xfun = @(p) max(0,x1-p(2)).^p(1)
                                                # bs1 = lsqnonlin(Jfun1, bs10)
                                                bs1 = np.linalg.lstsq(Jfun1, bs10)
                                                # bs1 = np.polyfit(Jfun1, bs10)
                                                if opt_subq.resmodel:
                                                    yhat1 = Jfun1[bs1] + y1
                                    else:
                                        if opt_subq.model == 'ksmooth':
                                            x1 = Xhsub[rec, i]
                                            y1 = np.squeeze(qXhfsub[rec, j, i])
                                            yphatm[k, j, :, i] = ks_regress.ks_regress(y1, x1,
                                                                                       opt_subq.xp, opt_subq.bdw, optk)

                                            # NOTE: We kernel-smooth to a fixed set of prediction nodes,
                                            # then interpolate these later with constant extrapolation. Re-use of the
                                            # raw input data is generally too slow. if ~isempty(opt_subq.resmodel);
                                            # yhat1 = ks_regress(y1,x1,x1,opt_subq.bdw,optk); end
                                            if not len(opt_subq.resmodel) == 0:
                                                yhat1 = Xplin.Xplin(x1, opt_subq.xp, 1) * np.squeeze(yphatm[k, j, :, i])
                            if opt_subq.model != 'ksmooth':
                                bs[k, j, :, i] = bs1
                            if opt_subq.resmodel:
                                x1res = Xhsub[rec, i]
                                if fractional_subdelta == 1:
                                    y1res = np.abs(y1 - yhat1) / yhat1
                                else:
                                    y1res = np.abs(y1 - yhat1)
                            if opt_subq.resmodel == 'polynomial':
                                rec2 = np.arange(0, nrec)
                                if opt_subq.logXres == 1:
                                    rec2 = rec2[np.argwhere(x1res[rec2] > 0)[:, 0]]
                                if opt_subq.logYres == 1:
                                    rec2 = rec2[np.argwhere(y1res[rec2] > 0)[:, 0]]
                                nrec2 = len(rec2)
                                x1res = x1res[rec2]
                                y1res = y1res[rec2]
                                if opt_subq.logXres == 1:
                                    x1rest = np.log(x1res)
                                else:
                                    x1rest = x1res
                                if opt_subq.logYres == 1:
                                    y1rest = np.log(y1res)
                                else:
                                    y1rest = y1res
                                if opt_subq.nparbres == 3:
                                    Xres1 = np.concatenate((np.ones((nrec2, 1)), x1rest, x1rest ** 2), axis=0)
                                else:
                                    Xres1 = np.concatenate((np.ones((nrec2, 1)), x1rest), axis=0)
                                bsres1 = np.linalg.solve(Xres1, y1rest)
                            else:
                                if opt_subq.resmodel == 'gaussian':
                                    bsres100 = np.std(x1res)
                                    Xfun = lambda p=None: np.exp(- (0.5 * p[0] ** 2) * x1res ** 2)
                                    # optJ = struct('Y', y1res, 'Xfun', Xfun, 'nonnegb', opt_subq.nonnegbres)
                                    # opto = struct('optJ', optJ)
                                    optJ.Y = y1res
                                    optJ.Xfun = Xfun
                                    optJ.nonnegb = opt_subq.nonnegbres
                                    opto.optJ = optJ
                                    # preshat = optimize_theta(bsres100, opt_subq.presL, opt_subq.presH,
                                    #                          J_fit_partially_linear_model, opt_subq.nrres, opto)
                                    # Jmins1, outJ = J_fit_partially_linear_model(preshat, optJ)
                                    preshat = []
                                    Jmins1, outJ = [[], []]
                                    bsres1 = np.concatenate((outJ.b, preshat), axis=0)
                                    Jminsres[k, j, i] = Jmins1
                                else:
                                    if opt_subq.resmodel == 'powerlaw':
                                        bsres100 = np.linalg.solve(x1res, y1res)
                                        bsres10 = np.concatenate((bsres100, 1), axis=0)
                                        Jfun2 = lambda p=None: p[0] * x1res ** p[1] - y1res
                                        # bsres1 = lsqnonlin(Jfun2, bsres10)
                                        bsres1 = np.linalg.lstsq(Jfun2, bsres10)
                                        # bsres1 = np.polyfit(Jfun2, bsres10)
                                    else:
                                        if str(opt_subq.resmodel) == str('ksmooth') == 1:
                                            ypreshatm[k, j, :, i] = np.amax(0,
                                                                            ks_regress.ks_regress(y1res, x1res,
                                                                                                  opt_subq.xpres,
                                                                                                  opt_subq.bdwres,
                                                                                                  optkres))
                                            # NOTE:
                                            # We kernel-smooth to a fixed set of prediction nodes, then interpolate
                                            # these later with constant extrapolation. Re-use of the raw input data
                                            # is generally too slow.
                            if opt_subq.resmodel and opt_subq.resmodel != 'ksmooth':
                                bsres[k, j, :, i] = bsres1
                            if opt_subq.doplot == 1:
                                if opt_subq.bymonth == 1:
                                    Xpcmin1 = np.amin(Xpc(np.argwhere(monthp == k)[:, 0], i))
                                    Xpcmax1 = np.amax(Xpc(np.argwhere(monthp == k)[:, 0], i))
                                else:
                                    Xpcmin1 = np.amin(Xpc[:, i])
                                    Xpcmax1 = np.amax(Xpc[:, i])
                                xp1 = np.transpose(np.linspace(Xpcmin1, Xpcmax1, 200))
                                if opt_subq.model == 'polynomial':
                                    if opt_subq.logX == 1:
                                        xp1t = np.log(xp1)
                                    else:
                                        xp1t = xp1
                                    if opt_subq.nparb == 4:
                                        yphat1t = bs1[0] + bs1[1] * xp1t + bs1[2] * xp1t ** 2 + bs1[3] * xp1t ** 3
                                    else:
                                        if opt_subq.nparb == 3:
                                            yphat1t = bs1[0] + bs1[1] * xp1t + bs1[2] * xp1t ** 2
                                        else:
                                            yphat1t = bs1[0] + bs1[1] * xp1t
                                    if opt_subq.logY == 1:
                                        yphat1 = np.exp(yphat1t)
                                    else:
                                        yphat1 = yphat1t
                                else:
                                    if opt_subq.model == 'tanh':
                                        if opt_subq.nparb == 3:
                                            yphat1 = bs1[0] + bs1[1] * np.tanh(bs1[2] * xp1)
                                        else:
                                            yphat1 = bs1[0] * np.tanh(bs1[1] * xp1)
                                    else:
                                        if opt_subq.model == 'powerlaw':
                                            if opt_subq.nparb == 2:
                                                yphat1 = bs1[0] * xp1 ** bs1[1]
                                            if opt_subq.nparb == 3:
                                                yphat1 = np.amax(0, bs1[0] * xp1 ** bs1[1] + bs1[2])
                                        else:
                                            if opt_subq.model == 'ksmooth':
                                                yphat1 = Xplin.Xplin(xp1, opt_subq.xp, 1) * np.squeeze(
                                                    yphatm[k, j, :, i])
                                if opt_subq.resmodel == 'polynomial':
                                    if opt_subq.logXres == 1:
                                        xp1t = np.log(xp1)
                                    else:
                                        xp1t = xp1
                                    if opt_subq.nparbres == 3:
                                        ypreshat1t = bsres1[0] + bsres1[1] * xp1t + bsres1[2] * xp1t ** 2
                                    else:
                                        ypreshat1t = bsres1[0] + bsres1[1] * xp1t
                                    if opt_subq.logYres == 1:
                                        ypreshat1 = np.exp(ypreshat1t)
                                    else:
                                        ypreshat1 = ypreshat1t
                                else:
                                    if opt_subq.resmodel == 'gaussian':
                                        ypreshat1 = bsres1[0] * np.exp(- (0.5 * bsres1[1] ** 2) * xp1 ** 2)
                                    else:
                                        if opt_subq.resmodel == 'powerlaw':
                                            ypreshat1 = bsres1[0] * xp1 ** bsres1[1]
                                        else:
                                            if opt_subq.resmodel == 'ksmooth':
                                                ypreshat1 = Xplin.Xplin(xp1, opt_subq.xpres, 1) * np.squeeze(
                                                    ypreshatm[k, j, :, i])
                                plt.figure(ifig0 + k)
                                plt.subplot(4, 6, j)
                                plt.plot(Xhsub[rec, i], np.squeeze(qXhfsub[rec, j, i]), 'k.', xp1, yphat1, 'r-',
                                         'LineWidth', 2)
                                if opt_subq.resmodel:
                                    plt.figure(ifig1 + k)
                                    plt.subplot(4, 6, j)
                                    plt.plot(x1res, y1res, 'k.', xp1, ypreshat1, 'r-', 'LineWidth', 2)
                    if (opt_subq.bymonth == 1) and (verbose > 0):
                        print(np.array(['Done ', str(k), ' of 12 months']))
                if opt_subq.model != 'ksmooth':
                    out.bs = bs
                if opt_subq.resmodel and not opt_subq.resmodel != 'ksmooth':
                    out.bsres = bsres
                if opt_subq.model == 'powerlaw' or opt_subq.model == 'tanh':
                    out.Jmins = Jmins

                del x1, y1, Xreg1
                if verbose > 0:
                    print('Done fitting regression models for subdelta quantiles')
                    # Next calculate the baseline variability and prepare for main loop over projection periods
                    # First we need to establish the time stamps for the subdelta projection data:
        if (len(yearpminv) == 0) or (len(yearpmaxv) == 0):
            nyrssub = yearmaxsub - yearminsub + 1
            yearpminv = np.arange(np.amin(yearp), np.amax(yearp) + nyrssub, nyrssub)
            yearpmaxv = np.arange(np.amin(yearp) + nyrssub - 1, np.amax(yearp) + nyrssub, nyrssub)
            # If not matching the subdelta variability by most-similar days, then we need to define
            # a set of multiannual year periods for which the subdelta variability will be repeated over
            # (usually most convenient to divide the prediciton period into 10 or 20-year periods).
            # If not defined by input these periods are defined here.
            # Although it should not be necessary to split the full prediction in match_subdelta_deltamean=1,
            # tests show that the calculation is almost 2x faster if it is split up (with identical final results).
        selp = np.argwhere(yearpmaxv <= max(yearp))[:, 0]
        yearpminv = yearpminv[selp]
        yearpmaxv = yearpmaxv[selp]
        out.yearpminv = yearpminv
        out.yearpmaxv = yearpmaxv
        dtdhf1 = np.amin(tdhf - int(np.floor(tdhf)))
        if delta_timescale == 'daily':
            difftdhf = 1 / 24
            tdpf = np.transpose((np.arange(tdpfmin + dtdhf1, tdpfmax + 1 + difftdhf, difftdhf)))
            if verbose > 0:
                print('Assuming regular hourly increments for subdelta variability')
        elif (np.amax(np.diff(tdhf)) - np.amin(np.diff(tdhf))) < 0.001:
            difftdhf = np.mean(np.diff(tdhf))
            tdpf = np.transpose((np.arange(tdpfmin + dtdhf1, tdpfmax + 1 + difftdhf, difftdhf)))
            if verbose > 0:
                print(np.array(
                    ['Assuming regular increments ', str(difftdhf, '%3.2f'), ' day(s) for subdelta variability']))
        else:
            # In this case we allow some irregularity in the subdelta time increments (as in e.g. xCO2 hindcast data);
            # we repeat the subdelta time increments relative to the start of each prediction block.
            if verbose > 0:
                print(
                    'Building projection subdelta timestamps by repeating increments '
                    'from beginning of subdelta reference period')
            m = 0
            for i in np.arange(0, len(yearpminv)).reshape(-1):
                if i == 0:
                    tdpf = (date.toordinal(yearpminv[i], 1, 1) + (tdhfsub - date.toordinal(yearminsub, 1, 1)))
                else:
                    tdpf = np.concatenate(tdpf,
                                          (date.toordinal(yearpminv[i], 1, 1) + (
                                                  tdhfsub - date.toordinal(yearminsub, 1, 1))))
                # tdpf[np.arange(m + 1, m + ntfsub+1)] = datenum(yearpminv[i],1,1) + (tdhfsub - datenum(yearminsub,1,1))
                # tdpf[np.arange(m + 1, m + ntfsub + 1)] = (date.toordinal(yearpminv[i], 1, 1) +
                #                                                (tdhfsub - date.toordinal(yearminsub, 1, 1)))
                m = m + ntfsub
            tdpf = tdpf[:]
        ntpf = len(tdpf)
        # yearpf,monthpf,__ = datevec(tdpf)
        yearpf = (tdpf + datetime(1990, 1, 1)).year
        monthpf = (tdpf + datetime(1990, 1, 1)).month
        daypf = (tdpf + datetime(1990, 1, 1)).day
        # yrdaypf = tdpf - datenum(yearpf,1,1)
        yrdaypf = tdpf - date.toordinal(yearpf, 1, 1)
        Xpcf = np.nan * np.ones((ntpf, ns))
        if verbose > 0:
            print('Calculating baseline variability to which subdelta variation will be added')
        if remove_deltascale_var > 0:
            Xpci = np.interp(tdp, Xpc, tdpf, 'linear', 'extrap')
            # Note: Here we allow linear extrapolation, consistent with definition of subdaily variations above.
            # Note: use of Xplin.m (Xpci = Xplin.Xplin(tdpf,tdp,1)*Xpc) demands too much memory.
        else:
            Xpci = Xpcf
            for i in np.arange(0, ntp).reshape(-1):
                if delta_timescale == 'daily':
                    self1 = np.argwhere(int(np.floor(tdpf)) == int(np.floor(tdp[i])))[:, 0]
                else:
                    if delta_timescale == 'monthly':
                        self1 = np.argwhere(np.logical_and(yearpf == yearp[i], monthpf == monthp[i]))[:, 0]
                    else:
                        if delta_timescale == 'yearly':
                            self1 = np.argwhere(yearpf == yearp[i])[:, 0]
                Xpci[self1, :] = np.ones((len(self1), 1)) * Xpc[i, :]
        if verbose > 0:
            print('Done baseline variability')
        if match_subdelta_deltamean == 1:
            selhsub = np.nan * np.ones((ntp, nblocks))
        if match_subdelta_hourly == 1:
            selhsubhour = np.nan * np.ones((ntp, nblocks))
        correction_factor = np.ones((ntp, ns))
        # dtdhsub = tdh(seltsub) - datenum(yearminsub,1,1)
        dtdhsub = tdh[seltsub] - date.toordinal(yearminsub, 1, 1)
        # Main loop over projection periods for which finescale projections are required
        for j in np.arange(0, len(yearpminv)).reshape(-1):
            sel1 = np.argwhere(np.logical_and(yearp >= yearpminv[j], yearp <= yearpmaxv[j]))[:, 0]
            self1 = np.argwhere(np.logical_and(yearpf >= yearpminv[j], yearpf <= yearpmaxv[j]))[:, 0]
            nt1 = len(sel1)
            ntf1 = len(self1)
            if match_subdelta_deltamean == 0:
                if delta_timescale == 'daily' and allow_block_lag == 1:
                    # Identify starting index for the repeated hindcast variability
                    # dtdp1 = tdp(sel1(1)) - datenum(yearpminv(j),1,1)
                    dtdp1 = tdp(sel1[1]) - date.toordinal(yearpminv[j], 1, 1)
                    # Starting index within repeated hindcast series (deltascale)
                    ihsub1 = np.argwhere(abs(dtdhsub - dtdp1) <= tol_dtclimp)[:, 0]
                else:
                    ihsub1 = 1
                # Identify repeated hindcast series to use for this block
                selsub1 = np.arange(ihsub1, np.amin(ihsub1 - 1 + nt1, ntsub) + 1)
                if correct_subquantiles == 0:
                    # Note: This is not necesssary if correct_subquantiles=1 because in this case only
                    # the deltascale indices selsub1 are used
                    if str(delta_timescale) == str('daily') == 1:
                        ihfsub1 = np.argwhere(int(np.floor(tdhfsub)) == int(np.floor(tdhsub(selsub1[0]))))[:, 0][0]
                        ihfsub2 = np.argwhere(int(np.floor(tdhfsub)) == int(np.floor(tdhsub(selsub1[-1]))))[:, 0][-1]
                    else:
                        if str(delta_timescale) == str('monthly') == 1:
                            ihfsub1 = np.argwhere(np.logical_and(yearhfsub == yearhsub[selsub1[0]],
                                                                 monthhfsub == monthhsub[selsub1[0]]))[:, 0][0]
                            ihfsub2 = np.argwhere(np.logical_and(yearhfsub == yearhsub[selsub1[-1]],
                                                                 monthhfsub == monthhsub[selsub1[-1]]))[:, 0][-1]
                        else:
                            if str(delta_timescale) == str('yearly') == 1:
                                ihfsub1 = np.argwhere(yearhfsub == yearhsub[selsub1[0]])[:, 0][0]
                                ihfsub2 = np.argwhere(yearhfsub == yearhsub[selsub1[-1]])[:, 0][-1]
                    dXhf1 = dXhfsub[np.arange(ihfsub1, ihfsub2), :]
                    if delta_timescale != 'daily' and len(dXhf1) == (ntf1 - 1):
                        dXhf1 = np.concatenate((dXhf1, dXhf1[-1, :]), axis=0)
                        if verbose > 0:
                            print(np.array(
                                ['Repeating last entry of reference period subdelta variability to apply to years ',
                                 str(yearpminv[j]), ' to ', str(yearpmaxv[j])]))
                # Add extra day if necessary to account for possible different no. leap years per decade
                # (e.g. only 2 in subdelta period 2010-2019)
                if delta_timescale == 'daily' and (ihsub1 - 1 + nt1 > ntsub):
                    if (ihsub1 - 1 + nt1) == (ntsub + 1):
                        if verbose > 0:
                            print(np.array(['Repeating last day of subdelta reference period to apply to years ',
                                            str(yearpminv[j]), ' to ', str(yearpmaxv[j])]))
                        selsub1 = np.concatenate((selsub1, ntsub), axis=0)
                        # The use of 'selsub1' avoids having to append several variables ('Xhsub' etc.) for use below.
                        if correct_subquantiles == 0:
                            dXhf1 = np.concatenate((dXhf1, dXhf1[np.arange(-24, -1), :]), axis=0)
                    else:
                        raise Exception(np.array(
                            ['Problem finding enough days of data from subdelta reference period to apply to years ',
                             str(yearpminv[j]), ' to ', str(yearpmaxv[j])]))
                if correct_substd == 1:
                    if verbose > 0:
                        print('Calculating correction factors for subdelta variance')
                    if opt_subq.model != 'ksmooth':
                        bs1 = np.nan * np.ones((opt_subq.nparb, ns))
                    correction_factor1 = np.nan * np.ones((nt1, ns))
                    for k in np.arange(0, ndivt_subq).reshape(-1):
                        if opt_subq.bymonth == 1:
                            rec1 = np.argwhere(monthp[sel1] == k)[:, 0]
                        else:
                            rec1 = np.arange(0, nt1)
                        nrec1 = len(rec1)
                        if opt_subq.model == 'polynomial':
                            bs1[:, :] = bs[k, :, :]
                            if opt_subq.logX == 1:
                                x1 = np.log(Xhsub[selsub1[rec1], :])
                                xp1 = np.log(Xpc[sel1[rec1], :])
                            else:
                                x1 = Xhsub[selsub1[rec1], :]
                                xp1 = Xpc[sel1[rec1], :]
                            if opt_subq.logY == 1:
                                if opt_subq.nparb == 3:
                                    correction_factor11 = np.exp(
                                        np.multiply((np.ones((nrec1, 1)) * bs1[1, :]), (xp1 - x1)) + np.multiply(
                                            (np.ones((nrec1, 1)) * bs1[2, :]), (xp1 ** 2 - x1 ** 2)))
                                else:
                                    correction_factor11 = np.exp(
                                        np.multiply((np.ones((nrec1, 1)) * bs1[1, :]), (xp1 - x1)))
                            else:
                                correction_factor11 = (np.ones((nrec1, 1)) * bs1[0, :] + np.multiply(
                                    (np.ones((nrec1, 1)) * bs1[1, :]), xp1) + np.multiply(
                                    (np.ones((nrec1, 1)) * bs1[2, :]), xp1 ** 2)) / (
                                                              np.ones((nrec1, 1)) * bs1[0, :] + np.multiply(
                                                          (np.ones((nrec1, 1)) * bs1[1, :]), x1) + np.multiply(
                                                          (np.ones((nrec1, 1)) * bs1[2, :]), x1 ** 2))
                            if opt_subq.logX == 1:
                                correction_factor11[xp1 == - np.inf] = 0
                                correction_factor11[x1 == - np.inf] = 1
                                correction_factor11[np.argwhere(np.logical_and(xp1 == - np.inf,
                                                                               x1 == - np.inf))[:, 0]] = 1
                            correction_factor1[rec1, :] = np.amax(0, correction_factor11)
                        else:
                            if opt_subq.model == 'ksmooth':
                                yphat1 = np.nan * np.ones((nrec1, ns))
                                yphat2 = yphat1
                                yphat11 = np.nan * np.ones((len(opt_subq.xp), 1))
                                for l in np.arange(0, ns).reshape(-1):
                                    yphat11 = yphatm[k, :, l]
                                    yphat1[:, l] = Xplin[Xpc[sel1[rec1], l], opt_subq.xp, 1] * yphat11
                                    yphat2[:, l] = Xplin[Xhsub[selsub1[rec1], l], opt_subq.xp, 1] * yphat11
                                correction_factor1[rec1, :] = np.amax(0, yphat1 / yphat2)
                    if verbose > 0:
                        print('Done correction factors for subdelta variance')
                else:
                    correction_factor1 = np.ones((nt1, ns))
                correction_factor[sel1, :] = correction_factor1
            if verbose > 0:
                print(np.array(
                    ['Adding subdelta variability to projections for period ', str(yearpminv[j]), ' through ',
                     str(yearpmaxv[j])]))
            if loop_over_projection_times == 0 and delta_timescale == 'daily':
                if fractional_subdelta == 1:
                    Xpcf[self1, :] = minXsub + np.multiply(
                        np.kron((np.multiply((Xpc[sel1, :] - minXsub), correction_factor1)), np.ones((24, 1))), dXhf1)
                    if remove_deltascale_var == 1:
                        Xpcf[self1, :] = Xpcf[self1, :] + (Xpci[self1, :] - np.kron(Xpc[sel1, :], np.ones((24, 1))))
                    #   Check: This should recover Xhf in case where Xpc=Xh:
                    #   If remove_deltascale_var=0:
                    #       Xpcf = minXsub + (Xh-minXsub)*(Xhfc1-minXsub)/(Xh-minXsub)
                    #            = minXsub + (Xhf-minXsub)                               (since Xhfc1=Xhf in this case)
                    #            = Xhf
                    #   If remove_deltascale_var=1:
                    #       Xpcf = minXsub + (Xh-minXsub)*(Xhfc1-minXsub)/(Xh-minXsub) + (Xhi-Xh)
                    #            = minXsub + (Xhf-Xhi+Xh-minXsub) + (Xhi-Xh)     (since Xhfc1=Xhf-Xhi+Xh in this case)
                    #            = Xhf
                else:
                    Xpcf[self1, :] = Xpci[self1, :] + np.multiply(np.kron(correction_factor1,
                                                                          np.ones((24, 1))), dXhf1)
                    # Check: This should recover Xhf in case where Xpc=Xh:
                    #   If remove_deltascale_var=0:
                    #       Xpcf = Xh + (Xhfc1-Xh)                                (since Xhi=Xh in this case)
                    #            = Xh + (Xhf-Xh)                                  (since Xhfc1=Xhf in this case)
                    #            = Xhf
                    #   If remove_deltascale_var=1:
                    #       Xpcf = Xhi + (Xhfc1-Xh)
                    #            = Xhi + (Xhf-Xhi+Xh-Xh)                          (since Xhfc1=Xhf-Xhi+Xh in this case)
                    #            = Xhf
            else:
                for iblock in np.arange(0, nblocks).reshape(-1):
                    if subdelta_by_blocks == 1:
                        selblock = np.argwhere(np.logical_and(np.logical_and(latp >= latL_blocks1[iblock],
                                                                             latp < latH_blocks1[iblock]),
                                                              np.logical_and(lonp >= lonL_blocks1[iblock],
                                                                             lonp < lonH_blocks1[iblock])))[:, 0]
                        pass
                    else:
                        selblock = np.arange(0, ns)
                    nselblock = len(selblock)
                    countf = 0
                    for i in np.arange(0, nt1).reshape(-1):
                        if delta_timescale == 'daily':
                            self2 = self1[int(np.floor(tdpf[self1])) == int(np.floor(tdp[sel1[i]]))]
                        elif delta_timescale == 'monthly':
                            self2 = self1[np.argwhere(np.logical_and(yearpf[self1] == yearp[sel1[i]],
                                                                     monthpf[self1]) == monthp[sel1[i]])[:, 0]]
                        elif delta_timescale == 'yearly':
                            self2 = self1[np.argwhere(yearpf[self1] == yearp[sel1[i]])]
                        nself2 = len(self2)
                        if match_subdelta_deltamean == 1:
                            diff1 = np.sum((Xhsub[:, selblock] - Xpc[sel1[i], selblock]) ** 2, 1)
                            selhsub[sel1[i], iblock] = np.argwhere(diff1 == np.amin(diff1))[:, 0][0]
                            if delta_timescale == 'daily':
                                selhfsub1 = np.argwhere(
                                    int(np.floor(tdhfsub)) == int(np.floor(tdhsub[selhsub[sel1[i], iblock]])))[:, 0]
                                # Time indices of selected subdelta variability from dXhfsub
                            elif delta_timescale == 'monthly':
                                selhfsub1 = np.argwhere(np.logical_and(
                                    yearhfsub == yearhsub[selhsub[sel1[i], iblock]],
                                    monthhfsub == monthhsub[selhsub[sel1[i], iblock]]))[:, 0]
                            elif delta_timescale == 'yearly':
                                selhfsub1 = np.argwhere(yearhfsub == yearhsub[selhsub[sel1[i], iblock]])[:, 0]
                            dXhfsub1 = dXhfsub[selhfsub1, selblock]
                            if correct_substd == 1:
                                if opt_subq.bymonth == 1:
                                    m = monthp[sel1[i]]
                                else:
                                    m = 1
                                if opt_subq.model != 'ksmooth':
                                    bs1 = np.nan * np.ones((opt_subq.nparb, nselblock))
                                if opt_subq.model == 'polynomial':
                                    bs1[:, :] = bs[m, :, selblock]
                                    if opt_subq.logX == 1:
                                        x1 = np.log(Xhsub[selhsub[sel1[i]], selblock])
                                        xp1 = np.log(Xpc[sel1[i], selblock])
                                    else:
                                        x1 = Xhsub[selhsub[sel1[i]], selblock]
                                        xp1 = Xpc[sel1[i], selblock]
                                    if opt_subq.logY == 1:
                                        if opt_subq.nparb == 3:
                                            correction_factor1 = np.exp(
                                                np.multiply(bs1[1, selblock], (xp1 - x1)) + np.multiply(
                                                    bs1[2, selblock], (xp1 ** 2 - x1 ** 2)))
                                        else:
                                            correction_factor1 = np.exp(np.multiply(bs1[1, selblock], (xp1 - x1)))
                                    else:
                                        correction_factor1 = np.amax(0, (
                                                bs1[0, selblock] + np.multiply(bs1[1, selblock], xp1) + np.multiply(
                                            bs1[2, selblock], xp1 ** 2)) / (bs1[0, selblock] + np.multiply(
                                            bs1[1, selblock], x1) + np.multiply(bs1[2, selblock], x1 ** 2)))
                                    if opt_subq.logX == 1:
                                        correction_factor1[xp1 == - np.inf] = 0
                                        correction_factor1[x1 == - np.inf] = 1
                                        correction_factor1[xp1 == np.logical_and[- np.inf, x1] == - np.inf] = 1
                                elif opt_subq.model == 'ksmooth':
                                    yphat1 = np.nan * np.ones((1, nselblock))
                                    yphat2 = yphat1
                                    yphat11 = np.nan * np.ones((len(opt_subq.xp), 1))
                                    for l in np.arange(0, nselblock).reshape(-1):
                                        yphat11 = yphatm[m, :, selblock[l]]
                                        yphat1[l] = Xplin[Xpc[sel1[i], selblock[l]], opt_subq.xp, 1] * yphat11
                                        yphat2[l] = Xplin[Xhsub[
                                            selhsub[sel1[i]], selblock[l]], opt_subq.xp, 1] * yphat11
                                    correction_factor1[rec1, :] = np.amax(0, yphat1 / yphat2)
                            else:
                                correction_factor1 = np.ones((1, nselblock))
                            correction_factor[sel1[i], selblock] = correction_factor1

                            # NOTE: The correction factors derived with match_subdelta_deltamean=1
                            # should be milder (closer to 1) than those derived with match_subdelta_deltamean=0.

                            rec1 = np.transpose((np.arange(0, 24)))
                            self2last = np.amin(self2) - 1
                            if match_subdelta_hourly == 1 and self2last > 0:
                                if fractional_subdelta == 1:
                                    Xpcf1 = minXsub + np.multiply((np.ones((nself2, 1)) * (
                                        np.multiply((Xpci[self2[0], selblock] - minXsub), correction_factor1))),
                                                                  dXhfsub1)
                                    if remove_deltascale_var == 1:
                                        Xpcf1 = Xpcf1 + (Xpci[self2, selblock] - np.ones((nself2, 1)) * Xpc[
                                            sel1[i], selblock])
                                else:
                                    Xpcf1 = np.ones((nself2, 1)) * Xpci[self2[0], selblock] + np.multiply(
                                        (np.ones((nself2, 1)) * correction_factor1), dXhfsub1)
                                diff2 = np.sum((Xpcf1 - np.ones((nself2, 1)) * Xpcf[self2last, selblock]) ** 2, 1)
                                selhsubhour[sel1[i], iblock] = np.argwhere(diff2 == np.amin(diff2))[:, 0]
                                if selhsubhour(sel1(i), iblock) > 1:
                                    rec1 = np.transpose(np.array([np.arange(selhsubhour[sel1[i], iblock], 24),
                                                                  np.arange(23, 24 -
                                                                          selhsubhour[sel1[i], iblock] - 1,
                                                                            - 1)]))
                                    # Rearrange the hourly subdeltas to best match the end of the previous day,
                                    # and 'reflect' indices off the end of the 24 hour interval in order to fill the
                                    # 24 hours without introducing a jump.
                            if fractional_subdelta == 1:
                                Xpcf[self2, selblock] = minXsub + np.multiply((np.ones((nself2, 1)) * (
                                    np.multiply((Xpc[sel1[i], selblock] - minXsub), correction_factor1))),
                                                                              dXhfsub1[rec1, :])
                                if remove_deltascale_var == 1:
                                    Xpcf[self2, selblock] = Xpcf[self2, selblock] + (
                                            Xpci[self2, selblock] - np.ones((nself2, 1)) * Xpc[sel1[i], selblock])
                            else:
                                Xpcf[self2, selblock] = Xpci[self2, selblock] + np.multiply(
                                    (np.ones((nself2, 1)) * correction_factor1), dXhfsub1[rec1, :])
                        else:
                            if correct_subquantiles == 1:
                                if opt_subq.model == 'polynomial':
                                    if opt_subq.logX == 1:
                                        x1 = np.log(Xhsub[selsub1[i], selblock])
                                        xp1 = np.log(Xpc[sel1[i], selblock])
                                    else:
                                        x1 = Xhsub[selsub1[i], selblock]
                                        xp1 = Xpc[sel1(i), selblock]
                                else:
                                    x1 = Xhsub[selsub1[i], selblock]
                                    xp1 = Xpc[sel1[i], selblock]
                                if opt_subq.resmodel == 'polynomial':
                                    if opt_subq.logXres == 1:
                                        x1res = np.log(Xhsub[selsub1[i], selblock])
                                        xp1res = np.log(Xpc[sel1[i], selblock])
                                    else:
                                        x1res = Xhsub[selsub1[i], selblock]
                                        xp1res = Xpc[sel1[i], selblock]
                                else:
                                    x1res = Xhsub[selsub1[i], selblock]
                                    xp1res = Xpc[sel1[i], selblock]
                                q1 = np.squeeze(qXhfsub[selsub1[i], :, selblock])
                                Isort1 = np.squeeze(IsortXhfsub[selsub1[i], :, selblock])
                                if nselblock == 1:
                                    q1 = np.transpose(q1)
                                    Isort1 = np.transpose(Isort1)
                                if opt_subq.logY == 1:
                                    q1 = np.log(q1)
                                if opt_subq.model != 'ksmooth':
                                    bs1 = np.nan * np.ones((opt_subq.nparb, nselblock))
                                if opt_subq.resmodel and opt_subq.resmodel != 'ksmooth':
                                    bsres1 = np.nan * np.ones((opt_subq.nparbres, nselblock))
                                if opt_subq.bymonth == 1:
                                    m = monthp[sel1[i]]
                                else:
                                    m = 1
                                for k in np.arange(0, nqsub).reshape(-1):
                                    if opt_subq.model != 'ksmooth':
                                        bs1[:, :] = bs[m, k, :, selblock]
                                    if opt_subq.model == 'polynomial':
                                        if opt_subq.nparb == 4:
                                            yphat1 = (bs1[0, :] + np.multiply(bs1[1, :], xp1) +
                                                      np.multiply(bs1[2, :], xp1 ** 2) +
                                                      np.multiply(bs1[3, :], xp1 ** 3))
                                            yhat1 = (bs1[0, :] + np.multiply(bs1[1, :], x1) +
                                                     np.multiply(bs1[2, :], x1 ** 2) +
                                                     np.multiply(bs1[3, :], x1 ** 3))
                                        else:
                                            if opt_subq.nparb == 3:
                                                yphat1 = (bs1[0, :] + np.multiply(bs1[1, :], xp1) +
                                                          np.multiply(bs1[2, :], xp1 ** 2))
                                                yhat1 = (bs1[0, :] + np.multiply(bs1[1, :], x1) +
                                                         np.multiply(bs1[2, :], x1 ** 2))
                                            else:
                                                yphat1 = bs1[0, :] + np.multiply(bs1[1, :], xp1)
                                                yhat1 = bs1[0, :] + np.multiply(bs1[1, :], x1)
                                    else:
                                        if opt_subq.model == 'tanh':
                                            if opt_subq.nparb == 3:
                                                yphat1 = (bs1[0, :] +
                                                          np.multiply(bs1[1, :],
                                                                      np.tanh(np.multiply(bs1[2, :], xp1))))
                                                yhat1 = (bs1[0, :] +
                                                         np.multiply(bs1[1, :],
                                                                     np.tanh(np.multiply(bs1[2, :], x1))))
                                            else:
                                                yphat1 = np.multiply(bs1[0, :],
                                                                     np.tanh(np.multiply(bs1[1, :], xp1)))
                                                yhat1 = np.multiply(bs1[0, :],
                                                                    np.tanh(np.multiply(bs1[1, :], x1)))
                                        else:
                                            if opt_subq.model == 'powerlaw':
                                                if opt_subq.nparb == 2:
                                                    yphat1 = np.multiply(bs1[0, :], xp1 ** bs1[1, :])
                                                    yhat1 = np.multiply(bs1[0, :], x1 ** bs1[1, :])
                                                else:
                                                    if opt_subq.nparb == 3:
                                                        yphat1 = max(0, bs1[0, :] * (xp1 ** bs1[1, :]) + bs1[2, :])
                                                        yhat1 = max(0, bs1[0, :] * (x1 ** bs1[1, :]) + bs1[2, :])
                                                        # stop
                                            else:
                                                if opt_subq.model == 'ksmooth':
                                                    yphat1 = np.nan * np.ones((1, nselblock))
                                                    yhat1 = yphat1
                                                    yphat11 = np.nan * np.ones((len(opt_subq.xp), 1))
                                                    for l in np.arange(0, nselblock).reshape(-1):
                                                        yphat11 = yphatm[m, k, :, selblock[l]]
                                                        yphat1[l] = Xplin[xp1[l], opt_subq.xp, 1] * yphat11
                                                        yhat1[l] = Xplin[x1[l], opt_subq.xp, 1] * yphat11
                                    if opt_subq.resmodel:
                                        if opt_subq.resmodel != 'ksmooth':
                                            bsres1[:, :] = bsres[m, k, :, selblock]
                                        if opt_subq.resmodel == 'polynomial':
                                            if opt_subq.nparbres == 3:
                                                ypreshat1 = bsres1[0, :] + np.multiply(bsres1[1, :],
                                                                                       xp1res) + np.multiply(
                                                    bsres1[2, :], xp1res ** 2)
                                                yreshat1 = bsres1[0, :] + np.multiply(bsres1[1, :],
                                                                                      x1res) + np.multiply(bsres1[2, :],
                                                                                                           x1res ** 2)
                                            else:
                                                ypreshat1 = bsres1[0, :] + np.multiply(bsres1[1, :], xp1res)
                                                yreshat1 = bsres1[0, :] + np.multiply(bsres1[1, :], x1res)
                                            if opt_subq.logYres == 1:
                                                ypreshat1 = np.exp(ypreshat1)
                                                yreshat1 = np.exp(yreshat1)
                                            fac1 = ypreshat1 / yreshat1
                                        else:
                                            if opt_subq.resmodel == 'gaussian':
                                                ypreshat1 = np.multiply(bsres1[0, :], np.exp(
                                                    np.multiply(- (0.5 * bsres1[1, :] ** 2), xp1res ** 2)))
                                                yreshat1 = np.multiply(bsres1[0, :], np.exp(
                                                    np.multiply(- (0.5 * bsres1[1, :] ** 2), x1res ** 2)))
                                                fac1 = ypreshat1 / yreshat1
                                            else:
                                                if opt_subq.resmodel == 'powerlaw':
                                                    fac1 = (xp1res / x1res) ** bsres1[1, :]
                                                else:
                                                    if opt_subq.resmodel == 'ksmooth':
                                                        ypm1 = np.nan * np.ones((2, nselblock))
                                                        for l in np.arange(0, nselblock).reshape(-1):
                                                            ypm1[:, l] = Xplin[np.array([[xp1res[l]], [
                                                                x1res[l]]]), opt_subq.xpres, 1] * np.squeeze(
                                                                ypreshatm[m, k, :, selblock[l]])
                                                        fac1 = ypm1[0, :] / ypm1[1, :]
                                                    else:
                                                        raise Exception(np.array([
                                                            'Inflation factor fac1 not yet coded for ',
                                                            'opt_subq.resmodel = ',
                                                            opt_subq.resmodel]))
                                    else:
                                        fac1 = np.ones((1, nselblock))
                                    if fractional_subdelta == 1:
                                        q1[k, :] = (minXsub +
                                                    np.multiply(np.multiply((yphat1 - minXsub), fac1),
                                                                [q1[k, :] - minXsub]) / (yhat1 - minXsub))
                                    else:
                                        q1[k, :] = yphat1 + np.multiply(fac1, (q1[k, :] - yhat1))
                                if opt_subq.logY == 1:
                                    q1 = np.exp(q1)
                                for k in np.arange(0, nselblock).reshape(-1):
                                    if remove_deltascale_var == 1:
                                        Xpcf[self2[Isort1[:, k]], selblock[k]] = Xpci[self2[Isort1[:, k]], selblock[
                                            k]] + q1[:, k]
                                    else:
                                        Xpcf[self2[Isort1[:, k]], selblock[k]] = q1[:, k]
                            else:
                                if fractional_subdelta == 1:
                                    Xpcf[self2, :] = minXsub + np.multiply((np.ones((nself2, 1)) * (
                                        np.multiply((Xpc[sel1[i], :] - minXsub), correction_factor[sel1[i], :]))),
                                                                           dXhf1[np.arange(countf,
                                                                                           countf + nself2), :])
                                    if remove_deltascale_var == 1:
                                        Xpcf[self2, :] = Xpcf[self2, :] + (
                                                Xpci[self2, :] - np.ones((nself2, 1)) * Xpc[sel1[i], :])
                                else:
                                    Xpcf[self2, :] = Xpci[self2, :] + np.multiply(
                                        (np.ones((nself2, 1)) * correction_factor[sel1[i], :]),
                                        dXhf1[np.arange(countf + 1, countf + nself2), :])
                        countf = countf + nself2
            if verbose > 0:
                print('Done adding subdelta variability to projections')
        if dtdhf1 == 0:
            Xpcf[-1, :] = Xpcf[-2, :]
        # If hourly data are 00:00, 01:00, ..., 23:00 then the tdpf set above will result in an extra
        # time at 00:00 on ceil(max(tdp))+1 --- this can be filled by repeating the last datum.
        if correct_substd == 1:
            out.correction_factor = correction_factor
        # Impose lower/upper limits if provided
        if XpcfL != [] or XpcfH != []:
            if verbose > 0:
                print('Imposing lower/upper limits on projections')
            Xpcfo = Xpcf
            if XpcfL != []:
                Xpcf = np.amax(XpcfL, Xpcf)
            if XpcfH != []:
                Xpcf = np.amin(XpcfH, Xpcf)
            if Xpcfo != Xpcf:
                out.Xpcfo = Xpcfo
            del Xpcfo
            # Note we do not subselect for not-NaN (too slow), so this will set any NaN values to XpcL
        if use_XpcfH_SWI > 0:
            optSWIf = optSWI
            optSWIf.calc_dav = 0
            optSWIf.year = yearpf
            if use_SWI_hav == 1:
                optSWIf.calc_hav = 1
            if verbose > 0:
                print('Calculating maximum hourly values of shortwave irradiance')
            out.XpcfH_SWI = np.nan * np.ones((len(tdpf), ns))
            for i in np.arange(0, ns).reshape(-1):
                if use_SWI_hav == 1:
                    __, outSWI = calc_surfaceSWI.calc_surfaceSWI(latp[i], lonp[i], yrdaypf, optSWIf)
                    out.XpcfH_SWI[:, i] = outSWI.Q0hav
                else:
                    out.XpcfH_SWI[:, i] = calc_surfaceSWI.calc_surfaceSWI(latp[i], lonp[i], yrdaypf, optSWIf)
                if np.mod(i, 100) == 0 and verbose > 0:
                    print(np.array(['Done ', str(i), ' out of ', str(ns), ' time series']))
            if verbose > 0:
                print('Done maximum hourly values of shortwave irradiance')
            if use_XpcfH_SWI == 1:
                if verbose > 0:
                    print('Capping all hourly projections of shortwave irradiance '
                          'with instantaneous clear-sky maximum values')
                Xpcf = np.amin(out.XpcfH_SWI, Xpcf)

                # Warning: This is can lead to widespread correction of values near sunrise/sunset due to use of
                # hourly-average irradiance in ERA5

            else:
                if use_XpcfH_SWI == 2:
                    if verbose > 0:
                        print('Calculating daily maximum values of projected shortwave irradiance')
                    Xpcf_dmax = np.nan * Xpcf
                    tdu = np.unique(int(np.floor(tdpf)))
                    for i in np.arange(0, len(tdu)).reshape(-1):
                        self1 = np.argwhere(int(np.floor(tdpf)) == tdu[i])[:, 0]
                        Xpcf_dmax[self1, :] = np.ones((len(self1), 1)) * np.amax(Xpcf[self1, :])
                    if verbose > 0:
                        print('Done daily maximum values of projected shortwave irradiance')
                    if verbose > 0:
                        print(
                            'Capping selected hourly projections of shortwave irradiance ',
                            'with instantaneous clear-sky maximum values')
                        # Note this subselection can be slow
                    rec = np.argwhere(Xpcf[:] >= frcrit_capSWI * Xpcf_dmax[:])[:, 0]
                    Xpcf[rec] = np.amin(out.XpcfH_SWI[rec], Xpcf[rec])
        if verbose > 0:
            print('Done imposing lower/upper limits on finescale projections (if any)')
        # if verbose > 0:
        #     print('Done finescale projections')
        #     toc
        if recalculate_Xpc == 1:
            out.Xpco2 = Xpc
            if delta_timescale == 'daily':
                if verbose > 0:
                    print('Recalculating daily averages from sub-daily projections')
                for i in np.arange(0, ntp).reshape(-1):
                    Xpc[i, :] = np.mean(Xpcf[int(np.floor(tdpf)) == int(np.floor(tdp[i])), :])
                if verbose > 0:
                    print('Done recalculated daily averages')
            elif delta_timescale == 'monthly':
                if verbose > 0:
                    print('Recalculating monthly averages from sub-monthly projections')
                for i in np.arange(0, ntp).reshape(-1):
                    Xpc[i, :] = np.mean(Xpcf[yearpf == np.logical_and(yearp[i], monthpf) == monthp[i], :])
                if verbose > 0:
                    print('Done recalculated monthly averages')
                elif str(delta_timescale) == str('yearly') == 1:
                    if verbose > 0:
                        print('Recalculating yearly averages from sub-monthly projections')
                    for i in np.arange(0, ntp).reshape(-1):
                        Xpc[i, :] = np.mean(Xpcf[yearpf == yearp[i], :])
                    if verbose > 0:
                        print('Done recalculated yearly averages')
        out.Xpcf = Xpcf
        out.tdpf = tdpf
        if match_subdelta_deltamean == 1:
            out.selhsub = selhsub
        if match_subdelta_hourly == 1:
            out.selhsubhour = selhsubhour

    return Xpc, out
