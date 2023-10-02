
function [Xpc,out] = fn_biascorrect_projections(Xp,tdp,Xh,tdh,yearminref,yearmaxref,method,opt)
%
%function [Xpc,out] = fn_biascorrect_projections(Xp,tdp,Xh,tdh,yearminref,yearmaxref,method,opt)
%
%Computes bias-corrected projections (Xpc [ntimes*nseries]) given:
%   Xp  = input projection time series [ntimes*nseries]
%   tdp = input projection date numbers [ntimes,1]
%   Xh  = hindcast time series [ntimesh*nseries]
%   tdh = hindcast date numbers [ntimes,1]
%   yearminref,yearmaxref = reference period min/max year
%   method = bias correction method:
%            0 = delta-change (e.g. Kay and Butenschon, 2018)
%            1 = empirical quantile mapping (e.g. 'QM' method in Cannon et al., 2015)
%            2 = repeated hindcast correction (RHC, e.g. Willems and Vrac, 2011)
%
%Methods (0,1) start with the projection time series data (Xp) and aim to correct the statistical
%variability level for consistency with the hindcast time series, either by correcting climatological
%(e.g. monthly) mean values (for method=1, 'delta-change'), or by correcting the quantiles of variability
%(for method=1, 'empirical quantile mapping'). Method 2 takes a different approach, starting with the
%hindcast data, repeating it in blocks to fill the projection time points, and then perturbating these
%repeated data using climatic trends in quantiles (e.g. decadal-scale changes) that are calculated from
%the projection data. See text within function below for further advice about choice of method.
%
%By default, reference climatologies and quantiles are defined on a monthly basis (opt.use_month = 1).
%For other resolutions, set opt.use_month = 0 and provide a set of central climatological times opt.tclimc 
%as yeardays within range 0-365 (e.g. tclimc=0.5:365.5 for daily delta-change). We can also specify a
%tolerance of +/- tol_dtclim (default = 0) about the central times tclimc, to enable more
%robust calculations of  climatologies and quantile corrections. For example if we want to include
%(December,January,February) data when calculating means/quantiles for January data we set opt.use_month = 1
%and opt.tol_dtclim = 1
%
%Input scalar opt.fractional determines whether a fractional (1) or absolute (0, default) correction is
%applied to each series.
%
%For method=2 (Repeated Hindcast Correction) we also need to input the target blocks of years
%for which projections are required e.g. opt.yearpminv = 2000:10:2090; opt.yearpmaxv = 2009:10:2099.
%Hindcast data from the reference period (yearminref to yearmaxref inclusive) will be inserted to fill
%these blocks, and then corrected using calculated trends in quantiles (differences between the
%projection blocks calculated from the projection model data). The last day of the reference period
%will be duplicated to fill gaps resulting from different numbers of leap years per decade.
%
%The function will also add 'subdelta variability' (variability on a timescale shorter than that
%of the delta changes, i.e. the temporal resolution of the projection data), for example to achieve
%hourly projections given daily averages from an ESM. For this we must set opt.use_subdelta = 1 (def=0)
%and also supply fine-resolution hindcast time series data, with corresponding datenumbers 
%(opt.Xhf [ntimesf*nseries], opt.tdhf [ntimesf*1]).
%
%The subdelta variability is imported from the hindcast as a repeated multiannual chunk specified by
%(opt.yearminsub,opt.yearmaxsub), and is used to fill data for prediction chunks specified by
%(opt.yearpminv,opt.yearpminv), e.g. opt.yearpminv = 2000:10:2090; opt.yearpmaxv = 2009:10:2099.
%The subdelta variability is imposed either as additive corrections (opt.fractional_subdelta=0, default)
%or as multiplicative corrections (opt.fractional_subdelta=1) to the deltascale averages, assuming the 
%deltascale averages are either monthly or daily averages. This variability may be further corrected to 
%improve consistency with the projected deltascale average values, either by inflating/deflating the
%subdelta variance (opt.correct_substd=1) or by correcting the subdelta quantiles (opt.correct_subquantiles=1).
%By default, the subdelta variability is defined wrt a linear interpolation from the deltascale
%averages to the finescale time points (opt.remove_deltascale_var=1, default), but in some cases it may make more
%sense to define it wrt the deltascale averages as piecewise constant values (e.g. for imposing hourly
%variability on daily average shortwave irradiance as piecewise constant values we set opt.remove_deltascale_var=0). 
%
%If method=1 (Empirical Quantile Mapping), it is likely that the subdelta variability will need to be
%adjusted for consistency with the projected deltascale averages (opt.correct_substd=1 or opt.correct_subquantiles=1).
%In such cases, an options substructure opt.opt_subq is used to specify modelling options for correcting
%this subdelta variability. By default, a simple two-parameter linear regression model is a applied to
%model the standard deviation or deviation quantiles as a function of deltascale mean values, but several
%other options are available (consult code for details). If correcting deviation quantiles (opt.correct_subquantiles=1)
%it is also possible to model the residual variance of the deviation quantiles wrt the fitted model values,
%using an additional 'residual model' (see code for details).
%
%As an alternative to correcting the subdelta variability, the variability can be selected from the
%hindcast deltascale time period that best matches each of the projected deltascale means (opt.match_subdelta_deltamean=1).
%However, this approach can run into problems because i) sufficiently-similar deltascale periods might not be found
%in the hindcast dataset, when considering all spatial locations together (to preserve realistic spatial variability), or
%ii) the subselection procedure may generate large discontinuities between deltascale periods at the subdelta time scale
%(e.g. between 2300 and 0000 the next day).
%
%If method=2 (Repeated Hindcast Correction), it is more likely that the subdelta variability can be imposed
%without further correction, because the the projected changes between decades for each quantile should be
%relatively small (this is one major advantage of the RHC method).
%
%Corrected projections are capped at lower and/or upper limits if given as input (opt.XpcL,XpcH,XpcfL,XpcfH).
%Scalars can be used here as well as arrays of size(Xp) (e.g. input opt.XpcL = 0 for nonnegative variables).
%In addition, for surface shortwave irradiance, an upper limit can be calculated (and imposed) internally using astronomical
%formulae for the clear-sky surface shortwave irradiance (see code for details).
%
%By default, deltascale averages are recalculated from the finescale (subdelta) time series (opt.recalculate_Xpc=1, default).
%This can take some time so might be switched off for testing purposes.
%
%Uses:  fnperiodic_distm.m (calculate separations for periodic variables such as month/day-of-year)
%       Xplin.m (simple linear interpolation with optional constant extrapolation)
%       calc_surfaceSWI.m (calculates maximum daily/hourly surface shortwave irradiance, for use as upper limits on projected values)
%       ks_regress.m (used for fitting kernel-smoothing regression models if correcting subdelta variance/quantiles)
%       J_fit_partially_linear_model.m (used for fitting partially-linear regression models if correcting subdelta variance/quantiles)
%       optimize_theta.m (used for fitting partially-linear regression models if correcting subdelta variance/quantiles)
%
%%Phil Wallhead 10/05/2023


%Some General Advice About Choice Of Method
%------------------------------------------
%
%In general we would recommend either method=1 (Empirical Quantile Mapping) or method=2 (Repeated Hindcast Correction)
%as these will more thoroughly correct the level of variability in the projections. The simple delta-change approach
%(method=0) may clearly under-correct the projections in cases where the main temporal variability signal is not
%seasonal (e.g. with surface wind speeds), but might yet be optimal in cases where: i) seasonal variability dominates, and
%ii) the hindcast data are too few to allow robust estimation of quantiles (but does allow reasonable mean estimates).
%
%Between methods (1,2), method 1 (Empirical Quantile Mapping) might be preferred in cases where: 
%i) Subdelta variability is not a major concern; 
%ii) It is preferable to have projections that are closer to the original projections in character;
%iii) Artifacts due to interpolation of original projection data (leading to spatial oversmoothing),
%potential loss of multivariate correlation structure due to the univariate implementation, and
%potential distortion of climatic trends due to quantile correction, are not such of a concern. 
%Note that more elaborate implementations of quantile mapping have been developed to address the issues in iii) separately,
%though to our knowledge not yet to address all together.
%iv) It is desired to preserve the contribution of internal variability to the uncertainty in climatic change 
%between decades.
%
%Alternatively, method 2 (Repeated Hindcast Correction) might be preferred in cases where: 
%i) Subdelta variability is significant and strongly dependent on deltascale mean values. RHC is preferable
%over EQM is such cases because the perturbations to deltascale mean values (e.g. monthly means, or daily means)
%that are imposed by RHC due to interdecadal/climatic changes are often smaller than those imposed by EQM to correct 
%the statistical variability level of the original projection data. Therefore it is more likely that the subdelta
%variability can be imported from the hindcast dataset with no corrections (or more minor corrections) to account
%for changes in the deltascale mean values.
%ii) It is preferable to have projections that are closer to the hindcast data in character. The RHC method
%puts less 'faith' in the model used to produce the original projection data, in the sense that only certain
%aspects of the original model data --- i.e. the climatic/decadal-scale variability --- is used in the final
%corrected projections.
%iii) Artifacts due to spatial interpolation, inaccuracy in the multivariate (cross)-correlation structure of the
%original projection data, and loss of cross-correlations or distortion of climatic trends due to the simple EQM 
%approach provided herein, are a major cause for concern. The RHC method preserves the full multivariate spatio-temporal
%variability and (cross) correlation structure of the hindcast data, at the cost of introducing some moderate temporal
%discontinuity between repetition blocks (in our experience not very noticeable), and some decadal-scale periodicity in the
%corrected projections that may be considered artifactual/spurious.
%iv) Loss of uncertainty due to varying internal variability between (multi)decadal repetition block is not a major
%concern. This may be the case when the corrected projections are to be used for comparing two (multi)decadal periods
%(present-day vs. future). In such cases, the fact that the internal year-to-year variation is identical in both periods
%may actually be an advantage, allowing a more precise estimate of climatic change (i.e. closer to what would be obtained
%by averaging over many realizations of the internal variability).


[ntp,ns] = size(Xp);
out = [];
if nargin<8; opt = []; end
if isfield(opt,'verbose')==1; verbose = opt.verbose; else verbose = 1; end
if isfield(opt,'seasonal')==1; seasonal = opt.seasonal; else seasonal = 1; end
if isfield(opt,'use_month')==1; use_month = opt.use_month; else use_month = 1; end
%1 to use month-of-year (1-12) rather than yrday (day-of-year 0-365) to compute climatological statistics and deltas.
%%%if isfield(opt,'dyrday')==1; dyrday = opt.dyrday; else dyrday = 7; end
if isfield(opt,'tclimc')==1; tclimc = opt.tclimc; else tclimc = 1:12; end
%Climatological time centres used to calculate climatological (means,deltas) and quantiles (def = 1:12 for monthly delta change with use_month=1).
if isfield(opt,'ndays_year')==1; ndays_year = opt.ndays_year; else ndays_year = 365.25; end
%Best approach here seems to be to use ndays_year = 365.25; if we use ndays_year = 365 then we
%will combine data from Dec 31st and Jan 1st when the preceding year is a leap year.
if isfield(opt,'tol_dtclim')==1; tol_dtclim = opt.tol_dtclim; else tol_dtclim = 0; end
%Tolerance level for climatological time used in computing climatological statistics (def = 0 for simple monthly delta change).
% Using tol_dtclim>0 may allow a more robust definition of climatology, e.g. if the no. of reference years (yearmaxref-yearminref+1) is limited.
% E.g. if aiming for a daily climatology (tclimc=(0.5:365.5)), it is probably a good idea to allow e.g. tol_dtclim = 5 (i.e. +/- 5 days tolerance)
% in order to allow robust climatological statistics and consequent deltas.
if isfield(opt,'tol_dtclimp')==1; tol_dtclimp = opt.tol_dtclimp; else tol_dtclimp = 0; end
%Tolerance level for climatological time used in correcting projections (def = 0 for simple monthly delta change).
% This may be different to tol_dtclim because tol_dtclim may be expanded to allow more robust climatological statistics. 
if isfield(opt,'fractional')==1; fractional = opt.fractional; else fractional = 0; end
if isfield(opt,'minX')==1; minX = opt.minX; else minX = []; end
%Reference minimum value to avoid potential division by ~0 in fractional methods.
if isfield(opt,'allow_block_lag')==1; allow_block_lag = opt.allow_block_lag; else allow_block_lag = 1; end

if (size(Xh,2)>0 && ns==1)
    %In case of multiple hindcast series and one projection series, apply the single projection series to all hindcast climatologies
    ns = size(Xh,2); Xp = Xp*ones(1,ns);
end

if (method==0)
    if isfield(opt,'correct_iavar')==1; correct_iavar = opt.correct_iavar; else correct_iavar = 0; end
    if isfield(opt,'dyear_iavar')==1; dyear_iavar = opt.dyear_iavar; else dyear_iavar = 10; end
    if isfield(opt,'use_coef_iavar')==1; use_coef_iavar = opt.use_coef_iavar; else use_coef_iavar = 1; end
    %use_coef_iavar=1: Correct the interannual coefficient of variation rather than standard deviation, only for fractional deltas (def)
    %use_coef_iavar=2: Correct the interannual coefficient of variation rather than standard deviation for all deltas
end

if (method==1)
    if isfield(opt,'qsel')==1; qsel = opt.qsel; else qsel = []; end
    qsel = qsel(:);
    %Optional set of fixed quantiles to interpolate to (e.g. 0.005:0.01:0.995) --- otherwise the sorted data are used.
    %Experience to date suggests that it is usually not helpful to interpolate to fixed quantiles --- better to leave qsel empty.
end

if (method==2)
    if isfield(opt,'qsel')==1; qsel = opt.qsel; else qsel = 0.05:0.1:0.95; end
    qsel = qsel(:);
end


if isfield(opt,'XpcL')==1; XpcL = opt.XpcL; else XpcL = []; end %Optional lower limits for Xpc
if isfield(opt,'XpcH')==1; XpcH = opt.XpcH; else XpcH = []; end %Optional upper limits for Xpc
if isfield(opt,'use_XpcH_SWI')==1; use_XpcH_SWI = opt.use_XpcH_SWI; else use_XpcH_SWI = 0; end 
%If opt.use_XpcH_SWI=1, maximum daily-average values of shortwave irradiance are calculated using calc_surfaceSWI.m
%and these are used to impose an upper limit on Xpc
if (use_XpcH_SWI==1)
    if isfield(opt,'latp')==1; latp = opt.latp; else latp = []; end %Latitudes for each time series (degN)
    if isfield(opt,'lonp')==1; lonp = opt.lonp; else lonp = []; end %Longitudes for each time series (degE)
    if isfield(opt,'optSWI')==1; optSWI = opt.optSWI; else optSWI = []; end %Options for maximum (clear-sky) SWI calculation.
    if (isempty(optSWI))
        optSWI = struct('f_model',1,'decl_model',1,'use_eqtime',1,'Q_model',0,'cloud_model',0,'ndays_year',365.2422,'use_tday_lag',1);
        %Default options for calc_surfaceSWI.m: These have been found to give best performance (minimum overshoot)
        %in tests using hourly ERA5 data for ROHO800 model (western Norway).
    end
    if isfield(opt,'use_SWI_hav')==1; use_SWI_hav = opt.use_SWI_hav; else use_SWI_hav = 1; end 
    %1 (def) to use hourly averages rather than instantaneous values in calculation of hourly data (consistent with ERA5).
end

if isfield(opt,'yearpminv')==1; yearpminv = opt.yearpminv; else yearpminv = []; end
if isfield(opt,'yearpmaxv')==1; yearpmaxv = opt.yearpmaxv; else yearpmaxv = []; end
%Optional input of set of prediction year periods over which the subdelta variability will be resampled
%and/or deltascale variability will be repeated (if method=2).

if isfield(opt,'legacy')==1; legacy = opt.legacy; else legacy = 0; end
if (legacy==1)
    allow_block_lag = 0;
end

if isfield(opt,'use_subdelta')==1; use_subdelta = opt.use_subdelta; else use_subdelta = 0; end
%1 to use subdelta variability from input finescale hindcast (opt.Xhf,opt,tdhf) to correct the projections
if (use_subdelta==1)
    if isfield(opt,'tdpfmin')==1; tdpfmin = opt.tdpfmin; else tdpfmin = floor(min(tdp)); end
    if isfield(opt,'tdpfmax')==1; tdpfmax = opt.tdpfmax; else tdpfmax = ceil(max(tdp)); end
    if isfield(opt,'Xhf')==1; Xhf = opt.Xhf; else Xhf = []; end
    if isfield(opt,'tdhf')==1; tdhf = opt.tdhf; else tdhf = []; end
    if isfield(opt,'yearminsub')==1; yearminsub = opt.yearminsub; else yearminsub = yearminref; end
    if isfield(opt,'yearmaxsub')==1; yearmaxsub = opt.yearmaxsub; else yearmaxsub = yearmaxref; end
    %Optional input of minimum/maximum years (yearminsub,yearmaxsub) to sample subdelta variability from
    %-- by default set equal to (yearminsub,yearmaxsub).
    if isfield(opt,'fractional_subdelta')==1; fractional_subdelta = opt.fractional_subdelta; else fractional_subdelta = 0; end
    %fractional_subdelta = 1 to treat the subdelta variability as a fractional perturbation (better for non-negative variables)
    if isfield(opt,'minXsub')==1; minXsub = opt.minXsub; else minXsub = 0; end
    %Reference minimum value to avoid potential division by ~0 in fractional subdelta methods.
    if isfield(opt,'loop_over_projection_times')==1; loop_over_projection_times = opt.loop_over_projection_times; else loop_over_projection_times = 0; end
    %Using matlab's kron.m below gives identical results to the loop calculation, but it ~800x faster
    if isfield(opt,'correct_substd')==1; correct_substd = opt.correct_substd; else correct_substd = 0; end
    if isfield(opt,'correct_subquantiles')==1; correct_subquantiles = opt.correct_subquantiles; else correct_subquantiles = 0; end
    %opt.correct_subquantiles = 1 to model the subdelta-scale quantiles (e.g. hourly) as functions of delta-scale means (e.g. daily means)
    if isfield(opt,'remove_deltascale_var')==1; remove_deltascale_var = opt.remove_deltascale_var; else remove_deltascale_var = 1; end

    if (correct_substd==1||correct_subquantiles==1)
        if isfield(opt,'opt_subq')==1; opt_subq = opt.opt_subq; else opt_subq = []; end
        %Options structure to specify the models for subdelta quantiles
        %Defaults:
        if (isempty(opt_subq)==1)
            if verbose>0; disp('Using default options for models to correct subdelta quantiles'); end
            opt_subq = struct('model','polynomial','logX',0,'logY',0,'nparb',2,'bymonth',1,'resmodel',[],'doplot',1);
        end
        if ~isfield(opt_subq,'logX'); opt_subq.logX = 0; end
        if ~isfield(opt_subq,'logY'); opt_subq.logY = 0; end
        if ~isfield(opt_subq,'nonnegb'); opt_subq.nonnegb = 0; end
        if ~isfield(opt_subq,'nr'); opt_subq.nr = 1; end
        if ~isfield(opt_subq,'logXres'); opt_subq.logXres = 0; end
        if ~isfield(opt_subq,'logYres'); opt_subq.logYres = 0; end
        if ~isfield(opt_subq,'nonnegbres'); opt_subq.nonnegbres = 0; end
        if ~isfield(opt_subq,'nrres'); opt_subq.nrres = 1; end
        
        if (correct_subquantiles==1)
            if ~strcmp(opt_subq.model,'polynomial'); opt_subq.logX = 0; opt_subq.logY = 0; end
            if ~strcmp(opt_subq.resmodel,'polynomial'); opt_subq.logXres = 0; opt_subq.logYres = 0; end
            %Log transforms only applicable to polynomial models
            if strcmp(opt_subq.resmodel,'gaussian')==1; opt_subq.nparbres = 2; end
        end
    end

    if isfield(opt,'match_subdelta_deltamean')==1; match_subdelta_deltamean = opt.match_subdelta_deltamean; else match_subdelta_deltamean = 0; end
    %1 to select the subdelta variability from the hindcast day/month that best matches the projected deltascale mean
    if (match_subdelta_deltamean==1 || correct_subquantiles==1); loop_over_projection_times = 1; end
    if isfield(opt,'match_subdelta_hourly')==1; match_subdelta_hourly = opt.match_subdelta_hourly; else match_subdelta_hourly = 1; end
    if match_subdelta_deltamean==0; match_subdelta_hourly = 0; end
    if isfield(opt,'XpcfL')==1; XpcfL = opt.XpcfL; else XpcfL = []; end %Optional lower limits for Xpcf
    if isfield(opt,'XpcfH')==1; XpcfH = opt.XpcfH; else XpcfH = []; end %Optional upper limits for Xpcf 
    if isfield(opt,'use_XpcfH_SWI')==1; use_XpcfH_SWI = opt.use_XpcfH_SWI; else use_XpcfH_SWI = 0; end
    %If opt.use_XpcfH_SWI>0, maximum hourly values of shortwave irradiance are calculated using calc_surfaceSWI.m
    %and these are used to impose an upper limit on Xpcf. Then if:
    %opt.use_XpcfH_SWI = 1: All hourly projections are limited by the maximum (clear-sky) values.
    %opt.use_XpcfH_SWI = 2: Only hourly projections >=opt.frcrit_capSWI*daily maximum are limited by the clear-sky values.
    if (use_XpcfH_SWI==2)
        if isfield(opt,'frcrit_capSWI')==1; frcrit_capSWI = opt.frcrit_capSWI; else frcrit_capSWI = 0.25; end
    end
    if isfield(opt,'recalculate_Xpc')==1; recalculate_Xpc = opt.recalculate_Xpc; else recalculate_Xpc = 1; end
    %opt.recalculate_Xpc = 1: Recalculate the daily average projections using the hourly projections after limits imposed.
    
    if isfield(opt,'subdelta_by_blocks')==1; subdelta_by_blocks = opt.subdelta_by_blocks; else subdelta_by_blocks = 0; end
    %1 to divide subdelta variability calculation into latitude-longitude blocks, to obtain better matches to hindcast daily means (match_subdelta_deltamean=1).
    %This is probably not a good idea unless the target variable has very limited spatial correlation between grid cells,
    %since the blockwise calculation will not preserve spatial correlations across block boundaries.
    if (subdelta_by_blocks==1)
        if isfield(opt,'latp')==1; latp = opt.latp; else latp = []; end %Latitudes for each time series (degN)
        if isfield(opt,'lonp')==1; lonp = opt.lonp; else lonp = []; end %Longitudes for each time series (degE)
        if isfield(opt,'latL_blocks')==1; latL_blocks = opt.latL_blocks; else latL_blocks = []; end %Lower latitude boundaries of blocks to divide calculation by (one row).
        if isfield(opt,'latH_blocks')==1; latH_blocks = opt.latH_blocks; else latH_blocks = []; end %Upper latitude boundaries of blocks to divide calculation by (one row).
        if isfield(opt,'lonL_blocks')==1; lonL_blocks = opt.lonL_blocks; else lonL_blocks = []; end %Lower longitude boundaries of blocks to divide calculation by (one row).
        if isfield(opt,'lonH_blocks')==1; lonH_blocks = opt.lonH_blocks; else lonH_blocks = []; end %Upper longitude boundaries of blocks to divide calculation by (one row).
        latL_blocks1 = kron(latL_blocks(:)',ones(1,length(lonL_blocks))); %Lower latitude boundaries for all blocks.
        latH_blocks1 = kron(latH_blocks(:)',ones(1,length(lonH_blocks))); %Upper latitude boundaries for all blocks.
        lonL_blocks1 = kron(ones(1,length(latL_blocks)),lonL_blocks(:)'); %Lower longitude boundaries for all blocks.
        lonH_blocks1 = kron(ones(1,length(latH_blocks)),lonH_blocks(:)'); %Upper longitude boundaries for all blocks.
        nblocks = length(latL_blocks1);
    else
        nblocks = 1;
    end
end



[yearp,monthp,~] = datevec(tdp);
yrdayp = tdp - datenum(yearp,1,1);
if seasonal==1; ntclimc = length(tclimc); else ntclimc = 1; end

if isempty(minX)
    minX = min([Xp(:); Xh(:)]) - 1*(max([Xp(:); Xh(:)])-min([Xp(:); Xh(:)]));
    %Set reference minimum to minimum minus the range --- this will avoid potential division by ~0,
    %and also avoids inaccuracy due to small variations very far from zero causing correction factors
    %very close to 1 (e.g. for pressure).
    %NOTE: minX is only used for correction when fractional = 1.
    %NOTE: Since the ratios are defined by e.g. (qXh-minX)/(qXp-minX), the above minX will restrict
    %the possible range of ratio values between 0.5 and 2.
end
out.minX = minX;

if method==0; methodstr = 'delta-change'; end
if method==1; methodstr = 'empirical quantile mapping'; end
if method==2; methodstr = 'repeated hindcast correction'; end
if verbose>0; disp(['Calculating bias-corrected projections using ',methodstr,' method']);tic; end



if (method==0 || method==1)
    %Calculate projection climatologies / reference quantiles
    nseltclim = NaN*ones(ntclimc,1);
    if (method==0)
        Xpclim = NaN*ones(ntclimc,ns);
        if correct_iavar==1; Xpstdclim = Xpclim; end
    end
    if method==1; qXpref = cell(ntclimc,1); end
    for i=1:ntclimc %Loop over climatological time divisions (e.g. month-of-year)
        if (seasonal==1)
            %Calculate temporal separations in a periodic sense using fnperiodic_distm.m
            if use_month==1; dist1 = fnperiodic_distm(monthp,tclimc(i),12); else dist1 = fnperiodic_distm(yrdayp,tclimc(i),ndays_year); end
% REDO            selt = find(dist1<=tol_dtclim); %Choose all data within tolerance tol_dtclim of climatological time centre tclimc(i)
        else
            selt = 1:ntp;
        end
        
% REDO        seltclim = selt(yearp(selt)>=yearminref & yearp(selt)<=yearmaxref);
        nseltclim(i) = length(seltclim);
        
        if (method==0) %Delta-change
            Xpclim(i,:) = mean(Xp(seltclim,:)); %Projection climatology during reference period
            
            if (correct_iavar==1) %Calculate the interannual variability of the input projection over the baseline period (Xstd0)
                X1 = [ones(nseltclim(i),1) year(seltclim)];
                b1 = X1\Xp(seltclim,:); %Linear regression using 'backslash'
                Xphat1 = X1*b1;
                r1 = Xp(seltclim,:) - Xphat1; %Residuals defined by linear regression
                Xpstdclim(i,:) = std(r1);
                if (use_coef_iavar==1 && fractional==1); Xpstdclim(i,:) = std(r1./Xphat1); end
                if use_coef_iavar==2; Xpstdclim(i,:) = std(r1./Xphat1); end
            end
            
        elseif (method==1) %Quantile Mapping
            if ~isempty(qsel); qXpref{i} = quantile(Xp(seltclim,:),qsel); else qXpref{i} = sort(Xp(seltclim,:)); end
            %Interpolate to fixed set of quantiles if qsel is specified, otherwise use the sorted values.
        end
    end
    out.nseltclim = nseltclim;
    if (method==0)
        out.Xpclim = Xpclim;
        if (correct_iavar==1); out.Xpstdclim = Xpstdclim; end
    end
    if method==1; out.qXpref = qXpref; end
end



if (method==2)
    %Calculate projection model trends in quantiles between repeated hindcast periods
    nperiods = length(yearpminv);
    nq = length(qsel);
    qXp = NaN*ones(ntclimc,nperiods,nq,ns); 
    for i=1:ntclimc %Loop over climatological time divisions (e.g. month-of-year)
        if (seasonal==1)
            %Calculate temporal separations in a periodic sense using fnperiodic_distm.m
            if use_month==1; dist1 = fnperiodic_distm(monthp,tclimc(i),12); else dist1 = fnperiodic_distm(yrdayp,tclimc(i),ndays_year); end
% REDO            selt0 = find(dist1<=tol_dtclim); %Choose all data within tolerance tol_dtclim of climatological time centre tclimc(i)
        else
            selt0 = 1:ntp;
        end
        
        for j=1:nperiods %Loop over repeated hindcast periods
% REDO            selt = selt0(yearp(selt0)>=yearpminv(j) & yearp(selt0)<=yearpmaxv(j));
            qXp(i,j,1:nq,1:ns) = quantile(Xp(selt,:),qsel);
        end
    end
    out.qXp = qXp;
end




%%Apply the delta changes / quantile corrections to generate bias-corrected projections.
[yearh,monthh,~] = datevec(tdh);
if use_month==0; yrdayh = tdh - datenum(yearh,1,1); end

% REDO seltclimh0 = find(yearh>=yearminref & yearh<=yearmaxref);
nthref = length(seltclimh0);

if all(diff(yearh)==1)
    delta_timescale = 'yearly';
elseif all(diff(monthh)==1|diff(monthh)==-11)
    delta_timescale = 'monthly';
elseif all(diff(tdh)==1)   
    delta_timescale = 'daily';
else
    delta_timescale = [];
end    
    


%Preallocation for hindcast statistic arrays
nseltclimh = NaN*ones(ntclimc,1); Xpc = NaN*ones(ntp,ns);
if (method==0)
    Xhclim = NaN*ones(ntclimc,ns); 
    if correct_iavar==1; Xhstdclim = Xhclim; end
elseif (method==1)
    qXhref = cell(ntclimc,1);
elseif (method==2)
    Fhref = NaN*ones(nthref,ns);    
end

%First we calculate the climatology (or reference quantiles) for the hindcast model.
%For simple delta change, the corrected projections (Xpc) are also calculated within the same loop.
for i=1:ntclimc
    if (seasonal==1)
        if (use_month==1)
            dist1 = fnperiodic_distm(monthh(seltclimh0),tclimc(i),12);
            dist1p = fnperiodic_distm(monthp,tclimc(i),12);
        else
            dist1 = fnperiodic_distm(yrdayh(seltclimh0),tclimc(i),ndays_year);
            dist1p = fnperiodic_distm(yrdayp,tclimc(i),ndays_year);
        end
% REDO        seltclimh = seltclimh0(dist1<=tol_dtclim); %Choose all data within tolerance tol_dtclim of climatological time centre tclimc(i)
% REDO        seltp = find(dist1p<=tol_dtclimp); %Choose all data within tolerance tol_dtclimp of climatological time centre tclimc(i)
    else
        seltclimh = seltclimh0;
        seltp = 1:ntp;
    end
    nseltclimh(i) = length(seltclimh);

    if (method==0) %Delta-change
        Xhclim(i,:) = mean(Xh(seltclimh,:)); %Hindcast climatology during reference period
        nseltp = length(seltp);

        %Here we calculate simple delta-change projections with no correction of interannual variance or quantiles.
        if (fractional==1)
            Xpc(seltp,:) = minX + (Xp(seltp,:)-minX).*(ones(nseltp,1)*Xhclim(i,:)-minX)./(ones(nseltp,1)*Xpclim(i,:)-minX);
        else
            Xpc(seltp,:) = ones(nseltp,1)*Xhclim(i,:) + Xp(seltp,:) - ones(nseltp,1)*Xpclim(i,:);
        end

        if (correct_iavar==1) %Calculate the interannual standard deviation of the hindcast over the baseline period (Xstd0h)
            Xh1 = [ones(nseltclimh(i),1) yearh(seltclimh)];
            bh1 = Xh1\Xh(seltclimh,:);
            Xhhat1 = Xh1*bh1;
            rh1 = Xh(seltclimh,:) - Xhhat1; %Residuals defined by linear regression
            Xhstdclim(i,:) = std(rh1);
            if (use_coef_iavar==1 && fractional==1); Xhstdclim(i,:) = std(rh1./Xhhat1); end
            if use_coef_iavar==2; Xhstdclim(i,:) = std(rh1./Xhhat1); end
        end

    elseif (method==1) %Quantile correction
        if ~isempty(qsel)
            qXhref{i} = quantile(Xh(seltclimh,:),qsel);
        else
            if nseltclimh(i)<nseltclim(i) %If fewer reference data from hindcast, repeat hindcast data
                seltclimhc = [seltclimh; seltclimh(1:(nseltclim(i)-nseltclimh(i)))];
            elseif nseltclimh(i)>nseltclim(i) %If fewer reference data from projection, subset hindcast data
                seltclimhc = seltclimh(1:nseltclim(i));
            else
                seltclimhc = seltclimh;
            end
            qXhref{i} = sort(Xh(seltclimhc,:));
        end

    elseif (method==2) %Repeated Hindcast Correction
        %Here we need to record the cdf value F for each time point in the repeated hindcast series
% REDO        if seasonal==1; selthref1 = find(dist1<=tol_dtclimp); else selthref1 = 1:nthref; end
        selth1 = seltclimh0(selthref1);
        nh1 = length(selthref1);

        for j=1:nh1
            for k=1:ns
                Fhref(selthref1(j),k) = (sum(Xh(seltclimh,k)<Xh(selth1(j),k))+0.5)/nseltclimh(i);
                %NOTE: This is consistent with how quantiles are calculated within quantile.m.
                %      There is no need for interpolation because the values Xh(seltclimh1(j),k) will
                %      always correspond exactly to one of the values in Xh(seltclimh,k).
                %NOTE: For leap years, the final year-day (yrday=365.5) will be compared to a subset
                %      seltclimh that is equally large as for the other year days, as long as the
                %      tolerance tol_dtclim>0.
            end
        end
    end
end
out.nseltclimh = nseltclimh;
if method==0; out.Xhclim = Xhclim; end
if method==1; out.qXhref = qXhref; end
if method==2; out.Fhref = Fhref; end


%Next we apply the correction factors to the projections (if not applied already for simple delta change)
if (method==0 && correct_iavar==1) %Delta-change projections corrected for interannual variance.
    %Here we inflate the residuals from linear regression of ndivdec decadal/multiannual periods,
    %using the ratio of standard deviations of linear trend residuals during the baseline period.
    yearpmin = min(yearp); yearpmax = max(yearp);
    yearpspan = yearpmax - yearpmin;
    ndivdec = round(yearpspan/dyear_iavar);
    Rstdclim = Xhstdclim./Xpstdclim;

    Xpco = Xpc; %Delta-change projections without correction of interannual variance
    for i=1:ntclimc
        if (seasonal==1)
            if use_month==1; dist1p = fnperiodic_distm(monthp,tclimc(i),12); else dist1p = fnperiodic_distm(yrdayp,tclimc(i),ndays_year); end
% REDO            seltp = find(dist1p<=tol_dtclimp); %Choose all data within tolerance tol_dtclimp of climatological time centre tclimc(i)
        else
            seltp = 1:ntp;
        end

        for j=1:ndivdec %Loop over decadal divisions to apply linear detrending and adjust residual interannual variance.
            if (j==ndivdec)
                seltp2 = seltp(yearp(seltp)>=yearpmin+(j-1)*dyear_iavar);
            else
% REDO                seltp2 = seltp(yearp(seltp)>=yearpmin+(j-1)*dyear_iavar & yearp(seltp)<yearpmin+j*dyear_iavar);
            end
            nseltp2 = length(seltp2);
            Xp1 = [ones(nseltp2,1) tdp(seltp2)];
            bp1 = Xp1\Xp(seltp2,:);
            Xphat1 = Xp1*bp1;
            rp1 = Xpco(seltp2,:) - Xphat1; %Residuals defined by linear regression
            Xpc(seltp2,:) = Xphat1 + (ones(nseltp2,1)*Rstdclim(i,:)).*rp1; %Inflate the interannual variance using ratio of stds
            if (use_coef_iavar==1 && fractional==1); Xpc(seltp2,:) = Xphat1.*(1 + Rstdclim(i,:).*rp1./Xphat1); end
            if use_coef_iavar==2; Xpc(seltp2,:) = Xphat1.*(1 + Rstdclim(i,:).*rp1./Xphat1); end
        end
    end
    out.Xhstdclim = Xhstdclim; out.Rstdclim = Rstdclim;
    
elseif (method==1) %Quantile correction
    %Here we adjust the projections using factors derived by interpolating the ratios/differences
    %of quantiles during the reference period. This is an implementation of simple empirical
    %quantile mapping (QM, see Cannon et al., 2015) which does not preserve climatic trends,
    %and which potentially lead to unrealistic distortion of climate change signals (e.g. Cannon et al., 2015).
    %However, such distortion effects should be minimal as long as the projected climatic trends are weak
    %relative to the full standard deviation of ESM variability during the reference period.
    Rqref = cell(ntclimc,1);
    for i=1:ntclimc
        if (seasonal==1)
            if use_month==1; dist1p = fnperiodic_distm(monthp,tclimc(i),12); else dist1p = fnperiodic_distm(yrdayp,tclimc(i),ndays_year); end
% REDO            seltp = find(dist1p<=tol_dtclimp); %Choose all data within tolerance tol_dtclimp of climatological time centre tclimc(i)
        else
            seltp = 1:ntp;
        end
        qXpref1 = qXpref{i}; qXhref1 = qXhref{i};
        Rqref1 = (qXhref1 - minX) ./ (qXpref1 - minX); %Use minX to stabilize the ratio, see note below.
        Rqref{i} = Rqref1;
        
        for j=1:ns
            qXpref11 = qXpref1(:,j);
            
            if (fractional==1)
                Rqref11 = Rqref1(:,j);
                if (min(diff(qXpref11))==0) %Replace any duplicated nodal values with average for linear interpolation
                    qXpref11o = qXpref11; Rqref11o = Rqref11;
                    qXpref11 = unique(qXpref11); Rqref11 = NaN*qXpref11;
                    for k=1:length(qXpref11); Rqref11(k) = mean(Rqref11o(qXpref11o==qXpref11(k))); end
                end
                Xpc(seltp,j) = minX + (Xp(seltp,j)-minX) .* (Xplin(Xp(seltp,j),qXpref11,1)*Rqref11);
                %Apply quantile correction, interpolating between correction factors.
                %NOTE: We extrapolate a constant (1 passed to Xplin.m), meaning that the correction factors for values beyond the
                %      range of variability in the reference period are assigned the correction factor for the
                %      lowest or highest value during the reference period.
                %NOTE: A stability constant (minX) is used to stabilize the calculation.
                %      The calculation is not very sensitive to the value of this constant, as long as it is large
                %      enough to avoid ~0 values in the denominator of Rq011.
            else
                dqref11 = qXhref1(:,j) - qXpref11;
                if (min(diff(qXpref11))==0) %Replace any duplicated nodal values with average for linear interpolation
                    qXpref11o = qXpref11; dqref11o = dqref11;
                    qXpref11 = unique(qXpref11); dqref11 = NaN*qXpref11;
                    for k=1:length(qXpref11); dqref11(k) = mean(dqref11o(qXpref11o==qXpref11(k))); end
                end
                Xpc(seltp,j) = Xp(seltp,j) + Xplin(Xp(seltp,j),qXpref11,1)*dqref11;
                %Apply quantile correction, interpolating between additive corrections (with constant extrapolation)
            end
            %NOTE: Neither of these quantile corrections are strictly trend-preserving, since the corrections 
            %      depend on the original values of the projections; hence e.g. the correction to the median will
            %      vary (slightly) over time as the median value increases/decreases.
        end
    end
    out.qXhref = qXhref; out.Rqref = Rqref;
    
elseif (method==2) %Repeated hindcast correction

    iref = find(yearpminv==yearminref & yearpmaxv==yearmaxref);
    dtdhref = tdh(seltclimh0) - datenum(yearminref,1,1); %Days in hindcast since start of reference period

    dqXp = NaN*qXp;
    for i=1:nperiods %Loop over repeated time blocks
% REDO        sel1 = find(yearp>=yearpminv(i) & yearp<=yearpmaxv(i));
        nt1 = length(sel1);

        if (strcmp(delta_timescale,'daily')==1 && allow_block_lag==1)
            dtdp1 = tdp(sel1(1)) - datenum(yearpminv(i),1,1); %Days since start of projection block i
% REDO            ihref1 = find(abs(dtdhref-dtdp1)<=tol_dtclimp); %Starting index within repeated hindcast series
            %For example if the projection time series only runs from Jan 2nd of the projection block,
            %we should start the repeated hindcast block also from Jan 2nd.
        else
            ihref1 = 1;
        end
        selth1 = seltclimh0(ihref1:min(nt1+(ihref1-1),nthref)); %Indices within Xh, tdh etc. of the repeated hindcast data for this block
        Fhref1 = Fhref(ihref1:min(nt1+(ihref1-1),nthref),:); %Quantiles (cdf values) of the repeated hindcast data for this block

        if (i==iref && 0==1)
            Xpc(sel1,:) = Xh(seltclimh0,:);
        else

            for j=1:ntclimc
                if (seasonal==1)
                    if use_month==1; dist1 = fnperiodic_distm(monthh(selth1),tclimc(j),12); else dist1 = fnperiodic_distm(yrdayh(selth1),tclimc(j),ndays_year); end
                    if (legacy==1)
% REDO                        selthref1 = find(dist1<=tol_dtclim);
                    else
% REDO                        selthref1 = find(dist1<=tol_dtclimp); %Choose all data within tolerance tol_dtclimp of climatological time centre tclimc(i)
                        %Note: This should pick up only one time index: tol_dtclimp is only to allow for small
                        %      discrepancies e.g. daily average timestamps not located exactly at midday.
                    end
                else
                    selthref1 = 1:length(selth1);
                end
                seltp1 = (sel1(1)-1)+selthref1;

                for l=1:ns
                    if (fractional==1)
                        dqXp1 = (squeeze(qXp(j,i,:,l))-minX)./(squeeze(qXp(j,iref,:,l))-minX);
                        Xpc(seltp1,l) = minX + (Xh(selth1(selthref1),l)-minX) .* Xplin(Fhref1(selthref1,l),qsel(:),1)*dqXp1(:);
                    else
                        dqXp1 = squeeze(qXp(j,i,:,l)) - squeeze(qXp(j,iref,:,l));
                        Xpc(seltp1,l) = Xh(selth1(selthref1),l) + Xplin(Fhref1(selthref1,l),qsel(:),1)*dqXp1(:);
                    end
                    %Interpolate correction linearly from the fixed quantiles (e.g. 0.005:0.01:0.995) to the
                    %quantile/CDF value represented by the reference hindcast datum, with constant extrapolation.
                    dqXp(j,i,:,l) = dqXp1; %Record perturbation factors for analysis purposes
                end
            end
            if (strcmp(delta_timescale,'daily')==1 && isnan(Xpc(sel1(end),1)))
                if verbose>0; disp(['Repeating last day of hindcast reference period to apply to years ',int2str(yearpminv(i)),' to ',int2str(yearpmaxv(i))]); end
                Xpc(sel1(end),:) = Xpc(sel1(end-1),:);
            elseif nt1>(nthref+1)
                error(['Not yet coded re-use of hindcast data for nt1 = ',int2str(nt1),' and nthref = ',int2str(nthref)])
            end
        end
        if verbose>0; disp(['Done repeated hindcast corrected projections for years ',int2str(yearpminv(i)),' through ',int2str(yearpmaxv(i))]); tic; end
    end
    out.dqXp = dqXp;
end



%Impose lower/upper limits if provided
if (~isempty(XpcL)||~isempty(XpcH))
    Xpco = Xpc; %Can be a useful diagnostic to check that capping is not too severe.
    if ~isempty(XpcL); Xpc = max(XpcL,Xpc); end
    if ~isempty(XpcH); Xpc = min(XpcH,Xpc); end
    if ~isequaln(Xpco,Xpc); out.Xpco = Xpco; end
    clear Xpco
end

if (use_XpcH_SWI==1)
    if verbose>0; disp('Calculating maximum daily-average values of shortwave irradiance'); end
    optSWI.calc_dav = 1; optSWI.year = yearp;
    [~,outSWI] = calc_surfaceSWI(latp,lonp,yrdayp,optSWI);
    if verbose>0; disp('Done maximum daily-average values of shortwave irradiance'); end
    out.XpcH_SWI = outSWI.Q0dav;
    Xpc = min(out.XpcH_SWI,Xpc); %Cap the corrected daily-average projections with clear-sky maximum values
end

if verbose>0; disp(['Done bias-corrected projections using ',methodstr,' method']); toc; end





%%%If required, add 'subdelta' variability at the extra temporal resolution of the hindcast (e.g. hourly)
if (use_subdelta==1 && ~isempty(Xhf) && ~isempty(tdhf))
    if verbose>0; disp('Calculating finescale projections with subdelta variability'); tic; end
    
    if (correct_subquantiles==1)
        if strcmp(delta_timescale,'daily')
            nqsub = 24;
        elseif strcmp(delta_timescale,'monthly')
            nqsub = 28;
        else
            error('Subdelta variability quantiles not yet defined if hindcast data opt.Xh neither daily nor monthly')
        end
    end
    
% REDO    seltsub = find(yearh>=yearminsub & yearh<=yearmaxsub);
    Xhsub = Xh(seltsub,:); tdhsub = tdh(seltsub); %Hindcast output during subdelta period
    [yearhsub,monthhsub,dayhsub] = datevec(tdhsub);
    ntsub = length(tdhsub);

    [yearhf,monthhf,dayhf] = datevec(tdhf);
% REDO    seltfsub = find(yearhf>=yearminsub & yearhf<=yearmaxsub);
    Xhfsub = Xhf(seltfsub,:); tdhfsub = tdhf(seltfsub); %Finescale hindcast output during subdelta period
    [yearhfsub,monthhfsub,~] = datevec(tdhfsub);
    ntfsub = length(tdhfsub);

    if (remove_deltascale_var==1)
        Xhi = interp1(tdh,Xh,tdhf,'linear','extrap'); %Deltascale variation interpolated to fine scale
        dXhf = Xhf - Xhi;
        clear Xhi
        %Note: Here we allow linear extrapolation, consistent with projection baseline variability (Xpci) calculated below.
    end
    
    if correct_subquantiles==0; dXhfsub = NaN*Xhfsub; end
    if correct_substd==1; dXhfstdsub = NaN*Xhsub; dXhfminsub = NaN*Xhsub; dXhfmaxsub = NaN*Xhsub; end
    if correct_subquantiles==1; qXhfsub = NaN*ones(ntsub,nqsub,ns); IsortXhfsub = qXhfsub; end
    m = 0;
    for i=1:ntsub
        if (strcmp(delta_timescale,'daily')==1) %Hourly variability about daily delta scale
            sel1 = find(yearhf==yearhsub(i) & monthhf==monthhsub(i) & dayhf==dayhsub(i));
        elseif (strcmp(delta_timescale,'monthly')==1) %Hourly/daily variability about monthly delta scale
            sel1 = find(yearhf==yearhsub(i) & monthhf==monthhsub(i));
        elseif (strcmp(delta_timescale,'yearly')==1) %Hourly/daily variability about yearly delta scale
            sel1 = find(yearhf==yearhsub(i));
        end
        sel2 = m+1:m+length(sel1);

        if (remove_deltascale_var==1)
            if (correct_subquantiles==1)
                Xhfc1 = dXhf(sel1,:);
                %These are the finescale anomalies about the deltascale variation (day-to-day, or month-to-month).
            else
                Xhfc1 = dXhf(sel1,:) + Xhsub(i,:);
                %This is the finescale variation with the signal from deltascale variation (day-to-day, or month-to-month) removed.
            end
        else
            Xhfc1 = Xhf(sel1,:);
        end

        if (correct_subquantiles==1)
            [qXhfsub(i,:,:),IsortXhfsub(i,:,:)] = sort(Xhfc1);
        else
            if (fractional_subdelta==1)
                dXhfsub(sel2,:) = (Xhfc1-minXsub) ./ (ones(length(sel1),1)*Xhsub(i,:)-minXsub);
                %Check: mean(dXhfsub(sel2,:)) %This should be = 1 for all time series IFF remove_deltascale_var=0 and the hindcast daily/monthly averages were calculated by simple averages of all hourly/daily data within each day/month.
            else
                dXhfsub(sel2,:) = Xhfc1 - Xhsub(i,:);
                %Check: mean(dXhfsub(sel2,:)) %This should be = 0 for all time series IFF remove_deltascale_var=0 and the hindcast daily/monthly averages were calculated by simple averages of all hourly/daily data within each day/month.
            end
        end

        %Record subdelta statistics
        if (correct_substd==1)
            dXhfstdsub(i,:) = std(dXhfsub(sel2,:)); dXhfminsub(i,:) = min(dXhfsub(sel2,:)); dXhfmaxsub(i,:) = max(dXhfsub(sel2,:));
        end
        m = m + length(sel1);
    end
    if correct_substd==1; out.dXhfstdsub = dXhfstdsub; out.dXhfminsub = dXhfminsub; out.dXhfmaxsub = dXhfmaxsub; out.dXhfrangesub = dXhfmaxsub-dXhfminsub; end
    if correct_subquantiles==1; out.IsortXhfsub = IsortXhfsub; end

    if (correct_substd==1)
        if opt_subq.bymonth==1; ndivt_subq = 12; str1 = ' by month-of-year'; else ndivt_subq = 1; str1 = []; end
        if verbose>0; disp(['Fitting ',opt_subq.model,' models',str1,' to correct subdelta standard deviation as function of delta-scale mean']); end
        
        if strcmp(opt_subq.model,'ksmooth')
            yphatm = NaN*ones(ndivt_subq,length(opt_subq.xp),ns); 
            optk = struct('kfunc',opt_subq.kfunc,'lr_order',opt_subq.lr_order);
            if isfield(opt_subq,'ndclosest'); optk.ndclosest = opt_subq.ndclosest; end
        else
            bs = NaN*ones(ndivt_subq,opt_subq.nparb,ns);
        end

        ifig0 = 130;
        for k=1:ndivt_subq %Loop over months (if opt_subq.bymonth=1)
            if opt_subq.bymonth==1; reck = find(monthhsub==k); else reck = 1:ntsub; end
            for i=1:ns
                rec = reck(~isnan(Xhsub(reck,i)) & ~isnan(dXhfstdsub(reck,i)));

                if (strcmp(opt_subq.model,'polynomial')==1)
                    if opt_subq.logX==1; rec = rec(Xhsub(rec,i)>0); end
                    if opt_subq.logY==1; rec = rec(dXhfstdsub(rec,i)>0); end
                    nrec = length(rec);
                    if opt_subq.logX==1; x1 = log(Xhsub(rec,i)); else x1 = Xhsub(rec,i); end
                    if opt_subq.logY==1; y1 = log(dXhfstdsub(rec,i)); else y1 = dXhfstdsub(rec,i); end
                    if opt_subq.nparb==3; Xreg1 = [ones(nrec,1) x1 x1.^2]; else Xreg1 = [ones(nrec,1) x1]; end
                    %Fit model for subdelta variability std as function of daily mean values
                    %[bs,bsse,Cbs,bsint,statss] = fit_m(y1,Xreg1,struct('doplotset',1,'selxshow',1)); stop;
                    bs1 = Xreg1\y1; %Linear regression using 'backslash'
                elseif (strcmp(opt_subq.model,'ksmooth')==1)
                    x1 = Xhsub(rec,i); y1 = dXhfstdsub(rec,i);
                    yphatm(k,:,i) = max(0, ks_regress(y1,x1,opt_subq.xp,opt_subq.bdw,optk));
                    %NOTE: We kernel-smooth to a fixed set of prediction nodes, then interpolate these later with
                    %      constant extrapolation. Re-use of the raw input data is generally too slow.
                end

                if ~strcmp(opt_subq.model,'ksmooth'); bs(k,:,i) = bs1; end

                if (opt_subq.doplot==1 && i==1)
                    if (strcmp(opt_subq.model,'polynomial')==1)
                        if opt_subq.bymonth==1; Xpcmin1 = min(Xpc(monthp==k,i)); Xpcmax1 = max(Xpc(monthp==k,i)); else Xpcmin1 = min(Xpc(:,i)); Xpcmax1 = max(Xpc(:,i)); end
                        xp1 = linspace(Xpcmin1,Xpcmax1,200)';
                        if opt_subq.logX==1; xp1t = log(xp1); else xp1t = xp1; end
                        if (opt_subq.nparb==4)
                            yphat1t = bs1(1)+bs1(2)*xp1t+bs1(3)*xp1t.^2+bs1(4)*xp1t.^3;
                        elseif (opt_subq.nparb==3)
                            yphat1t = bs1(1)+bs1(2)*xp1t+bs1(3)*xp1t.^2;
                        else
                            yphat1t = bs1(1)+bs1(2)*xp1t;
                        end
                        if opt_subq.logY==1; yphat1 = exp(yphat1t); else yphat1 = yphat1t; end
                    elseif (strcmp(opt_subq.model,'ksmooth')==1)
                        xp1 = opt_subq.xp; yphat1 = squeeze(yphatm(k,:,i))';
                    end

                    figure(ifig0); 
                    subplot(3,4,k);plot(Xhsub(rec,i),dXhfstdsub(rec,i),'k.',xp1,yphat1,'r-','LineWidth',2)
                    xlabel(opt_subq.xstr,'FontSize',11)
                    ylabel(opt_subq.ystr,'FontSize',11)
                    title(opt_subq.tstr{k},'FontWeight','bold','FontSize',11)
                    if isfield(opt_subq,'xlim'); set(gca,'XLim',opt_subq.xlim); end
                end
            end
            if (opt_subq.bymonth==1 && verbose>0); disp(['Done ',int2str(k),' of 12 months']); end
        end
        if ~strcmp(opt_subq.model,'ksmooth'); out.bs = bs; end

        clear x1 y1 Xreg1
        if verbose>0; disp('Done fitting regression models for subdelta standard deviation'); end

    elseif (correct_subquantiles==1)

        if opt_subq.bymonth==1; ndivt_subq = 12; str1 = ' by month-of-year'; else ndivt_subq = 1; str1 = []; end
        if verbose>0; disp(['Fitting ',opt_subq.model,' models',str1,' to correct subdelta quantiles as function of delta-scale mean']); end
        
        if strcmp(opt_subq.model,'ksmooth')
            yphatm = NaN*ones(ndivt_subq,nqsub,length(opt_subq.xp),ns); 
            optk = struct('kfunc',opt_subq.kfunc,'lr_order',opt_subq.lr_order);
            if isfield(opt_subq,'ndclosest'); optk.ndclosest = opt_subq.ndclosest; end
        else
            bs = NaN*ones(ndivt_subq,nqsub,opt_subq.nparb,ns); Jmins = NaN*ones(ndivt_subq,nqsub,ns);
        end
        if strcmp(opt_subq.resmodel,'ksmooth')
            ypreshatm = NaN*ones(ndivt_subq,nqsub,length(opt_subq.xpres),ns); 
            optkres = struct('kfunc',opt_subq.kfuncres,'lr_order',opt_subq.lr_orderres);
            if isfield(opt_subq,'ndclosestres'); optkres.ndclosest = opt_subq.ndclosestres; end
        elseif ~isempty(opt_subq.resmodel) 
            bsres = NaN*ones(ndivt_subq,nqsub,opt_subq.nparbres,ns); Jminsres = NaN*ones(ndivt_subq,nqsub,ns); 
        end

        ifig0 = 150; ifig1 = 170;
        if ns>1; opt_subq.doplot = 0; end

        for k=1:ndivt_subq %Loop over months (if opt_subq.bymonth=1)
            if opt_subq.bymonth==1; reck = find(monthhsub==k); else reck = 1:ntsub; end

            %if k>1; opt_subq.doplot=0; end
            if (opt_subq.doplot==1)
                figure(ifig0+k); clf;
                if ~isempty(opt_subq.resmodel); figure(ifig1+k); clf; end
            end

            for j=1:nqsub %Loop over quantiles of subdelta variability
                for i=1:ns %Loop over spatial locations

                    rec = reck;
                    nrec = length(rec);
                    if (strcmp(opt_subq.model,'polynomial')==1)
                        if opt_subq.logX==1; rec = rec(Xhsub(rec,i)>0); end
                        if opt_subq.logY==1; rec = rec(squeeze(qXhfsub(rec,j,i)>0)); end
                        nrec = length(rec);
                        if opt_subq.logX==1; x1 = log(Xhsub(rec,i)); else x1 = Xhsub(rec,i); end
                        if opt_subq.logY==1; y1 = squeeze(log(qXhfsub(rec,j,i))); else y1 = squeeze(qXhfsub(rec,j,i)); end
                        if (opt_subq.nparb==4)
                            Xreg1 = [ones(nrec,1) x1 x1.^2 x1.^3]; 
                        elseif (opt_subq.nparb==3)
                            Xreg1 = [ones(nrec,1) x1 x1.^2]; 
                        else 
                            Xreg1 = [ones(nrec,1) x1]; 
                        end
                        %Fit model for subdelta quantiles as function of daily mean values
                        %[bs1,bs1se,Cbs1,bs1int,statss1] = fit_m(y1,Xreg1,struct('doplotset',0,'selxshow',1));
                        bs1 = Xreg1\y1; %Linear regression using 'backslash'
                        if ~isempty(opt_subq.resmodel); yhat1 = Xreg1*bs1; end
                    elseif (strcmp(opt_subq.model,'tanh')==1)
                        x1 = Xhsub(rec,i); y1 = squeeze(qXhfsub(rec,j,i));
                        bs100 = 1/max(x1);
                        if opt_subq.nparb==3; Xfun = @(p) [ones(nrec,1) tanh(p(1)*x1)]; else Xfun = @(p) tanh(p(1)*x1); end
                        optJ = struct('Y',y1,'Xfun',Xfun,'nonnegb',opt_subq.nonnegb); opto = struct('optJ',optJ);
                        phat = optimize_theta(bs100,opt_subq.pL,opt_subq.pH,@J_fit_partially_linear_model,opt_subq.nr,opto);                        
                        [Jmins1,outJ] = J_fit_partially_linear_model(phat,optJ);
                        bs1 = [outJ.b; phat];
                        if ~isempty(opt_subq.resmodel); yhat1 = outJ.Ym; end
                        Jmins(k,j,i) = Jmins1;
                    elseif (strcmp(opt_subq.model,'powerlaw')==1)
                        x1 = Xhsub(rec,i); y1 = squeeze(qXhfsub(rec,j,i));
                        if (opt_subq.nparb==2)
                            bs100 = x1\y1;
                            bs10 = [log(bs100); 1];
                            Xfun = @(p) x1.^p(1);
                            optJ = struct('Y',y1,'Xfun',Xfun,'nonnegb',opt_subq.nonnegb); opto = struct('optJ',optJ);
                            phat = optimize_theta(bs10(2),opt_subq.pL,opt_subq.pH,@J_fit_partially_linear_model,opt_subq.nr,opto);
                            [Jmins1,outJ] = J_fit_partially_linear_model(phat,optJ);
                            bs1 = [outJ.b phat];
                            if ~isempty(opt_subq.resmodel); yhat1 = outJ.Ym; end
                            Jmins(k,j,i) = Jmins1;
                        elseif (opt_subq.nparb==3)
                            bs100 = x1\y1;
                            bs10 = [log(bs100); 1; 0];
                            Jfun1 = @(p) max(0,p(1)*x1.^p(2)+p(3)) - y1;
                            %Xfun = @(p) max(0,x1-p(2)).^p(1);
                            bs1 = lsqnonlin(Jfun1,bs10);
                            if ~isempty(opt_subq.resmodel); yhat1 = Jfun1(bs1) + y1; end
                        end
                    elseif (strcmp(opt_subq.model,'ksmooth')==1)
                        x1 = Xhsub(rec,i); y1 = squeeze(qXhfsub(rec,j,i));
                        yphatm(k,j,:,i) = ks_regress(y1,x1,opt_subq.xp,opt_subq.bdw,optk);
                        %NOTE: We kernel-smooth to a fixed set of prediction nodes, then interpolate these later with
                        %      constant extrapolation. Re-use of the raw input data is generally too slow.
                        %if ~isempty(opt_subq.resmodel); yhat1 = ks_regress(y1,x1,x1,opt_subq.bdw,optk); end
                        if ~isempty(opt_subq.resmodel); yhat1 = Xplin(x1,opt_subq.xp,1)*squeeze(yphatm(k,j,:,i)); end
                    end

                    if ~strcmp(opt_subq.model,'ksmooth'); bs(k,j,:,i) = bs1; end

                    if ~isempty(opt_subq.resmodel)
                        x1res = Xhsub(rec,i);
                        if fractional_subdelta==1; y1res = abs(y1-yhat1)./yhat1; else y1res = abs(y1-yhat1); end
                    end
                    if (strcmp(opt_subq.resmodel,'polynomial')==1)
                        rec2 = 1:nrec;
                        if opt_subq.logXres==1; rec2 = rec2(x1res(rec2)>0); end
                        if opt_subq.logYres==1; rec2 = rec2(y1res(rec2)>0); end
                        nrec2 = length(rec2);
                        x1res = x1res(rec2); y1res = y1res(rec2);
                        if opt_subq.logXres==1; x1rest = log(x1res); else x1rest = x1res; end
                        if opt_subq.logYres==1; y1rest = log(y1res); else y1rest = y1res; end
                        if (opt_subq.nparbres==3)
                            Xres1 = [ones(nrec2,1) x1rest x1rest.^2];
                        else
                            Xres1 = [ones(nrec2,1) x1rest];
                        end
                        bsres1 = Xres1\y1rest;
                    elseif (strcmp(opt_subq.resmodel,'gaussian')==1)
                        bsres100 = std(x1res);
                        Xfun = @(p) exp(-(0.5*p(1)^2)*x1res.^2); %p(1) here is 1/(standard deviation)
                        optJ = struct('Y',y1res,'Xfun',Xfun,'nonnegb',opt_subq.nonnegbres); opto = struct('optJ',optJ);
                        preshat = optimize_theta(bsres100,opt_subq.presL,opt_subq.presH,@J_fit_partially_linear_model,opt_subq.nrres,opto);
                        [Jmins1,outJ] = J_fit_partially_linear_model(preshat,optJ);
                        bsres1 = [outJ.b preshat];
                        Jminsres(k,j,i) = Jmins1;
                    elseif (strcmp(opt_subq.resmodel,'powerlaw')==1)
                        bsres100 = x1res\y1res;
                        bsres10 = [bsres100; 1];
                        Jfun2 = @(p) p(1)*x1res.^p(2) - y1res;
                        bsres1 = lsqnonlin(Jfun2,bsres10);
                    elseif (strcmp(opt_subq.resmodel,'ksmooth')==1)
                        ypreshatm(k,j,:,i) = max(0, ks_regress(y1res,x1res,opt_subq.xpres,opt_subq.bdwres,optkres));
                        %NOTE: We kernel-smooth to a fixed set of prediction nodes, then interpolate these later with
                        %      constant extrapolation. Re-use of the raw input data is generally too slow.
                    end

                    if (~isempty(opt_subq.resmodel) && ~strcmp(opt_subq.resmodel,'ksmooth')); bsres(k,j,:,i) = bsres1; end

                    if (opt_subq.doplot==1)
                        if opt_subq.bymonth==1; Xpcmin1 = min(Xpc(monthp==k,i)); Xpcmax1 = max(Xpc(monthp==k,i)); else Xpcmin1 = min(Xpc(:,i)); Xpcmax1 = max(Xpc(:,i)); end
                        xp1 = linspace(Xpcmin1,Xpcmax1,200)';
                        
                        if (strcmp(opt_subq.model,'polynomial')==1)
                            if opt_subq.logX==1; xp1t = log(xp1); else xp1t = xp1; end
                            if (opt_subq.nparb==4)
                                yphat1t = bs1(1)+bs1(2)*xp1t+bs1(3)*xp1t.^2+bs1(4)*xp1t.^3;
                            elseif (opt_subq.nparb==3)
                                yphat1t = bs1(1)+bs1(2)*xp1t+bs1(3)*xp1t.^2; 
                            else 
                                yphat1t = bs1(1)+bs1(2)*xp1t; 
                            end
                            if opt_subq.logY==1; yphat1 = exp(yphat1t); else yphat1 = yphat1t; end
                        elseif (strcmp(opt_subq.model,'tanh')==1)
                            if opt_subq.nparb==3; yphat1 = bs1(1)+bs1(2)*tanh(bs1(3)*xp1); else yphat1 = bs1(1)*tanh(bs1(2)*xp1); end
                        elseif (strcmp(opt_subq.model,'powerlaw')==1)
                            if opt_subq.nparb==2; yphat1 = bs1(1)*xp1.^bs1(2); end
                            if opt_subq.nparb==3; yphat1 = max(0,bs1(1)*xp1.^bs1(2)+bs1(3)); end
                        elseif (strcmp(opt_subq.model,'ksmooth')==1)
                            yphat1 = Xplin(xp1,opt_subq.xp,1)*squeeze(yphatm(k,j,:,i));
                        end

                        if (strcmp(opt_subq.resmodel,'polynomial')==1)
                            if opt_subq.logXres==1; xp1t = log(xp1); else xp1t = xp1; end
                            if (opt_subq.nparbres==3)
                                ypreshat1t = bsres1(1) + bsres1(2)*xp1t + bsres1(3)*xp1t.^2;
                            else
                                ypreshat1t = bsres1(1) + bsres1(2)*xp1t;
                            end
                            if opt_subq.logYres==1; ypreshat1 = exp(ypreshat1t); else ypreshat1 = ypreshat1t; end
                        elseif (strcmp(opt_subq.resmodel,'gaussian')==1)
                            ypreshat1 = bsres1(1)*exp(-(0.5*bsres1(2)^2)*xp1.^2);
                        elseif (strcmp(opt_subq.resmodel,'powerlaw')==1)
                            ypreshat1 = bsres1(1)*xp1.^bsres1(2);
                        elseif (strcmp(opt_subq.resmodel,'ksmooth')==1)
                            ypreshat1 = Xplin(xp1,opt_subq.xpres,1)*squeeze(ypreshatm(k,j,:,i));
                        end

                        figure(ifig0+k); subplot(4,6,j);plot(Xhsub(rec,i),squeeze(qXhfsub(rec,j,i)),'k.',xp1,yphat1,'r-','LineWidth',2)
                        if ~isempty(opt_subq.resmodel); figure(ifig1+k); subplot(4,6,j);plot(x1res,y1res,'k.',xp1,ypreshat1,'r-','LineWidth',2); end
                    end
                end
            end
            if (opt_subq.bymonth==1 && verbose>0); disp(['Done ',int2str(k),' of 12 months']); end
        end
        if ~strcmp(opt_subq.model,'ksmooth'); out.bs = bs; end
        if (~isempty(opt_subq.resmodel) && ~strcmp(opt_subq.resmodel,'ksmooth')); out.bsres = bsres; end
        if any(strcmp(opt_subq.model,{'powerlaw','tanh'})); out.Jmins = Jmins; end
        %if opt_subq.doplot==1; stop; end

        clear x1 y1 Xreg1
        if verbose>0; disp('Done fitting regression models for subdelta quantiles'); end
    end

    
    %Next calculate the baseline variability and prepare for main loop over projection periods
    %First we need to establish the time stamps for the subdelta projection data:
    if (isempty(yearpminv) || isempty(yearpmaxv))
        nyrssub = yearmaxsub-yearminsub+1;
        yearpminv = min(yearp):nyrssub:max(yearp); yearpmaxv = min(yearp)+nyrssub-1:nyrssub:max(yearp);
        %If not matching the subdelta variability by most-similar days, then we need to define
        %a set of multiannual year periods for which the subdelta variability will be repeated over
        %(usually most convenient to divide the prediciton period into 10 or 20-year periods).
        %If not defined by input these periods are defined here.
        %Although it should not be necessary to split the full prediction in match_subdelta_deltamean=1,
        %tests show that the calculation is almost 2x faster if it is split up (with identical final results).
    end
% REDO    selp = find(yearpmaxv<=max(yearp)); yearpminv = yearpminv(selp); yearpmaxv = yearpmaxv(selp);
    out.yearpminv = yearpminv; out.yearpmaxv = yearpmaxv;
    
    dtdhf1 = min(tdhf-floor(tdhf)); %Initial increment after 00:00
    if (strcmp(delta_timescale,'daily')==1)
        difftdhf = 1/24;
        tdpf = (tdpfmin+dtdhf1:difftdhf:tdpfmax+1)';
        if verbose>0;
            disp('Assuming regular hourly increments for subdelta variability');
        end %This is to avoid problems in computation of daily averages due to potentially-limited precision
    elseif ((max(diff(tdhf))-min(diff(tdhf)))<1e-3)
        difftdhf = mean(diff(tdhf)); %Finescale time increment (in days, assumed constant)
        tdpf = (tdpfmin+dtdhf1:difftdhf:tdpfmax+1)';
        if verbose>0; disp(['Assuming regular increments ',num2str(difftdhf,'%3.2f'),' day(s) for subdelta variability']); end 
    else
        %In this case we allow some irregularity in the subdelta time increments (as in e.g. xCO2 hindcast data);
        %we repeat the subdelta time increments relative to the start of each prediction block.
        if verbose>0; disp('Building projection subdelta timestamps by repeating increments from beginning of subdelta reference period'); end 
        m = 0;
        for i=1:length(yearpminv)
            tdpf(m+1:m+ntfsub) = datenum(yearpminv(i),1,1) + (tdhfsub-datenum(yearminsub,1,1));
            m = m + ntfsub;
        end
        tdpf = tdpf(:);
    end
    ntpf = length(tdpf);
    [yearpf,monthpf,~] = datevec(tdpf);
    yrdaypf = tdpf - datenum(yearpf,1,1);
    Xpcf = NaN*ones(ntpf,ns); %Preallocate final 'finescale' projection time series, with subdelta variability added

    
    if verbose>0; disp('Calculating baseline variability to which subdelta variation will be added'); end
    if (remove_deltascale_var>0)
        Xpci = interp1(tdp,Xpc,tdpf,'linear','extrap');
        %Note: Here we allow linear extrapolation, consistent with definition of subdaily variations above.
        %Note: use of Xplin.m (Xpci = Xplin(tdpf,tdp,1)*Xpc) demands too much memory.
    else
        Xpci = Xpcf;
        for i=1:ntp
            if (strcmp(delta_timescale,'daily')==1)
                self1 = find(floor(tdpf)==floor(tdp(i)));
            elseif (strcmp(delta_timescale,'monthly')==1)
                self1 = find(yearpf==yearp(i) & monthpf==monthp(i));
            elseif (strcmp(delta_timescale,'yearly')==1)
                self1 = find(yearpf==yearp(i));
            end
            Xpci(self1,:) = ones(length(self1),1)*Xpc(i,:);
        end
    end
    if verbose>0; disp('Done baseline variability'); end

    if match_subdelta_deltamean==1; selhsub = NaN*ones(ntp,nblocks); end %Indices of selected most-similar days from the daily hindcast time series
    if match_subdelta_hourly==1; selhsubhour = NaN*ones(ntp,nblocks); end %Indices of selected first-hours that minimize jumps between days.
    correction_factor = ones(ntp,ns);

    dtdhsub = tdh(seltsub) - datenum(yearminsub,1,1); %Days in hindcast since start of subdelta reference period

    %Main loop over projection periods for which finescale projections are required
    for j=1:length(yearpminv)
% REDO        sel1 = find(yearp>=yearpminv(j) & yearp<=yearpmaxv(j));
% REDO        self1 = find(yearpf>=yearpminv(j) & yearpf<=yearpmaxv(j));
        nt1 = length(sel1); %No. of projection times in current time chunk
        ntf1 = length(self1); %No. of finescale projection times in current time chunk

        if (match_subdelta_deltamean==0) %If not matching by day, identify the hindcast indices to repeat, and calculate correcion factors (if not correcting quantiles).
            
            if (strcmp(delta_timescale,'daily')==1 && allow_block_lag==1)
                %Identify starting index for the repeated hindcast variability
                dtdp1 = tdp(sel1(1)) - datenum(yearpminv(j),1,1); %Days since start of projection block j
% REDO                ihsub1 = find(abs(dtdhsub-dtdp1)<=tol_dtclimp); %Starting index within repeated hindcast series (deltascale)
            else
                ihsub1 = 1;
            end

            %Identify repeated hindcast series to use for this block
            selsub1 = ihsub1:min(ihsub1-1+nt1,ntsub);
            if (correct_subquantiles==0) %Identify indices for finescale data and repeated finescale data for this block
                %Note: This is not necesssary if correct_subquantiles=1 because in this case only the deltascale indices selsub1 are used
                if (strcmp(delta_timescale,'daily')==1)
                    ihfsub1 = find(floor(tdhfsub)==floor(tdhsub(selsub1(1))),1,'first');
                    ihfsub2 = find(floor(tdhfsub)==floor(tdhsub(selsub1(end))),1,'last');
                elseif (strcmp(delta_timescale,'monthly')==1)
                    ihfsub1 = find(yearhfsub==yearhsub(selsub1(1)) & monthhfsub==monthhsub(selsub1(1)),1,'first');
                    ihfsub2 = find(yearhfsub==yearhsub(selsub1(end)) & monthhfsub==monthhsub(selsub1(end)),1,'last');  
                elseif (strcmp(delta_timescale,'yearly')==1)
                    ihfsub1 = find(yearhfsub==yearhsub(selsub1(1)),1,'first');
                    ihfsub2 = find(yearhfsub==yearhsub(selsub1(end)),1,'last');                      
                end
                dXhf1 = dXhfsub(ihfsub1:ihfsub2,:);

                if (strcmp(delta_timescale,'daily')==0 && length(dXhf1(:,1))==(ntf1-1))
                    dXhf1 = [dXhf1; dXhf1(end,:)]; %#ok<AGROW> %Repeat last entry if necessary due to varying no. leap years per block.
                    if verbose>0; disp(['Repeating last entry of reference period subdelta variability to apply to years ',int2str(yearpminv(j)),' to ',int2str(yearpmaxv(j))]); end
                end
            end

            %Add extra day if necessary to account for possible different no. leap years per decade (e.g. only 2 in subdelta period 2010-2019)
            if (strcmp(delta_timescale,'daily')==1 && (ihsub1-1+nt1>ntsub))
                if (ihsub1-1+nt1==ntsub+1)
                    if verbose>0; disp(['Repeating last day of subdelta reference period to apply to years ',int2str(yearpminv(j)),' to ',int2str(yearpmaxv(j))]); end
                    selsub1 = [selsub1 ntsub]; %#ok<AGROW> %Repeat last day of subdelta period
                    %The use of 'selsub1' avoids having to append several variables ('Xhsub' etc.) for use below.
                    if correct_subquantiles==0; dXhf1 = [dXhf1; dXhf1(end-23:end,:)]; end %#ok<AGROW>
                else
                    error(['Problem finding enough days of data from subdelta reference period to apply to years ',int2str(yearpminv(j)),' to ',int2str(yearpmaxv(j))])
                end
            end

            if (correct_substd==1) %We can also precalculate correction factors for the subdelta standard deviation
                if verbose>0; disp('Calculating correction factors for subdelta variance'); end

                if ~strcmp(opt_subq.model,'ksmooth'); bs1 = NaN*ones(opt_subq.nparb,ns); end
                correction_factor1 = NaN*ones(nt1,ns);
                for k=1:ndivt_subq %Loop over months (if opt_subq.bymonth=1)
                    if opt_subq.bymonth==1; rec1 = find(monthp(sel1)==k); else rec1 = 1:nt1; end
                    nrec1 = length(rec1);

                    if strcmp(opt_subq.model,'polynomial')
                        bs1(:,:) = bs(k,:,:);
                        if opt_subq.logX==1; x1 = log(Xhsub(selsub1(rec1),:)); xp1 = log(Xpc(sel1(rec1),:)); else x1 = Xhsub(selsub1(rec1),:); xp1 = Xpc(sel1(rec1),:); end
                        if (opt_subq.logY==1)
                            if (opt_subq.nparb==3)
                                correction_factor11 = exp((ones(nrec1,1)*bs1(2,:)).*(xp1-x1) + (ones(nrec1,1)*bs1(3,:)).*(xp1.^2-x1.^2));
                            else
                                correction_factor11 = exp((ones(nrec1,1)*bs1(2,:)).*(xp1-x1));
                            end
                        else
                            correction_factor11 = (ones(nrec1,1)*bs1(1,:)+(ones(nrec1,1)*bs1(2,:)).*xp1+(ones(nrec1,1)*bs1(3,:)).*xp1.^2) ./ ...
                                (ones(nrec1,1)*bs1(1,:)+(ones(nrec1,1)*bs1(2,:)).*x1+(ones(nrec1,1)*bs1(3,:)).*x1.^2);
                        end
                        if (opt_subq.logX==1)
                            correction_factor11(xp1(:)==-Inf) = 0; %Zero daily mean in projection - set correction_factor = 0 (=> no variability in projection)
                            correction_factor11(x1(:)==-Inf) = 1; %Zero daily mean in hindcast - set correction_factor = 1 (=> no variability in projection)
                            correction_factor11(xp1(:)==-Inf & x1(:)==-Inf) = 1; %Zero daily mean in projection and hindcast - set correction_factor = 1 (=> no variability in projection)
                        end
                        correction_factor1(rec1,:) = max(0, correction_factor11);

                    elseif strcmp(opt_subq.model,'ksmooth')
                        yphat1 = NaN*ones(nrec1,ns); yphat2 = yphat1; yphat11 = NaN*ones(length(opt_subq.xp),1);
                        for l=1:ns
                            yphat11(:) = yphatm(k,:,l);
                            yphat1(:,l) = Xplin(Xpc(sel1(rec1),l),opt_subq.xp,1)*yphat11;
                            yphat2(:,l) = Xplin(Xhsub(selsub1(rec1),l),opt_subq.xp,1)*yphat11;
                        end
                        correction_factor1(rec1,:) = max(0, yphat1./yphat2);
                    end
                end

                if verbose>0; disp('Done correction factors for subdelta variance'); end
            else
                correction_factor1 = ones(nt1,ns);
            end
            correction_factor(sel1,:) = correction_factor1;
        end

        if verbose>0; disp(['Adding subdelta variability to projections for period ',int2str(yearpminv(j)),' through ',int2str(yearpmaxv(j))]); end

        if (loop_over_projection_times==0 && strcmp(delta_timescale,'daily')==1)
            if (fractional_subdelta==1)
                Xpcf(self1,:) = minXsub +  kron(((Xpc(sel1,:)-minXsub).*correction_factor1),ones(24,1)) .* dXhf1;
                if remove_deltascale_var==1;  Xpcf(self1,:) = Xpcf(self1,:) + (Xpci(self1,:)-kron(Xpc(sel1,:),ones(24,1))); end %Because dXhf1 in this case lacks the deltascale variation
                %Check: This should recover Xhf in case where Xpc=Xh:
                %   If remove_deltascale_var=0:
                %       Xpcf = minXsub + (Xh-minXsub)*(Xhfc1-minXsub)/(Xh-minXsub)
                %            = minXsub + (Xhf-minXsub)                               (since Xhfc1=Xhf in this case)
                %            = Xhf
                %
                %   If remove_deltascale_var=1:
                %       Xpcf = minXsub + (Xh-minXsub)*(Xhfc1-minXsub)/(Xh-minXsub) + (Xhi-Xh)
                %            = minXsub + (Xhf-Xhi+Xh-minXsub) + (Xhi-Xh)             (since Xhfc1=Xhf-Xhi+Xh in this case)
                %            = Xhf
            else
                Xpcf(self1,:) = Xpci(self1,:) + kron(correction_factor1,ones(24,1)) .* dXhf1;
                %Check: This should recover Xhf in case where Xpc=Xh:
                %   If remove_deltascale_var=0:
                %       Xpcf = Xh + (Xhfc1-Xh)                                       (since Xhi=Xh in this case) 
                %            = Xh + (Xhf-Xh)                                         (since Xhfc1=Xhf in this case)
                %            = Xhf
                %
                %   If remove_deltascale_var=1:
                %       Xpcf = Xhi + (Xhfc1-Xh)
                %            = Xhi + (Xhf-Xhi+Xh-Xh)                                 (since Xhfc1=Xhf-Xhi+Xh in this case)
                %            = Xhf
            end

        else

            for iblock=1:nblocks %Optional loop over spatial blocks e.g. to limit memory consumption
                if (subdelta_by_blocks==1)
% REDO                    selblock = find(latp>=latL_blocks1(iblock) & latp<latH_blocks1(iblock) & lonp>=lonL_blocks1(iblock) & lonp<lonH_blocks1(iblock));
                else
                    selblock = 1:ns;
                end
                nselblock = length(selblock);

                countf = 0;
                for i=1:nt1 %Loop over projection times (deltascale)
                    if (strcmp(delta_timescale,'daily')==1)
                        self2 = self1(floor(tdpf(self1))==floor(tdp(sel1(i)))); %Times to fill in finescale projection time series
                    elseif (strcmp(delta_timescale,'monthly')==1)
                        self2 = self1(yearpf(self1)==yearp(sel1(i)) & monthpf(self1)==monthp(sel1(i)));
                    elseif (strcmp(delta_timescale,'yearly')==1)
                        self2 = self1(yearpf(self1)==yearp(sel1(i)));
                    end
                    nself2 = length(self2);

                    if (match_subdelta_deltamean==1)
                        diff1 = sum((Xhsub(:,selblock) - Xpc(sel1(i),selblock)).^2,2);
                        selhsub(sel1(i),iblock) = find(diff1==min(diff1),1,'first'); %Choose the 'most-similar' day/month from hindcast wrt euclidean distance of daily/monthly averages
                        if (strcmp(delta_timescale,'daily')==1)
                            selhfsub1 = find(floor(tdhfsub)==floor(tdhsub(selhsub(sel1(i),iblock))));
                            %Time indices of selected subdelta variability from dXhfsub
                        elseif (strcmp(delta_timescale,'monthly')==1)
                            selhfsub1 = find(yearhfsub==yearhsub(selhsub(sel1(i),iblock)) & monthhfsub==monthhsub(selhsub(sel1(i),iblock)));
                        elseif (strcmp(delta_timescale,'yearly')==1)
                            selhfsub1 = find(yearhfsub==yearhsub(selhsub(sel1(i),iblock)));
                        end
                        dXhfsub1 = dXhfsub(selhfsub1,selblock);

                        if (correct_substd==1) %Correct the subdelta standard deviation for consistency with the deltascale mean
                            if opt_subq.bymonth==1; m = monthp(sel1(i)); else m = 1; end

                            if ~strcmp(opt_subq.model,'ksmooth'); bs1 = NaN*ones(opt_subq.nparb,nselblock); end
                            if (strcmp(opt_subq.model,'polynomial')==1)
                                bs1(:,:) = bs(m,:,selblock);
                                if opt_subq.logX==1; x1 = log(Xhsub(selhsub(sel1(i)),selblock)); xp1 = log(Xpc(sel1(i),selblock)); else x1 = Xhsub(selhsub(sel1(i)),selblock); xp1 = Xpc(sel1(i),selblock); end
                                if (opt_subq.logY==1)
                                    if (opt_subq.nparb==3)
                                        correction_factor1 = exp(bs1(2,selblock).*(xp1-x1) + bs1(3,selblock).*(xp1.^2-x1.^2));
                                    else
                                        correction_factor1 = exp(bs1(2,selblock).*(xp1-x1));
                                    end
                                else
                                    correction_factor1 = max(0, (bs1(1,selblock)+bs1(2,selblock).*xp1+bs1(3,selblock).*xp1.^2) ./ ...
                                        (bs1(1,selblock)+bs1(2,selblock).*x1+bs1(3,selblock).*x1.^2));
                                end
                                if (opt_subq.logX==1)
                                    correction_factor1(xp1==-Inf) = 0; %Zero daily mean in projection - set correction_factor = 0 (=> no variability in projection)
                                    correction_factor1(x1==-Inf) = 1; %Zero daily mean in hindcast - set correction_factor = 1 (=> no variability in projection)
                                    correction_factor1(xp1==-Inf & x1==-Inf) = 1; %Zero daily mean in projection and hindcast - set correction_factor = 1 (=> no variability in projection)
                                end    
                            elseif (strcmp(opt_subq.model,'ksmooth')==1)
                                yphat1 = NaN*ones(1,nselblock); yphat2 = yphat1; yphat11 = NaN*ones(length(opt_subq.xp),1);
                                for l=1:nselblock
                                    yphat11(:) = yphatm(m,:,selblock(l));
                                    yphat1(l) = Xplin(Xpc(sel1(i),selblock(l)),opt_subq.xp,1)*yphat11;
                                    yphat2(l) = Xplin(Xhsub(selhsub(sel1(i)),selblock(l)),opt_subq.xp,1)*yphat11;
                                end
                                correction_factor1(rec1,:) = max(0, yphat1./yphat2);
                            end

                        else
                            correction_factor1 = ones(1,nselblock);
                        end
                        correction_factor(sel1(i),selblock) = correction_factor1;
                        %NOTE: The correction factors derived with match_subdelta_deltamean=1 should be milder (closer to 1)
                        %      than those derived with match_subdelta_deltamean=0.
                        
                        rec1 = (1:24)';
                        self2last = min(self2)-1;
                        if (match_subdelta_hourly==1 && self2last>0) %Modify rec1 to minimize discontinuities
                            if (fractional_subdelta==1)
                                Xpcf1 = minXsub + (ones(nself2,1)*((Xpci(self2(1),selblock)-minXsub).*correction_factor1)).*dXhfsub1;
                                if remove_deltascale_var==1; Xpcf1 = Xpcf1 + (Xpci(self2,selblock)-ones(nself2,1)*Xpc(sel1(i),selblock)); end %Because dXhfsub1 in this case lacks the deltascale variation
                            else
                                Xpcf1 = ones(nself2,1)*Xpci(self2(1),selblock) + (ones(nself2,1)*correction_factor1).*dXhfsub1;
                            end

                            diff2 = sum((Xpcf1 - ones(nself2,1)*Xpcf(self2last,selblock)).^2,2);
                            selhsubhour(sel1(i),iblock) = find(diff2==min(diff2),1,'first'); %Choose the 'most-similar' hour that would minimize jump between days
                            if (selhsubhour(sel1(i),iblock)>1)
                                rec1 = [selhsubhour(sel1(i),iblock):24 23:-1:24-(selhsubhour(sel1(i),iblock)-1)]';
                                %Rearrange the hourly subdeltas to best match the end of the previous day, and 'reflect' indices off
                                %the end of the 24 hour interval in order to fill the 24 hours without introducing a jump.
                            end
                        end
                        if (fractional_subdelta==1)
                            Xpcf(self2,selblock) = minXsub + (ones(nself2,1)*((Xpc(sel1(i),selblock)-minXsub).*correction_factor1)).*dXhfsub1(rec1,:);
                            if remove_deltascale_var==1; Xpcf(self2,selblock) = Xpcf(self2,selblock) + (Xpci(self2,selblock)-ones(nself2,1)*Xpc(sel1(i),selblock)); end %Because dXhfsub1 in this case lacks the deltascale variation
                        else
                            Xpcf(self2,selblock) = Xpci(self2,selblock) + (ones(nself2,1)*correction_factor1).*dXhfsub1(rec1,:);
                        end

                    else

                        if (correct_subquantiles==1)
                            if (strcmp(opt_subq.model,'polynomial')==1)
                                if opt_subq.logX==1; x1 = log(Xhsub(selsub1(i),selblock)); xp1 = log(Xpc(sel1(i),selblock)); else x1 = Xhsub(selsub1(i),selblock); xp1 = Xpc(sel1(i),selblock); end
                            else
                                x1 = Xhsub(selsub1(i),selblock); xp1 = Xpc(sel1(i),selblock);
                            end
                            if (strcmp(opt_subq.resmodel,'polynomial')==1)
                                if opt_subq.logXres==1; x1res = log(Xhsub(selsub1(i),selblock)); xp1res = log(Xpc(sel1(i),selblock)); else x1res = Xhsub(selsub1(i),selblock); xp1res = Xpc(sel1(i),selblock); end
                            else
                                x1res = Xhsub(selsub1(i),selblock); xp1res = Xpc(sel1(i),selblock);
                            end

                            q1 = squeeze(qXhfsub(selsub1(i),:,selblock));
                            Isort1 = squeeze(IsortXhfsub(selsub1(i),:,selblock));
                            if nselblock==1; q1 = q1'; Isort1 = Isort1'; end
                            if opt_subq.logY==1; q1 = log(q1); end

                            if ~strcmp(opt_subq.model,'ksmooth'); bs1 = NaN*ones(opt_subq.nparb,nselblock); end
                            if (~isempty(opt_subq.resmodel) && ~strcmp(opt_subq.resmodel,'ksmooth')); bsres1 = NaN*ones(opt_subq.nparbres,nselblock); end
                            if opt_subq.bymonth==1; m = monthp(sel1(i)); else m = 1; end
                            for k=1:nqsub
                                if ~strcmp(opt_subq.model,'ksmooth'); bs1(:,:) = bs(m,k,:,selblock); end

                                if (strcmp(opt_subq.model,'polynomial')==1)
                                    if (opt_subq.nparb==4)
                                        yphat1 = bs1(1,:) + bs1(2,:).*xp1 + bs1(3,:).*xp1.^2 + bs1(4,:).*xp1.^3;
                                        yhat1 = bs1(1,:) + bs1(2,:).*x1 + bs1(3,:).*x1.^2 + bs1(4,:).*x1.^3;
                                    elseif (opt_subq.nparb==3)
                                        yphat1 = bs1(1,:) + bs1(2,:).*xp1 + bs1(3,:).*xp1.^2;
                                        yhat1 = bs1(1,:) + bs1(2,:).*x1 + bs1(3,:).*x1.^2;
                                    else
                                        yphat1 = bs1(1,:) + bs1(2,:).*xp1;
                                        yhat1 = bs1(1,:) + bs1(2,:).*x1;
                                    end
                                elseif (strcmp(opt_subq.model,'tanh')==1)
                                    if (opt_subq.nparb==3)
                                        yphat1 = bs1(1,:)+bs1(2,:).*tanh(bs1(3,:).*xp1);
                                        yhat1 = bs1(1,:)+bs1(2,:).*tanh(bs1(3,:).*x1);
                                    else
                                        yphat1 = bs1(1,:).*tanh(bs1(2,:).*xp1);
                                        yhat1 = bs1(1,:).*tanh(bs1(2,:).*x1);
                                    end
                                elseif (strcmp(opt_subq.model,'powerlaw')==1)
                                    if (opt_subq.nparb==2)
                                        yphat1 = bs1(1,:).*xp1.^bs1(2,:);
                                        yhat1 = bs1(1,:).*x1.^bs1(2,:);
                                    elseif (opt_subq.nparb==3)
                                        %yphat1 = max(0,bs1(1,:)*xp1.^bs1(2,:)+bs1(3,:));
                                        %yhat1 = max(0,bs1(1,:)*x1.^bs1(2,:)+bs1(3,:));
                                        stop %WIP!!!
                                    end
                                elseif strcmp(opt_subq.model,'ksmooth')
                                    yphat1 = NaN*ones(1,nselblock); yhat1 = yphat1; yphat11 = NaN*ones(length(opt_subq.xp),1);
                                    for l=1:nselblock
                                        yphat11(:) = yphatm(m,k,:,selblock(l));
                                        yphat1(l) = Xplin(xp1(l),opt_subq.xp,1)*yphat11;
                                        yhat1(l) = Xplin(x1(l),opt_subq.xp,1)*yphat11;
                                    end
                                end

                                if (~isempty(opt_subq.resmodel))
                                    if ~strcmp(opt_subq.resmodel,'ksmooth'); bsres1(:,:) = bsres(m,k,:,selblock); end
                                    if strcmp(opt_subq.resmodel,'polynomial')
                                        if (opt_subq.nparbres==3)
                                            ypreshat1 = bsres1(1,:) + bsres1(2,:).*xp1res + bsres1(3,:).*xp1res.^2;
                                            yreshat1 = bsres1(1,:) + bsres1(2,:).*x1res + bsres1(3,:).*x1res.^2;
                                        else
                                            ypreshat1 = bsres1(1,:) + bsres1(2,:).*xp1res;
                                            yreshat1 = bsres1(1,:) + bsres1(2,:).*x1res;
                                        end
                                        if opt_subq.logYres==1; ypreshat1 = exp(ypreshat1); yreshat1 = exp(yreshat1); end
                                        fac1 = ypreshat1./yreshat1;
                                    elseif (strcmp(opt_subq.resmodel,'gaussian')==1)
                                        ypreshat1 = bsres1(1,:).*exp(-(0.5*bsres1(2,:).^2).*xp1res.^2);
                                        yreshat1 = bsres1(1,:).*exp(-(0.5*bsres1(2,:).^2).*x1res.^2);
                                        fac1 = ypreshat1./yreshat1;
                                    elseif strcmp(opt_subq.resmodel,'powerlaw')
                                        fac1 = (xp1res./x1res) .^ bsres1(2,:);
                                    elseif strcmp(opt_subq.resmodel,'ksmooth')
                                        ypm1 = NaN*ones(2,nselblock);
                                        for l=1:nselblock
                                            ypm1(:,l) = Xplin([xp1res(l);x1res(l)],opt_subq.xpres,1)*squeeze(ypreshatm(m,k,:,selblock(l)));
                                        end
                                        fac1 = ypm1(1,:)./ypm1(2,:);
                                    else
                                        error(['Inflation factor fac1 not yet coded for opt_subq.resmodel = ',opt_subq.resmodel])
                                    end
                                else
                                    fac1 = ones(1,nselblock);
                                end

                                if (fractional_subdelta==1)
                                    q1(k,:) = minXsub + (yphat1-minXsub).*fac1.*(q1(k,:)-minXsub)./(yhat1-minXsub);
                                else
                                    q1(k,:) =  yphat1 + fac1.*(q1(k,:)-yhat1);
                                end
                            end
                            if opt_subq.logY==1; q1 = exp(q1); end

                            for k=1:nselblock
                                if (remove_deltascale_var==1)
                                    Xpcf(self2(Isort1(:,k)),selblock(k)) = Xpci(self2(Isort1(:,k)),selblock(k)) + q1(:,k);
                                else
                                    Xpcf(self2(Isort1(:,k)),selblock(k)) = q1(:,k); %Quantile-corrected finescale projections
                                end
                            end

                        else %Simply add the hindcast variability with possible corrections for subdelta standard deviation
                            if (fractional_subdelta==1)
                                Xpcf(self2,:) = minXsub + (ones(nself2,1)*((Xpc(sel1(i),:)-minXsub).*correction_factor(sel1(i),:))).*dXhf1(countf+1:countf+nself2,:);
                                if remove_deltascale_var==1; Xpcf(self2,:) = Xpcf(self2,:) + (Xpci(self2,:)-ones(nself2,1)*Xpc(sel1(i),:)); end %Because dXhf1 in this case lacks the deltascale variation
                            else
                                Xpcf(self2,:) = Xpci(self2,:) + (ones(nself2,1)*correction_factor(sel1(i),:)).*dXhf1(countf+1:countf+nself2,:);
                            end
                        end
                    end
                    countf = countf + nself2;
                end
            end
        end
        if verbose>0; disp('Done adding subdelta variability to projections'); end
    end
    if dtdhf1==0; Xpcf(end,:) = Xpcf(end-1,:); end
    %If hourly data are 00:00, 01:00, ..., 23:00 then the tdpf set above will result in an extra
    %time at 00:00 on ceil(max(tdp))+1 --- this can be filled by repeating the last datum.
    if correct_substd==1; out.correction_factor = correction_factor; end

    
    %Impose lower/upper limits if provided
    if (~isempty(XpcfL)||~isempty(XpcfH))
        if verbose>0; disp('Imposing lower/upper limits on projections'); end
        Xpcfo = Xpcf; %Can be a useful diagnostic to check that capping is not too severe.
        if ~isempty(XpcfL); Xpcf = max(XpcfL,Xpcf); end
        if ~isempty(XpcfH); Xpcf = min(XpcfH,Xpcf); end
        if ~isequaln(Xpcfo,Xpcf); out.Xpcfo = Xpcfo; end
        clear Xpcfo
        %Note we do not subselect for not-NaN (too slow), so this will set any NaN values to XpcL
    end
    
    if (use_XpcfH_SWI>0)
        optSWIf = optSWI; optSWIf.calc_dav = 0; optSWIf.year = yearpf;
        if use_SWI_hav==1; optSWIf.calc_hav = 1; end
        if verbose>0; disp('Calculating maximum hourly values of shortwave irradiance'); end
        out.XpcfH_SWI = NaN*ones(length(tdpf),ns);
        for i=1:ns %We have to loop over the time series here, else we run out of memory.
            if (use_SWI_hav==1)
                [~,outSWI] = calc_surfaceSWI(latp(i),lonp(i),yrdaypf,optSWIf);
                out.XpcfH_SWI(:,i) = outSWI.Q0hav;
            else
                out.XpcfH_SWI(:,i) = calc_surfaceSWI(latp(i),lonp(i),yrdaypf,optSWIf);
            end
            if (mod(i,100)==0 && verbose>0); disp(['Done ',int2str(i),' out of ',int2str(ns),' time series']); end
        end
        if verbose>0; disp('Done maximum hourly values of shortwave irradiance'); end
        
        if (use_XpcfH_SWI==1)
            if verbose>0; disp('Capping all hourly projections of shortwave irradiance with instantaneous clear-sky maximum values'); end
            Xpcf = min(out.XpcfH_SWI,Xpcf); %Cap the corrected hourly projections with instantaneous clear-sky maximum values
            %Warning: This is can lead to widespread correction of values near sunrise/sunset due to use of hourly-average irradiance in ERA5
            
        elseif (use_XpcfH_SWI==2)
            if verbose>0; disp('Calculating daily maximum values of projected shortwave irradiance'); end
            Xpcf_dmax = NaN*Xpcf;
            tdu = unique(floor(tdpf));
            for i=1:length(tdu)
                self1 = find(floor(tdpf)==tdu(i));
                Xpcf_dmax(self1,:) = ones(length(self1),1)*max(Xpcf(self1,:));
            end   
            if verbose>0; disp('Done daily maximum values of projected shortwave irradiance'); end
            if verbose>0; disp('Capping selected hourly projections of shortwave irradiance with instantaneous clear-sky maximum values'); end
% REDO            rec = find(Xpcf(:)>=frcrit_capSWI*Xpcf_dmax(:)); %Note this subselection can be slow
            Xpcf(rec) = min(out.XpcfH_SWI(rec),Xpcf(rec));
        end
    end
    if verbose>0; disp('Done imposing lower/upper limits on finescale projections (if any)'); end
    if verbose>0; disp('Done finescale projections'); toc; end
    
    if (recalculate_Xpc==1)
        out.Xpco2 = Xpc;
        if (strcmp(delta_timescale,'daily')==1) %Hourly variability about daily delta scale
            if verbose>0; disp('Recalculating daily averages from sub-daily projections'); end
            for i=1:ntp; Xpc(i,:) = mean(Xpcf(floor(tdpf)==floor(tdp(i)),:)); end
            if verbose>0; disp('Done recalculated daily averages'); end
        elseif (strcmp(delta_timescale,'monthly')==1) %Hourly/daily variability about monthly delta scale
            if verbose>0; disp('Recalculating monthly averages from sub-monthly projections'); end
            for i=1:ntp; Xpc(i,:) = mean(Xpcf(yearpf==yearp(i) & monthpf==monthp(i),:)); end
            if verbose>0; disp('Done recalculated monthly averages'); end
        elseif (strcmp(delta_timescale,'yearly')==1) %Hourly/daily/monthly variability about yearly delta scale
            if verbose>0; disp('Recalculating yearly averages from sub-monthly projections'); end
            for i=1:ntp; Xpc(i,:) = mean(Xpcf(yearpf==yearp(i),:)); end
            if verbose>0; disp('Done recalculated yearly averages'); end
        end
    end
    
    out.Xpcf = Xpcf; out.tdpf = tdpf;
    if match_subdelta_deltamean==1; out.selhsub = selhsub; end
    if match_subdelta_hourly==1; out.selhsubhour = selhsubhour; end
end



