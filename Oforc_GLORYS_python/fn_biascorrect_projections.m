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

[ntp,ns] = size(Xp);
out = [];
if nargin<8;
    opt = [];
end
if isfield(opt,'verbose')==1;
    verbose = opt.verbose;
else
    verbose = 1;
end
if isfield(opt,'seasonal')==1;
    seasonal = opt.seasonal;
else seasonal = 1;
end
if isfield(opt,'use_month')==1;
    use_month = opt.use_month;
else
    use_month = 1;
end
%1 to use month-of-year (1-12) rather than yrday (day-of-year 0-365) to compute climatological statistics and deltas.
%%%if isfield(opt,'dyrday')==1; dyrday = opt.dyrday; else dyrday = 7; end
if isfield(opt,'tclimc')==1;
    tclimc = opt.tclimc;
else
    tclimc = 1:12;
end
%Climatological time centres used to calculate climatological (means,deltas) and quantiles (def = 1:12 for monthly delta change with use_month=1).
if isfield(opt,'ndays_year')==1;
    ndays_year = opt.ndays_year;
else
    ndays_year = 365.25;
end
%Best approach here seems to be to use ndays_year = 365.25; if we use ndays_year = 365 then we
%will combine data from Dec 31st and Jan 1st when the preceding year is a leap year.
if isfield(opt,'tol_dtclim')==1;
    tol_dtclim = opt.tol_dtclim;
else
    tol_dtclim = 0;
end
%Tolerance level for climatological time used in computing climatological statistics (def = 0 for simple monthly delta change).
% Using tol_dtclim>0 may allow a more robust definition of climatology, e.g. if the no. of reference years (yearmaxref-yearminref+1) is limited.
% E.g. if aiming for a daily climatology (tclimc=(0.5:365.5)), it is probably a good idea to allow e.g. tol_dtclim = 5 (i.e. +/- 5 days tolerance)
% in order to allow robust climatological statistics and consequent deltas.
if isfield(opt,'tol_dtclimp')==1;
    tol_dtclimp = opt.tol_dtclimp;
else
    tol_dtclimp = 0;
end
%Tolerance level for climatological time used in correcting projections (def = 0 for simple monthly delta change).
% This may be different to tol_dtclim because tol_dtclim may be expanded to allow more robust climatological statistics.
if isfield(opt,'fractional')==1;
    fractional = opt.fractional;
else
    fractional = 0;
end
if isfield(opt,'minX')==1;
    minX = opt.minX;
else
    minX = [];
end
%Reference minimum value to avoid potential division by ~0 in fractional methods.
if isfield(opt,'allow_block_lag')==1;
    allow_block_lag = opt.allow_block_lag;
else
    allow_block_lag = 1;
end

if (size(Xh,2)>0 && ns==1)
    %In case of multiple hindcast series and one projection series, apply the single projection series to all hindcast climatologies
     size(Xh,2);
    Xp = Xp*ones(1,ns);
end

if (method==2)
    if isfield(opt,'qsel')==1; qsel = opt.qsel; else qsel = 0.05:0.1:0.95; end
    qsel = qsel(:);
end

if isfield(opt,'XpcL')==1;
    XpcL = opt.XpcL;
else
    XpcL = [];
end %Optional lower limits for Xpc
if isfield(opt,'XpcH')==1;
    XpcH = opt.XpcH;
else
    XpcH = [];
end %Optional upper limits for Xpc
if isfield(opt,'use_XpcH_SWI')==1;
    use_XpcH_SWI = opt.use_XpcH_SWI;
else
    use_XpcH_SWI = 0;
end
%If opt.use_XpcH_SWI=1, maximum daily-average values of shortwave irradiance are calculated using calc_surfaceSWI.m
%and these are used to impose an upper limit on Xpc
if (use_XpcH_SWI==1)
    if isfield(opt,'latp')==1;
        latp = opt.latp;
    else
        latp = [];
    end %Latitudes for each time series (degN)
    if isfield(opt,'lonp')==1;
        lonp = opt.lonp;
    else
        lonp = [];
    end %Longitudes for each time series (degE)

    if isfield(opt,'optSWI')==1;
        optSWI = opt.optSWI;
    else
        optSWI = [];
    end %Options for maximum (clear-sky) SWI calculation.
    if (isempty(optSWI))
        optSWI = struct('f_model',1,'decl_model',1,'use_eqtime',1,'Q_model',0,'cloud_model',0,'ndays_year',365.2422,'use_tday_lag',1);
        %Default options for calc_surfaceSWI.m: These have been found to give best performance (minimum overshoot)
        %in tests using hourly ERA5 data for ROHO800 model (western Norway).
    end
    if isfield(opt,'use_SWI_hav')==1;
        use_SWI_hav = opt.use_SWI_hav;
    else
        use_SWI_hav = 1;
    end
    %1 (def) to use hourly averages rather than instantaneous values in calculation of hourly data (consistent with ERA5).
end

if isfield(opt,'yearpminv')==1;
    yearpminv = opt.yearpminv;
else
    yearpminv = [];
end
if isfield(opt,'yearpmaxv')==1;
    yearpmaxv = opt.yearpmaxv;
else
    yearpmaxv = [];
end
%Optional input of set of prediction year periods over which the subdelta variability will be resampled
%and/or deltascale variability will be repeated (if method=2).

if isfield(opt,'legacy')==1;
    legacy = opt.legacy;
else
    legacy = 0;
end
if (legacy==1)
    allow_block_lag = 0;
end

if isfield(opt,'use_subdelta')==1;
    use_subdelta = opt.use_subdelta;
else
    use_subdelta = 0;
end
%1 to use subdelta variability from input finescale hindcast (opt.Xhf,opt,tdhf) to correct the projections
if (use_subdelta==1)
    if isfield(opt,'tdpfmin')==1;
        tdpfmin = opt.tdpfmin;
    else
        tdpfmin = floor(min(tdp));
    end
    if isfield(opt,'tdpfmax')==1;
        tdpfmax = opt.tdpfmax;
    else
        tdpfmax = ceil(max(tdp));
    end
    if isfield(opt,'Xhf')==1;
        Xhf = opt.Xhf;
    else
        Xhf = [];
    end
    if isfield(opt,'tdhf')==1;
        tdhf = opt.tdhf;
    else
        tdhf = [];
    end
    if isfield(opt,'yearminsub')==1;
        yearminsub = opt.yearminsub;
    else
        yearminsub = yearminref;
    end
    if isfield(opt,'yearmaxsub')==1;
        yearmaxsub = opt.yearmaxsub;
    else
        yearmaxsub = yearmaxref;
    end
    %Optional input of minimum/maximum years (yearminsub,yearmaxsub) to sample subdelta variability from
    %-- by default set equal to (yearminsub,yearmaxsub).
    if isfield(opt,'fractional_subdelta')==1;
        fractional_subdelta = opt.fractional_subdelta;
    else
        fractional_subdelta = 0;
    end
    %fractional_subdelta = 1 to treat the subdelta variability as a fractional perturbation (better for non-negative variables)
    if isfield(opt,'minXsub')==1;
        minXsub = opt.minXsub;
    else minXsub = 0;
    end
    %Reference minimum value to avoid potential division by ~0 in fractional subdelta methods.
    if isfield(opt,'loop_over_projection_times')==1;
        loop_over_projection_times = opt.loop_over_projection_times;
    else
        loop_over_projection_times = 0;
    end
    %Using matlab's kron.m below gives identical results to the loop calculation, but it ~800x faster
    if isfield(opt,'remove_deltascale_var')==1;
        remove_deltascale_var = opt.remove_deltascale_var;
    else
        remove_deltascale_var = 1;
    end

    if isfield(opt,'XpcfL')==1;
        XpcfL = opt.XpcfL;
    else
        XpcfL = [];
    end %Optional lower limits for Xpcf
    if isfield(opt,'XpcfH')==1;
        XpcfH = opt.XpcfH;
    else
        XpcfH = [];
    end %Optional upper limits for Xpcf
    if isfield(opt,'use_XpcfH_SWI')==1;
        use_XpcfH_SWI = opt.use_XpcfH_SWI;
    else
        use_XpcfH_SWI = 0;
    end
    %If opt.use_XpcfH_SWI>0, maximum hourly values of shortwave irradiance are calculated using calc_surfaceSWI.m
    %and these are used to impose an upper limit on Xpcf. Then if:
    %opt.use_XpcfH_SWI = 1: All hourly projections are limited by the maximum (clear-sky) values.
    %opt.use_XpcfH_SWI = 2: Only hourly projections >=opt.frcrit_capSWI*daily maximum are limited by the clear-sky values.

    if (use_XpcfH_SWI==2)
        if isfield(opt,'frcrit_capSWI')==1;
            frcrit_capSWI = opt.frcrit_capSWI;
        else
            frcrit_capSWI = 0.25;
        end
    end
    if isfield(opt,'recalculate_Xpc')==1;
        recalculate_Xpc = opt.recalculate_Xpc;
    else
        recalculate_Xpc = 1;
    end
    %opt.recalculate_Xpc = 1: Recalculate the daily average projections using the hourly projections after limits imposed.

    if isfield(opt,'subdelta_by_blocks')==1;
        subdelta_by_blocks = opt.subdelta_by_blocks;
    else
        subdelta_by_blocks = 0;
    end
    %1 to divide subdelta variability calculation into latitude-longitude blocks, to obtain better matches to hindcast daily means (match_subdelta_deltamean=1).
    %This is probably not a good idea unless the target variable has very limited spatial correlation between grid cells,
    %since the blockwise calculation will not preserve spatial correlations across block boundaries.
    if (subdelta_by_blocks==1)
        if isfield(opt,'latp')==1;
            latp = opt.latp;
        else
            latp = [];
        end %Latitudes for each time series (degN)
        if isfield(opt,'lonp')==1;
            lonp = opt.lonp;
        else
            lonp = [];
        end %Longitudes for each time series (degE)

        if isfield(opt,'latL_blocks')==1;
            latL_blocks = opt.latL_blocks;
        else
            latL_blocks = [];
        end %Lower latitude boundaries of blocks to divide calculation by (one row).
        if isfield(opt,'latH_blocks')==1;
            latH_blocks = opt.latH_blocks;
        else
            latH_blocks = [];
        end %Upper latitude boundaries of blocks to divide calculation by (one row).
        if isfield(opt,'lonL_blocks')==1;
            lonL_blocks = opt.lonL_blocks;
        else
            lonL_blocks = [];
        end %Lower longitude boundaries of blocks to divide calculation by (one row).
        if isfield(opt,'lonH_blocks')==1;
            lonH_blocks = opt.lonH_blocks;
        else
            lonH_blocks = [];
        end %Upper longitude boundaries of blocks to divide calculation by (one row).

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
if seasonal==1;
    ntclimc = length(tclimc);
else
    ntclimc = 1;
end

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

if method==2;
    methodstr = 'repeated hindcast correction';
end
if verbose>0;
    disp(['Calculating bias-corrected projections using ',methodstr,' method']);
    tic;
end

if (method==2)
    %Calculate projection model trends in quantiles between repeated hindcast periods
    nperiods = length(yearpminv);
    nq = length(qsel);
    qXp = NaN*ones(ntclimc, nperiods, nq, ns);
    for i=1:ntclimc %Loop over climatological time divisions (e.g. month-of-year)
        if (seasonal==1)
            % Calculate temporal separations in a periodic sense using fnperiodic_distm.m
            if use_month==1;
                dist1 = fnperiodic_distm(monthp,tclimc(i),12);
            else
                dist1 = fnperiodic_distm(yrdayp,tclimc(i),ndays_year);
            end
            selt0 = find(dist1<=tol_dtclim); % Choose all data within tolerance tol_dtclim of climatological time centre tclimc(i)
        else
            selt0 = 1:ntp;
        end

        % Loop over repeated hindcast periods
        for j=1:nperiods
            selt = selt0(yearp(selt0)>=yearpminv(j) & yearp(selt0)<=yearpmaxv(j));
            qXp(i,j,1:nq,1:ns) = quantile(Xp(selt,:),qsel);
        end
    end
    out.qXp = qXp;
end

%%Apply the delta changes / quantile corrections to generate bias-corrected projections.
[yearh,monthh,~] = datevec(tdh);
if use_month==0;
    yrdayh = tdh - datenum(yearh,1,1);
end

seltclimh0 = find(yearh>=yearminref & yearh<=yearmaxref);
nthref = length(seltclimh0);

if all(diff(tdh)==1)
    delta_timescale = 'daily';
elseif all(diff(monthh)==1|diff(monthh)==-11)
    delta_timescale = 'monthly';
end

%Preallocation for hindcast statistic arrays
nseltclimh = NaN*ones(ntclimc,1);
Xpc = NaN*ones(ntp,ns);
Fhref = NaN*ones(nthref,ns);

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
        seltclimh = seltclimh0(dist1<=tol_dtclim); %Choose all data within tolerance tol_dtclim of climatological time centre tclimc(i)
        seltp = find(dist1p<=tol_dtclimp); %Choose all data within tolerance tol_dtclimp of climatological time centre tclimc(i)
    else
        seltclimh = seltclimh0;
        seltp = 1:ntp;
    end
    nseltclimh(i) = length(seltclimh);

    %Here we need to record the cdf value F for each time point in the repeated hindcast series
    if seasonal==1;
        selthref1 = find(dist1<=tol_dtclimp);
    else
        selthref1 = 1:nthref;
    end
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
out.nseltclimh = nseltclimh;
%if method==0; out.Xhclim = Xhclim; end
%if method==1; out.qXhref = qXhref; end
if method==2;
    out.Fhref = Fhref;
end


%%Next we apply the correction factors to the projections (if not applied already for simple delta change)

    iref = find(yearpminv==yearminref & yearpmaxv==yearmaxref);
    dtdhref = tdh(seltclimh0) - datenum(yearminref,1,1); %Days in hindcast since start of reference period

    dqXp = NaN*qXp;
    for i=1:nperiods %Loop over repeated time blocks
        sel1 = find(yearp>=yearpminv(i) & yearp<=yearpmaxv(i));
        nt1 = length(sel1);

        if (allow_block_lag==1)
            dtdp1 = tdp(sel1(1)) - datenum(yearpminv(i),1,1); %Days since start of projection block i
            ihref1 = find(abs(dtdhref-dtdp1)<=tol_dtclimp); %Starting index within repeated hindcast series
            %For example if the projection time series only runs from Jan 2nd of the projection block,
            %we should start the repeated hindcast block also from Jan 2nd.
        else
            ihref1 = 1;
        end
        % Indices within Xh, tdh etc. of the repeated hindcast data for this block
        selth1 = seltclimh0(ihref1:min(nt1+(ihref1-1),nthref));
        % Quantiles (cdf values) of the repeated hindcast data for this block
        Fhref1 = Fhref(ihref1:min(nt1+(ihref1-1),nthref),:);

        if (i==iref && 0==1)
            Xpc(sel1,:) = Xh(seltclimh0,:);
        else

            for j=1:ntclimc
                if (seasonal==1)
                    if use_month==1;
                        dist1 = fnperiodic_distm(monthh(selth1),tclimc(j),12);
                    else
                        dist1 = fnperiodic_distm(yrdayh(selth1),tclimc(j),ndays_year);
                    end
                    if (legacy==1)
                        selthref1 = find(dist1<=tol_dtclim);
                    else
                        selthref1 = find(dist1<=tol_dtclimp); %Choose all data within tolerance tol_dtclimp of climatological time centre tclimc(i)
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
                if verbose>0;
                    disp(['Repeating last day of hindcast reference period to apply to years ',int2str(yearpminv(i)),' to ',int2str(yearpmaxv(i))]);
                end
                Xpc(sel1(end),:) = Xpc(sel1(end-1),:);
            elseif nt1>(nthref+1)
                error(['Not yet coded re-use of hindcast data for nt1 = ',int2str(nt1),' and nthref = ',int2str(nthref)])
            end
        end
        if verbose>0;
            disp(['Done repeated hindcast corrected projections for years ',int2str(yearpminv(i)),' through ',int2str(yearpmaxv(i))]);
            tic;
        end
    end
    out.dqXp = dqXp;

%Impose lower/upper limits if provided
if (~isempty(XpcL)||~isempty(XpcH))
    Xpco = Xpc; %Can be a useful diagnostic to check that capping is not too severe.
    if ~isempty(XpcL);
        Xpc = max(XpcL,Xpc);
    end
    if ~isempty(XpcH);
        Xpc = min(XpcH,Xpc);
    end
    if ~isequaln(Xpco,Xpc);
        out.Xpco = Xpco;
    end
    clear Xpco
end

if (use_XpcH_SWI==1)
    if verbose>0;
        disp('Calculating maximum daily-average values of shortwave irradiance');
    end
    optSWI.calc_dav = 1;
    optSWI.year = yearp;
    [~,outSWI] = calc_surfaceSWI(latp,lonp,yrdayp,optSWI);
    if verbose>0;
        disp('Done maximum daily-average values of shortwave irradiance');
    end
    out.XpcH_SWI = outSWI.Q0dav;
    Xpc = min(out.XpcH_SWI,Xpc); %Cap the corrected daily-average projections with clear-sky maximum values
end

if verbose>0;
    disp(['Done bias-corrected projections using ',methodstr,' method']);
    toc;
end

%%%If required, add 'subdelta' variability at the extra temporal resolution of the hindcast (e.g. hourly)
if (use_subdelta==1 && ~isempty(Xhf) && ~isempty(tdhf))
    if verbose>0;
        disp('Calculating finescale projections with subdelta variability');
        tic;
    end

    if all(diff(tdh)==1)
        delta_timescale = 'daily';
        nqsub = 24;
    elseif all(diff(monthh)==1|diff(monthh)==-11)
        delta_timescale = 'monthly';
        nqsub = 28;
    else
        error('Subdelta variability quantiles not yet defined if hindcast data opt.Xh neither daily nor monthly')
    end

    seltsub = find(yearh>=yearminsub & yearh<=yearmaxsub);
    Xhsub = Xh(seltsub,:); tdhsub = tdh(seltsub); %Hindcast output during subdelta period
    [yearhsub,monthhsub,dayhsub] = datevec(tdhsub);
    ntsub = length(tdhsub);

    [yearhf,monthhf,dayhf] = datevec(tdhf);
    seltfsub = find(yearhf>=yearminsub & yearhf<=yearmaxsub);
    Xhfsub = Xhf(seltfsub,:); tdhfsub = tdhf(seltfsub); %Finescale hindcast output during subdelta period
    [yearhfsub,monthhfsub,~] = datevec(tdhfsub);

    if (remove_deltascale_var==1)
        Xhi = interp1(tdh,Xh,tdhf,'linear','extrap'); %Deltascale variation interpolated to fine scale
        dXhf = Xhf - Xhi;
        clear Xhi
        %Note: Here we allow linear extrapolation, consistent with projection baseline variability (Xpci) calculated below.
    end
    dXhfsub = NaN*Xhfsub;

    m = 0;
    for i=1:ntsub
        % Hourly variability about daily delta scale
        if (strcmp(delta_timescale,'daily')==1)
            sel1 = find(yearhf==yearhsub(i) & monthhf==monthhsub(i) & dayhf==dayhsub(i));
        % Hourly/daily variability about monthly delta scale
        elseif (strcmp(delta_timescale,'monthly')==1)
            sel1 = find(yearhf==yearhsub(i) & monthhf==monthhsub(i));
        end
        sel2 = m+1:m+length(sel1);

        if (remove_deltascale_var==1)
            Xhfc1 = dXhf(sel1,:) + Xhsub(i,:);
            %This is the finescale variation with the signal from deltascale variation (day-to-day, or month-to-month) removed.
        else
            Xhfc1 = Xhf(sel1,:);
        end

        if (fractional_subdelta==1)
            dXhfsub(sel2,:) = (Xhfc1-minXsub) ./ (ones(length(sel1),1)*Xhsub(i,:)-minXsub);
            %Check: mean(dXhfsub(sel2,:)) %This should be = 1 for all time series IFF remove_deltascale_var=0 and the hindcast daily/monthly averages were calculated by simple averages of all hourly/daily data within each day/month.
        else
            dXhfsub(sel2,:) = Xhfc1 - Xhsub(i,:);
            %Check: mean(dXhfsub(sel2,:)) %This should be = 0 for all time series IFF remove_deltascale_var=0 and the hindcast daily/monthly averages were calculated by simple averages of all hourly/daily data within each day/month.
        end

        m = m + length(sel1);
    end

    %Next calculate the baseline variability and prepare for main loop over projection periods
    dtdhf1 = min(tdhf-floor(tdhf)); %Initial increment after 00:00
    difftdhf = tdhf(2)-tdhf(1); %Finescale time increment (in days, assumed constant) --- This is unsafe for hourly finescale due to potentially-limited precision
    if strcmp(delta_timescale,'daily')==1;
        difftdhf = 1/24;
    end
    tdpf = (tdpfmin+dtdhf1:difftdhf:tdpfmax+1)';
    ntpf = length(tdpf);
    %Note: We set the hourly projection timestamps for consistency with the hourly hindcast timestamps.

    [yearpf,monthpf,~] = datevec(tdpf);
    yrdaypf = tdpf - datenum(yearpf,1,1);
    Xpcf = NaN*ones(ntpf,ns); %Preallocate final 'finescale' projection time series, with subdelta variability added

    if verbose>0;
        disp('Calculating baseline variability to which subdelta variation will be added');
    end
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
            end
            Xpci(self1,:) = ones(length(self1),1)*Xpc(i,:);
        end
    end
    if verbose>0;
        disp('Done baseline variability');
    end

    if (isempty(yearpminv) || isempty(yearpmaxv))
        nyrssub = yearmaxsub-yearminsub+1;
        yearpminv = min(yearp):nyrssub:max(yearp);
        yearpmaxv = min(yearp)+nyrssub-1:nyrssub:max(yearp);

        %If not matching the subdelta variability by most-similar days, then we need to define
        %a set of multiannual year periods for which the subdelta variability will be repeated over
        %(usually most convenient to divide the prediciton period into 10 or 20-year periods).
        %If not defined by input these periods are defined here.
        %Although it should not be necessary to split the full prediction in match_subdelta_deltamean=1,
        %tests show that the calculation is almost 2x faster if it is split up (with identical final results).
    end
    selp = find(yearpmaxv<=max(yearpf)); yearpminv = yearpminv(selp); yearpmaxv = yearpmaxv(selp);
    out.yearpminv = yearpminv; out.yearpmaxv = yearpmaxv;
    correction_factor = ones(ntp,ns);

    dtdhsub = tdh(seltsub) - datenum(yearminsub,1,1); %Days in hindcast since start of subdelta reference period

    %Main loop over projection periods for which finescale projections are required
    for j=1:length(yearpminv)
        sel1 = find(yearp>=yearpminv(j) & yearp<=yearpmaxv(j));
        self1 = find(yearpf>=yearpminv(j) & yearpf<=yearpmaxv(j));
        % No. of projection times in current time chunk
        nt1 = length(sel1);
        % No. of finescale projection times in current time chunk
        ntf1 = length(self1);

        if (allow_block_lag==1)
            %Identify starting index for the repeated hindcast variability

            % Days since start of projection block j
            dtdp1 = tdp(sel1(1)) - datenum(yearpminv(j),1,1);

            % Starting index within repeated hindcast series (deltascale)
            ihsub1 = find(abs(dtdhsub-dtdp1)<=tol_dtclimp);
        else
            ihsub1 = 1;
        end

        %Identify repeated hindcast series to use for this block
        selsub1 = ihsub1:min(ihsub1-1+nt1,ntsub);
        %Note: This is not necesssary if correct_subquantiles=1 because in this case
        %      only the coarsescale indices selsub1 are needed
        if (strcmp(delta_timescale,'daily')==1)
            ihfsub1 = find(floor(tdhfsub)==floor(tdhsub(selsub1(1))),1,'first');
            ihfsub2 = find(floor(tdhfsub)==floor(tdhsub(selsub1(end))),1,'last');
        elseif (strcmp(delta_timescale,'monthly')==1)
            ihfsub1 = find(yearhfsub==yearhsub(selsub1(1)) & monthhfsub==monthhsub(selsub1(1)),1,'first');
            ihfsub2 = find(yearhfsub==yearhsub(selsub1(end)) & monthhfsub==monthhsub(selsub1(end)),1,'last');
        end
        dXhf1 = dXhfsub(ihfsub1:ihfsub2,:);

        if (strcmp(delta_timescale,'monthly')==1 && length(dXhf1(:,1))==(ntf1-1))
            % Repeat last day if necessary due to varying no. leap years per block.
            dXhf1 = [dXhf1; dXhf1(end,:)];

            if verbose>0;
                disp(['Repeating last day of subdelta reference period to apply to years ',int2str(yearpminv(j)),' to ',int2str(yearpmaxv(j))]);
            end
        end

        %Add extra day if necessary to account for possible different no. leap years per decade (e.g. only 2 in subdelta period 2010-2019)
        if (strcmp(delta_timescale,'daily')==1 && (ihsub1-1+nt1>ntsub))
            if (ihsub1-1+nt1==ntsub+1)
                if verbose>0;
                    disp(['Repeating last day of subdelta reference period to apply to years ',int2str(yearpminv(j)),' to ',int2str(yearpmaxv(j))]);
                end
                selsub1 = [selsub1 ntsub]; %#ok<AGROW> %Repeat last day of subdelta period
                %The use of 'selsub1' avoids having to append several variables ('Xhsub' etc.) for use below.
                dXhf1 = [dXhf1; dXhf1(end-23:end,:)];
            else
                error(['Problem finding enough days of data from subdelta reference period to apply to years ',int2str(yearpminv(j)),' to ',int2str(yearpmaxv(j))])
            end
        end
        correction_factor1 = ones(nt1,ns);
        correction_factor(sel1,:) = correction_factor1;

        if verbose>0;
            disp(['Adding subdelta variability to projections for period ',int2str(yearpminv(j)),' through ',int2str(yearpmaxv(j))]);
        end

        if (loop_over_projection_times==0 && strcmp(delta_timescale,'daily')==1)
            if (fractional_subdelta==1)
                Xpcf(self1,:) = minXsub +  kron(((Xpc(sel1,:)-minXsub).*correction_factor1),ones(24,1)) .* dXhf1;
                if remove_deltascale_var==1;
                    Xpcf(self1,:) = Xpcf(self1,:) + (Xpci(self1,:)-kron(Xpc(sel1,:),ones(24,1)));
                end % Because dXhf1 in this case lacks the deltascale variation

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

            for iblock=1:nblocks %Optional loop over blocks e.g. to limit memory consumption
                if (subdelta_by_blocks==1)
                    selblock = find(latp>=latL_blocks1(iblock) & latp<latH_blocks1(iblock) & lonp>=lonL_blocks1(iblock) & lonp<lonH_blocks1(iblock));
                else
                    selblock = 1:ns;
                end
                nselblock = length(selblock);

                countf = 0;
                for i=1:nt1 %Loop over projection times
                    if (strcmp(delta_timescale,'daily')==1)
                        self2 = self1(floor(tdpf(self1))==floor(tdp(sel1(i)))); %Times to fill in finescale projection time series
                    elseif (strcmp(delta_timescale,'monthly')==1)
                        self2 = self1(yearpf(self1)==yearp(sel1(i)) & monthpf(self1)==monthp(sel1(i)));
                    end
                    nself2 = length(self2);

                    if (fractional_subdelta==1)
                        Xpcf(self2,:) = minXsub + (ones(nself2,1)*((Xpc(sel1(i),:)-minXsub).*correction_factor(sel1(i),:))).*dXhf1(countf+1:countf+nself2,:);
                        if remove_deltascale_var==1; Xpcf(self2,:) = Xpcf(self2,:) + (Xpci(self2,:)-ones(nself2,1)*Xpc(sel1(i),:)); end %Because dXhf1 in this case lacks the deltascale variation
                    else
                        Xpcf(self2,:) = Xpci(self2,:) + (ones(nself2,1)*correction_factor(sel1(i),:)).*dXhf1(countf+1:countf+nself2,:);
                    end
                    countf = countf + nself2;

                end
            end
        end
        if verbose>0;
            disp('Done adding subdelta variability to projections');
        end
    end
    if dtdhf1==0;
        Xpcf(end,:) = Xpcf(end-1,:);
    end
    %If hourly data are 00:00, 01:00, ..., 23:00 then the tdpf set above will result in an extra
    %time at 00:00 on ceil(max(tdp))+1 --- this can be filled by repeating the last datum.

    %Impose lower/upper limits if provided
    if (~isempty(XpcfL)||~isempty(XpcfH))
        if verbose>0;
            disp('Imposing lower/upper limits on projections');
        end
        Xpcfo = Xpcf; %Can be a useful diagnostic to check that capping is not too severe.
        if ~isempty(XpcfL);
            Xpcf = max(XpcfL,Xpcf);
        end
        if ~isempty(XpcfH);
            Xpcf = min(XpcfH,Xpcf);
        end
        if ~isequaln(Xpcfo,Xpcf);
            out.Xpcfo = Xpcfo;
        end
        clear Xpcfo
        %Note we do not subselect for not-NaN (too slow), so this will set any NaN values to XpcL
    end

    if (use_XpcfH_SWI>0)
        optSWIf = optSWI; optSWIf.calc_dav = 0; optSWIf.year = yearpf;
        if use_SWI_hav==1;
            optSWIf.calc_hav = 1;
        end
        if verbose>0;
            disp('Calculating maximum hourly values of shortwave irradiance');
        end
        out.XpcfH_SWI = NaN*ones(length(tdpf),ns);
        for i=1:ns %We have to loop over the time series here, else we run out of memory.
            if (use_SWI_hav==1)
                [~,outSWI] = calc_surfaceSWI(latp(i),lonp(i),yrdaypf,optSWIf);
                out.XpcfH_SWI(:,i) = outSWI.Q0hav;
            else
                out.XpcfH_SWI(:,i) = calc_surfaceSWI(latp(i),lonp(i),yrdaypf,optSWIf);
            end
            if (mod(i,100)==0 && verbose>0);
                disp(['Done ',int2str(i),' out of ',int2str(ns),' time series']);
            end
        end
        if verbose>0;
            disp('Done maximum hourly values of shortwave irradiance');
        end

        if (use_XpcfH_SWI==1)
            if verbose>0;
                disp('Capping all hourly projections of shortwave irradiance with instantaneous clear-sky maximum values');
            end
            Xpcf = min(out.XpcfH_SWI,Xpcf); %Cap the corrected hourly projections with instantaneous clear-sky maximum values
            %Warning: This is can lead to widespread correction of values near sunrise/sunset due to use of hourly-average irradiance in ERA5

        elseif (use_XpcfH_SWI==2)
            if verbose>0;
                disp('Calculating daily maximum values of projected shortwave irradiance');
            end
            Xpcf_dmax = NaN*Xpcf;
            tdu = unique(floor(tdpf));
            for i=1:length(tdu)
                self1 = find(floor(tdpf)==tdu(i));
                Xpcf_dmax(self1,:) = ones(length(self1),1)*max(Xpcf(self1,:));
            end
            if verbose>0;
                disp('Done daily maximum values of projected shortwave irradiance');
            end
            if verbose>0;
                disp('Capping selected hourly projections of shortwave irradiance with instantaneous clear-sky maximum values');
            end
            rec = find(Xpcf(:)>=frcrit_capSWI*Xpcf_dmax(:)); %Note this subselection can be slow
            Xpcf(rec) = min(out.XpcfH_SWI(rec),Xpcf(rec));
        end
    end
    if verbose>0;
        disp('Done imposing lower/upper limits on finescale projections (if any)');
    end
    if verbose>0;
        disp('Done finescale projections');
        toc;
    end

    if (recalculate_Xpc==1)
        out.Xpco2 = Xpc;
        if (strcmp(delta_timescale,'daily')==1) %Hourly variability about daily delta scale
            if verbose>0;
                disp('Recalculating daily averages from sub-daily projections');
            end
            for i=1:ntp;
                Xpc(i,:) = mean(Xpcf(floor(tdpf)==floor(tdp(i)),:));
            end
            if verbose>0;
                disp('Done recalculated daily averages');
            end
        elseif (strcmp(delta_timescale,'monthly')==1) %Hourly/daily variability about monthly delta scale
            if verbose>0;
                disp('Recalculating monthly averages from sub-monthly projections');
            end
            for i=1:ntp;
                Xpc(i,:) = mean(Xpcf(yearpf==yearp(i) & monthpf==monthp(i),:));
            end
            if verbose>0;
                disp('Done recalculated monthly averages');
            end
        end
    end

    out.Xpcf = Xpcf; out.tdpf = tdpf;
end