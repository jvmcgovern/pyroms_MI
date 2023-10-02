
%NOTE: This is a simplified script to demonstrate/test calculation of bias-corrected projections using
%      the matlab function fn_biascorrect_projections.m, applied to a test dataset of NORESM hindcast
%      and projection data covering the North Sea.
%      For an example of a real script used to make bias-corrected datasets and various analysis plots
%      see example_scripts/make_model_biascorrected_forcings.m.


varstrp = 'no3'; %Variable name for projection data (raw and bias-corrected)

varstrh = varstrp; if strcmp(varstrp,'no3'); varstrh = 'no3no2'; end %Hindcast variable name (could be different)


%% Read in the raw projection data from netcdf.
fnamep = ['../testdata/NorthSea/',varstrp,'_Omon_NorESM2-MM_historical_ssp585_r1i1p1f1_gr_198001-210012_NorthSea.nc'];
Xp = ncread(fnamep,varstrp);
disp(['Loaded ',varstrp,' from: ',fnamep])
yearp1 = (1980:2099)'; %Sequence of years covered by projection data  
%Projection data times are in a strange unit -- easier to just redefine the times rather than read from file:
yearp = kron(yearp1,ones(12,1));
monthp = kron(ones(length(yearp1),1),(1:12)');
dayp = 15*ones(length(yearp1)*12,1);
tdp = datenum(yearp,monthp,dayp);
tyrp = fndatenum_to_decyr(tdp);
ntp = length(tdp);


%% Read in the hindcast data from netcdf.
fnameh = ['../testdata/NorthSea/NORESMOCv1p2_40N_1980_2020_',varstrh,'_reanal_NorthSea.nc'];
Xh = ncread(fnameh,varstrh);
disp(['Loaded ',varstrh,' from: ',fnameh])
%No need to reread the spatial coordinate variables --- these are the same as in the projection data
tdh = datenum(1970,1,1) + ncread(fnameh,'time'); %NORESM reanalysis 'time' variable in units 'days since 1970-01-01'
[yearh,monthh,dayh] = datevec(tdh);
tyrh = fndatenum_to_decyr(tdh);
nth = length(tdh);
%We take the spatial coordinates from the hindcast file --- the coordinates in the projection file
%are the same except that the range 0 to 360 is used for longitude rather than -180 to 180 as in the hindcast file. 
lat1 = ncread(fnameh,'lat'); lon1 = ncread(fnameh,'lon'); %2D arrays, since NORESM data not on regular lat-lon grid
z1 = ncread(fnameh,'depth'); %Nominal depths at midpoints of depth bands
nz1 = length(z1);


%% Reorganize 4D arrays (Xp,Xh) into 2D arrays [ntp,nth x ns] for use in fn_biascorrect_projections.m.
nx1 = size(Xp,1); ny1 = size(Xp,2);
ns0 = nz1*nx1*ny1; %Total number of grid cells before parsing out the land mask
Xhc = Xh; Xpc = Xp;
Xh = NaN*ones(nth,ns0); Xp = NaN*ones(ntp,ns0);
lat = NaN*ones(1,ns0); lon = lat; z = lat; q = 0;
for i=1:nz1
    for j=1:ny1
        for k=1:nx1
            q = q+1;
            Xh(1:nth,q) = squeeze(Xhc(k,j,i,1:nth));
            Xp(1:ntp,q) = squeeze(Xpc(k,j,i,1:ntp));
            lat(q) = lat1(k,j); lon(q) = lon1(k,j);
            z(q) = z1(i);
        end
    end
end
clear Xhc Xpc


%% Parse out the land mask, recording wet grid selection indices (needed later for writing to netcdf).
selan = find(~isnan(Xh(1,:)) & ~isnan(Xp(1,:)));
%Note: In this example both hindcast and projection data have the exact same land mask (check: sum(isnan(Xh(1,:))~=isnan(Xp(1,:)))=0).
%      If this is not the case it may be desirable to fill with nearest grid cells for data missing
%      in projection and/or hindcast datasets (filling land masks where one of (projection,hindcast) is available).
Xh = Xh(:,selan); Xp = Xp(:,selan); 
lat = lat(selan); lon = lon(selan); z = z(selan);
ns = length(selan); %Total number of grid cells after parsing out the land mask


%% Set options for fn_biascorrect_projections.m.
%Note: these comprise a few basic 'essential' option parameters (Xp,tdp,Xh,tdh,yearminref,yearmaxref,method)
%      plus various extra option parameters passed via the structure variable 'opt'
yearminref = 2000; yearmaxref = 2019; %Define reference (training) period (inclusive year limits)
method = 0; %0 for Delta Change (DC), e.g. Kay and Butenschon (2018)
            %1 for Empirical Quantile Mapping (EQM), e.g. Cannon et al. (2015)
            %2 (def) for Repeated Hindcast Correction (RHC), e.g. Willems and Vrac (2011)
            %Can also input as string ('DC','EQM,'RHC')
            %See guidance text within fn_biascorrect_projections.m for a descriptive and discussion of the different methods.
opt = struct('seasonal',1,'use_month',1,'fractional',1,'tol_dtclim',0,'qsel',0.05:0.1:0.95,...
    'yearpminv',1980:20:2080,'yearpmaxv',1999:20:2099); %See fn_biascorrect_projections.m for description of option parameters
%Note: The reference parameter choices were based on some limited testing of different value combinations
%      with the aim of optimizing the skill in reproducing monthly quantiles during the validation period (1980-1999),
%      considering *the whole pan-Arctic domain* i.e. at all latitudes >40N. These values might not be optimal
%      for the North Sea, and might not necessarily produce most-plausible results for many decades in the future.
%      A general finding seems to be, however, that the RHC method (2) is a lot less sensitive to these choices
%      of parameter values, and might therefore be considered a more robust method.
if method==0; mstr1 = 'DC'; end
if method==1; mstr1 = 'EQM'; end
if method==2; mstr1 = 'RHC'; end


%% Compute quantile-corrected projections (Xpc) covering 120 years (ntp=1440) and ns=7605 wet grid cells.
[Xpc,out] = fn_biascorrect_projections(Xp,tdp,Xh,tdh,yearminref,yearmaxref,method,opt); 
%DC method takes ~0.3 seconds
%EQM method takes ~4.2 seconds (~6.3 seconds with opt.use_quantile = 0)
%RHC method takes ~17 seconds (~23 seconds with opt.use_quantile = 0)
%Times recorded on Dell Precision 5530 Intel(R) Core(TM) i7-8850H CPU @ 2.60GHz, 32GB RAM, running Matlab R2019b


%% Compare with reference solution
fnamepcr = ['../testdata/NorthSea/reference_solutions/',varstrh,'_NorESM2-NORESMreanalv1_1980-2099_NorthSea_',mstr1,'_v1_ref.nc'];
Xpcrc = ncread(fnamepcr,varstrh);

Xpcr = NaN*ones(ntp,ns0); q = 0;
for i=1:nz1
    for j=1:ny1
        for k=1:nx1
            q = q+1;
            Xpcr(:,q) = squeeze(Xpcrc(k,j,i,:));
        end
    end
end
clear Xpcrc
Xpcr = Xpcr(:,selan);
disp(['Maximum discrepancy with reference solution = ',num2str(max(max(abs(Xpcr-Xpc)))),' mol/m3']) %Note: typical rounding errors due to data compression ~1e-9 mol/m3 (negligible) 


%% Plot the solution!
%It is recommended to check/analyse the solution with several different plots, of which
%the quantile-quantile plots using the validation hindcast dataset (1980-1999) are especially
%recommended and may form the basis of a choise between different bias correction methods.
%See example_scripts/make_model_biascorrected_forcings.m for examples.
%Here we limit ourselves to perhaps the most basic analysis which is to compare the hindcast,
%projection, and bias-corrected projection time series at a selected point.
doplot_timeseries = 1;
if (doplot_timeseries==1)
    tfontsize = 11;
    lfontsize = 13;
    lfontweight = 'bold';
    interpreter1 = 'tex';
    showtitle = 1;
    varstr1 = varstrh;
    tolstr = [', tol=',int2str(opt.tol_dtclim),' month(s)'];
    if method==0; qstr1 = []; else qstr1 = [', nq=',int2str(length(opt.qsel))]; end
    if method==0; cstr1 = 'climatology'; else cstr1 = 'quantile'; end

    %latp1 = 59.24; lonp1 = 2.92; %N. North Sea: both EQM and RHC reference solutions look reasonable here
    latp1 = 54.88; lonp1 = 2.85; %S. North Sea: reference EQM solution looks suspect due to low and constant winter no3no2 --- solution looks better with tol_dtclim = 3
    dist1 = fn_spherical_distm(lat,lon,latp1,lonp1);
    sels = find(dist1==min(dist1)); sels = sels(z(sels)==5);
    if ns==1; sels = 1; end

    figure(40); clf;
    show_Xp = 1;
    show_Xpc = 1;
    Xh_on_top = 1;
    tyrminshow = 1980;
    
%    selh = find(yearh>=tyrminshow); % REDO
%    selp = find(yearp>=tyrminshow); % REDO
    zstr1 = [int2str(round(z(sels))),'m'];
    tstr1 = [varstr1,' monthly averages at (',num2str(lat(sels),'%3.1f'),'N,',num2str(lon(sels),'%3.1f'),'E,',zstr1,') from NORESM reanalysis v1 (black)'];
    
    if show_Xp==1; plot(tyrp(selp),Xp(selp,sels),'r-'); tstr1 = [tstr1,'; NORESM2-SSP585 (red)']; end
    hold on; plot(tyrh(selh),Xh(selh,sels),'k-')
    if show_Xpc==1; hold on; plot(tyrp(selp),Xpc(selp,sels),'c-'); tstr1 = [tstr1,';\newlineNORESM2-SSP585 with monthly ',cstr1,' correction to NORESM reanalysis v1 (cyan)']; end
    if Xh_on_top==1; hold on; plot(tyrh(selh),Xh(selh,sels),'k-'); end
    tstr1 = [tstr1,' (',mstr1,tolstr,qstr1,')'];
    axis tight
    
    if showtitle==1; title(tstr1,'FontSize',tfontsize,'FontWeight','bold','Interpreter',interpreter1); end
    xlabel('Year','FontSize',lfontsize,'FontWeight',lfontweight,'Interpreter',interpreter1);
    ylabel([varstr1,' [mol/m3]'],'FontSize',lfontsize,'FontWeight',lfontweight,'Interpreter',interpreter1);
    box on
 
    %Useful post-processing to remove x10 label in y ticks:
    %ax = gca; ax.YAxis.Exponent = 0;
end


%% Write to NetCDF if req'd
write_nc = 1;
if (write_nc==1)
    verstr = ['_',mstr1,'_v1'];
    
    yearwmin = 1980; yearwmax = 2099; %Limits defining subset of years for writing to netcdf
    modelstr1 = 'NorESM2-NORESMreanalv1';
    fnamepc = [varstrh,'_',modelstr1,'_',int2str(yearwmin),'-',int2str(yearwmax),'_NorthSea',verstr,'.nc'];

    %Subset of time steps for writing
%    sel1 = find(yearp>=yearwmin & yearp<=yearwmax); % REDO
    np1 = length(sel1);

    %Many different ways of writing to netcdf. Here we use the format from the hindcast file
    %modify the number of time steps to fit the bias-corrected projections.
    finfo = ncinfo(fnameh);
    finfo.Filename = fnamepc;
    varstrt = 'time';
    finfo.Dimensions(strcmp({finfo.Dimensions.Name},varstrt)).Length = np1;
    finfo.Variables(strcmp({finfo.Variables.Name},varstrt)).Dimensions.Length = np1;
    dimnames1 = {finfo.Variables(strcmp({finfo.Variables.Name},varstrh)).Dimensions.Name};
    finfo.Variables(strcmp({finfo.Variables.Name},varstrh)).Dimensions(strcmp(dimnames1,varstrt)).Length = np1;
    
    %Add history variable as global attribute. This is important to document the options used for the bias correction.
    sourcehstr = 'NORESM reanalysis v1';
    history1 = [date,': Bias-corrected projections calculated using make_testdata_biascorrected_forcings.m with method=',int2str(method),' (',mstr1,')',...
        ', yearminref=',int2str(yearminref),', yearmaxref=',int2str(yearmaxref),': ',...
        newline,'             ',sourcehstr,' file(s): ',fnameh,';',newline,'             ESM file(s): ',fnamep,'.'];
    history1 = [history1,newline,'options used for fn_biascorrect_projections.m:',newline,evalc('disp(opt)')];
    finfo.Attributes(strcmp({finfo.Attributes.Name},'history')).Value = history1;
    
    ncwriteschema(fnamepc,finfo);
    
    timep = tdp(sel1) - datenum(1970,1,1); %NORESM reanalysis 'time' variable in units 'days since 1970-01-01'
    ncwrite(fnamepc,varstrt,timep);
    
    nccopy_vars({'depth','lat','lon'},fnameh,fnamepc,struct('write_atts',0,'newfile',0));
    ncwrite(fnamepc,'depth',z1)
    
    ncwriteatt(fnamepc,varstrh,'long_name',varstrh);

    disp(['Writing bias-corrected projections to: ',fnamepc]); tic;
    for j=1:np1
        X1c = NaN*ones(1,ns0);
        X1c(:,selan) = Xpc(sel1(j),:);
        X1c = reshape(X1c,nx1,ny1,nz1);
        ncwrite(fnamepc,varstrh,X1c,[1 1 1 sel1(j)]);
        if mod(j,100)==0; disp(['Done writing ',int2str(j),' of out ',int2str(np1),' time steps to netcdf']); end
    end
    clear X1c
    disp(['Done writing bias-corrected projections to ',fnamepc]); toc;
    
    %NOTE: Despite the fact that finfo.Variables.DeflateLevel = 5 for all variables, it seem ncwrite.m does not write
    %      compressed netcdf data. I therefore compress the netcdf afterwards using NCO (-> typically ~10x smaller file size):
    %      ncks -L 5 <fnamepc> tmp.nc
    %      mv tmp.nc <fnamepc>
end
