model = 'roho800'

oceanmod = 1
small_screen = 1

model_proj = 'noresm'
scenstr1 = 'ssp585'
use_daily = 1
interp_to_gregorian = 1
if use_daily==0
    interp_to_gregorian = 0
end
test_method = 0

seasonal = 1
use_month = 0
if use_daily==1
    use_month = 0
end
correct_iavar = 0
#Hindcast year range for file names
yearmin0 = 2000
yearmax0 = 2019
tclimc = 0.5:365.5

method = 2 #1 for Empirical Quantile Mapping (EQM), 2 for Repeated Hindcast Correction (RHC)
yearminref = 2010 yearmaxref = 2019 #Reference period to define climatologies / quantile correction factors
legacy = 0

use_subdelta = 1
remove_deltascale_var = 1
correct_substd = 0
correct_subquantiles = 1
if method==2
    correct_subquantiles = 0
end
match_subdelta_deltamean = 0
match_subdelta_hourly = 0
use_calc_surfaceSWI = 1

if (oceanmod==1)
    dir1 = 'D:\USERS\Phil\ROHO800\ROHO800_input\projections\'
else
    dir0 = 'C:\Data\ROHO800\ROHO800_input\projections\'
    dir1 = dir0
end

##Grid point subselection indices (used only for atmospheric forcings)
latsel = 1:25 lonsel = 1:37 #Whole ERA5 domain for ROHO800

nlatsel = length(latsel) nlonsel = length(lonsel)
ns = nlatsel*nlonsel

##Atmospheric forcing variables
#varstr = 'Tair'
#varstr = 'Pair'
#varstr = 'Qair'
#varstr = 'cloud'
#varstr = 'Uwind'
#varstr = 'Vwind'
#varstr = 'swrad'
#varstr = 'lwrad_down'
#varstr = 'rain'

##Boundary condition variables
#varstr = 'temp'
#varstr = 'salt'
#varstr = 'zeta'
#varstr = 'ubar'
#varstr = 'vbar'
#varstr = 'u'
#varstr = 'v'

#varstr = 'O3_c'
#varstr = 'O3_TA'
#varstr = 'O2_o'
#varstr = 'N1_p'
#varstr = 'N3_n'
varstr = 'N5_s'

#bcstrv = '_west'
#bcstrv = '_north'
#bcstrv = '_east'
#bcstrv = '_south'
bcstrv = {'_west','_north','_east','_south'}
#bcstrv = {'_north','_east','_south'}
if ~iscell(bcstrv)
    bcstrv = {bcstrv}
end

do_ocean_bcs = 0
dir0 = 'D:\USERS\Phil\ROHO800\ROHO800_input\ERA5\concatenated\'
source0str = 'ERA5 hindcast'
if any(strcmp(varstr,{'temp','salt','zeta','ubar','vbar','u','v','O3_c','O3_TA','O2_o','N1_p','N3_n','N5_s'}))
    do_ocean_bcs = 1
    dir0 = 'D:\USERS\Phil\ROHO800\ROHO800_input\old_vertgrid_pre_v2e\'
    if any(strcmp(varstr,{'temp','salt','zeta','ubar','vbar','u','v'}))
        source0str = 'GLORYS12 reanalysis'
    else
        source0str = 'NORESM reanalysis v1'
    end
end
if (do_ocean_bcs==1 && test_method==1)
    yearminref = 2014
    yearmaxref = 2019
end
if do_ocean_bcs==1
    nbcstr = length(bcstrv)
else
    nbcstr = 1
end

for ibcstr=1:nbcstr
if do_ocean_bcs==1
    bcstr = bcstrv{ibcstr}
end

XpcL = [np.nan, np.nan, np.nan, np.nan,  np.nan,  np.nan,  np.nan,  0]
XpcH = [np.nan, np.nan, np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan]
XpcfL = [np.nan, np.nan, np.nan, np.nan,  np.nan,  np.nan,  np.nan,  0]
XpcfH = [np.nan, np.nan, np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan]
minX = [np.nan, np.nan, np.nan, np.nan,  np.nan,  np.nan,  np.nan,  np.nan]
minXsub = [np.nan, np.nan, np.nan, np.nan,  np.nan,  0,  np.nan,  -1e-4]
atmvars =    ['T2M', 'pmer', 'Q', 'U10M', 'V10M', 'SSR', 'STRD', 'TP']
tol_dtclim = [50, 50, 50, 5, 5, 100, 100, 100]
fractional = [0, 0, 0, 0, 0, 0, 0, 0]
fractional_subdelta = [1, 1, 1, 1, 1, 1, 1, 1]
seasonal = [1, 1, 1, 0, 0, 1, 1, 1]
atmtstr = ['time', 'time', 'time', 'time', 'time', 'time', 'time', 'time']
frcrit_capSWI = 0.0

atmvars =   ['temp', 'salt', 'zeta','ubar','vbar','u','v', 'DIC', 'TALK', 'NO3', 'PO4', 'Si', 'O2', 'NH4', 'FER']
use_daily = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
use_month = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
tol_dtclim = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
fractional = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
fractional_subdelta = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
seasonal = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
atmtstr = ['temp_time', 'salt_time', 'zeta_time', 'v2d_time', 'v2d_time', 'v3d_time', 'v3d_time', 'dic_time', 'talk_time', 'no3_time', 'po4_time', 'si_time', 'o2_time', 'nh4_time', 'fer_time']

tstrv = ['January','February','March','April','May','June','July','August','September','October','November','December']
if use_daily==1
    scstr2 = 'daily'
    if seasonal==1
        tolstr = [', tol=',int2str(tol_dtclim),' day(s)']
    else
        tolstr = ', non-seasonal'
    end
else
    scstr2 = 'monthly'
    if seasonal==1
        tolstr = [', tol=',int2str(tol_dtclim),' month(s)']
    else
        tolstr = ', non-seasonal'
    end

end

if (use_daily==1)
    fname11 = 'NorESM2-MM.atm.day_HIST_SSP585_2000-2100_ROHO800.nc'
    fname1 = [dir1,fname11]
    if any(strcmp(varstr1,{'U','V'}))
        X1 = squeeze(ncread(fname1,varstr1,[lonsel(1) latsel(1) 1 1],[nlonsel nlatsel 1 Inf]))
    else
        X1 = squeeze(ncread(fname1,varstr1,[lonsel(1) latsel(1) 1],[nlonsel nlatsel Inf]))
    end
    if strcmp(varstr,'Tair')
        X1 = X1 - 273.15
    end
    if strcmp(varstr,'rain')
        # Convert [m/s] to [kg/m2/s] using density of rainwater = 1000 kg/m3
        X1 = X1*1000
    end

    time1 = ncread(fname1,'time')
    lat11 = squeeze(ncread(fname1,'lat',latsel(1),nlatsel))
    lon11 = squeeze(ncread(fname1,'lon',lonsel(1),nlonsel))

    td11 = (datenum(2001,1,1):datenum(2001,12,31))' # An example year with 365 days
    [~,month11,day11] = datevec(td11)
    year1 = kron((2000:2099)',ones(365,1)) # Repeating 365-day year
    month1 = kron(ones(100,1),month11)
    day1 = kron(ones(100,1),day11)
    year1 = year1(2:end)
    month1 = month1(2:end)
    day1 = day1(2:end) #For some reason we are missing Jan 1st 2000.
    td1 = datenum(year1,month1,day1) #This gives the datenumbers of the original data, in agreement with 'cdo sinfon'.
    #Note that they are skipping instances of February 29th, hence max(diff(td1))=2.
    disp(['Loaded ',varstr1,' from: ',fname1])

    fname01 = ['roho800_',varstr,'_',int2str(yearmin0),'_',int2str(yearmax0),'_daymean.nc']
    fname0 = [dir0,'daymean\',fname01]
    X0 = squeeze(ncread(fname0,varstrc,[lonsel(1) latsel(1) 1],[nlonsel nlatsel Inf]))
    time0 = ncread(fname0,varstrt)
    td0 = datenum(1948,1,1) + time0
    [year0,month0,day0] = datevec(td0)
    yrday0 = td0 - datenum(year0,1,1)
    tyr0 = fndatenum_to_decyr(td0)
    disp(['Loaded ',varstrc,' from: ',fname0])

    dtdmean = unique(td0-floor(td0)) #This is the shift of daymean timestamps wrt 00:00
    #For cdo, timestamps of daily means use the mean of the timestamps for all hourly data used,
    #which is the mean of 0000,0100,...,2300 in the case of non-flux data and 0030,0130,...,2330
    #in the case of flux data. We correct the projection time stamps for consistency with this:
    td1 = td1 + dtdmean #NORESM daily mean timestamps were at 0000 of each day.
    yrday1 = td1 - datenum(year1,1,1)

else
    if (do_ocean_bcs==1)
        fname11 = 'NorESM2-MM.ocn.mon_HIST_SSP585_1990-2100_ROHO800.nc'
        fname1 = [dir1,fname11]
        if any(strcmp(varstr,{'temp','salt','zeta','ubar','vbar','u','v'}))
            fname01 = 'roho800_bry_GLORYS_20070115_to_20200115_N25_fix.nc'
        else
            fname01 = 'roho800_v2ax_bry_NORESMOCv1p2reanal_20070115_to_20211215_closest.nc'
        end
        fname0 = [dir0,fname01]
        if any(strcmp(varstr,{'zeta','ubar','vbar','u','v'}))
            if method~=2 stop end
            #For the circulation variables we assume no climatic change, but in order to preserve the
            #eddy variability we do not force with climatology. Instead we force with repeated (and
            #uncorrected) blocks of hindcast variability. This is achieved by applying the RHC method
            #(method=2) with constant (unity-valued) projection data.
            X1 = 1
        else
            X1 = ncread(fname1,varstr1)
        end
        if any(strcmp(varstr,{'O3_c','O3_TA','O2_o','N1_p','N3_n','N5_s'}))
            X1 = X1*1000
        end #Convert [mol/m3] to [mmol/m3]
        lat1 = ncread(fname1,'latitude') lon1 = ncread(fname1,'longitude')
        depth_bnds = ncread(fname1,'lev_bnds') z11 = mean(depth_bnds)'
        nz11 = length(z11)

        X0 = ncread(fname0,varstrc)
        fnameg = [dir1,'../Grid/ROHO800_grid_fix5_unsmooth1.nc']
        Lats0 = ncread(fnameg,'lat_rho') Lons0 = ncread(fnameg,'lon_rho')
        H0 = ncread(fnameg,'h') M0 = ncread(fnameg,'mask_rho')

        year11 = (1990:2100)'
        td00 = datenum(1948,1,1)
        timeunit = 1/86400
    else
        fname1 = [dir1,'NorESM2-MM.atm.day_HIST_SSP585_2000-2100_ROHO800_monmean.nc']
        fname0 = [dir0,'roho800_',varstr,'_',int2str(yearmin0),'_',int2str(yearmax0),'_monmean.nc']
        if any(strcmp(varstr1,{'U','V'}))
            X1 = squeeze(ncread(fname1,varstr1,[lonsel(1) latsel(1) 1 1],[nlonsel nlatsel 1 Inf]))
        else
            X1 = squeeze(ncread(fname1,varstr1,[lonsel(1) latsel(1) 1],[nlonsel nlatsel Inf]))
        end
        if strcmp(varstr,'Tair')
            X1 = X1 - 273.15
        end
        X0 = squeeze(ncread(fname0,varstrc,[lonsel(1) latsel(1) 1],[nlonsel nlatsel Inf]))
        lat11 = squeeze(ncread(fname1,'lat',latsel(1),nlatsel))
        lon11 = squeeze(ncread(fname1,'lon',lonsel(1),nlonsel))
        year11 = (2000:2099)' td00 = datenum(1948,1,1) timeunit = 1
    end

    time1 = ncread(fname1,'time')
    year1 = kron(year11,ones(12,1))
    month1 = kron(ones(length(year11),1),(1:12)')
    day1 = 15*ones(length(year11)*12,1)
    td1 = datenum(year1,month1,day1)
    disp(['Loaded ',varstr1,' from: ',fname1])

    time0 = ncread(fname0,varstrt)
    td0 = td00 + time0*timeunit
    [year0,month0,day0] = datevec(td0)
    tyr0 = fndatenum_to_decyr(td0)
    disp(['Loaded ',varstrc,' from: ',fname0])
end

if (use_subdelta==1)
    fname0f1 = ['roho800_',varstr,'_',int2str(yearmin0),'_',int2str(yearmax0),'.nc']
    fname0f = [dir0,fname0f1]
    X0f = squeeze(ncread(fname0f,varstrc,[lonsel(1) latsel(1) 1],[nlonsel nlatsel Inf])) #'Finescale' hindcast data
    time0f = ncread(fname0f,varstrt)
    td0f = datenum(1948,1,1) + time0f
    [year0f,month0f,day0f,hour0f,minutes0f] = datevec(td0f)
    hour0fc = hour0f + minutes0f/60 #Decimal hours (may be 0.5,1.5,2.5,... etc. for flux variables)
    yrday0f = td0f - datenum(year0f,1,1)
    tyr0f = fndatenum_to_decyr(td0f)
    disp(['Loaded ',varstrc,' from: ',fname0f])
end

if (do_ocean_bcs==1)
    ##
    #Match each boundary condition location in X0 to a nearest grid point in X1
    #This should be consistent with the approach in make_model_inputs.m
    if strcmp(bcstr,'_south')
        lat01 = Lats0(:,1)
        lon01 = Lons0(:,1)
        h0 = H0(:,1)
        m0 = M0(:,1)
    end
    if strcmp(bcstr,'_north')
        lat01 = Lats0(:,end)
        lon01 = Lons0(:,end)
        h0 = H0(:,end)
        m0 = M0(:,end)
    end
    if strcmp(bcstr,'_west')
        lat01 = Lats0(1,:)'
        lon01 = Lons0(1,:)'
        h0 = H0(1,:)'
        m0 = M0(1,:)'
    end
    if strcmp(bcstr,'_east')
        lat01 = Lats0(end,:)'
        lon01 = Lons0(end,:)'
        h0 = H0(end,:)'
        m0 = M0(end,:)'
    end
    [Tcline,theta_s,theta_b,Vtransform,Vstretching,N,hc] = setROMSvgrid('roho800','_v2bg')

    if any(strcmp(varstr,{'zeta','ubar','vbar'}))
        N = 1
    end
    sz0 = size(X0)
    n0 = sz0(1)
    nt0 = length(td0)
    nt1 = length(td1)
    X11 = squeeze(X1(:,:,1,1)) isnan11 = isnan(X11)
    [nx1,ny1] = size(X11)
    ns = N*n0

    #Loop over n0 target profiles: interpolate closest NORESM grid point data over depth, store in X10
    interpolate_ESM = 1
    use_lin = 1
    if any(strcmp(varstr,{'zeta','ubar','vbar','u','v'}))
        interpolate_ESM = 0
    end
    disp('Interpolating NORESM projections to the ROHO800 boundary grid points')
    extrap_constz = 1     #1 to extrapolate a constant in z direction (def = 1, safer)
    mindist10 = NaN*ones(n0,1)
    X10 = NaN*ones(n0,N,nt1)
    z10 = NaN*ones(n0,N)
    if (interpolate_ESM==1)
        Fv = cell(1,nz11)
        for k=1:nz11
            #Horizontal linear interpolation, set triangulation for this NORESM depth level
            v = squeeze(X1(:,:,k,1))
            sel1k = find(~isnan(v))
            if (length(sel1k)>2) #Need at least 3 points for scattered interpolation
                x = lon1(sel1k)
                y = lat1(sel1k)
                v = v(sel1k)
                Fv{k} = scatteredInterpolant(x,y,v,'linear','none')
            end
        end
    end
    for j=1:n0
        if (m0(j)==1)
            #Calculate depth levels
            z01 = roms_zdepth(h0(j),hc,theta_s,theta_b,N,Vtransform,Vstretching)
            z01 = -1*z01    #ROMS depth levels are negative need to flip sign for use in Xplin.m below
            z10(j,:) = z01'

            if (interpolate_ESM==1)
                dist10 = fn_spherical_distm(lat1(:),lon1(:),lat01(j),lon01(j))
                dist10m = reshape(dist10,nx1,ny1)
                mindist10(j) = min(dist10(isnan11(:)==0))
                selxy1 = find(dist10==mindist10(j) & isnan11(:)==0) selxy1 = selxy1(1)
                [selx1,sely1] = find(dist10m==mindist10(j) & isnan11==0) selx1 = selx1(1) sely1 = sely1(1)

                #Linearly interpolate to ROMS depths levels (z01) for point j
                X11_nn = squeeze(X1(selx1,sely1,:,:)) #Nearest-neighbour NORESM data on NORESM vertical grid

                if (use_lin==1)
                    #Get linear scattered interpolants, looping over depth (shallow-to-deep) and time
                    X11_lin = NaN*X11_nn
                    xi = lon01(j)
                    yi = lat01(j)
                    for k=1:nz11
                        #Horizontal linear interpolation, set triangulation for this NORESM depth level
                        v = squeeze(X1(:,:,k,1))
                        sel1k = find(~isnan(v))
                        if (length(sel1k)>2) #Need at least 3 points for scattered interpolation
                            F = Fv{k}
                            for l=1:nt1
                                X111 = squeeze(X1(:,:,k,l))
                                F.Values = X111(sel1k)
                                X11_lin(k,l) = F(xi,yi)
                            end
                        end
                    end
                end
                X11 = X11_nn
                if use_lin==1
                    sel_lin = find(~isnan(X11_lin))
                    X11(sel_lin) = X11_lin(sel_lin)
                end

                selz = find(~isnan(X11(:,1)))
                X10(j,:,:) = Xplin(z01,z11(selz),extrap_constz)*X11(selz,:)

            else
                X10(j,:,:) = ones(N,nt1)
            end
        end
        if mod(j,10)==0
            disp(['Done ',int2str(j),' of ',int2str(n0),' boundary grid points'])
        end
    end
    disp('Done interpolation')

    #Rearrange (X1,X10) into 2D arrays (nt x nseries)
    X0c = X0
    X0 = NaN*ones(nt0,ns)
    X1 = NaN*ones(nt1,ns)
    lat1 = NaN*ones(1,ns)
    lon1 = lat1
    z1 = lat1
    for i=1:N
        for j=1:n0
            if N>1
                X0(:,(i-1)*n0+j) = squeeze(X0c(j,i,:))
            end
            X1(:,(i-1)*n0+j) = squeeze(X10(j,i,:))
            lat1((i-1)*n0+j) = lat01(j)
            lon1((i-1)*n0+j) = lon01(j)
            z1((i-1)*n0+j) = z10(j,i)
        end
    end
    if N==1
        X0 = X0c'
    end
    clear X0c X10

else
    if (ns>1)
        lat1 = NaN*ones(1,ns)
        lon1 = lat1
        X1c = NaN*ones(size(X1,3),ns)
        X0c = NaN*ones(size(X0,3),ns)
        if use_subdelta==1
            X0fc = NaN*ones(size(X0f,3),ns)
        end
        for i=1:nlonsel
            for j=1:nlatsel
                lat1((i-1)*nlatsel+j) = lat11(j)
                lon1((i-1)*nlatsel+j) = lon11(i)
                X1c(:,(i-1)*nlatsel+j) = squeeze(X1(i,j,:))
                X0c(:,(i-1)*nlatsel+j) = squeeze(X0(i,j,:))
                if use_subdelta==1
                    X0fc(:,(i-1)*nlatsel+j) = squeeze(X0f(i,j,:))
                end
            end
        end
        X1 = X1c
        X0 = X0c
        clear X1c X0c
        if use_subdelta==1
            X0f = X0fc
            clear X0fc
        end
    else
        lat1 = lat11
        lon1 = lon11
    end
end

if (use_daily==1 && interp_to_gregorian==1)
    td1o = td1
    X1o = X1
    td1 = (datenum(2000,1,2)+dtdmean:datenum(2099,12,31)+dtdmean)'
    [year1,month1,day1] = datevec(td1)
    yrday1 = td1 - datenum(year1,1,1)
    tyr1 = fndatenum_to_decyr(td1)

    X1 = interp1(td1o,X1o,td1) #This will fill the missing Feb 29ths.
end

#Correct bogus negative values
if (strcmp(varstr,'rain')==1)
    Xmin = 1e-9 #See min(X0(X0(:)>1e-18)) and min(X0f(X0f(:)>1e-18))
    X0 = max(Xmin,X0)
    if use_subdelta==1
        X0f = max(Xmin,X0f)
    end
end

## Check formula for maximum (clear-sky) irradiance, if needed
test_maxswrad_formulae = 0
if (strcmp(varstr,'swrad')==1 && test_maxswrad_formulae==1)
    sel0f = 1:length(year0f)
    #sel0f = find(year0f>=2007)

    sels = 1:ns
    #sels = find(lat1==60.5 & lon1==7)

    #optSWI = struct('ROMS_model',1) #Overshoot ~ 32.6 W/m2
    optSWI = struct('f_model',1,'decl_model',1,'use_eqtime',1,'Q_model',0,'cloud_model',0,...
    'ndays_year',365.2422,'use_tday_lag',1,'year',year0f(sel0f))
    # This is the 'best' model, minimizing overshoot (~6.5 W/m2 at 60.5N,7E ~21.4 W/m2 over full ROHO800 domain)
    # Note: switching off the Equation-of-time correction (use_eqtime=0) exacerbates overshoot (~14.3 W/m2 at 60.5N,7E ~32.8 W/m2 over full ROHO800 domain)
    # Note: switching off leap year correction (use_tday_lag=0) makes no significant difference (~6.5 W/m2 at 60.5N,7E ~20.0 W/m2 over full ROHO800 domain)

    maxswrad0f = calc_surfaceSWI(lat1(sels),lon1(sels),yrday0f(sel0f),optSWI)

    max(max(X0f(sel0f,sels) - maxswrad0f)) #Overshoot (should ideally be <0)
    #figureplot(yrday0f,maxswrad0f,'r.',yrday0f,X0f,'k.')
    #figureplot(hour0f,maxswrad0f,'r.',hour0f,X0f,'k.')
    stop
end

## Apply delta change and quantile correction using function fn_biascorrect_projections.m
tdp = td1
Xp = X1
[yearp,monthp,dayp] = datevec(tdp)
yrdayp = tdp - datenum(yearp,1,1)
tyrp = fndatenum_to_decyr(tdp)
np = length(tyrp)
opt = struct('seasonal',seasonal,'use_month',use_month,'use_subdelta',use_subdelta,'fractional',fractional,...
    'correct_iavar',correct_iavar,'dyear_iavar',10,'XpcL',XpcL,'XpcH',XpcH,'legacy',legacy)
if fractional==1
    opt.minX = minX
end
if use_daily==1
    opt.tclimc = tclimc
    opt.tol_dtclim = tol_dtclim
    opt.tol_dtclimp = 1/24
end
if (strcmp(varstr,'swrad')==1)
    opt.use_XpcH_SWI = use_calc_surfaceSWI
    opt.latp = lat1 opt.lonp = lon1
    opt.optSWI = struct('f_model',1,'decl_model',1,'use_eqtime',1,'Q_model',0,'cloud_model',0,...
        'ndays_year',365.2422,'use_tday_lag',1)
    opt.use_SWI_hav = 1
end

#opt.yearpminv = 2000:20:2080 opt.yearpmaxv = 2019:20:2099
opt.yearpminv = 2000:10:2090
opt.yearpmaxv = 2009:10:2099
if (method==2)
    opt.yearpminv = 2000:10:2090
    opt.yearpmaxv = 2009:10:2099
    #if test_method==1 opt.yearpminv = 2008:6:2026 opt.yearpmaxv = 2013:6:2031 end
end

if (use_subdelta==1) #Extra options needed if adding subdelta variability
    #opt.tdpfmax = datenum(2039,12,31)
    opt.Xhf = X0f opt.tdhf = td0f
    opt.remove_deltascale_var = remove_deltascale_var
    opt.correct_substd = correct_substd
    opt.correct_subquantiles = correct_subquantiles
    if (correct_substd==1||correct_subquantiles==1)
        opt.opt_subq = opt_subq
    end

    opt.match_subdelta_deltamean = match_subdelta_deltamean
    if match_subdelta_deltamean==0
        opt.yearminsub = yearminref
        opt.yearmaxsub = yearmaxref
    end
    opt.match_subdelta_hourly = match_subdelta_hourly

    opt.fractional_subdelta = fractional_subdelta
    if fractional_subdelta==1
        opt.minXsub = minXsub
    end
    opt.XpcfL = XpcfL opt.XpcfH = XpcfH

    if (strcmp(varstr,'swrad')==1)
        opt.use_XpcfH_SWI = 2*use_calc_surfaceSWI
        opt.frcrit_capSWI = frcrit_capSWI
    end
    # This can be time-consuming -- switch off for development purposes
    opt.recalculate_Xpc = 1
    # opt.recalculate_Xpc = 0
end

opt2 = opt
opt2.tol_dtclim = tol_dtclim
#opt2.qsel = 0.005:0.01:0.995 #This doesn't really help for method=1
if (method==2)
    if use_daily==1
        opt2.qsel = 0.005:0.01:0.995
    else
        opt2.qsel = 0.05:0.1:0.95
    end
end

###Compute quantile-corrected projections
[Xpc2,out2] = fn_biascorrect_projections(Xp,tdp,X0,td0,yearminref,yearmaxref,method,opt2)

if (use_subdelta==1)
    tdpf = out2.tdpf
    [yearpf,monthpf,daypf,hourpf,minutespf] = datevec(tdpf)
    hourpfc = hourpf + minutespf/60
    yrdaypf = tdpf - datenum(yearpf,1,1)
    tyrpf = fndatenum_to_decyr(tdpf)
    if doXpc==1
        Xpcf = out.Xpcf
    end
    Xpcf2 = out2.Xpcf
    npf = length(tdpf)

end

## Write to NetCDF if req'd
write_nc = 1
if (write_nc==1)
    ##
    verstr = '_v3'

    if (do_ocean_bcs==1)
        sel1 = find(yearp>=2000 & yearp<=2099)
        if any(strcmp(varstr,{'temp','salt','zeta','ubar','vbar','u','v'}))
            modelstr1 = 'NorESM2-GLORYS12reanal'
        else
            modelstr1 = 'NorESM2-NorESMreanal1'
        end
        fnamep1 = ['roho800_bry_',modelstr1,'_',scenstr1,'_2000-2099',verstr,'.nc']
        fname0f = fname0
        fname0str = fname01
        varstr1c = varstrc
        use_ncwriteschema = 0
        time_1day = 86400
        tdpf = tdp
        yearpf = yearp
    else
        sel1 = 1:npf
        fnamep1 = ['roho800_',varstr,'_NorESM2-ERA5_',scenstr1,'_',int2str(min(yearpf)),'-',int2str(max(yearpf)),verstr,'.nc']
        fname0str = [fname0f1,', ',fname01]
        varstr1c = varstr
        use_ncwriteschema = 1
        time_1day = 1
    end
    np1 = length(sel1)
    fnamep = [dir1,fnamep1]
    disp(['Writing bias-corrected projections to: ',fnamep1])
    tic

    if (use_ncwriteschema==1)
        finfo = ncinfo(fname0c) #varnamev = {finfo.Variables.Name}
        finfo.Filename = fnamep
        finfo.Dimensions(strcmp({finfo.Dimensions.Name},varstrt)).Length = np1
        finfo.Variables(strcmp({finfo.Variables.Name},varstrt)).Dimensions.Length = np1
        dimnames1 = {finfo.Variables(strcmp({finfo.Variables.Name},varstr1c)).Dimensions.Name}
        finfo.Variables(strcmp({finfo.Variables.Name},varstr1c)).Dimensions(strcmp(dimnames1,varstrt)).Length = np1

        history1 = finfo.Attributes(strcmp({finfo.Attributes.Name},'history')).Value
        history1 = [date,': Bias-corrected projections calculated using make_roho_biascorrected_forcings.m with method=',int2str(method),': ',...
            newline,'             ',source0str,' file(s): ',fname0str,'',newline,'             ESM file: ',fname11,'.',newline,history1]
        history1 = [history1,newline,'options used for fn_biascorrect_projections.m:',newline,evalc('disp(opt2)')]
        if isfield(opt2,'opt_subq')
            history1 = [history1,newline,'where opt_subq = ',newline,evalc('disp(opt_subq)'),newline,'(see opt2c_',varstr,'_ROHO800',verstr,'.mat)']
        end
        finfo.Attributes(strcmp({finfo.Attributes.Name},'history')).Value = history1

        ncwriteschema(fnamep,finfo)
    end

    timepf = time_1day*(tdpf(sel1) - datenum(1948,1,1))
    ncwrite(fnamep,varstrt,timepf)
    if (do_ocean_bcs==1 && strcmp(varstr,'zeta')==1)
        nccopy_vars({'theta_b','hc','theta_s','Tcline','s_rho','lat_psi','lon_psi','angle','h','lat_rho','lon_rho'},fname0,fnamep,struct('write_atts',0,'newfile',0))
    elseif (do_ocean_bcs==0)
        ncwrite(fnamep,'lat',lat11)
        ncwrite(fnamep,'lon',lon11)
    end
    yearpminv = 2000:10:2090
    yearpmaxv = 2009:10:2099
    nchunks = length(yearpminv)
    for i=1:nchunks
        self1 = find(yearpf(sel1)>=yearpminv(i) & yearpf(sel1)<=yearpmaxv(i))
        npf1 = length(self1)

        if (do_ocean_bcs==1)
            X1fc = NaN*ones(n0,N,npf1)
            for j=1:npf1
                X1fc(1:n0,1:N,j) = reshape(Xpc2(sel1(self1(j)),:),n0,N)
            end
        else
            X1fc = NaN*ones(nlonsel,nlatsel,npf1)
            for j=1:npf1
                X1fc(1:nlonsel,1:nlatsel,j) = reshape(Xpcf2(sel1(self1(j)),:),nlatsel,nlonsel)'
            end
        end

        if (do_ocean_bcs==1 && N==1)
            ncwrite(fnamep,varstr1c,squeeze(X1fc),[1 self1(1)])
        else
            ncwrite(fnamep,varstr1c,X1fc,[1 1 self1(1)])
        end
        disp(['Done writing ',varstr,' bias-corrected projections for years ',int2str(yearpminv(i)),' to ',int2str(yearpmaxv(i))])
    end
    disp(['Done writing bias-corrected projections to ',fnamep1])
    toc
    clear X1fc
end
end