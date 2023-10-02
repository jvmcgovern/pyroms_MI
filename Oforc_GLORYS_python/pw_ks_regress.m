

function [yp,out] = ks_regress(y,x,xp,bdw,opt)

%function [yp,out] = ks_regress(y,x,xp,bdw,opt)
%
%Kernel smoothing regression of y over m dimensions from x to xp
%Parameter 'bdw' [ni*m], [1*m], or scalar, describing smoothing bandwidths 
%for ni interpolation points and m smoothing dimensions.
%Optional parameters (opt.param):
%'kfunc' scalar defining form of kernel function:
%   kfunc = 1: e^(-s) kernel (default)
%   kfunc = 2: Gaussian kernel
%   kfunc = 3: sinc kernel
%   kfunc = 4: box-car running mean
%   kfunc = 5: Cauchy/Lorentzian kernel
%'obserrv' [nd*1] observational error variances to weight the data (optional)
%'rad' [1*m] allowed radii of influence (default = Inf, but this can be slow)
%
%Can also output weight matrix 'wt' for replicated calculations,
%and the indices 'recs' of the data used for the interpolation.
%Reference: Hastie et al., 2009.  The Elements of Statistical Learning.
%
%Validated against npreg.R from package "np" for local constant and local linear regressions
%See testks_regress.m for examples.
%
%Uses:  conformsizes.m (included in this file)
%       fnperiodic_distm.m
%
%%Phil Wallhead 10/02/2019




[ni,m] = size(xp);
nd = length(y);
if nd==0; disp('No data for kernel smoothing!'); stop; end

if nargin<5; opt = []; end
if isfield(opt,'verbose')==1; verbose = opt.verbose; else verbose = 0; end
if isfield(opt,'nblocks')==1; nblocks = opt.nblocks; else nblocks = []; end
    %Can specify here the number of blocks into which to divide the prediction set
if isfield(opt,'subset_wrt_xp')==1; subset_wrt_xp = opt.subset_wrt_xp; else subset_wrt_xp = (m<3); end
    %1 to subset the predictions with respect to the range of prediction x values xp
if isfield(opt,'one_by_one')==1; one_by_one = opt.one_by_one; else one_by_one = 0; end
    %1 to make predictions one by one (not usually fastest)
    if one_by_one==1; nblocks = ni; subset_wrt_xp = 0; end
if isfield(opt,'nimax')==1; nimax = opt.nimax; else nimax = Inf; end
    if (ni>nimax)
        nblocks = ceil(ni/nimax);
        if verbose>0; disp(['ks_regress: Prediction set (ni = ',int2str(ni),') exceeds maximum size (nimax = ',int2str(nimax),'); setting nblocks = ',int2str(nblocks)]); end
    end
    if (isempty(nblocks))
        if (nd>1e3 && ni>1e3)
            nblocks = 10;
            if verbose>0; disp(['ks_regress: Dividing prediction set into ',int2str(nblocks),' blocks to limit matrix sizes']); end
        else
            nblocks = 1;
        end
    end
if isfield(opt,'kfunc')==1; kfunc = opt.kfunc; else kfunc = 1; end 
    %Choice of smoothing kernel (def = exponential (Laplacian))
if isfield(opt,'one_sided')==1; one_sided = opt.one_sided; else one_sided = zeros(1,m); end
    %Set one_sided(i) = -1 to only consider influence of data with values x(:,i) <= xp(i)
    %    one_sided(i) = 1 to only consider influence of data with values x(:,i) >= xp(i)
if isfield(opt,'obserrv')==1; obserrv = opt.obserrv; else obserrv = []; end %Optional observational error variances for weighting
if isfield(opt,'yprior')==1; yprior = opt.yprior; else yprior = []; end %Optional prior values at the prediction points to stabilize estimates
    if (length(yprior)==1 && ni>1); yprior = yprior*ones(ni,1); end
if isfield(opt,'priorerrv')==1; priorerrv = opt.priorerrv; else priorerrv = ones(ni,1); end %Optional prior error variances for weighting
    if (length(priorerrv)==1 && ni>1); priorerrv = priorerrv*ones(ni,1); end
    if (isempty(obserrv) && ~isempty(priorerrv)); obserrv = ones(nd,1); end
if isfield(opt,'use_rad')==1; use_rad = opt.use_rad; else use_rad = 0; end
    %1 to use buffer vector "rad" to limit initial data set such that nds<=ndsmax (def = 0)
    if (use_rad==1)
        if isfield(opt,'variable_rad')==1; variable_rad = opt.variable_rad; else variable_rad = 0; end
        %1 to vary rad iteratively in order to limit initial data set such that ndsmin<=nds<=ndsmax (def = 0)
        if isfield(opt,'rad')==1; rad = opt.rad; else rad = []; end
        %maximum "radius" of allowed influence of data points (default = 10 * bdw)
        %Beware: restricting rad can lead to irregularities
        %Note: rad must be in units of the input coordinates (e.g. in degrees latitude/longitude if aslatlon==1)
        %Note: user inputs target or default rad/radfac; rad will be changed if nds<ndsmin or nds>ndsmax
        if isfield(opt,'radfac')==1; radfac = opt.radfac; else radfac = []; end
        if (isempty(rad)==1)
            if (isempty(radfac)==1)
                if verbose>0; disp('Assuming default radii of influence rad = 10 * bdw (radfac = 10)'); end
                radfac = 10;
                rad = radfac*bdw;
            else
                rad = radfac*bdw;
            end
        end
        if isfield(opt,'radfacmin')==1; radfacmin = opt.radfacmin; else radfacmin = 2; end
        %Minimum constraint on rad: radmin = radfacmin*bdw
        if isfield(opt,'radfacmax')==1; radfacmax = opt.radfacmax; else radfacmax = []; end
        %Maximum constraint on rad: radmax = radfacmax*bdw
        if isempty(radfacmax)==0; rad = radfacmax*bdw; else radfacmax = Inf; end
        if isfield(opt,'ndsmin')==1; ndsmin = opt.ndsmin; else ndsmin = min(nd,20); end
        %Increase ndsmin to avoid irregular smoothing in data-sparse regions
        if isfield(opt,'ndsmax')==1; ndsmax = opt.ndsmax; else ndsmax = 5000; end
        %Target maximum number of data (to limit matrix size)
        if isfield(opt,'incradmax')==1; incradmax = opt.incradmax; else incradmax = 1e2; end
    end
if isfield(opt,'ndclosest')==1; ndclosest = opt.ndclosest; else ndclosest = []; end
    %Maximum number of closest data (nearest neighbours) to use for smoothing (only for ni1 = 1)
    if ~isempty(ndclosest); nblocks = ni; subset_wrt_xp = 0; end %Loop over prediction points one-by-one if ndclosest is restricted
if isfield(opt,'lr_order')==1; lr_order = opt.lr_order; else lr_order = 0; end
    %Order of local polynomial regression (0 for local-constant Nadaraya-Watson method)
    %Can be vector [1 * m]
if isfield(opt,'xc')==1; xc = opt.xc; else xc = []; end
    %Additional linear covariates (an alternative to adding to x and setting extra bdws = Inf)  
if isfield(opt,'xpc')==1; xpc = opt.xpc; else xpc = []; end
    %Additional linear prediction covariates (an alternative to adding to xp and setting extra bdws = Inf) 
if isfield(opt,'lr_orderc')==1; lr_orderc = opt.lr_orderc; else lr_orderc = 0; end
    %Order of local polynomial regression over additional covariates  
    %Can be vector [1 * mc]
if isfield(opt,'period')==1; period = opt.period; else period = NaN*ones(1,m); end    
    %Can specify looping period such that dx_ij = mod(x_i-x_j,period)
if isfield(opt,'aslatlon')==1; aslatlon = opt.aslatlon; else aslatlon = []; end  
    %[1*2] vector specifying indices in (x,xp) of (latitude,longitude) for spherical distance treatment
    %Note: the assumed spherical distance bandwidth in km is bdw(aslatlon(1)),]; bdw(aslatlon(2)) is not used
    if (length(aslatlon)==1) 
        if aslatlon==1; aslatlon = [1 2]; end %Quick input aslatlon = 1 will assume aslatlon = [1 2]
        if aslatlon==0; aslatlon = []; end %Quick input aslatlon = 0 will assume aslatlon = []
    end
    if (isempty(aslatlon)==0 && isfield(opt,'rad')==0)
        if (verbose>0)
            disp('Warning: if interpolating over lat/lon, the allowed radii rad should be input (since cannot be inferred from bdw)')
            disp('Assuming rad = Inf for lat/lon')
        end
        rad(aslatlon) = Inf;
    end
if isfield(opt,'calcunc')==1; calcunc = opt.calcunc; else calcunc = 0; end
    %1 to calculate standard errors and confidence/prediction intervals (mpL,mpH)/(ypL,ypH)
if isfield(opt,'variable_res_var')==1; variable_res_var = opt.variable_res_var; else variable_res_var = 1; end
    %1 (def) to allow variable residual variance V = V(x), otherwise V = V0
if isfield(opt,'nneg')==1; nneg = opt.nneg; else nneg = 0; end
if isfield(opt,'vbdw')==1; vbdw = opt.vbdw; else vbdw = bdw; end
if isfield(opt,'vlr_order')==1; vlr_order = opt.vlr_order; else vlr_order = lr_order; end
if isfield(opt,'alph')==1; alph = opt.alph; else alph = 0.05; end
if isfield(opt,'check_evenness')==1; check_evenness = opt.check_evenness; else check_evenness = 0; end
    %1 to check evenness of the weight matrix (for dominant data, def = 0)
if isfield(opt,'domfrac')==1; domfrac = opt.domfrac; else domfrac = 0.95; end
if isfield(opt,'x_coverage')==1; x_coverage = opt.x_coverage; else x_coverage = []; end
    %Optional "coverage variable" size nd*1
    %When supplied, out.xp_coverage gives the number of different integer values of x_coverage used for each prediction point
if isfield(opt,'coverage_type')==1; coverage_type = opt.coverage_type; else coverage_type = 0; end
    %0 for simple coverage: range of values of x_coverage for data contributing to yp at each point
    %1 for number of unique values of x_coverage
if isfield(opt,'fr_wt_coverage')==1; fr_wt_coverage = opt.fr_wt_coverage; else fr_wt_coverage = 0; end
    %Minimum relative weight a data point must have to be included in the coverage calculation 
if isfield(opt,'max_sdist_coverage')==1; max_sdist_coverage = opt.max_sdist_coverage; else max_sdist_coverage = 2; end
    %Maximum scaled distance from the prediction point that a data point must have to be included in the coverage calculation
    
%Plotting options
if isfield(opt,'doplot')==1; doplot = opt.doplot; else doplot = 0; end
if isfield(opt,'ifig')==1; ifig = opt.ifig; else ifig = 17; end
if isfield(opt,'nrows')==1; nrows = opt.nrows; else nrows = 1; end
if isfield(opt,'ncols')==1; ncols = opt.ncols; else ncols = 1; end
if isfield(opt,'isub')==1; isub = opt.isub; else isub = 1; end
if isfield(opt,'plottype')==1; plottype = opt.plottype; else plottype = 1; end
if isfield(opt,'tstr')==1; tstr = opt.tstr; else tstr = []; end

%Options for plots against covariate (e.g. time)
if isfield(opt,'doplotc')==1; doplotc = opt.doplotc; else doplotc = 0; end
if isempty(xc)==1; doplotc = 0; end
if (doplotc==1)
    if isfield(opt,'max_sdist_plotc')==1; max_sdist_plotc = opt.max_sdist_plotc; else max_sdist_plotc = 2; end
    %Maximum scaled distance from the prediction point for a datum to be included in the plot vs. covariate
    if isfield(opt,'selpp')==1; selpp = opt.selpp; else selpp = 1; end
    %If smoothing blocks of > 1 prediction points, selpp specifies which prediction point will be plotted
    if isfield(opt,'xcpp')==1; xcpp = opt.xcpp; else xcpp = [min(xc); max(xc)]; end
    %Vector of covariate values at which to predict the response in the plot
    ncpp = length(xcpp);
    if isfield(opt,'xcoff')==1; xcoff = opt.xcoff; else xcoff = 0; end
    %Offset to apply to covariate x axis in plots
    if isfield(opt,'ifigc')==1; ifigc = opt.ifigc; else ifigc = 18; end
    if isfield(opt,'kcstep')==1; kcstep = opt.kcstep; else kcstep = Inf; end
    %Set kcstep>1 to only plot for every (kcstep)^th valid prediction point
    if isfield(opt,'xpp')==1; xpp = opt.xpp; else xpp = []; end
    %Alternatively provide a set of points at which to make a plot
    npp = length(xpp(:,1));
    if isfield(opt,'nrowsc')==1; nrowsc = opt.nrowsc; else nrowsc = 5; end
    if isfield(opt,'ncolsc')==1; ncolsc = opt.ncolsc; else ncolsc = 3; end
    if isfield(opt,'ycol')==1; ycol = opt.ycol; else ycol = 'k'; end
    if isfield(opt,'ypcol')==1; ypcol = opt.ypcol; else ypcol = 'r'; end
    if isfield(opt,'ymarkersize')==1; ymarkersize = opt.ymarkersize; else ymarkersize = 8; end
    if isfield(opt,'linewidth')==1; linewidth = opt.linewidth; else linewidth = 2; end
    if isfield(opt,'lfontsize')==1; lfontsize = opt.lfontsize; else lfontsize = 10; end
    if isfield(opt,'xstr')==1; xstr = opt.xstr; else xstr = 't'; end
    if isfield(opt,'ystr')==1; ystr = opt.ystr; else ystr = 'y'; end
    if isfield(opt,'xlimc')==1; xlimc = opt.xlimc; else xlimc = []; end
    if isfield(opt,'ylimc')==1; ylimc = opt.ylimc; else ylimc = []; end
end



bdw = conformsizes([ni m],bdw);
if (use_rad==1)
    if length(rad)==1; rad = rad*ones(1,m); end
end
if length(lr_order)==1; lr_order = lr_order*ones(1,m); end
if length(one_sided)==1; one_sided = one_sided*ones(1,m); end
mc = 0; isubc = 1; kc = 0;
if (isempty(xc)==0) 
    mc = length(xc(1,:));
    if isempty(xpc); xpc = zeros(ni,mc); end %Additional prediction covariates default to zero
    if length(xpc(:,1))==1; xpc = ones(ni,1)*xpc; end 
    if length(lr_orderc)==1; lr_orderc = lr_orderc*ones(1,mc); end
    if (doplotc==1) 
        figure(ifigc); clf;
        nsubc = nrowsc*ncolsc;
    end
else
    doplotc = 0; nsubc = 0;
end
%No. linear parameters fitted locally
nb1 = sum(lr_order);
nbc = sum(lr_orderc);
nb = 1 + nb1 + nbc;  %Note this does not yet allow for interaction terms


%Divide the prediction set into blocks if req'd (to limit matrix size)
if (nblocks>1)
    if subset_wrt_xp==0; bblock = [0 round((ni/nblocks:ni/nblocks:ni))]; end
    if (subset_wrt_xp==1)
        nblocks1 = round(nblocks^(1/m)); nblocks = nblocks1^m;
        minxp = min(xp); maxxp = max(xp);
        dxp = (maxxp - minxp)/nblocks1;
        xpL1 = NaN*ones(nblocks1,m); xpH1 = NaN*ones(nblocks1,m);
        for i=1:m; xpL1(:,i) = (minxp(i):dxp(i):maxxp(i)-dxp(i))'; xpH1(:,i) = (minxp(i)+dxp(i):dxp(i):maxxp(i))'; end
        if (m==1)
            xpL = xpL1; xpH = xpH1;
        elseif (m==2)
            [X,Y] = ndgrid(xpL1(:,1),xpL1(:,2)); xpL = [X(:) Y(:)]; clear X Y
            [X,Y] = ndgrid(xpH1(:,1),xpH1(:,2)); xpH = [X(:) Y(:)]; clear X Y
            %[xpL xpH]
            %stop
        else
            error('Block division calculation not yet coded for m>2')
        end
    end
    
    I1 = NaN*ones(ni,1); I2 = cell(1,nblocks);
    xsm = I2; radm = NaN*ones(nblocks,m);
    ypf = I1; ndsf = I1; 
    if check_evenness==1; ndsdomf = I1; end
    Kf = I2; wtf = I2; recsf = I2;
    if isempty(x_coverage)==0; xp_coveragef = zeros(ni,1); end
    if nb>1; bf = NaN*ones(nb,ni); end
    if (calcunc==1)
        rhatf = I2; varhatf = I2; varpf = I1;
        mpsef = I1; mpLf = I1; mpHf = I1;
        ypsef = I1; ypLf = I1; ypHf = I1;
    end
end


%Add prior values as pseudodata at the prediction points to stabilize estimates if req'd
if (~isempty(yprior))
    y = [y; yprior];
    x = [x; xp];
    if ~isempty(obserrv); obserrv = [obserrv; priorerrv]; end
    if (~isempty(xc))
        error('Cannot use additional covariates if adding prior values as pseudodata')
    end
end
            
            

%Master loop of nblocks----------------------------------------------------
for iblock = 1:nblocks
    if (nblocks>1)
        if subset_wrt_xp==0; recp1 = bblock(iblock)+1:bblock(iblock+1); end
        if (subset_wrt_xp==1)
            recp1 = 1:ni;
            for i=1:m
% REDO                recp1 = recp1(xp(recp1,i)>=xpL(iblock,i) & xp(recp1,i)<=xpH(iblock,i));
            end
        end
        
        xp1 = xp(recp1,:); bdw1 = bdw(recp1,:);
        if isempty(xpc)==0; xp1c = xpc(recp1,:); end  %Predictions values of the covariate coords (default to zero)
    else
        xp1 = xp; bdw1 = bdw;
        if isempty(xpc)==0; xp1c = xpc; end           %Predictions values of the covariate coords (default to zero)
    end
    ni1 = length(xp1(:,1)); if use_rad==1; rad1 = rad; end
   
    
    if (ni1>0)
        %Establish subset recs of influential data
        recs = find(~isnan(y) & sum(isnan(x),2)==0);
        %Adapt for one-sided kernels
        xp1min = min(xp1,[],1); xp1max = max(xp1,[],1);
        for i=1:m
            if (one_sided(m)==-1)
% REDO                recs = recs(x(recs,i)<=xp1max(i));
            elseif (one_sided(m)==1)
% REDO                recs = recs(x(recs,i)>=xp1min(i));
            end
        end
        nds = length(recs);

        
        if (use_rad==1) %Use vector of allowed buffers (rad) to limit initial data matrix size
            recs0 = recs; nrecs0 = nds;
            incrad = 0; flag_exit = 0;

            if (nrecs0<ndsmin)                  %Full data set too small
                if verbose>1; disp(['ks_regress: Full data set smaller than target minimum size (nds = ',int2str(nds),'), hence using full set']); end
            else
                rad1 = rad;
% % REDO vvv               while ((nds<ndsmin||nds>ndsmax) && incrad<=incradmax && flag_exit==0)
                while ((nds<ndsmin||nds>ndsmax) && flag_exit==0)
                    recs = recs0;
                    for i=1:m
                        if (~isnan(period(i)))
                            sels = fnperiodic_limit(x(recs,i),xp1min(i)-rad1(i),xp1max(i)+rad1(i),period(i));
                            recs = recs(sels);
                        else
                            if (one_sided(m)==-1)
% REDO                                recs = recs(x(recs,i)>=xp1min(i)-rad1(i) & x(recs,i)<=xp1max(i));
                            elseif (one_sided(m)==1)
% REDO                                recs = recs(x(recs,i)>=xp1min(i) & x(recs,i)<=xp1max(i)+rad1(i));
                            else
% REDO                                recs = recs(x(recs,i)>=xp1min(i)-rad1(i) & x(recs,i)<=xp1max(i)+rad1(i));
                            end
                        end
                        %Exclude data with any one coordinate beyond allowed distance from xp1min, xp1max
                    end
                    nds = length(recs);
                    
                    if (variable_rad==0)
                        flag_exit = 1;
                    elseif (nds<ndsmin) %Try to increase nds by increasing buffer zones
                        if (any(2*rad1>radfacmax*min(bdw1,[],1)))
                            if verbose>1; disp(['ks_regress: Data set smaller that target minimum size (nds = ',int2str(nds),') but rad cannot be increased from ',num2str(rad1),' without exceeding (radfacmax=',num2str(radfacmax),')*bandwidth for one or more coords']); end
                            flag_exit = 1;
                        else
                            rad1 = 2*rad1; incrad = incrad + 1;
                        end
                    elseif (nds>ndsmax) %Try to reduce nds by reducing buffer zones
                        if (any(0.9*rad1<radfacmin*max(bdw1,[],1)))
                            if verbose>1; disp(['ks_regress: Data set exceeds target maximum size (nds = ',int2str(nds),') but rad cannot be reduced from ',num2str(rad1),' without dropping below (radfacmin=',num2str(radfacmin),')*bandwidth for one or more coords']); end
                            flag_exit = 1;
                        else
                            rad1 = 0.9*rad1; incrad = incrad + 1;
                        end
                    end 
                    if (incrad==incradmax && verbose>1); disp(['ks_regress: Reached maximum subset adjustment iterations incradmax = ',int2str(incradmax),', nds = ',int2str(nds)]); end
                end
                if (incrad>0 && nds<nrecs0 && verbose>1); disp(['ks_regress: Data set reduced from ',int2str(nrecs0),' to ',int2str(nds),'; rad = ',num2str(rad1),'; incrad = ',int2str(incrad)]); end
            end
        end
        
        if nds==0; if verbose>1; disp('ks_regress: No data within allowed radius for kernel smoothing'); end; yp = NaN*ones(ni,1); end
        
        
        if (nds>0)
            xs = x(recs,:); ys = y(recs);
            if isempty(xc)==0; xcs = xc(recs,:); end
            if nblocks>1; xsm{iblock} = xs; if use_rad==1; radm(iblock,:) = rad1; end; end
            if isempty(x_coverage)==0; x_coverages = x_coverage(recs); end
            
            if (nds>1e3 && ni1>1e3)
                if verbose>1; disp(['ks_regress: Warning: working with large matrices (ni1 = ',int2str(ni1),', nds = ',int2str(nds),'), could be slow!']); end
            end
            
            
            %Calculate total squared separation matrix (ssm)
            ssm = zeros(ni1,nds); iscartesian = ones(1,m);
            
            %First treat spherical distances using (lat,lon) inputs
            if (isempty(aslatlon)==0)
                iscartesian(aslatlon) = 0;
                pi180 = pi/180; R = 6378.137;
                lat1 = xp1(:,aslatlon(1))*pi180; lon1 = xp1(:,aslatlon(2))*pi180;
                lat2 = xs(:,aslatlon(1))*pi180; lon2 = xs(:,aslatlon(2))*pi180;
                coslat1coslat2 = cos(lat1)*cos(lat2)';
                dlat = lat1*ones(1,nds)-ones(ni1,1)*lat2';
                dlon = lon1*ones(1,nds)-ones(ni1,1)*lon2';
                distm1 = R*2*asin(sqrt(sin(dlat/2).^2+coslat1coslat2.*sin(dlon/2).^2));
                %Haversine formula for distance in km, poached from code from m_lldist.m by Rich Pawlowicz (validated against same)
                if (one_sided(aslatlon(1))==-1)
                    distm1(dlat<0) = Inf; %Exclude where data > model
                elseif (one_sided(aslatlon(1))==1)
                    distm1(dlat>0) = Inf; %Exclude where data < model
                end
                if (one_sided(aslatlon(2))==-1)
                    distm1(dlon<0) = Inf; %Exclude where data > model
                elseif (one_sided(aslatlon(1))==1)
                    distm1(dlon>0) = Inf; %Exclude where data < model
                end
                ssm = ssm + distm1.^2./(bdw1(:,aslatlon(1)).^2*ones(1,nds));
            end
            
            %Loop over coordinates, adding scaled squared distance if cartesian
            for i=1:m
                if (iscartesian(i)==1)
                    if (isnan(period(i)))
                        distm1 = xp1(:,i)*ones(1,nds)-ones(ni1,1)*xs(:,i)';
                        if (one_sided(i)==-1)
                            distm1(distm1<0) = Inf; %Exclude where data > model
                        elseif (one_sided(i)==1)
                            distm1(distm1>0) = Inf; %Exclude where data < model
                        end
                    end
                    if (~isnan(period(i)))
                        distm1 = fnperiodic_distm(xp1(:,i),xs(:,i),period(i));
                    end
                    ssm = ssm + distm1.^2./(bdw1(:,i).^2*ones(1,nds));
                end
            end

            
            %Calculate kernel matrix from sum of squares matrix ssm
            if kfunc==1; rhom = sqrt(ssm); K = exp(-rhom); end
            if kfunc==2; rhom = ssm/2; K = exp(-rhom); end
            if kfunc==3; rhom = sqrt(ssm); K = sinc(rhom); end
% REDO            if kfunc==4; K = abs(ssm)<=1; end
            if kfunc==5; K = 1./(1+ssm); end
            
            
            %Correct for artifacts due to finite numerical precision
            if (kfunc==1||kfunc==2)
                minrho = min(rhom,[],2); sel = find(minrho>700); nsel = length(sel);
                %Rows where min(rho)>700 require special treatment to allow for finite numerical precision
                if (nsel>0)
                    if verbose>1; disp('Some rows of weight matrix modified to allow for finite numerical precision'); end
                    if (kfunc==1)
                        disp('Bad numerics for kfunc = 1, consider increasing bdw')
                        stop
                    end
                    if (kfunc==2)
                        for i=1:nsel
                            closest = find(ssm(sel(i),:)==min(ssm(sel(i),:)));
                            K(sel(i),:) = zeros(1,nds); K(sel(i),closest) = 1/length(closest);
                            %Usually, this means that the closest datum takes all the weight
                        end
                    end
                end
            end
            
            
            %Additional weighting for observational error
            if ~isempty(obserrv); K=K./(ones(ni1,1)*reshape(obserrv(recs),1,nds)); end

            
            %Limit to closest data if req'd (only if smoothing to one prediction point)
            if (ni1==1 && isempty(ndclosest)==0)
                Ks = sort(K,'descend');
                Kmin = Ks(min(length(Ks),ndclosest));
% REDO                selclosest = find(K>=Kmin);
                K = K(selclosest);
                recs = recs(selclosest); nds = length(recs);
                xs = xs(selclosest,:); ys = ys(selclosest);
                if isempty(xc)==0; xcs = xcs(selclosest,:); end
                if nblocks>1; xsm{iblock} = xs; end
                if isempty(x_coverage)==0; x_coverages = x_coverages(selclosest); end
            end
            
            
            %Normalize kernel so that rows sum to 1
            K = K./(sum(K,2)*ones(1,nds));
            

            %Check evenness of weight distribution
            if (check_evenness==1)
                Ks = sort(K,2,'descend'); Ks = cumsum(Ks,2);
                ndsdom = NaN*ones(ni1,1);
                for i=1:ni1
                    if (sum(~isnan(Ks(i,:)))>0)
                        ndsdom(i) = find(Ks(i,:)>domfrac,1,'first');
                        %This is the number of data comprising >= domfrac of the weight
% REDO vvv                       if (ndsdom(i)<=2 && ndsdom(i)<=0.01*nds && verbose>1)
                        if (verbose>1)
                            disp(['ks_regress: Warning: kernel dominated by only ',int2str(ndsdom(i)),' data (= ',num2str(100*ndsdom(i)/nds),'% of all data)'])
                        end
                    end
                end
            end
            
            
            
            %Calculate weight matrix for linear prediction yp = wt*ys
            if (nb==1)                        %Nadaraya-Watson, local constant regression (possibly with additional linear covariates)
                wt = K;
            end
            if (nb>1)                         %Local linear regression
                wt = NaN*ones(ni1,nds); b = NaN*ones(nb,ni1);
                for i=1:ni1
                    X1 = [ones(nds,1) NaN*ones(nds,nb-1)]; Xp1 = [1 zeros(1,nb-1)]; q = 1;
                    if (nb1>0) 
                        for j=1:m
                            for k=1:lr_order(j); X1(:,q+1) = (xs(:,j) - xp1(i,j)*ones(nds,1)).^k; q = q + 1; end
                        end
                    end
                    if (nbc>0)
                        for j=1:mc
                            for k=1:lr_orderc(j); X1(:,q+1) = (xcs(:,j) - xp1c(i,j)*ones(nds,1)).^k; q = q + 1; end
                        end
                        %if lr_orderc>0; X1(:,q+1:q+mc) = xcs - ones(nds,1)*xp1c(i,:); end
                    end
                    W1 = diag(K(i,:));
                    R1 = (X1'*W1*X1)\(X1'*W1);
                    wt(i,:) = Xp1*R1; b(:,i) = R1*ys;
                    %size(R1)
                    %i ni1]
                end
            end
            
            

            
            %If x_coverage supplied, calculate xp_coverage for these prediction points
            if (isempty(x_coverage)==0)
                xp_coverage = NaN*ones(ni1,1);
                if (coverage_type==0) %Simple coverage: range of x_coverages for data contributing to yp
                    for i=1:ni1
                        sel = find(wt(i,:)>fr_wt_coverage*mean(wt(i,:)) & ssm(i,:)<max_sdist_coverage^2);
                        if isempty(sel)==0; xp_coverage(i) = range(x_coverages(sel)); else xp_coverage(i) = 0; end
                    end 
                end
                if (coverage_type==1)
                    for i=1:ni1
                        sel = find(wt(i,:)>fr_wt_coverage*mean(wt(i,:)) & ssm(i,:)<max_sdist_coverage^2);
                        xp_coverage(i) = length(unique(floor(x_coverages(sel)))); %#ok<FNDSB>
                        %Note: this will return 0 if sel if empty
                    end
                end
            end
            
            
            
            %Calculate smoothing estimates
            yp = wt*ys;
            if nneg==1; yp = max(0,yp); end
            
            
            %Calculate uncertainties if req'd (takes more time)
            if (calcunc==1)
                %Assuming: y = m(x) + V(x)*eps
                
                %First calculate residuals by smoothing to the data
                optv = opt; optv.lr_order = vlr_order; optv.calcunc = 0;
                optv.nblocks = 1; %NOTE: To get an accurate weight matrix, cannot divide into blocks
                [yhat,out2] = ks_regress(ys,xs,xs,vbdw,optv);
                Khat = out2.K;
                rhat = ys - yhat; rhat2 = rhat.^2;
                
                if (variable_res_var==1)
                    %Estimate V(x) using the normalized weighted RSS (see rlpoly_Stata.pdf, p10)
                    wthat = Khat./(sum(Khat,2)*ones(1,nds));
                    varhat = wthat*rhat2;
                    ndshatv = sum(wthat>0,2);
                    varhat = varhat.*ndshatv./max(1,(ndshatv-1));   %varhat cannot be bias-corrected where ndshat = 1
                    %Estimate V(xp1) similarly
                    wtvarp = K./(sum(K,2)*ones(1,nds));
                    varp = wtvarp*rhat2;
                    ndsvarpv = sum(wtvarp>0,2);
                    varp = varp.*ndsvarpv./max(1,(ndsvarpv-1));
                else
                    var0hat = sum(rhat2)/(nds-1);
                    varhat = var0hat*ones(nds,1);
                    ndsvarpv = nds*ones(ni1,1);
                    varp = var0hat*ones(ni1,1);
                end
                
                
                %mp = wt*ys => C[mp] = wt*C[ys]*wt' = wt*V(x)*wt' assuming uncorrelated errors
                Vhat = diag(varhat);
                Cmp = wt*Vhat*wt';
                mp = yp; mpvar = diag(Cmp); mpse = sqrt(mpvar);
                tfac = tinv(1-alph/2,ndsvarpv-1);
                mpL = mp - tfac.*mpse; mpH = mp + tfac.*mpse;
                ypse = sqrt(mpvar + varp);
                ypL = yp - tfac.*ypse; ypH = yp + tfac.*ypse;
                if (nneg==1)
                    mpL = max(0.,mpL); ypL = max(0.,ypL);
                end
                if (nb>1) %b = R1*ys => C[b] = R1*C[ys]*R1' = R1*V(x)*R1' assuming uncorrelated errors
                    %stop
                end
            end
            
            
            if (nblocks>1)
                ypf(recp1) = yp; ndsf(recp1)=nds; 
                if check_evenness==1; ndsdomf(recp1) = ndsdom; end
                Kf{iblock} = K; wtf{iblock} = wt; recsf{iblock} = recs;
                if isempty(x_coverage)==0; xp_coveragef(recp1) = xp_coverage; end
                if nb>1; bf(:,recp1) = b; end
                if (calcunc==1)
                    rhatf{iblock} = rhat; varhatf{iblock} = varhat; varpf(recp1) = varp;
                    mpsef(recp1) = mpse; mpLf(recp1) = mpL; mpHf(recp1) = mpH;
                    ypsef(recp1) = ypse; ypLf(recp1) = ypL; ypHf(recp1) = ypH;
                end
            end

% REDO VVV           if (doplotc==1 && isubc<=nsubc)
            if (doplotc==1)
                if (mod(kc,kcstep)==0 || min(sum((xpp-ones(npp,1)*xp1(selpp,:)).^2,2))==0) %Plot vs. covariate (e.g. time)
                    Xcpp = [ones(ncpp,1) NaN*ones(ncpp,nb-1)]; q = 1;  %Covariate prediction matrix
                    if lr_order>0; Xcpp(:,q+1:q+m) = zeros(ncpp,m); q = q + m; end
                    if lr_orderc>0; Xcpp(:,q+1:q+mc) = xcpp; end
                    ypp = Xcpp*b(:,selpp); %Predictions for plotting along covariate dimension
                    selcs = find(ssm(selpp,:)<max_sdist_plotc^2);
                    
                    subplot(nrowsc,ncolsc,isubc)
                    plot(xcs(selcs,1)+xcoff,ys(selcs),'.','Color',ycol,'MarkerSize',ymarkersize);
                    hold on; plot(xcpp+xcoff,ypp,'Color',ypcol,'LineStyle','-','LineWidth',linewidth)
                    xlabel(xstr,'FontSize',lfontsize); ylabel(ystr,'FontSize',lfontsize);
                    title(['x = ',num2str(xp1(selpp,1)),', y = ',num2str(xp1(selpp,2))],'FontSize',lfontsize)
                    if isempty(xlimc)==0; set(gca,'XLim',xlimc); end
                    if isempty(ylimc)==0; set(gca,'YLim',ylimc); end
                    isubc = isubc + 1;
                end
            end
            kc = kc + 1;
        end
    end%ni1>0
    if (mod(iblock,1e4)==0 && verbose>0); disp(['Done ',int2str(iblock),' of ',int2str(nblocks),' blocks']); end
end%Master loop over nblocks-----------------------------------------------

if (nblocks>1)
    yp = ypf; nds = ndsf; 
    if check_evenness==1; ndsdom = ndsdomf; end
    if isempty(x_coverage)==0; xp_coverage = xp_coveragef; end
    if nb>1; b = bf; end
    if (calcunc==1)
        rhat = rhatf; varp = varpf;
        mpse = mpsef; mpL = mpLf; mpH = mpHf;
        ypse = ypsef; ypL = ypLf; ypH = ypHf;
    end
else
    nds = nds*ones(ni,1);
end





%Plot if req'd
if (doplot==1)
    figure(ifig)
    subplot(nrows,ncols,isub)
    
    if (m==1)
        if (plottype==0)
            plot(x,y,'k.',xp,yp,'r-')
            if calcunc==1; hold on;plot(xp,yp-yppe,'r--',xp,yp-yppe,'r--'); end
            xlabel('x'); ylabel('y');
        end
        if (plottype==1)
            plot(y,x,'k.',yp,xp,'r-')
            if calcunc==1; hold on;plot(yp-yppe,xp,'r--',yp+yppe,xp,'r--'); end
            if (length(yp)==1)
                hold on;plot(yp,xp,'ro');
                if calcunc==1; hold on; errorbar_x(yp,xp,yppe,yppe,'o',struct('color','r')); end
            end
            xlabel('d'); ylabel('z')
            axis ij
            if isempty(tstr)==0; title(tstr,'FontSize',10); end
        end
    end
end





%Store auxiliary results in output structure
if (nargout>1 && max(nds)>0)
    out.nblocks = nblocks;
    if (nblocks>1) 
        out.xs = xsm; if use_rad==1; out.rad = radm; end
    else
        out.xs = xs; if use_rad==1; out.rad = rad1; end
    end
    out.K = K; out.nds = nds; out.wt = wt; out.recs = recs;
    if check_evenness==1; out.ndsdom = ndsdom; end
    if isempty(x_coverage)==0; out.xp_coverage = xp_coverage; end
    if nb>1; out.b = b; end
    if (calcunc==1)
        out.rhat = rhat; out.varhat = varhat; out.varp = varp;
        out.mpse = mpse; out.mpL = mpL; out.mpH = mpH;  %Confidence intervals
        out.ypse = ypse; out.ypL = ypL; out.ypH = ypH;  %Prediction intervals
    end
else
    out = [];
end






function [A1c,A2c,A3c,A4c] = conformsizes(nm,A1,A2,A3,A4)

%function [A1c,A2c,A3c,A4c] = conformsizes(nm,A1,A2,A3,A4)
%Conform up to 4 array sizes to [nm(1)*nm(2)] by repetition
%If nm is empty it is taken from size(A1)
%%Phil Wallhead 26/04/2014

if isempty(nm); nm = size(A1); end
n = nm(1); m = nm(2);

A1c = A1;
if length(A1)==1; A1c = A1*ones(n,m); end
if (length(A1(:,1))>1 && length(A1(1,:))==1); A1c = A1*ones(1,m); end
if (length(A1(:,1))==1 && length(A1(1,:))>1); A1c = ones(n,1)*A1; end

if (nargin>2)
    A2c = A2;
    if length(A2)==1; A2c = A2*ones(n,m); end
    if (length(A2(:,1))>1 && length(A2(1,:))==1); A2c = A2*ones(1,m); end
    if (length(A2(:,1))==1 && length(A2(1,:))>1); A2c = ones(n,1)*A2; end
end

if (nargin>3)
    A3c = A3;
    if length(A3)==1; A3c = A3*ones(n,m); end
    if (length(A3(:,1))>1 && length(A3(1,:))==1); A3c = A3*ones(1,m); end
    if (length(A3(:,1))==1 && length(A3(1,:))>1); A3c = ones(n,1)*A3; end
end

if (nargin>4)
    A4c = A4;
    if length(A4)==1; A4c = A4*ones(n,m); end
    if (length(A4(:,1))>1 && length(A4(1,:))==1); A4c = A4*ones(1,m); end
    if (length(A4(:,1))==1 && length(A4(1,:))>1); A4c = ones(n,1)*A4; end
end
