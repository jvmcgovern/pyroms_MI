

function [thetabest,out] = optimize_theta(thetai,thetamin,thetamax,Jfun,nr,opt)
      
%function [thetabest,out] = optimize_theta(thetai,thetamin,thetamax,Jfun,nr,opt)
%
%Minimize Jfun with respect to search parameters theta
%
%Uses:  npermutek (optional) for uniform exhaustive sampling
%       hessianest (optional) to calculate Hessian matrix
%
%%Phil Wallhead 23/12/2017


nparsh = length(thetamin);
sztheta = size(thetai);
thetai = thetai(:); thetamin = thetamin(:); thetamax = thetamax(:);

if nargin<5; nr = 1; end
if nargin<6; opt = []; end
if isfield(opt,'optJ')==1; optJ = opt.optJ; else optJ = []; end
if isfield(opt,'thetastart')==1; thetastart = opt.thetastart(:); else thetastart = []; end
if isfield(opt,'randomstart')==1; randomstart = opt.randomstart; else randomstart = 0; end
    %1 to randomize the start even for nr = 1
if isfield(opt,'fnoise')==1; fnoise = opt.fnoise; else fnoise = []; end
if isfield(opt,'maxfunevals')==1; maxfunevals = opt.maxfunevals; else maxfunevals = 2e4; end
if isfield(opt,'maxiter')==1; maxiter = opt.maxiter; else maxiter = 1e4; end    
if isfield(opt,'TolFun')==1; TolFun = opt.TolFun; else TolFun = 1e-6; end
if isfield(opt,'TolCon')==1; TolCon = opt.TolCon; else TolCon = 1e-6; end
if isfield(opt,'TolX')==1; TolX = opt.TolX; else TolX = 1e-6; end
if isfield(opt,'algorithm')==1; algorithm = opt.algorithm; else algorithm = 0; end
    %0 for fmincon/fminbnd (default)
    %1 to use snobfit.m
    %Could set default to snobfit.m if you don't have the Optimisation toolbox 
if isfield(opt,'usefminbnd')==1; usefminbnd = opt.usefminbnd; else usefminbnd = 1; end
    %1 (def) to force use of fminbnd for 1D problems (seems to be most efficient)
    %This does however require another routine to calculate the hessian if req'd
    if (nparsh==1 && usefminbnd==1); algorithm = 0; end    
if isfield(opt,'opta')==1; opta = opt.opta; else opta = []; end     
    %Additional options structure to be supplied to optimization algorithm
if isfield(opt,'Hcalc')==1; Hcalc = opt.Hcalc; else Hcalc = 0; end    
    %0 to switch off the hessian calculation (def)
    %1 to use fmincon (gives similar estimate to hessianest in most cases)
    %1.5 to use fminunc (gives similar estimate to hessianest in most cases)
    %2 to use hessianest from DERIVEST suite by John D'Errico (def~=0)
if isfield(opt,'divide_and_conquer')==1; divide_and_conquer = opt.divide_and_conquer; else divide_and_conquer = []; end
    %1 to divide the allowed range according to the starting value sampling design 
if isfield(opt,'rdesign')==1; rdesign = opt.rdesign; else rdesign = []; end
    %0 for uniform random sampling
    %1 for uniform exhaustive sampling
    %2 for uniform Latin-Hypercube sampling
    if (isempty(rdesign)==1)
        if nparsh==1; rdesign = 1; end    %Regular sampling default for 1D
        if nparsh>1; rdesign = 2; end     %LHS sampling default for (>1)D
    end
    if (isempty(divide_and_conquer))
        divide_and_conquer = 0;    
        if (nparsh==1 && isfield(opt,'rdesign')==0); divide_and_conquer = 1; end %By default, divide and conquer only for 1D case
    end
    if divide_and_conquer==1; rdesign = 1; end %Regular sampling for divide_and_conquer strategy
if isfield(opt,'alph')==1; alph = opt.alph; else alph = 0.05; end
if isfield(opt,'ndof')==1; ndof = opt.ndof; else ndof = Inf; end
if isfield(opt,'verbose')==1; verbose = opt.verbose; else verbose = 0; end
    


    
options = optimset('TolFun',TolFun,'TolCon',TolCon,'TolX',TolX,...
    'MaxFunEvals',maxfunevals,'MaxIter',maxiter);
%clear LASTN; LASTN = maxNumCompThreads(1); %#ok<NASGU>
if (algorithm==1)   %For snobfit, not clear how to make repeat restarts 
    nr = 1; if Hcalc==1; Hcalc = 1.5; end
end
if (algorithm==2)
    opta.Verbosity = verbose;
end


%Set design for random starting values
nr1 = nr;
if (nr>1 || randomstart==1)
    if rdesign==0; Xr = rand(nparsh,nr); end    %Uniform random sampling
    if (rdesign==1)                             %Regular sampling
        nr1 = nr.^(1/nparsh);
        if (abs(nr1-round(nr1))>1e-3)
            nr1 = floor(nr1); nr = nr1^nparsh;
            if verbose>0; disp(['Rounding nr1 down to nearest (1/nparsh)th root, corrected nr = ',int2str(nr)]); end
        end
        Xr = npermutek((0.5:nr1-0.5)/nr1,nparsh); Xr = Xr';
    end
    if rdesign==2; Xr = lhsdesign(nr,nparsh); Xr = Xr'; end
%    thetarange = thetamax - thetamin;
end
thetarange = thetamax - thetamin;





thetar = zeros(nparsh,nr); Jminr = zeros(1,nr); nfunevalsr = zeros(1,nr);
Jminbest = 1e8; thetabest = thetai; Hbest = NaN*ones(nparsh);
for i=1:nr%Loop over restarts----------------------------------------------
    
    if (nr>1 && isempty(thetastart)==1) %Use sampling design matrix
        thetai1 = thetamin + thetarange.*Xr(:,i);
    elseif (isempty(thetastart)==0)  %Use given first-guess starting value and add random noise
        thetai1 = thetastart.*(1 + fnoise*randn(nparsh,1));
        thetai1 = min(thetamax,thetai1); thetai1 = max(thetamin,thetai1);
    else
        thetai1 = thetai;
    end

    thetamin1 = thetamin; thetamax1 = thetamax;
    if (nr>1 && divide_and_conquer==1)
        dtheta = thetarange/(2*nr1);
        thetamin1 = thetai1 - dtheta; thetamax1 = thetai1 + dtheta;
        if verbose>1; disp('[thetamin1 thetamax1] = '); disp([thetamin1 thetamax1]); end
    end
    
    if (algorithm==0)
        if (nparsh>1 || usefminbnd==0)
            if (Hcalc==1)
                [theta,Jmin,exitflag,outa,lambda,grad,H] = fmincon(@(x)Jfun(x,optJ),...
                    thetai1,[],[],[],[],thetamin1,thetamax1,[],options); %#ok<ASGLU>
            else
                [theta,Jmin,exitflag,outa] = fmincon(@(x)Jfun(x,optJ),...
                    thetai1,[],[],[],[],thetamin1,thetamax1,[],options); %#ok<ASGLU>
            end
        else
            %fminbnd is more efficient for 1D problems
            [theta,Jmin,exitflag,outa] = fminbnd(@(x)Jfun(x,optJ),thetamin1,thetamax1);     %#ok<ASGLU>
        end
    end
    if (algorithm==1)
        [theta,Jmin,outa] = fnsnobfit(@(x)Jfun(x,optJ),[],[],thetamin1,thetamax1);
        theta = theta(:);
    end
    if (algorithm==2)
        opta.Generator = @(x) (thetamin1(:)' + rand(1,nparsh).*(thetamax1(:)'-thetamin1(:)'));
        [theta,Jmin,outa] = anneal(@(x)Jfun(x,optJ),thetai1(:)',opta);
        theta = theta(:);
    end
    if (algorithm==2.1)
        [theta,Jmin,outa] = sim_anl(@(x)Jfun(x,optJ),thetai1,thetamin1,thetamax1,TolFun);
    end
    
    if (verbose>0)
        vL = abs(theta-thetamin)./(thetamax-thetamin);
        vH = abs(theta-thetamax)./(thetamax-thetamin);
        if (min([vL; vH])<1e-3); disp('Warning: Solution near boundary'); end
    end
    
    if (nr>1)
        if (verbose>1)
            [i nr] %#ok<NOPRT>
            if (i>1)
                [Jmin Jminbest] %#ok<NOPRT>
                if (Jminbest-Jmin>0.1 && divide_and_conquer==0)
                    disp('Warning!!! LOCAL MINIMUM INFERRED ----------------------------')
                end
            else
                Jmin %#ok<NOPRT>
            end
            theta(1:min(9,nparsh))' %#ok<NOPRT>
        end
        thetar(:,i) = theta; Jminr(i) = Jmin; nfunevalsr(i) = outa.funcCount;
    end
    if (Jmin<Jminbest || i==1)
        Jminbest = Jmin; thetabest = theta; nfunevalsbest = outa.funcCount;
        if ((usefminbnd==0||nparsh>1) && Hcalc==1); Hbest = H; end
    end
end%nr---------------------------------------------------------------------





if (nr>1 && verbose>0 && divide_and_conquer==0)
    if (max(Jminr)-Jminbest>0.1)
        disp('Warning!!! LOCAL MINIMA WERE INFERRED ----------------------------')
        pause(0.5)
    end
end



if (Hcalc>0)
    vL = abs(thetabest-thetamin)./(thetamax-thetamin);
    vH = abs(thetabest-thetamax)./(thetamax-thetamin);
    
    if (min([vL; vH])<1e-3)
        if verbose>0; disp('Warning: Solution near boundary, Hessian-based SEs invalid (returning NaN)'); end
        Hbest = NaN*ones(nparsh);
        
    else
        
        if (Hcalc==1 && nparsh>1 && verbose>0); disp('Using fmincon to estimate Hessian'); end
        if (Hcalc==1.5 || (Hcalc==1 && nparsh==1 && usefminbnd==1))
            %Compute the hessian using fminunc started at thetabest
            if (Hcalc==1 && verbose>0)
                disp('Cannot use fmincon started at thetabest to estimate Hessian - it quits to early');
            end
            if verbose>0; disp('Using fminunc started at thetabest to estimate Hessian'); end
            [thetacheck,Jmincheck,exitflag,outa,grad,Hbest] = fminunc(@(x)Jfun(x,optJ),thetabest,options); %#ok<ASGLU>
            if (max(abs(thetacheck-thetabest))>2*1e-4 || abs(Jmincheck-Jminbest)>1e-3)
                if verbose>0; disp('fminunc gives different fit to fmincon/fminbnd: something strange happened!'); end
                [thetacheck thetabest] %#ok<NOPRT>
                [Jmincheck Jminbest] %#ok<NOPRT>
                stop
            end
        end
        if (Hcalc==2) %Use hessianest.m from John D'Errico's DERIVEST suite
            if verbose>0; disp('Using hessianest to estimate Hessian'); end
            Hbest = hessianest(@(x)Jfun(x,optJ),thetabest,thetamin,thetamax);
        end
    end
    
    Ctheta = (0.5*Hbest)\eye(nparsh);
    thetase = sqrt(diag(Ctheta));
    fac = tinv(1-alph/2,ndof);
    if (sztheta(1)>1) 
        thetase = thetase(:);
        thetaCIs = [thetabest-fac*thetase thetabest+fac*thetase];
    else
        thetase = thetase(:)'; 
        thetaCIs = [thetabest'-fac*thetase; thetabest'+fac*thetase];
    end
end


%Return thetabest same size as thetai
if sztheta(1)>1; thetabest = thetabest(:); else thetabest = thetabest(:)'; end


%Return other outputs
if (nargout>1)
    out = [];
    out.thetabest = thetabest; out.Jminbest = Jminbest; out.nfunevalsbest = nfunevalsbest;
    if Hcalc>0; out.Hbest = Hbest; out.Ctheta = Ctheta; out.thetase = thetase; out.thetaCIs = thetaCIs; end
    if (nr>1)
        out.thetar = thetar; out.Jminr = Jminr; out.nfunevalsr = nfunevalsr;
    end
end
