

function [J,out] = J_fit_partially_linear_model(theta,opt)


%function [J,out] = J_fit_partially_linear_model(theta,opt)
%
%Calculates J = -2log(likelihood) for partially linear model:
%
%           Y = X(theta)*b + eps
%
%where the linear matrix of covariates is a nonlinear function of nonlinear
%parameters theta, and is specified by an input function opt.Xfun.
%Normal errors eps are assumed.
%
%%Phil Wallhead 03/09/2021


if nargin==1; opt = []; end
if isfield(opt,'Y')==1; Y = opt.Y; else Y = []; end
if isfield(opt,'Xfun')==1; Xfun = opt.Xfun; else Xfun = []; end
if isfield(opt,'nonnegb')==1; nonnegb = opt.nonnegb; else nonnegb = 0; end
if isfield(opt,'Vs')==1; Vs = opt.Vs; else Vs = []; end %Input Vs to assume V_i = V * Vs_i
%if isfield(opt,'sumlogVs')==1; sumlogVs = opt.sumlogVs; else sumlogVs = []; end %Supplying precalculated sum(log(Vs)) doesn't seem to improve speed

n = length(Y);
X = feval(Xfun,theta);
if ~isempty(Vs)
    R2 = diag(sqrt(1./Vs)); Xs = R2*X; Ys = R2*Y;
    if nonnegb==1; b = lsqnonneg(Xs,Ys); else b = Xs\Ys; end
    Ym = X*b;
    VhatMLE = sum(((Y-Ym).^2)./Vs)/n;
    %if isempty(sumlogVs); sumlogVs = sum(log(Vs)); end
    sumlogVs = sum(log(Vs)); %Special case of log(det(C)), where C is the (scaled) covariance matrix (see e.g. Jccomplex.m)
    
    J = n + sumlogVs + n*log(VhatMLE); %Note that if Vs_i were a constant the sumlogVs would cancel with n*log(1/Vs)
else
    if nonnegb==1; b = lsqnonneg(X,Y); else b = X\Y; end
    Ym = X*b;
    VhatMLE = sum((Y-Ym).^2)/n;
    
    J = n + n*log(VhatMLE);
end


out = [];
if (nargout>1)
    out.X = X; out.b = b; out.Ym = Ym; out.VhatMLE = VhatMLE; out.SSQ = n*VhatMLE;
    out.Vhat = VhatMLE*n/(n-length(b));
end

