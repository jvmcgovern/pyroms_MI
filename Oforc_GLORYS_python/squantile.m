

function Y = squantile(X,P)

%function Y = squantile(X,P)
%
%Simple quantile function that gives results consistent with quantile.m.
%For convenience when quantile.m from stats toolbox is lacking.
%NOTE: X must be a vector for squantile.m (can be a matrix for quantile.m).
%
%X = vector or matrix of data [n1*n2]
%P = vector of quantiles [np*1] or [1*np].
%Y = output quantiles [np*n2], or [1*np] if X is vector and P is [1*np]    (same behaviour as quantile.m)
%
%Uses: Xplin.m
%
%%Phil Wallhead 26/09/2023

%Example tests of squantile.m vs. quantile.m:
%x = rand(1000,10); qsel = 0:0.01:1; y = squantile(x,qsel); yc = quantile(x,qsel); max(abs(y(:)-yc(:))) %<1e-15
%x = rand(1000,1); qsel = 0:0.01:1; y = squantile(x,qsel); yc = quantile(x,qsel); max(abs(y(:)-yc(:))) %<1e-15
%x = rand(1000,10); qsel = (0:0.01:1)'; y = squantile(x,qsel); yc = quantile(x,qsel); max(abs(y(:)-yc(:))) %<1e-15

[np1,np2] = size(P);
if (np1>1 && np2>1); error('Second input must be a scalar or a non-empty vector'); end
P = P(:); np = length(P);

[n1,n2] = size(X);
if n1==1; X = X(:); n2c = 1; else n2c = n2; end

Y = NaN*ones(np,n2c);
for i=1:n2c
    X1 = X(:,i);
    Xs1 = sort(X1(~isnan(X1))); ns1 = length(Xs1);
    if (ns1>0)
        Pn = (0.5:1:ns1-0.5)'/ns1; %Nodal quantiles corresponding to Xs
        Y(:,i) = Xplin(P,Pn,1)*Xs1; %Linear interpolation with constant extrapolation
    end
end

if (n2c==1 && np2>1); Y = Y'; end
