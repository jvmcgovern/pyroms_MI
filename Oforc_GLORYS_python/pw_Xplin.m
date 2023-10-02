
function X = Xplin(t,tn,extrapconst)

%function X = Xplin(t,tn,extrapconst=[0 0])
%
%Builds 1D linear interpolation matrix from tn to t.
%Nodes tn MUST all not-NaN and in increasing order.
%Input extrapconst = [1 1] to extrapolate constant beyond [first last] nodes
%(otherwise defaults to linear extrapolation, constant if only one node).
%Any missing target points (t=NaN) will result in NaN rows in X, unless
%there is only one node point, in which case a constant is extrapolated.
%
%Example:
% x = 1:9; y = rand(9,1); xp = 0:0.1:10; 
% tic;Xz0 = Xplin(xp,x,[0 1]); Xz1 = Xplin(xp,x,[1 0]); toc
% figure; plot(x,y,'ko',xp,Xz0*y,'k-',xp,Xz1*y,'r-')
%
%%Phil Wallhead 23/06/2021

t = t(:);
tn = tn(:);

if (sum(isnan(tn)) > 0)
    error('One or more node values tn = NaN in Xplin.m');
% elseif (min(diff(tn)) <= 0)
elseif le(min(diff(tn)), 0)
    error('Nodes tn MUST be in increasing order');
end

if nargin<3;
    extrapconst=[0 0];
end

if length(extrapconst)==1;
    extrapconst = extrapconst*ones(1,2);
end

m = length(tn);

if (m==1)
    
    X = ones(length(t),1);
    
elseif (m>1)
    
    X = NaN*ones(length(t),m);
    sel = find(~isnan(t));
    X(sel,:) = 0;

    for i=1:length(sel)
%        if (t(sel(i))<=tn(1))
        if le(t(sel(i)), tn(1))
            if (extrapconst(1)==1)
                X(sel(i),1) = 1.;
            else
                X(sel(i),1) = 1 - (t(sel(i))-tn(1))/(tn(2)-tn(1));
                X(sel(i),2) = 1 - X(sel(i),1);
            end
        elseif (t(sel(i))>tn(m))
            if (extrapconst(2)==1)
                X(sel(i),m) = 1.;
            else
                X(sel(i),m-1) = 1 - (t(sel(i))-tn(m-1))/(tn(m)-tn(m-1));
                X(sel(i),m) = 1 - X(sel(i),m-1);
            end
        else
            n1=find(tn<t(sel(i)),1,'last');
            X(sel(i),n1) = 1 - (t(sel(i))-tn(n1))/(tn(n1+1)-tn(n1));
            X(sel(i),n1+1) = 1 - X(sel(i),n1);
        end
    end
    
else
    X = [];
end
