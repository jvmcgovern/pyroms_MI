function distm = fnperiodic_distm(xp,x,period)

%function distm = fnperiodic_distm(xp,x,period)
%
%Calculate the distance matrix between xp (rows) and x (columns) over a periodic dimension.
%Works for lower-limit integers (e.g. days since 1st January, yrday = 0-364, period = 365),
%upper-limit integers (e.g. month of the year, month = 1-12, period = 12),
%and continuous variables (e.g. hours since start of year, period = 12*365).
%
%%Phil Wallhead 10/02/2019
%
%Example 1: distance from yrday = 364, where yrday = days since 1st January
%distm = fnperiodic_distm(364,0:364,365) %Assuming a non-leap year
%
%Example 2: distance from month = 12 (December), where month = month of the year (1-12)
%distm = fnperiodic_distm(12,1:12,12)
%
%Example 3: distance from hr = 364*12+11, where hr = hours since start of year
%distm = fnperiodic_distm(364*12+11,0:(365*12-1),365*12) %Assuming a non-leap year

n = length(x);
np = length(xp);

M1 = mod(xp(:)*ones(1,n),period);

M2 = mod(ones(np,1)*x(:)',period);
xmax = max(M1,M2);
xmin = min(M1,M2);
distm = min(xmax-xmin,xmin+period-xmax);

end