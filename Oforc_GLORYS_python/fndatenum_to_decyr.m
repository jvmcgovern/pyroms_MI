

function decyr = fndatenum_to_decyr(td,distort,nyrdays0,year0)

%function decyr = fndatenum_to_decyr(td,distort=0,nyrdays0=365.25,year0=2000)
%
%Converts datenumber (from datenum.m) to decimal years
%If distort = 1, time is distorted to keep exact agreement between year and floor(decyr)
%If distort = 0 (def), time is not distorted but a fixed no. of days in the year is
%assumed (nyrdays0) and years are computed with respect to a reference year year0
%
%See also fndecyr_to_datenum.m
%
%%Phil Wallhead 13/04/2017

if nargin<2; distort = 0; end
if nargin<3; nyrdays0 = 365.25; end
if nargin<4; year0 = 2000; end

[year,~,~] = datevec(td);
td0 = datenum(year,1,1);

if (distort==1)
    nyrdays = datenum(year+1,1,1) - td0;
    decyr = year + (td-td0)./nyrdays;
    %Time is distorted to guarantee that floor(decyr) = year
else
    decyr = year0 + (td-datenum(year0,1,1))/nyrdays0;
    %Here time is not distorted, but floor(decyr) may not always equal year
end