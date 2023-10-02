

function distm = fn_spherical_distm(lat1,lon1,lat2,lon2,radsin,R)

%function distm = fn_spherical_distm(lat1,lon1,lat2,lon2,radsin=0,R=6378.137km)
%
%Calculates matrix distm of spherical distances in km given vectors (lat1, lon1) of length n1 and
%vectors (lat2, lon2) of length n2. Input of degrees is assumed unless radsin = 1.
%
%%Phil Wallhead 16/03/2015

lat1 = lat1(:); lon1 = lon1(:); lat2 = lat2(:); lon2 = lon2(:);
n1 = length(lat1); n2 = length(lat2);
if nargin<5; radsin = 0; end
if nargin<6; R = 6378.137; end %Default radius of the Earth in km

if (radsin==1)
    lat1r = lat1; lon1r = lon1;
    lat2r = lat2; lon2r = lon2;
else
    pi180 = pi/180;
    lat1r = lat1*pi180; lon1r = lon1*pi180;
    lat2r = lat2*pi180; lon2r = lon2*pi180;
end

coslat1rcoslat2r = cos(lat1r)*cos(lat2r)';
dlatr = lat1r*ones(1,n2)-ones(n1,1)*lat2r';
dlonr = lon1r*ones(1,n2)-ones(n1,1)*lon2r';
distm = R*2*asin(sqrt(sin(dlatr/2).^2+coslat1rcoslat2r.*sin(dlonr/2).^2));