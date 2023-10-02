function [Q0,out] = calc_surfaceSWI(lat,lon,tday,opt)

%function [Q0,out] = calc_surfaceSWI(lat,lon,tday,opt)
%
%Calculate short-wave irradiance above sea surface [W/m2] as instantaneous (Q0)
%and daily-average values (out.Q0dav, if opt.calc_dav=1), given:
%lat = latitude [degrees]
%lon = longitude [degrees] (needed for hour angle calculation)
%tday = days since 1st January (0-364 or 0-365 for leap year, real-valued, UTC)
%opt.JD = Julian day (if needed, see below)
%
%%Phil Wallhead 29/01/2016

lat = lat(:)';
lon = lon(:)';
ns = length(lat);
tday = tday(:);
nt = length(tday);
if nargin<3; opt = []; end
if isfield(opt,'use_tday_lag')==1;
    use_tday_lag = opt.use_tday_lag;
else
    use_tday_lag = 0;
end
if isfield(opt,'year')==1;
    year = opt.year;
else
    year = [];
end
if isfield(opt,'f_model')==1;
    f_model = opt.f_model;
else
    f_model = 1;
end
if isfield(opt,'decl_model')==1;
    decl_model = opt.decl_model;
else
    decl_model = 1;
end
if isfield(opt,'use_eqtime')==1;
    use_eqtime = opt.use_eqtime;
else
    use_eqtime = 1;
end
if isfield(opt,'Q_model')==1;
    Q_model = opt.Q_model;
else
    Q_model = 0;
end
if isfield(opt,'cloud_model')==1;
    cloud_model = opt.cloud_model;
else
    cloud_model = 0;
end
if isfield(opt,'calc_hav')==1;
    calc_hav = opt.calc_hav;
else
    calc_hav = 0;
end
if isfield(opt,'calc_dav')==1;
    calc_dav = opt.calc_dav;
else
    calc_dav = 0;
end
if isfield(opt,'calc_Qs')==1;
    calc_Qs = opt.calc_Qs;
else
    calc_Qs = 0;
end
if isfield(opt,'JD')==1;
    JD = opt.JD;
else
    JD = [];
end
if isfield(opt,'e0')==1;
    e0 = opt.e0;
else
    e0 = zeros(nt,ns);
end %Vapour pressure (hPa=mbar)
if isfield(opt,'C')==1;
    C = opt.C;
else
    C = zeros(nt,ns);
end %Cloud cover (oktas, 0-1)
if isfield(opt,'ndays_year')==1;
    ndays_year = opt.ndays_year;
else
    ndays_year = 365.2422;
end
if isfield(opt,'ROMS_model')==1;
    ROMS_model = opt.ROMS_model;
else
    ROMS_model = 0;
end
if (ROMS_model==1)
    f_model = 0; decl_model = 0.1; use_eqtime = 0; Q_model = 0.1; cloud_model = 0.1;
    ndays_year = 365.2425;
end
if nargout>1;
    out = struct('f_model',f_model,'decl_model',decl_model,'Q_model',Q_model,'cloud_model',cloud_model,'ndays_year',ndays_year);
end
%Stepping parameters for computing daily averages (where analytical integration not possible).
nstep = 24;
dT = 1/nstep;

S0 = 1367;     %Solar constant (W/m^2)

if (use_tday_lag==1 && ~isempty(year))
    tdayc = tday - mod(year,4)/4;
else
    tdayc = tday;
end

%%Calculate the Sun-Earth distance factor (f)
if (f_model==0) %Assume constant Sun-Earth distance
    f = ones(nt,1);

elseif (f_model==1) %Use Liou et al (2002) eqn 2.2.9
               %Originally from Spencer (1971), see Brun and Pinardi (2007) Table A1
    an = [1.000110 0.034221 0.000719]';
    bn = [0.000000 0.001280 0.000077]';
    t = 2*pi*tdayc/ndays_year; %ndays_year=365 in original formula
    f = 0*t;
    for i=1:3
        n = i-1;
        f = f + an(i)*cos(n*t) + bn(i)*sin(n*t);
    end

elseif (f_model==2)
    g = (pi/180)*(357.528 + 0.9856003*(JD-2451545.0));
    f = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g); %Michalsky (1988) formula, see Brun and Pinardi (2007) Table A1

else
    error(['f_model = ',num2str(f_model),' not recognized'])
end

%%Calculate the declination (solar inclination) in radians
if (decl_model==0) %Use Brock (1981) eqn 1
    decl = 23.45*pi/180*sin(2*pi*(284+tdayc)/ndays_year); %ndays_year=365 in original formula

elseif (decl_model==0.1) %Use ROMS approx. (ana_srflx.h)
    yday = tdayc + 1; %ROMS uses yday ranging over 1-365 or 1-366, see function yearday in Utility/dateclock.F (yday = floor((275*1/9)) - 2*floor((1+9)/12) + 1 - 30)
    decl = 23.44*pi/180*cos((172.0-yday)*2.0*pi/ndays_year); %ndays_year=365.2425 in ROMS code

elseif (decl_model==1) %Use Liou et al (2002) eqn 2.2.10
    cn = [0.006918 -0.399912 -0.006758 -0.002697]';
    dn = [0.000000 0.0702570 0.0009070 0.0001480]';
    t = 2*pi*tdayc/ndays_year; %ndays_year=365 in original formula
    decl = 0*t;
    for i=1:4
        n = i-1;
        decl = decl + cn(i)*cos(n*t) + dn(i)*sin(n*t);
    end

else
    error(['decl_model = ',num2str(decl_model),' not recognized'])
end

latc = lat*pi/180;  %latitude in radians
H = acos(min(1,max(-1,-tan(decl).*tan(latc))));       %"Half-day" (see Liou et al., eqn 2.2.2)
%(see Liou02_Int_Geophys_Chapter2, eqn 2.2.2)
%Note: capping accounts for zero daylength during high-latitude winter
daylength = 24*H/pi;    %daylength in hours
hour = (tday - floor(tday))*24;  %Hours since 00:00 (real-valued)
fracyear = 2*pi*tdayc/ndays_year;

eqtime = 229.18*(0.000075+0.001868*cos(fracyear)-0.032077*sin(fracyear)-...
    0.014615*cos(2*fracyear)-0.040849*sin(2*fracyear));  %Equation of time correction [minutes]
%This is from NOAA formulae, see sunrisesunset.m, sunrisesunset2.m.

if (use_eqtime==1)
    h = (12.0-hour*ones(1,ns))*pi/12.0 - (ones(nt,1)*lon)*pi/180 - eqtime*pi/(60*12); %hour angle in radians (0-2*pi).
else
    h = (12.0-hour*ones(1,ns))*pi/12.0 - (ones(nt,1)*lon)*pi/180; %hour angle in radians (0-2*pi).
    %This is the 'negative' form used in ROMS, which is equivalent to positive form e.g. "15deg*(LST-12)", since cos(h) = cos(-h).
end
%Note: h is an [nt*ns] array.

if (calc_hav==1)
    dhv = (-29.5:1:29.5)/60 * -pi/12; %Increments in hour angle for +/- 30 minutes about prediction times
    nsteph = length(dhv);
    hm = h*ones(1,nsteph) + dhv;
end

 % Use SMS formula (Rosati and Miyakoda, 1988) as summarized in Byun and Pinardi (2007) Table 1
if (Q_model==0)
    tau = 0.7;
    A_a = 0.09;
    cos_theta = sin(decl*ones(1,ns)).*sin(ones(nt,1)*latc) + cos(ones(nt,1)*latc).*cos(decl*ones(1,ns)).*cos(h);
    SE = max(0., S0*(f*ones(1,ns)).*cos_theta);

    % This is needed to deal with potential zeros of cos_theta, which will produce NaN in QDir etc.
    rec = find(SE(:)>0);
    I1 = zeros(nt,ns); transfac = I1; QDir = I1; QDiff = I1;
    transfac(rec) = tau.^(1./cos_theta(rec));
    QDir(rec) = SE(rec).*transfac(rec);
    QDiff(rec) = 0.5*SE(rec).*(1-A_a-transfac(rec));
    Q0 = QDir + QDiff;

    % Calculate hourly averages centred on the prediction hours
    if (calc_hav==1)
        Q0hav = I1; transfacm = zeros(nt,nsteph); QDirm = zeros(nt,nsteph); QDiffm = zeros(nt,nsteph);
        for i=1:ns
            cos_thetam = (sin(decl)*sin(latc(i)))*ones(1,nsteph) + ...
                cos(latc(i))*(cos(decl)*ones(1,nsteph)).*cos(hm);
            SEm = max(0., S0*(f*ones(1,nsteph)).*cos_thetam);
            % This is needed to deal with potential zeros of cos_theta, which will produce NaN in QDir etc.
            rec = find(SEm(:)>0);
            transfacm(rec) = tau.^(1./cos_thetam(rec));
            QDirm(rec) = SEm(rec).*transfacm(rec);
            QDiffm(rec) = 0.5*SEm(rec).*(1-A_a-transfacm(rec));
            Q0m = QDirm + QDiffm;
            Q0hav(:,i) = mean(Q0m,2);
        end
    end
    if (calc_dav==1)
        Q0dav = I1;
        for i=1:ns
            % This is ~=0, so we do not need the subselection above.
            cos_thetam = (sin(decl)*sin(latc(i)))*ones(1,nstep) + ...
                (cos(latc(i))*cos(decl))*cos(2*pi*(dT/2:dT:1-dT/2));
            SEm = max(0., S0*(f*ones(1,nstep)).*cos_thetam);
            transfacm = tau.^(1./cos_thetam);
            QDirm = SEm.*transfacm;
            QDiffm = 0.5*SEm.*(1-A_a-transfacm);
            Q0m = QDirm + QDiffm;
            Q0dav(:,i) = mean(Q0m,2);
        end
    end

elseif (Q_model==0.1)
    % = 'zenith' in the ROMS code
    cos_theta = max(0, sin(decl*ones(1,ns)).*sin(ones(nt,1)*latc) + cos(ones(nt,1)*latc).*cos(decl*ones(1,ns)).*cos(h));
    Q0 = S0*(f*ones(1,ns)) .* cos_theta.^2 ./ (1.085*cos_theta + e0.*(2.7+cos_theta)*1e-3 + 0.1);
    % Zillman (1972) approx. (Eqn. 6 in Niemela et al., 2001) = ROMS formula before correcting for cloud cover.

    % Calculate hourly averages centred on the prediction hours
    if (calc_hav==1)
        Q0hav = zeros(nt,ns);
        for i=1:ns
            cos_thetam = max(0, (sin(decl)*sin(latc(i)))*ones(1,nsteph) + ...
                cos(latc(i))*(cos(decl)*ones(1,nsteph)).*cos(hm));
            Q0m = S0*(f*ones(1,nstep)) .* cos_thetam.^2 ./ (1.085*cos_thetam + (e0(:,i)*ones(1,nstep)).*(2.7+cos_thetam)*1e-3 + 0.1);
            Q0hav(:,i) = mean(Q0m,2);
        end
    end
    if (calc_dav==1)
        Q0dav = zeros(nt,ns);
        for i=1:ns
            cos_thetam = max(0, (sin(decl)*sin(latc(i)))*ones(1,nstep) + ...
                (cos(latc(i))*cos(decl))*cos(2*pi*(dT/2:dT:1-dT/2)));
            Q0m = S0*(f*ones(1,nstep)) .* cos_thetam.^2 ./ (1.085*cos_thetam + (e0(:,i)*ones(1,nstep)).*(2.7+cos_thetam)*1e-3 + 0.1);
            Q0dav(:,i) = mean(Q0m,2);
        end
    end

elseif (Q_model==0.2)
    cos_theta = max(0, sin(decl*ones(1,ns)).*sin(ones(nt,1)*latc) + cos(ones(nt,1)*latc).*cos(decl*ones(1,ns)).*cos(h));
    Q0 = S0*(f*ones(1,ns)) .* cos_theta.^2 ./ (1.2*cos_theta + e0.*(1.0+cos_theta)*1e-3 + 0.0455);
    % Shine (1984) improvement of Zillman (1972) model, aimed at reducing understimation in Arctic regions.

    % Calculate hourly averages centred on the prediction hours
    if (calc_hav==1)
        Q0hav = zeros(nt,ns);
        for i=1:ns
            cos_thetam = max(0, (sin(decl)*sin(latc(i)))*ones(1,nsteph) + ...
                cos(latc(i))*(cos(decl)*ones(1,nsteph)).*cos(hm));
            Q0m = S0*(f*ones(1,nstep)) .* cos_thetam.^2 ./ (1.2*cos_thetam + (e0(:,i)*ones(1,nstep)).*(1.0+cos_thetam)*1e-3 + 0.0455);
            Q0hav(:,i) = mean(Q0m,2);
        end
    end
    if (calc_dav==1)
        Q0dav = zeros(nt,ns);
        for i=1:ns
            cos_thetam = max(0, (sin(decl)*sin(latc(i)))*ones(1,nstep) + ...
                (cos(latc(i))*cos(decl))*cos(2*pi*(dT/2:dT:1-dT/2)));
            Q0m = S0*(f*ones(1,nstep)) .* cos_thetam.^2 ./ (1.2*cos_thetam + (e0(:,i)*ones(1,nstep)).*(1.0+cos_thetam)*1e-3 + 0.0455);
            Q0dav(:,i) = mean(Q0m,2);
        end
    end

else
    error(['Q_model = ',num2str(Q_model),' not recognized'])
end

if (calc_Qs==1)
    % Reed (1977) model for cloud transmission.
    if (cloud_model==0)
        % Solar altitude at noon, in degrees
        alph = 180/pi*asin(sin(decl*ones(1,ns)).*sin(ones(nt,1)*latc) + cos(ones(nt,1)*latc).*cos(decl*ones(1,ns)));
        out.Qs = Q0.*(1 - 0.62*C + 0.0019*alph);

    elseif (cloud_model==0.1) %Laevastu (1960) approx. based on mid-latitude oceans.
        out.Qs = Q0.*(1 - 0.6*C.^3);

    else
        error(['cloud_model = ',num2str(cloud_model),' not recognized'])
    end
end

if (nargout>1)
    out.f = f; out.decl = decl; out.daylength = daylength;
    if calc_hav==1; out.Q0hav = Q0hav; end
    if calc_dav==1; out.Q0dav = Q0dav; end
end