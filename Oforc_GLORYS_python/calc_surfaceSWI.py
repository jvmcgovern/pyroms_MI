import numpy as np


def calc_surfaceSWI(lat, lon, tday, opt):
    # function [Q0,out] = calc_surfaceSWI(lat,lon,tday,opt)

    # Calculate short-wave irradiance above sea surface [W/m2] as instantaneous (Q0)
    # and daily-average values (out.Q0dav, if opt.calc_dav=1), given:
    # lat = latitude [degrees]
    # lon = longitude [degrees] (needed for hour angle calculation)
    # tday = days since 1st January (0-364 or 0-365 for leap year, real-valued, UTC)
    # opt.JD = Julian day (if needed, see below)

    # Phil Wallhead 29/01/2016
    # Joe McGovern 31/08/2023

    lat = np.transpose(lat)
    lon = np.transpose(lon)
    ns = len(lat)
    tday = tday
    nt = len(tday)

    if hasattr(opt, 'use_tday_lag'):
        use_tday_lag = opt.use_tday_lag
    else:
        use_tday_lag = 0

    if hasattr(opt, 'year'):
        year = opt.year
    else:
        year = []

    if hasattr(opt, 'f_model'):
        f_model = opt.f_model
    else:
        f_model = 1

    if hasattr(opt, 'decl_model'):
        decl_model = opt.decl_model
    else:
        decl_model = 1

    if hasattr(opt, 'use_eqtime'):
        use_eqtime = opt.use_eqtime
    else:
        use_eqtime = 1

    if hasattr(opt, 'Q_model'):
        Q_model = opt.Q_model
    else:
        Q_model = 0

    if hasattr(opt, 'cloud_model'):
        cloud_model = opt.cloud_model
    else:
        cloud_model = 0

    if hasattr(opt, 'calc_hav'):
        calc_hav = opt.calc_hav
    else:
        calc_hav = 0

    if hasattr(opt, 'calc_dav'):
        calc_dav = opt.calc_dav
    else:
        calc_dav = 0

    if hasattr(opt, 'calc_Qs'):
        calc_Qs = opt.calc_Qs
    else:
        calc_Qs = 0

    if hasattr(opt, 'JD'):
        JD = opt.JD
    else:
        JD = []

    if hasattr(opt, 'e0'):
        e0 = opt.e0
    else:
        e0 = np.zeros((nt, ns))

    if hasattr(opt, 'C'):
        C = opt.C
    else:
        C = np.zeros((nt, ns))

    if hasattr(opt, 'ndays_year'):
        ndays_year = opt.ndays_year
    else:
        ndays_year = 365.2422

    if hasattr(opt, 'ROMS_model'):
        ROMS_model = opt.ROMS_model
    else:
        ROMS_model = 0

    if ROMS_model == 1:
        f_model = 0
        decl_model = 0.1
        use_eqtime = 0
        Q_model = 0.1
        cloud_model = 0.1
        ndays_year = 365.2425

    # out = struct('f_model', f_model, 'decl_model', decl_model, 'Q_model', Q_model,
    #              'cloud_model', cloud_model, 'ndays_year', ndays_year)
    class out:
        def __init__(self, f_model, decl_model, Q_model, cloud_model, ndays_year):
            self.f_model = f_model
            self.decl_model = decl_model
            self.Q_model = Q_model
            self.cloud_model = cloud_model
            self.ndays_year = ndays_year


    # Stepping parameters for computing daily averages (where analytical integration not possible)
    nstep = 24
    dT = 1 / nstep
    S0 = 1367

    if use_tday_lag == 1 and not len(year) == 0:
        tdayc = tday - np.mod(year, 4) / 4
    else:
        tdayc = tday

    # Calculate the Sun-Earth distance factor (f)
    if f_model == 0:
        f = np.ones((nt, 1))
    else:
        if f_model == 1:
            # Originally from Spencer (1971), see Brun and Pinardi (2007) Table A1
            an = np.transpose(np.array([1.00011, 0.034221, 0.000719]))
            bn = np.transpose(np.array([0.0, 0.00128, 7.7e-05]))
            t = 2 * np.pi * tdayc / ndays_year
            f = 0 * t
            for i in np.arange(1, 3 + 1).reshape(-1):
                n = i - 1
                f = f + an[i] * np.cos(n * t) + bn[i] * np.sin(n * t)
        else:
            if f_model == 2:
                g = (np.pi / 180) * (357.528 + 0.9856003 * (JD - 2451545.0))
                f = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
            else:
                raise Exception(np.array(['f_model = ', str(f_model), ' not recognized']))

    # Calculate the declination (solar inclination) in radians
    if decl_model == 0:
        decl = 23.45 * np.pi / 180 * np.sin(2 * np.pi * (284 + tdayc) / ndays_year)
    else:
        if decl_model == 0.1:
            yday = tdayc + 1
            decl = 23.44 * np.pi / 180 * np.cos((172.0 - yday) * 2.0 * np.pi / ndays_year)
        else:
            if decl_model == 1:
                cn = np.transpose(np.array([0.006918, - 0.399912, - 0.006758, - 0.002697]))
                dn = np.transpose(np.array([0.0, 0.070257, 0.000907, 0.000148]))
                t = 2 * np.pi * tdayc / ndays_year
                decl = 0 * t
                for i in np.arange(1, 4 + 1).reshape(-1):
                    n = i - 1
                    decl = decl + cn[i] * np.cos(n * t) + dn[i] * np.sin(n * t)
            else:
                raise Exception(np.array(['decl_model = ', str(decl_model), ' not recognized']))

    latc = lat * np.pi / 180

    H = np.arccos(np.amin(1, np.amax(- 1, np.multiply(- np.tan(decl), np.tan(latc)))))

    # (see Liou02_Int_Geophys_Chapter2, eqn 2.2.2)
    # Note: capping accounts for zero daylength during high-latitude winter
    daylength = 24 * H / np.pi

    hour = (tday - int(np.floor(tday))) * 24

    fracyear = 2 * np.pi * tdayc / ndays_year
    eqtime = 229.18 * (7.5e-05 + 0.001868 * np.cos(fracyear) - 0.032077 * np.sin(fracyear) - 0.014615 * np.cos(
        2 * fracyear) - 0.040849 * np.sin(2 * fracyear))

    # This is from NOAA formulae, see sunrisesunset.m, sunrisesunset2.m.

    if use_eqtime == 1:
        h = (12.0 - hour * np.ones((1, ns))) * np.pi / 12.0 - (
                np.ones((nt, 1)) * lon) * np.pi / 180 - eqtime * np.pi / (60 * 12)
    else:
        h = (12.0 - hour * np.ones((1, ns))) * np.pi / 12.0 - (np.ones((nt, 1)) * lon) * np.pi / 180
        # This is the 'negative' form used in ROMS, which is equivalent to positive form
        # e.g. "15deg*(LST-12)", since cos(h) = cos(-h).

    # Note: h is an [nt*ns] array.

    if calc_hav == 1:
        dhv = (np.arange(- 29.5, 29.5 + 1, 1)) / 60 * - np.pi / 12
        nsteph = len(dhv)
        hm = h * np.ones((1, nsteph)) + dhv

    if Q_model == 0:
        tau = 0.7
        A_a = 0.09
        cos_theta = np.multiply(np.sin(decl * np.ones((1, ns))), np.sin(np.ones((nt, 1)) * latc)) + np.multiply(
            np.multiply(np.cos(np.ones((nt, 1)) * latc), np.cos(decl * np.ones((1, ns)))), np.cos(h))
        SE = np.amax(0.0, np.multiply(S0 * (f * np.ones((1, ns))), cos_theta))
        rec = np.argwhere(SE > 0)
        I1 = np.zeros((nt, ns))
        transfac = I1
        QDir = I1
        QDiff = I1
        transfac[rec] = tau ** (1.0 / cos_theta(rec))
        QDir[rec] = np.multiply(SE[rec], transfac[rec])
        QDiff[rec] = np.multiply(0.5 * SE[rec], (1 - A_a - transfac[rec]))
        Q0 = QDir + QDiff
        if calc_hav == 1:
            Q0hav = I1
            transfacm = np.zeros((nt, nsteph))
            QDirm = np.zeros((nt, nsteph))
            QDiffm = np.zeros((nt, nsteph))
            for i in np.arange(1, ns + 1).reshape(-1):
                cos_thetam = (np.sin(decl) * np.sin(latc(i))) * np.ones((1, nsteph)) + np.multiply(
                    np.cos(latc(i)) * (np.cos(decl) * np.ones((1, nsteph))), np.cos(hm))
                SEm = np.amax(0.0, np.multiply(S0 * (f * np.ones((1, nsteph))), cos_thetam))
                rec = np.argwhere(SEm > 0)
                transfacm[rec] = tau ** (1.0 / cos_thetam(rec))
                QDirm[rec] = np.multiply(SEm[rec], transfacm[rec])
                QDiffm[rec] = np.multiply(0.5 * SEm[rec], (1 - A_a - transfacm[rec]))
                Q0m = QDirm + QDiffm
                Q0hav[:, i] = np.mean(Q0m, 2)
        if calc_dav == 1:
            Q0dav = I1
            for i in np.arange(1, ns + 1).reshape(-1):
                cos_thetam = (np.sin(decl) * np.sin(latc(i))) * np.ones((1, nstep)) + (
                        np.cos(latc(i)) * np.cos(decl)) * np.cos(
                    2 * np.pi * (np.arange(dT / 2, 1 - dT / 2 + dT, dT)))
                SEm = np.amax(0.0, np.multiply(S0 * (f * np.ones((1, nstep))), cos_thetam))
                transfacm = tau ** (1.0 / cos_thetam)
                QDirm = np.multiply(SEm, transfacm)
                QDiffm = np.multiply(0.5 * SEm, (1 - A_a - transfacm))
                Q0m = QDirm + QDiffm
                Q0dav[:, i] = np.mean(Q0m, 2)
    else:
        if Q_model == 0.1:
            cos_theta = np.amax(0, np.multiply(np.sin(decl * np.ones((1, ns))),
                                               np.sin(np.ones((nt, 1)) * latc)) + np.multiply(
                np.multiply(np.cos(np.ones((nt, 1)) * latc), np.cos(decl * np.ones((1, ns)))), np.cos(h)))
            Q0 = np.multiply(S0 * (f * np.ones((1, ns))), cos_theta ** 2) / (
                    1.085 * cos_theta + np.multiply(e0, (2.7 + cos_theta)) * 0.001 + 0.1)
            # Zillman (1972) approx. (Eqn. 6 in Niemela et al., 2001) = ROMS formula before correcting for cloud cover.
            if calc_hav == 1:
                Q0hav = np.zeros((nt, ns))
                for i in np.arange(1, ns + 1).reshape(-1):
                    cos_thetam = np.amax(0, (np.sin(decl) * np.sin(latc(i))) * np.ones((1, nsteph)) +
                                         np.multiply(np.cos(latc(i)) * (np.cos(decl) *
                                                                        np.ones((1, nsteph))), np.cos(hm)))
                    Q0m = np.multiply(S0 * (f * np.ones((1, nstep))), cos_thetam ** 2) / (
                            1.085 * cos_thetam + np.multiply((e0[:, i] * np.ones((1, nstep))),
                                                             (2.7 + cos_thetam)) * 0.001 + 0.1)
                    Q0hav[:, i] = np.mean(Q0m, 2)
            if calc_dav == 1:
                Q0dav = np.zeros((nt, ns))
                for i in np.arange(1, ns + 1).reshape(-1):
                    cos_thetam = np.amax(0, (np.sin(decl) * np.sin(latc(i))) * np.ones((1, nstep)) +
                                         (np.cos(latc(i)) * np.cos(decl)) * np.cos(2 * np.pi *
                                                                                   (np.arange(dT / 2, 1 - dT / 2 + dT,
                                                                                              dT))))
                    Q0m = np.multiply(S0 * (f * np.ones((1, nstep))), cos_thetam ** 2) / (
                            1.085 * cos_thetam + np.multiply((e0[:, i] * np.ones((1, nstep))),
                                                             (2.7 + cos_thetam)) * 0.001 + 0.1)
                    Q0dav[:, i] = np.mean(Q0m, 2)
        else:
            if Q_model == 0.2:
                cos_theta = np.amax(0, np.multiply(np.sin(decl * np.ones((1, ns))),
                                                   np.sin(np.ones((nt, 1)) * latc)) + np.multiply(
                    np.multiply(np.cos(np.ones((nt, 1)) * latc), np.cos(decl * np.ones((1, ns)))), np.cos(h)))
                Q0 = np.multiply(S0 * (f * np.ones((1, ns))), cos_theta ** 2) / (
                        1.2 * cos_theta + np.multiply(e0, (1.0 + cos_theta)) * 0.001 + 0.0455)
                # Shine (1984) improvement of Zillman (1972) model, aimed at reducing understimation in Arctic regions.
                if calc_hav == 1:
                    Q0hav = np.zeros((nt, ns))
                    for i in np.arange(1, ns + 1).reshape(-1):
                        cos_thetam = np.amax(0, (np.sin(decl) * np.sin(latc(i))) * np.ones((1, nsteph)) + np.multiply(
                            np.cos(latc(i)) * (np.cos(decl) * np.ones((1, nsteph))), np.cos(hm)))
                        Q0m = np.multiply(S0 * (f * np.ones((1, nstep))), cos_thetam ** 2) / (
                                1.2 * cos_thetam + np.multiply((e0[:, i] * np.ones((1, nstep))),
                                                               (1.0 + cos_thetam)) * 0.001 + 0.0455)
                        Q0hav[:, i] = np.mean(Q0m, 2)
                if calc_dav == 1:
                    Q0dav = np.zeros((nt, ns))
                    for i in np.arange(1, ns + 1).reshape(-1):
                        cos_thetam = np.amax(0, (np.sin(decl) * np.sin(latc(i))) * np.ones((1, nstep)) + (
                                np.cos(latc(i)) * np.cos(decl)) * np.cos(
                            2 * np.pi * (np.arange(dT / 2, 1 - dT / 2 + dT, dT))))
                        Q0m = np.multiply(S0 * (f * np.ones((1, nstep))), cos_thetam ** 2) / (
                                1.2 * cos_thetam + np.multiply((e0[:, i] * np.ones((1, nstep))),
                                                               (1.0 + cos_thetam)) * 0.001 + 0.0455)
                        Q0dav[:, i] = np.mean(Q0m, 2)
            else:
                raise Exception(np.array(['Q_model = ', str(Q_model), ' not recognized']))

    if calc_Qs == 1:
        if cloud_model == 0:
            alph = 180 / np.pi * np.arcsin(
                np.multiply(np.sin(decl * np.ones((1, ns))), np.sin(np.ones((nt, 1)) * latc)) + np.multiply(
                    np.cos(np.ones((nt, 1)) * latc), np.cos(decl * np.ones((1, ns)))))
            out.Qs = np.multiply(Q0, (1 - 0.62 * C + 0.0019 * alph))
        else:
            if cloud_model == 0.1:
                out.Qs = np.multiply(Q0, (1 - 0.6 * C ** 3))
            else:
                raise Exception(np.array(['cloud_model = ', str(cloud_model), ' not recognized']))

    out.f = f
    out.decl = decl
    out.daylength = daylength
    if calc_hav == 1:
        out.Q0hav = Q0hav
    if calc_dav == 1:
        out.Q0dav = Q0dav

    return Q0, out
