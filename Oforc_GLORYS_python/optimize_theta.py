import numpy as np


def optimize_theta(thetai, thetamin, thetamax, Jfun, nr, opt):
    # function [thetabest,out] = optimize_theta(thetai,thetamin,thetamax,Jfun,nr,opt)

    # Minimize Jfun with respect to search parameters theta

    # Uses:    npermutek (optional) for uniform exhaustive sampling
    #          hessianest (optional) to calculate Hessian matrix

    # Phil Wallhead 23/12/2017
    # Joe McGovern 26/09/2023

    nparsh = len(thetamin)
    sztheta = thetai.shape
    thetai = thetai
    thetamin = thetamin
    thetamax = thetamax
    if len(varargin) < 5:
        nr = 1

    if len(varargin) < 6:
        opt = []

    if hasattr(opt, 'optJ') == 1:
        optJ = opt.optJ
    else:
        optJ = []

    if hasattr(opt, 'thetastart') == 1:
        thetastart = opt.thetastart
    else:
        thetastart = []

    if hasattr(opt, 'randomstart') == 1:
        randomstart = opt.randomstart
    else:
        randomstart = 0

    # 1 to randomize the start even for nr = 1
    if hasattr(opt, 'fnoise') == 1:
        fnoise = opt.fnoise
    else:
        fnoise = []

    if hasattr(opt, 'maxfunevals') == 1:
        maxfunevals = opt.maxfunevals
    else:
        maxfunevals = 20000.0

    if hasattr(opt, 'maxiter') == 1:
        maxiter = opt.maxiter
    else:
        maxiter = 10000.0

    if hasattr(opt, 'TolFun') == 1:
        TolFun = opt.TolFun
    else:
        TolFun = 1e-06

    if hasattr(opt, 'TolCon') == 1:
        TolCon = opt.TolCon
    else:
        TolCon = 1e-06

    if hasattr(opt, 'TolX') == 1:
        TolX = opt.TolX
    else:
        TolX = 1e-06

    if hasattr(opt, 'algorithm') == 1:
        algorithm = opt.algorithm
    else:
        algorithm = 0

    # 0 for fmincon/fminbnd (default)
    # 1 to use snobfit.m
    # Could set default to snobfit.m if you don't have the Optimisation toolbox
    if hasattr(opt, 'usefminbnd') == 1:
        usefminbnd = opt.usefminbnd
    else:
        usefminbnd = 1

    # 1 (def) to force use of fminbnd for 1D problems (seems to be most efficient)
    # This does however require another routine to calculate the hessian if req'd
    if nparsh == 1 and usefminbnd == 1:
        algorithm = 0

    if hasattr(opt, 'opta') == 1:
        opta = opt.opta
    else:
        opta = []

    # Additional options structure to be supplied to optimization algorithm
    if hasattr(opt, 'Hcalc') == 1:
        Hcalc = opt.Hcalc
    else:
        Hcalc = 0

    # 0   to switch off the hessian calculation (def)
    # 1   to use fmincon (gives similar estimate to hessianest in most cases)
    # 1.5 to use fminunc (gives similar estimate to hessianest in most cases)
    # 2   to use hessianest from DERIVEST suite by John D'Errico (def~=0)

    if hasattr(opt, 'divide_and_conquer') == 1:
        divide_and_conquer = opt.divide_and_conquer
    else:
        divide_and_conquer = []

    # 1 to divide the allowed range according to the starting value sampling design
    if hasattr(opt, 'rdesign') == 1:
        rdesign = opt.rdesign
    else:
        rdesign = []

    # 0 for uniform random sampling
    # 1 for uniform exhaustive sampling
    # 2 for uniform Latin-Hypercube sampling
    if len(rdesign) == 0 == 1:
        if nparsh == 1:
            rdesign = 1
        if nparsh > 1:
            rdesign = 2

    if divide_and_conquer == []:
        divide_and_conquer = 0
        if nparsh == 1 and hasattr(opt, 'rdesign') == 0:
            divide_and_conquer = 1

    if divide_and_conquer == 1:
        rdesign = 1

    if hasattr(opt, 'alph') == 1:
        alph = opt.alph
    else:
        alph = 0.05

    if hasattr(opt, 'ndof') == 1:
        ndof = opt.ndof
    else:
        ndof = np.inf

    if hasattr(opt, 'verbose') == 1:
        verbose = opt.verbose
    else:
        verbose = 0

    options = optimset('TolFun', TolFun, 'TolCon', TolCon, 'TolX', TolX, 'MaxFunEvals', maxfunevals, 'MaxIter', maxiter)
    # clear LASTN; LASTN = maxNumCompThreads(1); ##ok<NASGU>
    if algorithm == 1:
        nr = 1
        if Hcalc == 1:
            Hcalc = 1.5

    if algorithm == 2:
        opta.Verbosity = verbose

    # Set design for random starting values
    nr1 = nr
    if nr > 1 or randomstart == 1:
        if rdesign == 0:
            Xr = np.random.rand(nparsh, nr)
        if rdesign == 1:
            nr1 = nr ** (1 / nparsh)
            if np.abs(nr1 - np.round(nr1)) > 0.001:
                nr1 = int(np.floor(nr1))
                nr = nr1 ** nparsh
                if verbose > 0:
                    print(np.array(['Rounding nr1 down to nearest (1/nparsh)th root, corrected nr = ', int2str(nr)]))
            Xr = npermutek((np.arange(0.5, nr1 - 0.5 + 1)) / nr1, nparsh)
            Xr = np.transpose(Xr)
        if rdesign == 2:
            Xr = lhsdesign(nr, nparsh)
            Xr = np.transpose(Xr)
        #    thetarange = thetamax - thetamin;

    thetarange = thetamax - thetamin
    thetar = np.zeros((nparsh, nr))
    Jminr = np.zeros((1, nr))
    nfunevalsr = np.zeros((1, nr))
    Jminbest = 100000000.0
    thetabest = thetai
    Hbest = np.nan * np.ones((nparsh, nparsh))
    for i in np.arange(1, nr + 1).reshape(-1):
        if nr > 1 and len(thetastart) == 0 == 1:
            thetai1 = thetamin + np.multiply(thetarange, Xr[:, i])
        else:
            if len(thetastart) == 0 == 0:
                thetai1 = np.multiply(thetastart, (1 + fnoise * np.random.randn(nparsh, 1)))
                thetai1 = np.amin(thetamax, thetai1)
                thetai1 = np.amax(thetamin, thetai1)
            else:
                thetai1 = thetai
        thetamin1 = thetamin
        thetamax1 = thetamax
        if nr > 1 and divide_and_conquer == 1:
            dtheta = thetarange / (2 * nr1)
            thetamin1 = thetai1 - dtheta
            thetamax1 = thetai1 + dtheta
            if verbose > 1:
                print('[thetamin1 thetamax1] = ')
                print(np.array([thetamin1, thetamax1]))
        if algorithm == 0:
            if nparsh > 1 or usefminbnd == 0:
                if Hcalc == 1:
                    theta, Jmin, exitflag, outa, lambda_, grad, H = fmincon(lambda x=None: Jfun(x, optJ), thetai1, [],
                                                                            [], [], [], thetamin1, thetamax1, [],
                                                                            options)
                else:
                    theta, Jmin, exitflag, outa = fmincon(lambda x=None: Jfun(x, optJ), thetai1, [], [], [], [],
                                                          thetamin1, thetamax1, [], options)
            else:
                # fminbnd is more efficient for 1D problems
                theta, Jmin, exitflag, outa = fminbnd(lambda x=None: Jfun(x, optJ), thetamin1, thetamax1)
        if algorithm == 1:
            theta, Jmin, outa = fnsnobfit(lambda x=None: Jfun(x, optJ), [], [], thetamin1, thetamax1)
            theta = theta
        if algorithm == 2:
            opta.Generator = lambda x=None: (np.transpose(thetamin1) + np.multiply(np.random.rand(1, nparsh), (
                        np.transpose(thetamax1) - np.transpose(thetamin1))))
            theta, Jmin, outa = anneal(lambda x=None: Jfun(x, optJ), np.transpose(thetai1), opta)
            theta = theta
        if algorithm == 2.1:
            theta, Jmin, outa = sim_anl(lambda x=None: Jfun(x, optJ), thetai1, thetamin1, thetamax1, TolFun)
        if verbose > 0:
            vL = np.abs(theta - thetamin) / (thetamax - thetamin)
            vH = np.abs(theta - thetamax) / (thetamax - thetamin)
            if np.amin(np.array([[vL], [vH]])) < 0.001:
                print('Warning: Solution near boundary')
        if nr > 1:
            if verbose > 1:
                np.array([i, nr])
                if i > 1:
                    np.array([Jmin, Jminbest])
                    if Jminbest - Jmin > 0.1 and divide_and_conquer == 0:
                        print('Warning!!! LOCAL MINIMUM INFERRED ----------------------------')
                else:
                    Jmin
                np.transpose(theta(np.arange(1, np.amin(9, nparsh) + 1)))
            thetar[:, i] = theta
            Jminr[i] = Jmin
            nfunevalsr[i] = outa.funcCount
        if Jmin < Jminbest or i == 1:
            Jminbest = Jmin
            thetabest = theta
            nfunevalsbest = outa.funcCount
            if (usefminbnd == 0 or nparsh > 1) and Hcalc == 1:
                Hbest = H

    if nr > 1 and verbose > 0 and divide_and_conquer == 0:
        if np.amax(Jminr) - Jminbest > 0.1:
            print('Warning!!! LOCAL MINIMA WERE INFERRED ----------------------------')
            pause(0.5)

    if Hcalc > 0:
        vL = np.abs(thetabest - thetamin) / (thetamax - thetamin)
        vH = np.abs(thetabest - thetamax) / (thetamax - thetamin)
        if np.amin(np.array([[vL], [vH]])) < 0.001:
            if verbose > 0:
                print('Warning: Solution near boundary, Hessian-based SEs invalid (returning NaN)')
            Hbest = np.nan * np.ones((nparsh, nparsh))
        else:
            if Hcalc == 1 and nparsh > 1 and verbose > 0:
                print('Using fmincon to estimate Hessian')
            if Hcalc == 1.5 or (Hcalc == 1 and nparsh == 1 and usefminbnd == 1):
                # Compute the hessian using fminunc started at thetabest
                if Hcalc == 1 and verbose > 0:
                    print('Cannot use fmincon started at thetabest to estimate Hessian - it quits to early')
                if verbose > 0:
                    print('Using fminunc started at thetabest to estimate Hessian')
                thetacheck, Jmincheck, exitflag, outa, grad, Hbest = fminunc(lambda x=None: Jfun(x, optJ), thetabest,
                                                                             options)
                if np.amax(np.abs(thetacheck - thetabest)) > 2 * 0.0001 or np.abs(Jmincheck - Jminbest) > 0.001:
                    if verbose > 0:
                        print('fminunc gives different fit to fmincon/fminbnd: something strange happened!')
                    np.array([thetacheck, thetabest])
                    np.array([Jmincheck, Jminbest])
                    stop
            if Hcalc == 2:
                if verbose > 0:
                    print('Using hessianest to estimate Hessian')
                Hbest = hessianest(lambda x=None: Jfun(x, optJ), thetabest, thetamin, thetamax)
        Ctheta = np.linalg.solve((0.5 * Hbest), np.eye(nparsh))
        thetase = np.sqrt(diag(Ctheta))
        fac = tinv(1 - alph / 2, ndof)
        if sztheta(1) > 1:
            thetase = thetase
            thetaCIs = np.array([thetabest - fac * thetase, thetabest + fac * thetase])
        else:
            thetase = np.transpose(thetase)
            thetaCIs = np.array([[np.transpose(thetabest) - fac * thetase], [np.transpose(thetabest) + fac * thetase]])

    # Return thetabest same size as thetai
    if sztheta(1) > 1:
        thetabest = thetabest
    else:
        thetabest = np.transpose(thetabest)

    # Return other outputs
    if nargout > 1:
        out = []
        out.thetabest = thetabest
        out.Jminbest = Jminbest
        out.nfunevalsbest = nfunevalsbest
        if Hcalc > 0:
            out.Hbest = Hbest
            out.Ctheta = Ctheta
            out.thetase = thetase
            out.thetaCIs = thetaCIs
        if nr > 1:
            out.thetar = thetar
            out.Jminr = Jminr
            out.nfunevalsr = nfunevalsr

    return thetabest, out
