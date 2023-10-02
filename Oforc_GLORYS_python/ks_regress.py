import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t
import fnperiodic_distm


def conformsizes(nm, A1, A2, A3, A4):
    # function [A1c,A2c,A3c,A4c] = conformsizes(nm,A1,A2,A3,A4)
    # Conform up to 4 array sizes to [nm(1)*nm(2)] by repetition
    # If nm is empty it is taken from size(A1)
    # Phil Wallhead 26/04/2014
    # Joe McGovern 31/08/2023

    if len(nm) == 0:
        nm = A1.shape

    n = nm[0]
    m = nm[1]
    A1c = A1

    if len(A1) == 1:
        A1c = A1 * np.ones((n, m))

    if len(A1[:, 0]) > 1 and len(A1[0, :]) == 1:
        A1c = A1 * np.ones((1, m))

    if len(A1[:, 0]) == 1 and len(A1[0, :]) > 1:
        A1c = np.ones((n, 1)) * A1

    if A2:
        A2c = A2
        if len(A2) == 1:
            A2c = A2 * np.ones((n, m))
        if len(A2[:, 0]) > 1 and len(A2[0, :]) == 1:
            A2c = A2 * np.ones((1, m))
    if len(A2[:, 0]) == 1 and len(A2[0, :]) > 1:
        A2c = np.ones((n, 1)) * A2

    if A3:
        A3c = A3
        if len(A3) == 1:
            A3c = A3 * np.ones((n, m))
        if len(A3[:, 0]) > 1 and len(A3[0, :]) == 1:
            A3c = A3 * np.ones((1, m))
    if len(A3[:, 0]) == 1 and len(A3[0, :]) > 1:
        A3c = np.ones((n, 1)) * A3

    if A4:
        A4c = A4
        if len(A4) == 1:
            A4c = A4 * np.ones((n, m))
        if len(A4[:, 0]) > 1 and len(A4[0, :]) == 1:
            A4c = A4 * np.ones((1, m))
    if len(A4[:, 0]) == 1 and len(A4[0, :]) > 1:
        A4c = np.ones((n, 1)) * A4

    return A1c, A2c, A3c, A4c


def ks_regress(y, x, xp, bdw, opt):
    # function [yp,out] = ks_regress(y,x,xp,bdw,opt)

    # Kernel smoothing regression of y over m dimensions from x to xp
    # Parameter 'bdw' [ni*m], [1*m], or scalar, describing smoothing bandwidths
    # for ni interpolation points and m smoothing dimensions.
    # Optional parameters (opt.param):
    # 'kfunc' scalar defining form of kernel function:
    #   kfunc = 1: e^(-s) kernel (default)
    #   kfunc = 2: Gaussian kernel
    #   kfunc = 3: sinc kernel
    #   kfunc = 4: box-car running mean
    #   kfunc = 5: Cauchy/Lorentzian kernel
    # 'obserrv' [nd*1] observational error variances to weight the data (optional)
    # 'rad' [1*m] allowed radii of influence (default = Inf, but this can be slow)

    # Can also output weight matrix 'wt' for replicated calculations,
    # and the indices 'recs' of the data used for the interpolation.
    # Reference: Hastie et al., 2009.  The Elements of Statistical Learning.

    # Validated against npreg.R from package "np" for local constant and local linear regressions
    # See testks_regress.m for examples.
    # Uses: conformsizes.m (included in this file)
    # fnperiodic_distm.m

    ## Phil Wallhead 10/02/2019
    ## Joe McGovern 31/08/2023 - conversion to python 3

    ni, m = xp.shape
    nd = len(y)

    if nd == 0:
        print('No data for kernel smoothing - stopping!')
        out = []
    else:

        # if len(varargin) < 5:
        #     opt = []

        if hasattr(opt, 'verbose'):
            verbose = opt.verbose
        else:
            verbose = 0

        if hasattr(opt, 'nblocks'):
            nblocks = opt.nblocks
        else:
            nblocks = []

        # Can specify here the number of blocks into which to divide the prediction set
        if hasattr(opt, 'subset_wrt_xp'):
            subset_wrt_xp = opt.subset_wrt_xp
        else:
            subset_wrt_xp = (m < 3)

        # 1 to subset the predictions with respect to the range of prediction x values xp
        if hasattr(opt, 'one_by_one'):
            one_by_one = opt.one_by_one
        else:
            one_by_one = 0

        # 1 to make predictions one by one (not usually fastest)
        if one_by_one == 1:
            nblocks = ni
            subset_wrt_xp = 0

        if hasattr(opt, 'nimax'):
            nimax = opt.nimax
        else:
            nimax = np.inf

        if ni > nimax:
            nblocks = np.ceil(ni / nimax)
            if verbose > 0:
                print(np.array(
                    ['ks_regress: Prediction set (ni = ', str(ni), ') exceeds maximum size (nimax = ', str(nimax),
                     '); setting nblocks = ', str(nblocks)]))

        if len(nblocks) == 0:
            if nd > 1000.0 and ni > 1000.0:
                nblocks = 10
                if verbose > 0:
                    print(np.array(
                        ['ks_regress: Dividing prediction set into ', str(nblocks), ' blocks to limit matrix sizes']))
            else:
                nblocks = 1

        if hasattr(opt, 'kfunc'):
            kfunc = opt.kfunc
        else:
            kfunc = 1

        # Choice of smoothing kernel (def = exponential (Laplacian))
        if hasattr(opt, 'one_sided'):
            one_sided = opt.one_sided
        else:
            one_sided = np.zeros((1, m))

        # Set one_sided[i] = -1 to only consider influence of data with values x(:,i) <= xp(i)
        #    one_sided[i] = 1 to only consider influence of data with values x(:,i) >= xp(i)
        if hasattr(opt, 'obserrv'):
            obserrv = opt.obserrv
        else:
            obserrv = []

        if hasattr(opt, 'yprior'):
            yprior = opt.yprior
        else:
            yprior = []

        if len(yprior) == 1 and ni > 1:
            yprior = yprior * np.ones((ni, 1))

        if hasattr(opt, 'priorerrv'):
            priorerrv = opt.priorerrv
        else:
            priorerrv = np.ones((ni, 1))

        if len(priorerrv) == 1 and ni > 1:
            priorerrv = priorerrv * np.ones((ni, 1))

        if len(obserrv) == 0 and not len(priorerrv) == 0:
            obserrv = np.ones((nd, 1))

        if hasattr(opt, 'use_rad'):
            use_rad = opt.use_rad
        else:
            use_rad = 0

        # 1 to use buffer vector "rad" to limit initial data set such that nds<=ndsmax (def = 0)
        if use_rad == 1:
            if hasattr(opt, 'variable_rad'):
                variable_rad = opt.variable_rad
            else:
                variable_rad = 0
            # 1 to vary rad iteratively in order to limit initial data set such that ndsmin<=nds<=ndsmax (def = 0)
            if hasattr(opt, 'rad'):
                rad = opt.rad
            else:
                rad = []
            # maximum "radius" of allowed influence of data points (default = 10 * bdw)
            # Beware: restricting rad can lead to irregularities
            # Note: rad must be in units of the input coordinates (e.g. in degrees latitude/longitude if aslatlon==1)
            # Note: user inputs target or default rad/radfac; rad will be changed if nds<ndsmin or nds>ndsmax
            if hasattr(opt, 'radfac'):
                radfac = opt.radfac
            else:
                radfac = []
            if len(rad) == 0 == 1:
                if len(radfac) == 0 == 1:
                    if verbose > 0:
                        print('Assuming default radii of influence rad = 10 * bdw (radfac = 10)')
                    radfac = 10
                    rad = radfac * bdw
                else:
                    rad = radfac * bdw
            if hasattr(opt, 'radfacmin'):
                radfacmin = opt.radfacmin
            else:
                radfacmin = 2
            # Minimum constraint on rad: radmin = radfacmin*bdw
            if hasattr(opt, 'radfacmax'):
                radfacmax = opt.radfacmax
            else:
                radfacmax = []
            # Maximum constraint on rad: radmax = radfacmax*bdw
            if len(radfacmax) == 0 == 0:
                rad = radfacmax * bdw
            else:
                radfacmax = np.inf
            if hasattr(opt, 'ndsmin'):
                ndsmin = opt.ndsmin
            else:
                ndsmin = np.amin(nd, 20)
            # Increase ndsmin to avoid irregular smoothing in data-sparse regions
            if hasattr(opt, 'ndsmax'):
                ndsmax = opt.ndsmax
            else:
                ndsmax = 5000
            # Target maximum number of data (to limit matrix size)
            if hasattr(opt, 'incradmax'):
                incradmax = opt.incradmax
            else:
                incradmax = 100.0

        if hasattr(opt, 'ndclosest'):
            ndclosest = opt.ndclosest
        else:
            ndclosest = []

        # Maximum number of closest data (nearest neighbours) to use for smoothing (only for ni1 = 1)
        if not len(ndclosest) == 0:
            nblocks = ni
            subset_wrt_xp = 0

        if hasattr(opt, 'lr_order'):
            lr_order = opt.lr_order
        else:
            lr_order = 0

        # Order of local polynomial regression (0 for local-constant Nadaraya-Watson method)
        # Can be vector [1 * m]
        if hasattr(opt, 'xc'):
            xc = opt.xc
        else:
            xc = []

        # Additional linear covariates (an alternative to adding to x and setting extra bdws = np.inf)
        if hasattr(opt, 'xpc'):
            xpc = opt.xpc
        else:
            xpc = []

        # Additional linear prediction covariates (an alternative to adding to xp and setting extra bdws = np.inf)
        if hasattr(opt, 'lr_orderc'):
            lr_orderc = opt.lr_orderc
        else:
            lr_orderc = 0

        # Order of local polynomial regression over additional covariates
        # Can be vector [1 * mc]
        if hasattr(opt, 'period'):
            period = opt.period
        else:
            period = np.nan * np.ones((1, m))

        # Can specify looping period such that dx_ij = mod(x_i-x_j,period)
        if hasattr(opt, 'aslatlon'):
            aslatlon = opt.aslatlon
        else:
            aslatlon = []

        # [1*2] vector specifying indices in (x,xp) of (latitude,longitude) for spherical distance treatment
        # Note: the assumed spherical distance bandwidth in km is bdw(aslatlon(1)),]; bdw(aslatlon(2)) is not used
        if len(aslatlon) == 1:
            if aslatlon == 1:
                aslatlon = np.array([1, 2])
            if aslatlon == 0:
                aslatlon = []

        if len(aslatlon) == 0 == 0 and hasattr(opt, 'rad') == 0:
            if verbose > 0:
                print('Warning: if interpolating over lat/lon, the allowed '
                      'radii rad should be input (since cannot be inferred from bdw)')
                print('Assuming rad = np.inf for lat/lon')
            rad[aslatlon] = np.inf

        if hasattr(opt, 'calcunc'):
            calcunc = opt.calcunc
        else:
            calcunc = 0

        # 1 to calculate standard errors and confidence/prediction intervals (mpL,mpH)/(ypL,ypH)
        if hasattr(opt, 'variable_res_var'):
            variable_res_var = opt.variable_res_var
        else:
            variable_res_var = 1

        # 1 (def) to allow variable residual variance V = V(x), otherwise V = V0
        if hasattr(opt, 'nneg'):
            nneg = opt.nneg
        else:
            nneg = 0

        if hasattr(opt, 'vbdw'):
            vbdw = opt.vbdw
        else:
            vbdw = bdw

        if hasattr(opt, 'vlr_order'):
            vlr_order = opt.vlr_order
        else:
            vlr_order = lr_order

        if hasattr(opt, 'alph'):
            alph = opt.alph
        else:
            alph = 0.05

        if hasattr(opt, 'check_evenness'):
            check_evenness = opt.check_evenness
        else:
            check_evenness = 0

        # 1 to check evenness of the weight matrix (for dominant data, def = 0)
        if hasattr(opt, 'domfrac'):
            domfrac = opt.domfrac
        else:
            domfrac = 0.95

        if hasattr(opt, 'x_coverage'):
            x_coverage = opt.x_coverage
        else:
            x_coverage = []

        # Optional "coverage variable" size nd*1
        # When supplied, out.xp_coverage gives the number of different integer values of
        # x_coverage used for each prediction point
        if hasattr(opt, 'coverage_type'):
            coverage_type = opt.coverage_type
        else:
            coverage_type = 0

        # 0 for simple coverage: range of values of x_coverage for data contributing to yp at each point
        # 1 for number of unique values of x_coverage
        if hasattr(opt, 'fr_wt_coverage'):
            fr_wt_coverage = opt.fr_wt_coverage
        else:
            fr_wt_coverage = 0

        # Minimum relative weight a data point must have to be included in the coverage calculation
        if hasattr(opt, 'max_sdist_coverage'):
            max_sdist_coverage = opt.max_sdist_coverage
        else:
            max_sdist_coverage = 2

        # Maximum scaled distance from the prediction point that a data point must
        # have to be included in the coverage calculation

        # Plotting options
        if hasattr(opt, 'doplot'):
            doplot = opt.doplot
        else:
            doplot = 0

        if hasattr(opt, 'ifig'):
            ifig = opt.ifig
        else:
            ifig = 17

        if hasattr(opt, 'nrows'):
            nrows = opt.nrows
        else:
            nrows = 1

        if hasattr(opt, 'ncols'):
            ncols = opt.ncols
        else:
            ncols = 1

        if hasattr(opt, 'isub'):
            isub = opt.isub
        else:
            isub = 1

        if hasattr(opt, 'plottype'):
            plottype = opt.plottype
        else:
            plottype = 1

        if hasattr(opt, 'tstr'):
            tstr = opt.tstr
        else:
            tstr = []

        # Options for plots against covariate (e.g. time)
        if hasattr(opt, 'doplotc'):
            doplotc = opt.doplotc
        else:
            doplotc = 0

        if len(xc) == 0 == 1:
            doplotc = 0

        if doplotc == 1:
            if hasattr(opt, 'max_sdist_plotc'):
                max_sdist_plotc = opt.max_sdist_plotc
            else:
                max_sdist_plotc = 2
            # Maximum scaled distance from the prediction point for a datum
            # to be included in the plot vs. covariate

            if hasattr(opt, 'selpp'):
                selpp = opt.selpp
            else:
                selpp = 1
            # If smoothing blocks of > 1 prediction points, selpp specifies
            # which prediction point will be plotted
            if hasattr(opt, 'xcpp'):
                xcpp = opt.xcpp
            else:
                xcpp = np.array([[np.amin(xc)], [np.amax(xc)]])

            # Vector of covariate values at which to predict the response in the plot
            ncpp = len(xcpp)
            if hasattr(opt, 'xcoff'):
                xcoff = opt.xcoff
            else:
                xcoff = 0
            # Offset to apply to covariate x axis in plots
            if hasattr(opt, 'ifigc'):
                ifigc = opt.ifigc
            else:
                ifigc = 18
            if hasattr(opt, 'kcstep'):
                kcstep = opt.kcstep
            else:
                kcstep = np.inf
            # Set kcstep>1 to only plot for every (kcstep)^th valid prediction point
            if hasattr(opt, 'xpp'):
                xpp = opt.xpp
            else:
                xpp = []
            # Alternatively provide a set of points at which to make a plot
            npp = len(xpp[:, 0])
            if hasattr(opt, 'nrowsc'):
                nrowsc = opt.nrowsc
            else:
                nrowsc = 5
            if hasattr(opt, 'ncolsc'):
                ncolsc = opt.ncolsc
            else:
                ncolsc = 3
            if hasattr(opt, 'ycol'):
                ycol = opt.ycol
            else:
                ycol = 'k'
            if hasattr(opt, 'ypcol'):
                ypcol = opt.ypcol
            else:
                ypcol = 'r'
            if hasattr(opt, 'ymarkersize'):
                ymarkersize = opt.ymarkersize
            else:
                ymarkersize = 8
            if hasattr(opt, 'linewidth'):
                linewidth = opt.linewidth
            else:
                linewidth = 2
            if hasattr(opt, 'lfontsize'):
                lfontsize = opt.lfontsize
            else:
                lfontsize = 10
            if hasattr(opt, 'xstr'):
                xstr = opt.xstr
            else:
                xstr = 't'
            if hasattr(opt, 'ystr'):
                ystr = opt.ystr
            else:
                ystr = 'y'
            if hasattr(opt, 'xlimc'):
                xlimc = opt.xlimc
            else:
                xlimc = []
            if hasattr(opt, 'ylimc'):
                ylimc = opt.ylimc
            else:
                ylimc = []

        bdw = conformsizes(np.array([ni, m]), bdw)
        if use_rad == 1:
            if len(rad) == 1:
                rad = rad * np.ones((1, m))

        if len(lr_order) == 1:
            lr_order = lr_order * np.ones((1, m))

        if len(one_sided) == 1:
            one_sided = one_sided * np.ones((1, m))

        mc = 0
        isubc = 1
        kc = 0
        if len(xc) == 0 == 0:
            mc = len(xc[0, :])
            if len(xpc) == 0:
                xpc = np.zeros((ni, mc))
            if len(xpc[:, 0]) == 1:
                xpc = np.ones((ni, 1)) * xpc
            if len(lr_orderc) == 1:
                lr_orderc = lr_orderc * np.ones((1, mc))
            if doplotc == 1:
                plt.figure(ifigc)

                nsubc = nrowsc * ncolsc
        else:
            doplotc = 0
            nsubc = 0

        # No. linear parameters fitted locally
        nb1 = sum(lr_order)
        nbc = sum(lr_orderc)
        nb = 1 + nb1 + nbc

        # Divide the prediction set into blocks if req'd (to limit matrix size)
        if nblocks > 1:
            if subset_wrt_xp == 0:
                bblock = np.array([0, np.round((np.arange(ni / nblocks, ni + ni / nblocks, ni / nblocks)))])
            if subset_wrt_xp == 1:
                nblocks1 = np.round(nblocks ** (1 / m))
                nblocks = nblocks1 ** m
                minxp = np.amin(xp)
                maxxp = np.amax(xp)
                dxp = (maxxp - minxp) / nblocks1
                xpL1 = np.nan * np.ones((nblocks1, m))
                xpH1 = np.nan * np.ones((nblocks1, m))
                for i in np.arange(1, m + 1).reshape(-1):
                    xpL1[:, i] = np.transpose((np.arange(minxp(i), maxxp(i) - dxp(i) + dxp(i), dxp(i))))
                    xpH1[:, i] = np.transpose((np.arange(minxp(i) + dxp(i), maxxp(i) + dxp(i), dxp(i))))
                if m == 1:
                    xpL = xpL1
                    xpH = xpH1
                else:
                    if m == 2:
                        X, Y = np.mgrid(xpL1[:, 0], xpL1[:, 1])
                        xpL = np.array([X, Y])
                        # clear('X','Y')
                        X, Y = np.mgrid(xpH1[:, 0], xpH1[:, 1])
                        xpH = np.array([X, Y])
                        # clear('X','Y')
                        # [xpL xpH]
                    # stop
                    else:
                        raise Exception('Block division calculation not yet coded for m>2')
            I1 = np.nan * np.ones((ni, 1))
            # I2 = cell(1,nblocks)
            I2 = np.zeros(1, nblocks)
            xsm = I2
            radm = np.nan * np.ones((nblocks, m))
            ypf = I1
            ndsf = I1
            if check_evenness == 1:
                ndsdomf = I1
            Kf = I2
            wtf = I2
            recsf = I2
            if len(x_coverage) == 0 == 0:
                xp_coveragef = np.zeros((ni, 1))
            if nb > 1:
                bf = np.nan * np.ones((nb, ni))
            if calcunc == 1:
                rhatf = I2
                varhatf = I2
                varpf = I1
                mpsef = I1
                mpLf = I1
                mpHf = I1
                ypsef = I1
                ypLf = I1
                ypHf = I1

        # Add prior values as pseudodata at the prediction points to stabilize estimates if req'd
        if not len(yprior) == 0:
            y = np.array([[y], [yprior]])
            x = np.array([[x], [xp]])
            if not len(obserrv) == 0:
                obserrv = np.array([[obserrv], [priorerrv]])
            if not len(xc) == 0:
                raise Exception('Cannot use additional covariates if adding prior values as pseudodata')

        # Master loop of nblocks----------------------------------------------------
        for iblock in np.arange(1, nblocks + 1).reshape(-1):
            if nblocks > 1:
                if subset_wrt_xp == 0:
                    recp1 = np.arange(bblock(iblock) + 1, bblock(iblock + 1) + 1)
                if subset_wrt_xp == 1:
                    recp1 = np.arange(1, ni + 1)
                    for i in np.arange(1, m + 1).reshape(-1):
                        recp1 = recp1[xp[recp1, i] >= xpL[iblock, i] & xp[recp1, i] <= xpH[iblock, i]]
                        pass
                xp1 = xp[recp1, :]
                bdw1 = bdw[recp1, :]
                if len(xpc) == 0 == 0:
                    xp1c = xpc[recp1, :]
            else:
                xp1 = xp
                bdw1 = bdw
                if len(xpc) == 0 == 0:
                    xp1c = xpc
            ni1 = len(xp1[:, 0])
            if use_rad == 1:
                rad1 = rad
            if ni1 > 0:
                # Establish subset recs of influential data
                recs = np.argwhere(np.logical_and(not np.isnan(y), np.sum(np.isnan(x), 2 - 1)) == 0)
                # Adapt for one-sided kernels
                xp1min = np.amin(xp1, [], 1)
                xp1max = np.amax(xp1, [], 1)
                for i in np.arange(1, m + 1).reshape(-1):
                    if one_sided[m] == - 1:
                        recs = recs[x[recs, i] <= xp1max[i]]
                        pass
                    else:
                        if one_sided[m] == 1:
                            recs = recs[x[recs, i] >= xp1min[i]]
                            pass
                nds = len(recs)
                if use_rad == 1:
                    recs0 = recs
                    nrecs0 = nds
                    incrad = 0
                    flag_exit = 0
                    if nrecs0 < ndsmin:
                        if verbose > 1:
                            print(np.array(['ks_regress: Full data set smaller than target minimum size (nds = ',
                                            str(nds), '), hence using full set']))
                    else:
                        rad1 = rad
                        while (nds < ndsmin or nds > ndsmax) and incrad <= incradmax and flag_exit == 0:

                            recs = recs0
                            for i in np.arange(1, m + 1).reshape(-1):
                                if not np.isnan(period[i]):
                                    sels = fnperiodic_limit(x[recs, i], xp1min[i] - rad1[i],
                                                            xp1max[i] + rad1[i], period[i])
                                    recs = recs(sels)
                                else:
                                    if one_sided[m] == - 1:
                                        recs = recs(x[recs, i] >= xp1min[i] - rad1[i] & x[recs, i] <= xp1max[i])
                                        pass
                                    else:
                                        if one_sided[m] == 1:
                                            recs = recs(x[recs, i] >= xp1min[i] & x[recs, i] <= xp1max[i] + rad1[i])
                                            pass
                                        else:
                                            recs = recs(
                                                x[recs, i] >= xp1min[i] - rad1[i] & x[recs, i] <= xp1max[i] + rad1[i])
                                            pass
                                # Exclude data with any one coordinate beyond allowed distance from xp1min, xp1max
                            nds = len(recs)
                            if variable_rad == 0:
                                flag_exit = 1
                            else:
                                if nds < ndsmin:
                                    if np.any(2 * rad1 > radfacmax * np.amin(bdw1, [], 1)):
                                        if verbose > 1:
                                            print(np.array(
                                                ['ks_regress: Data set smaller that target minimum size (nds = ',
                                                 str(nds), ') but rad cannot be increased from ', str(rad1),
                                                 ' without exceeding (radfacmax=', str(radfacmax),
                                                 ')*bandwidth for one or more coords']))
                                        flag_exit = 1
                                    else:
                                        rad1 = 2 * rad1
                                        incrad = incrad + 1
                                else:
                                    if nds > ndsmax:
                                        if np.any(0.9 * rad1 < radfacmin * np.amax(bdw1, [], 1)):
                                            if verbose > 1:
                                                print(np.array(
                                                    ['ks_regress: Data set exceeds target maximum size (nds = ',
                                                     str(nds), ') but rad cannot be reduced from ', str(rad1),
                                                     ' without dropping below (radfacmin=', str(radfacmin),
                                                     ')*bandwidth for one or more coords']))
                                            flag_exit = 1
                                        else:
                                            rad1 = 0.9 * rad1
                                            incrad = incrad + 1
                            if incrad == incradmax and verbose > 1:
                                print(np.array(['ks_regress: Reached maximum subset adjustment iterations incradmax = ',
                                                str(incradmax), ', nds = ', str(nds)]))

                        if incrad > 0 and nds < nrecs0 and verbose > 1:
                            print(np.array(
                                ['ks_regress: Data set reduced from ', str(nrecs0), ' to ', str(nds), '; rad = ',
                                 str(rad1), '; incrad = ', str(incrad)]))
                if nds == 0:
                    if verbose > 1:
                        print('ks_regress: No data within allowed radius for kernel smoothing')
                    yp = np.nan * np.ones((ni, 1))
                if nds > 0:
                    xs = x[recs, :]
                    ys = y[recs]
                    if len(xc) == 0 == 0:
                        xcs = xc[recs, :]
                    if nblocks > 1:
                        xsm[iblock] = xs
                        if use_rad == 1:
                            radm[iblock, :] = rad1
                    if len(x_coverage) == 0 == 0:
                        x_coverages = x_coverage[recs]
                    if nds > 1000.0 and ni1 > 1000.0:
                        if verbose > 1:
                            print(np.array(
                                ['ks_regress: Warning: working with large matrices (ni1 = ', str(ni1), ', nds = ',
                                 str(nds), '), could be slow!']))
                    # Calculate total squared separation matrix (ssm)
                    ssm = np.zeros((ni1, nds))
                    iscartesian = np.ones((1, m))
                    # First treat spherical distances using (lat,lon) inputs
                    if len(aslatlon) == 0 == 0:
                        iscartesian[aslatlon] = 0
                        pi180 = np.pi / 180
                        R = 6378.137
                        lat1 = xp1[:, aslatlon[0]] * pi180
                        lon1 = xp1[:, aslatlon[1]] * pi180
                        lat2 = xs[:, aslatlon[0]] * pi180
                        lon2 = xs[:, aslatlon[1]] * pi180
                        coslat1coslat2 = np.cos(lat1) * np.transpose(np.cos(lat2))
                        dlat = lat1 * np.ones((1, nds)) - np.ones((ni1, 1)) * np.transpose(lat2)
                        dlon = lon1 * np.ones((1, nds)) - np.ones((ni1, 1)) * np.transpose(lon2)
                        distm1 = R * 2 * np.arcsin(
                            np.sqrt(np.sin(dlat / 2) ** 2 + np.multiply(coslat1coslat2, np.sin(dlon / 2) ** 2)))
                        # Haversine formula for distance in km, poached from code from m_lldist.m by Rich Pawlowicz (validated against same)
                        if one_sided[aslatlon[0]] == - 1:
                            distm1[dlat < 0] = np.inf
                        else:
                            if one_sided[aslatlon[0]] == 1:
                                distm1[dlat > 0] = np.inf
                        if one_sided[aslatlon[1]] == - 1:
                            distm1[dlon < 0] = np.inf
                        else:
                            if one_sided[aslatlon[0]] == 1:
                                distm1[dlon > 0] = np.inf
                        ssm = ssm + distm1 ** 2.0 / (bdw1[:, aslatlon[0]] ** 2 * np.ones((1, nds)))
                    # Loop over coordinates, adding scaled squared distance if cartesian
                    for i in np.arange(1, m + 1).reshape(-1):
                        if iscartesian[i] == 1:
                            if np.isnan(period[i]):
                                distm1 = xp1[:, i] * np.ones((1, nds)) - np.ones((ni1, 1)) * np.transpose(xs[:, i])
                                if one_sided[i] == - 1:
                                    distm1[distm1 < 0] = np.inf
                                else:
                                    if one_sided[i] == 1:
                                        distm1[distm1 > 0] = np.inf
                            if not np.isnan(period[i]):
                                distm1 = fnperiodic_distm(xp1[:, i], xs[:, i], period[i])
                            ssm = ssm + distm1 ** 2.0 / (bdw1[:, i] ** 2 * np.ones((1, nds)))
                    # Calculate kernel matrix from sum of squares matrix ssm
                    if kfunc == 1:
                        rhom = np.sqrt(ssm)
                        K = np.exp(- rhom)
                    if kfunc == 2:
                        rhom = ssm / 2
                        K = np.exp(- rhom)
                    if kfunc == 3:
                        rhom = np.sqrt(ssm)
                        K = np.sinc(rhom)
                    if kfunc == 4:
                        K = abs(ssm) <= 1
                    if kfunc == 5:
                        K = 1.0 / (1 + ssm)
                    # Correct for artifacts due to finite numerical precision
                    if kfunc == 1 or kfunc == 2:
                        minrho = np.amin(rhom, [], 2)
                        sel = np.argwhere(minrho > 700)
                        nsel = len(sel)
                        # Rows where min(rho)>700 require special treatment to allow for finite numerical precision
                        if nsel > 0:
                            if verbose > 1:
                                print('Some rows of weight matrix modified to allow for finite numerical precision')
                            if kfunc == 1:
                                print('Bad numerics for kfunc = 1, consider increasing bdw')
                                # stop
                            if kfunc == 2:
                                for i in np.arange(1, nsel + 1).reshape(-1):
                                    closest = np.argwhere(ssm[sel[i], :] == np.amin(ssm[sel[i], :]))
                                    K[sel[i], :] = np.zeros((1, nds))
                                    K[sel[i], closest] = 1 / len(closest)
                                    # Usually, this means that the closest datum takes all the weight
                    # Additional weighting for observational error
                    if not len(obserrv) == 0:
                        K = K / (np.ones((ni1, 1)) * np.reshape(obserrv[recs], 1, nds))
                    # Limit to closest data if req'd (only if smoothing to one prediction point)
                    if ni1 == 1 and len(ndclosest) == 0 == 0:
                        # Ks = __builtint__.sorted(K,'descend')
                        Ks = np.sort(K)[::-1]
                        Kmin = Ks(np.amin(len(Ks), ndclosest))
                        selclosest = np.argwhere(K >= Kmin)
                        K = K[selclosest]
                        recs = recs[selclosest]
                        nds = len(recs)
                        xs = xs[selclosest, :]
                        ys = ys[selclosest]
                        if len(xc) == 0 == 0:
                            xcs = xcs[selclosest, :]
                        if nblocks > 1:
                            xsm[iblock] = xs
                        if len(x_coverage) == 0 == 0:
                            x_coverages = x_coverages(selclosest)
                    # Normalize kernel so that rows sum to 1
                    K = K / (np.sum(K, 2 - 1) * np.ones((1, nds)))
                    # Check evenness of weight distribution
                    if check_evenness == 1:
                        # Ks = __builtint__.sorted(K,2,'descend')
                        Ks = np.sort(K, 2)[::-1]
                        Ks = np.cumsum(Ks, 2)
                        ndsdom = np.nan * np.ones((ni1, 1))
                        for i in np.arange(1, ni1 + 1).np.reshape(-1):
                            if sum(not np.isnan(Ks[i, :])) > 0:
                                ndsdom[i] = np.argwhere(Ks[i, :] > domfrac, 1, 'first')
                                if ndsdom[i] <= 2 and ndsdom[i] <= 0.01 * nds and verbose > 1:
                                    print(np.array(['ks_regress: Warning: kernel dominated by only ',
                                                    str(ndsdom(i)), ' data (= ',
                                                    str(100 * ndsdom(i) / nds), '% of all data)']))
                    # Calculate weight matrix for linear prediction yp = wt*ys
                    if nb == 1:
                        wt = K
                    if nb > 1:
                        wt = np.nan * np.ones((ni1, nds))
                        b = np.nan * np.ones((nb, ni1))
                        for i in np.arange(1, ni1 + 1).np.reshape(-1):
                            X1 = np.array([np.ones((nds, 1)), np.nan * np.ones((nds, nb - 1))])
                            Xp1 = np.array([1, np.zeros((1, nb - 1))])
                            q = 1
                            if nb1 > 0:
                                for j in np.arange(1, m + 1).np.reshape(-1):
                                    for k in np.arange(1, lr_order(j) + 1).np.reshape(-1):
                                        X1[:, q + 1] = (xs[:, j] - xp1[i, j] * np.ones((nds, 1))) ** k
                                        q = q + 1
                            if nbc > 0:
                                for j in np.arange(1, mc + 1).np.reshape(-1):
                                    for k in np.arange(1, lr_orderc(j) + 1).np.reshape(-1):
                                        X1[:, q + 1] = (xcs[:, j] - xp1c[i, j] * np.ones((nds, 1))) ** k
                                        q = q + 1
                                # if lr_orderc>0; X1(:,q+1:q+mc) = xcs - ones(nds,1)*xp1c(i,:); end
                            W1 = np.diag(K[i, :])
                            R1 = np.linalg.solve((np.transpose(X1) * W1 * X1), (np.transpose(X1) * W1))
                            wt[i, :] = Xp1 * R1
                            b[:, i] = R1 * ys
                            # size(R1)
                    # i ni1]
                    # If x_coverage supplied, calculate xp_coverage for these prediction points
                    if len(x_coverage) == 0 == 0:
                        xp_coverage = np.nan * np.ones((ni1, 1))
                        if coverage_type == 0:
                            for i in np.arange(1, ni1 + 1).np.reshape(-1):
                                sel = np.argwhere(wt[i, :] > np.logical_and(fr_wt_coverage * np.mean(wt[i, :]),
                                                                            ssm[i, :]) < max_sdist_coverage ** 2)
                                if len(sel) == 0 == 0:
                                    xp_coverage[i] = range(x_coverages[sel])
                                else:
                                    xp_coverage[i] = 0
                        if coverage_type == 1:
                            for i in np.arange(1, ni1 + 1).np.reshape(-1):
                                sel = np.argwhere(wt[i, :] > np.logical_and(fr_wt_coverage * np.mean(wt[i, :]),
                                                                            ssm[i, :]) < max_sdist_coverage ** 2)
                                xp_coverage[i] = len(np.unique(int(np.floor(x_coverages[sel]))))
                                # Note: this will return 0 if sel if empty
                    # Calculate smoothing estimates
                    yp = wt * ys
                    if nneg == 1:
                        yp = np.amax(0, yp)
                    # Calculate uncertainties if req'd (takes more time)
                    if calcunc == 1:
                        # Assuming: y = m(x) + V(x)*eps
                        # First calculate residuals by smoothing to the data
                        optv = opt
                        optv.lr_order = vlr_order
                        optv.calcunc = 0
                        optv.nblocks = 1
                        yhat, out2 = ks_regress(ys, xs, xs, vbdw, optv)
                        Khat = out2.K
                        rhat = ys - yhat
                        rhat2 = rhat ** 2
                        if variable_res_var == 1:
                            # Estimate V(x) using the normalized weighted RSS (see rlpoly_Stata.pdf, p10)
                            wthat = Khat / (np.sum(Khat, 2 - 1) * np.ones((1, nds)))
                            varhat = wthat * rhat2
                            ndshatv = np.sum(wthat > 0, 2 - 1)
                            varhat = np.multiply(varhat, ndshatv) / np.amax(1, (ndshatv - 1))
                            # Estimate V(xp1) similarly
                            wtvarp = K / (np.sum(K, 2 - 1) * np.ones((1, nds)))
                            varp = wtvarp * rhat2
                            ndsvarpv = np.sum(wtvarp > 0, 2 - 1)
                            varp = np.multiply(varp, ndsvarpv) / np.amax(1, (ndsvarpv - 1))
                        else:
                            var0hat = sum(rhat2) / (nds - 1)
                            varhat = var0hat * np.ones((nds, 1))
                            ndsvarpv = nds * np.ones((ni1, 1))
                            varp = var0hat * np.ones((ni1, 1))
                        # mp = wt*ys => C[mp] = wt*C[ys]*wt' = wt*V(x)*wt' assuming uncorrelated errors
                        Vhat = np.diag(varhat)
                        Cmp = wt * Vhat * np.transpose(wt)
                        mp = yp
                        mpvar = np.diag(Cmp)
                        mpse = np.sqrt(mpvar)
                        # tfac = tinv(1 - alph / 2,ndsvarpv - 1)
                        tfac = t.ppf(1 - alph / 2, ndsvarpv - 1)
                        mpL = mp - np.multiply(tfac, mpse)
                        mpH = mp + np.multiply(tfac, mpse)
                        ypse = np.sqrt(mpvar + varp)
                        ypL = yp - np.multiply(tfac, ypse)
                        ypH = yp + np.multiply(tfac, ypse)
                        if nneg == 1:
                            mpL = np.amax(0.0, mpL)
                            ypL = np.amax(0.0, ypL)
                        if nb > 1:
                            # stop
                            pass
                    if nblocks > 1:
                        ypf[recp1] = yp
                        ndsf[recp1] = nds
                        if check_evenness == 1:
                            ndsdomf[recp1] = ndsdom
                        Kf[iblock] = K
                        wtf[iblock] = wt
                        recsf[iblock] = recs
                        if len(x_coverage) == 0 == 0:
                            xp_coveragef[recp1] = xp_coverage
                        if nb > 1:
                            bf[:, recp1] = b
                        if calcunc == 1:
                            rhatf[iblock] = rhat
                            varhatf[iblock] = varhat
                            varpf[recp1] = varp
                            mpsef[recp1] = mpse
                            mpLf[recp1] = mpL
                            mpHf[recp1] = mpH
                            ypsef[recp1] = ypse
                            ypLf[recp1] = ypL
                            ypHf[recp1] = ypH

                    if doplotc==1 and isubc<=nsubc:
                        if np.mod(kc, kcstep) == 0 or np.amin(
                                np.sum((xpp - np.ones((npp, 1)) * xp1[selpp, :]) ** 2, 2 - 1)) == 0:
                            Xcpp = np.array([np.ones((ncpp, 1)), np.nan * np.ones((ncpp, nb - 1))])
                            q = 1
                            if lr_order > 0:
                                Xcpp[:, np.arange[q + 1, q + m + 1]] = np.zeros((ncpp, m))
                                q = q + m
                            if lr_orderc > 0:
                                Xcpp[:, np.arange[q + 1, q + mc + 1]] = xcpp
                            ypp = Xcpp * b[:, selpp]
                            selcs = np.argwhere(ssm[selpp, :] < max_sdist_plotc ** 2)
                            plt.subplot(nrowsc, ncolsc, isubc)
                            plt.plot(xcs[selcs, 0] + xcoff, ys[selcs], '.', 'Color', ycol, 'MarkerSize', ymarkersize)
                            # hold('on')
                            plt.plot(xcpp + xcoff, ypp, 'Color', ypcol, 'LineStyle', '-', 'LineWidth', linewidth)
                            plt.xlabel(xstr, 'FontSize', lfontsize)
                            plt.ylabel(ystr, 'FontSize', lfontsize)
                            plt.title(np.array(['x = ', str(xp1(selpp, 1)), ', y = ', str(xp1(selpp, 2))]), 'FontSize',
                                      lfontsize)
                            if len(xlimc) == 0 == 0:
                                # set(gca,'XLim',xlimc)
                                plt.set_xlim(xlimc)
                            if len(ylimc) == 0 == 0:
                                # set(gca,'YLim',ylimc)
                                plt.set_ylim(ylimc)
                            isubc = isubc + 1
                    kc = kc + 1
            if np.mod(iblock, 10000.0) == 0 and verbose > 0:
                print(np.array(['Done ', str(iblock), ' of ', str(nblocks), ' blocks']))

        if nblocks > 1:
            yp = ypf
            nds = ndsf
            if check_evenness == 1:
                ndsdom = ndsdomf
            if len(x_coverage) == 0 == 0:
                xp_coverage = xp_coveragef
            if nb > 1:
                b = bf
            if calcunc == 1:
                rhat = rhatf
                varp = varpf
                mpse = mpsef
                mpL = mpLf
                mpH = mpHf
                ypse = ypsef
                ypL = ypLf
                ypH = ypHf
        else:
            nds = nds * np.ones((ni, 1))

        # #Plot if req'd
        # if doplot == 1:
        #     plt.figure(ifig)
        #     plt.subplot(nrows,ncols,isub)
        #     if m == 1:
        #         if plottype == 0:
        #             plt.plot(x, y, 'k.', xp, yp, 'r-')
        #             if calcunc == 1:
        #                 # hold('on')
        #                 plt.plot(xp, yp - yppe, 'r--', xp, yp - yppe,'r--')
        #             plt.xlabel('x')
        #             plt.ylabel('y')
        #         if plottype == 1:
        #             plt.plot(y,x,'k.',yp,xp,'r-')
        #             if calcunc == 1:
        #                 # hold('on')
        #                 plt.plot(yp - yppe,xp,'r--',yp + yppe,xp,'r--')
        #             if len(yp) == 1:
        #                 # hold('on')
        #                 plt.plot(yp,xp,'ro')
        #                 if calcunc == 1:
        #                     # hold('on')
        #                     errorbar_x(yp,xp,yppe,yppe,'o',struct('color','r'))
        #             plt.xlabel('d')
        #             plt.ylabel('z')
        #             plt.axis('ij')
        #             if len(tstr)==0 == 0:
        #                 plt.title(tstr,'FontSize',10)
        #
        # #Store auxiliary results in output structure
        # if nargout > 1 and np.amax(nds) > 0:
        #     out.nblocks = nblocks
        #     if nblocks > 1:
        #         out.xs = xsm
        #         if use_rad == 1:
        #             out.rad = radm
        #     else:
        #         out.xs = xs
        #         if use_rad == 1:
        #             out.rad = rad1
        #     out.K = K
        #     out.nds = nds
        #     out.wt = wt
        #     out.recs = recs
        #     if check_evenness == 1:
        #         out.ndsdom = ndsdom
        #     if len(x_coverage)==0 == 0:
        #         out.xp_coverage = xp_coverage
        #     if nb > 1:
        #         out.b = b
        #     if calcunc == 1:
        #         out.rhat = rhat
        #         out.varhat = varhat
        #         out.varp = varp
        #         out.mpse = mpse
        #         out.mpL = mpL
        #         out.mpH = mpH
        #         out.ypse = ypse
        #         out.ypL = ypL
        #         out.ypH = ypH
        # else:
        #     out = []

    return out
