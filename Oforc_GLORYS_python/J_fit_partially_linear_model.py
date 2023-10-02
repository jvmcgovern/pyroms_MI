import numpy as np
from scipy.optimize import nnls


def J_fit_partially_linear_model(theta, *args):
    # function [J,out] = J_fit_partially_linear_model(theta,opt)

    # Calculates J = -2log(likelihood) for partially linear model:

    #            Y = X(theta)*b + eps

    # where the linear matrix of covariates is a nonlinear function of nonlinear
    # parameters theta, and is specified by an input function opt.Xfun.
    # Normal errors eps are assumed.

    # Phil Wallhead 03/09/2021
    # Joe McGovern 26/09/2023 -- convert to python

    if len(args) == 0:
        opt = []
    else:
        opt = args[0]

        if hasattr(opt, 'Y') == 1:
            Y = opt.Y
        else:
            Y = []

        if hasattr(opt, 'Xfun') == 1:
            Xfun = opt.Xfun
        else:
            Xfun = []

        if hasattr(opt, 'nonnegb') == 1:
            nonnegb = opt.nonnegb
        else:
            nonnegb = 0

        if hasattr(opt, 'Vs') == 1:
            Vs = opt.Vs
        else:
            Vs = []

    # if hasattr(opt,'sumlogVs')==1; sumlogVs = opt.sumlogVs; else sumlogVs = []; end
    # Supplying precalculated sum(log(Vs)) doesn't seem to improve speed

    n = len(Y)
    X = feval(Xfun, theta)
    if not len(Vs) == 0:
        R2 = np.diag(np.sqrt(1.0 / Vs))
        Xs = R2 * X
        Ys = R2 * Y
        if nonnegb == 1:
            b = nnls(Xs, Ys)
        else:
            b = np.linalg.solve(Xs, Ys)
        Ym = X * b
        VhatMLE = sum(((Y - Ym) ** 2) / Vs) / n
        # if isempty(sumlogVs); sumlogVs = sum(log(Vs)); end
        sumlogVs = sum(np.log(Vs))
        J = n + sumlogVs + n * np.log(VhatMLE)
    else:
        if nonnegb == 1:
            b = nnls(X, Y)
        else:
            b = np.linalg.solve(X, Y)
        Ym = X * b
        VhatMLE = sum((Y - Ym) ** 2) / n
        J = n + n * np.log(VhatMLE)

    out = []
    if nargout > 1:
        out.X = X
        out.b = b
        out.Ym = Ym
        out.VhatMLE = VhatMLE
        out.SSQ = n * VhatMLE
        out.Vhat = VhatMLE * n / (n - len(b))

    return J, out
