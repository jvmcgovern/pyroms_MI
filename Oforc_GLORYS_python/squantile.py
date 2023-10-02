import numpy as np


def squantile(X=None, P=None):
    # function Y = squantile(X,P)

    # Simple quantile function that gives results consistent with quantile.m.
    # For convenience when quantile.m from stats toolbox is lacking.
    # NOTE: X must be a vector for squantile.m (can be a matrix for quantile.m).

    # X = vector or matrix of data [n1*n2]
    # P = vector of quantiles [np*1] or [1*np].
    # Y = output quantiles [np*n2], or [1*np] if X is vector and P is [1*np]    (same behaviour as quantile.m)

    # Uses: Xplin.m

    # Phil Wallhead 26/09/2023

    # Example tests of squantile.m vs. quantile.m:
    # x = rand(1000,10); qsel = 0:0.01:1; y = squantile(x,qsel); yc = quantile(x,qsel); max(abs(y(:)-yc(:))) #<1e-15
    # x = rand(1000,1); qsel = 0:0.01:1; y = squantile(x,qsel); yc = quantile(x,qsel); max(abs(y(:)-yc(:))) #<1e-15
    # x = rand(1000,10); qsel = (0:0.01:1)'; y = squantile(x,qsel); yc = quantile(x,qsel); max(abs(y(:)-yc(:))) #<1e-15

    np1, np2 = P.shape
    if np1 > 1 and np2 > 1:
        raise Exception('Second input must be a scalar or a non-empty vector')

    P = P
    np = len(P)
    n1, n2 = X.shape
    if n1 == 1:
        X = X
        n2c = 1
    else:
        n2c = n2

    Y = np.nan * np.ones((np, n2c))
    for i in np.arange(1, n2c + 1).reshape(-1):
        X1 = X[:, i]
        Xs1 = __builtint__.sorted(X1(not np.isnan(X1)))
        ns1 = len(Xs1)
        if ns1 > 0:
            Pn = np.transpose((np.arange(0.5, ns1 - 0.5 + 1, 1))) / ns1
            Y[:, i] = Xplin(P, Pn, 1) * Xs1

    if n2c == 1 and np2 > 1:
        Y = np.transpose(Y)

    return Y
