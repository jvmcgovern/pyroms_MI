import numpy as np


def Xplin(t, tn, extrapconst):
    # function X = Xplin(t,tn,extrapconst=[0 0])

    # Builds 1D linear interpolation matrix from tn to t.
    # Nodes tn MUST all not-NaN and in increasing order.
    # Input extrapconst = [1 1] to extrapolate constant beyond [first last] nodes
    # (otherwise defaults to linear extrapolation, constant if only one node).
    # Any missing target points (t=NaN) will result in NaN rows in X, unless
    # there is only one node point, in which case a constant is extrapolated.

    # Example:
    # x = 1:9; y = rand(9,1); xp = 0:0.1:10;
    # tic;Xz0 = Xplin(xp,x,[0 1]); Xz1 = Xplin(xp,x,[1 0]); toc
    # figure; plot(x,y,'ko',xp,Xz0*y,'k-',xp,Xz1*y,'r-')

    # Phil Wallhead 23/06/2021
    # Joe McGovern 31/08/2023 - converting to python 3

    if sum(np.isnan(tn)) > 0:
        raise Exception('One or more node values tn = NaN in Xplin.py')
    else:
        if np.amin(np.diff(tn)) <= 0:
            raise Exception('Nodes tn MUST be in increasing order')

    if not extrapconst:
        extrapconst = np.array([0, 0])

    # if extrapconst :
    if not isinstance(extrapconst, np.ndarray):
        extrapconst = extrapconst * np.ones((1, 2))

    m = len(tn) - 1
    if m == 1:
        X = np.ones((len(t), 1))
    else:
        if m > 1:
            X = np.nan * np.ones((len(t), m + 1))
            # sel = find(not np.isnan(t))
            sel = np.argwhere(t != np.nan)[:, 0]
            X[sel, :] = 0
            for i in np.arange(0, len(sel)).reshape(-1):
                if t[sel[i]] <= tn[0]:
                    if extrapconst[0] == 1:
                        X[sel[i], 0] = 1.0
                    else:
                        X[sel[i], 0] = 1 - (t[sel[i]] - tn[0]) / (tn[1] - tn[0])
                        X[sel[i], 1] = 1 - X(sel(i), 1)
                else:
                    if t[sel[i]] > tn[m]:
                        if extrapconst[0] == 1:
                            X[sel[i], m] = 1.0
                        else:
                            X[sel[i], m - 1] = 1 - (t[sel[i]] - tn[m - 1]) / (tn[m] - tn[m - 1])
                            X[sel[i], m] = 1 - X[sel[i], m - 1]
                    else:
                        n1 = np.argwhere(tn < t[sel[i]])[:, 0][-1]

                        X[sel[i], n1] = 1 - (t[sel[i]] - tn[n1]) / (tn[n1 + 1] - tn[n1])
                        X[sel[i], n1 + 1] = 1 - X[sel[i], n1]
        else:
            X = []

    return X
