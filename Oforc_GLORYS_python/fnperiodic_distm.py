import numpy as np


def fnperiodic_distm(xp, x, period):
    # function distm = fnperiodic_distm(xp,x,period)
    # Calculate the distance matrix between xp (rows) and x (columns) over a periodic dimension.
    # Works for lower-limit integers (e.g. days since 1st January, yrday = 0-364, period = 365),
    # upper-limit integers (e.g. month of the year, month = 1-12, period = 12),
    # and continuous variables (e.g. hours since start of year, period = 12*365).
    # #Phil Wallhead 10/02/2019
    # Example 1: distance from yrday = 364, where yrday = days since 1st January
    # distm = fnperiodic_distm(364,0:364,365) #Assuming a non-leap year
    # Example 2: distance from month = 12 (December), where month = month of the year (1-12)
    # distm = fnperiodic_distm(12,1:12,12)
    # Example 3: distance from hr = 364*12+11, where hr = hours since start of year
    # distm = fnperiodic_distm(364*12+11,0:(365*12-1),365*12) #Assuming a non-leap year

    if len(x.shape) == 0:
        lx = 1
    else:
        lx = x.shape[0]
    if len(xp.shape) == 0:
        lxp = 0
    else:
        lxp = xp.shape[0]
    M1 = np.mod(xp * np.ones((1, lx)), period)[0, :]
    M2 = np.mod(np.ones((lxp, 1)) * np.transpose(x), period)[:, 0]
    xmax = np.maximum(M1, M2)
    xmin = np.minimum(M1, M2)
    distm = np.minimum(xmax - xmin, xmin + period - xmax)
    return distm

