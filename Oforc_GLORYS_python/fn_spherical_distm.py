import numpy as np


def fn_spherical_distm(lat1, lon1, lat2, lon2, *args):
    # function distm = fn_spherical_distm(lat1,lon1,lat2,lon2,radsin=0,R=6378.137km)
    # Calculates matrix distm of spherical distances in km given vectors (lat1, lon1) of length n1 and
    # vectors (lat2, lon2) of length n2. Input of degrees is assumed unless radsin = 1.
    
    # Phil Wallhead 16/03/2015
    # Joe McGovern 26/09/2023 -- convert to python
    
    lat1 = lat1[:]
    lon1 = lon1[:]
    lat2 = lat2[:]
    lon2 = lon2[:]
    n1 = len(lat1)
    n2 = len(lat2)

    varargin = args
    if len(varargin) == 0:
        radsin = 0
        R = 6378.137
    
    elif len(varargin) == 1:
        radsin = args
        R = 6378.137

    else:
        radsin, R = args
    
    if radsin == 1:
        lat1r = lat1
        lon1r = lon1
        lat2r = lat2
        lon2r = lon2
    else:
        pi180 = np.pi / 180
        lat1r = lat1 * pi180
        lon1r = lon1 * pi180
        lat2r = lat2 * pi180
        lon2r = lon2 * pi180
    
    coslat1rcoslat2r = np.cos(lat1r) * np.transpose(np.cos(lat2r))
    dlatr = lat1r * np.ones((1, n2)) - np.ones((n1, 1)) * np.transpose(lat2r)
    dlonr = lon1r * np.ones((1, n2)) - np.ones((n1, 1)) * np.transpose(lon2r)
    distm = R * 2 * np.arcsin(np.sqrt(np.sin(dlatr / 2) ** 2 +
                                      np.multiply(coslat1rcoslat2r, np.sin(dlonr / 2) ** 2)))
    return distm
