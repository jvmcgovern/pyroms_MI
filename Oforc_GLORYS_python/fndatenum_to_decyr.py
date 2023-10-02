
from datetime import date, datetime

import numpy as np


def fndatenum_to_decyr(td, *args):
    # function decyr = fndatenum_to_decyr(td,distort=0,nyrdays0=365.25,year0=2000)

    # Converts datenumber (from datenum.m) to decimal years
    # If distort = 1, time is distorted to keep exact agreement between year and floor(decyr)
    # If distort = 0 (def), time is not distorted but a fixed no. of days in the year is
    # assumed (nyrdays0) and years are computed with respect to a reference year year0
    # See also fndecyr_to_datenum.m
    # Phil Wallhead 13/04/2017
    # Joe McGovern 26/09/2023

    varargin = args
    if len(varargin) == 0:
        distort = 0
        nyrdays0 = 365.25
        year0 = 2000

    elif len(varargin) == 1:
        distort = args
        nyrdays0 = 365.25
        year0 = 2000

    elif len(varargin) == 2:
        distort, nyrdays0 = args
        year0 = 2000

    else:
        distort, nyrdays0, year0 = args

    decyr = np.zeros_like(td)
    for ti in np.arange(0, len(td)):
        dfo = date.fromordinal(td)

        td0 = date.toordinal(date(dfo.year, 1, 1))
        if distort == 1:
            nyrdays = date.toordinal(date(dfo.year + 1, 1, 1)) - td0
            decyr[ti] = dfo.year + (td - td0) / nyrdays
            # Time is distorted to guarantee that floor(decyr) = year
        else:
            decyr[ti] = year0 + (td - date.toordinal(date(year0, 1, 1))) / nyrdays0
            # Here time is not distorted, but floor(decyr) may not always equal year

    return decyr
