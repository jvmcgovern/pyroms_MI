import glob
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import ScalarFormatter
from netCDF4 import Dataset as netcdf
from netCDF4 import num2date as n2d
from netCDF4 import date2index as d2i
from datetime import date, datetime
from scipy import interpolate
import numpy.matlib
import PyCO2SYS as pyco2
import time
import pickle
from pyhdf.SD import SD, SDC
import pandas as pd
from interp_Cgrid import *
import croco_vgrid as vgrd
import croco_glorys as glor
from math import cos, sin, asin, sqrt, radians

# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from scipy.interpolate import griddata
import scipy.io as scio
# from netCDF4 import date2index as d2i
# from netCDF4 import num2date as n2d
# from datetime import date, datetime
# from calendar import monthrange
# import cftime
# import numpy as np
# import numpy.ma as ma
# import croco_vgrid as vgrd
from progressbar import *

# import sys
# sys.path.insert(0,'')

floc = '/media/dskthree/atm_CO2_ssp/future/CMIP6GHGConcentrationProjections_1_2_1/'
atm119 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-IMAGE-ssp119-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm126 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm245 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm370 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm370_lowNTCF = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
                 'AerChemMIP_UoM-AIM-ssp370-lowNTCF-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm434 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-GCAM4-ssp434-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm460 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-GCAM4-ssp460-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm534_over = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
              'ScenarioMIP_UoM-REMIND-MAGPIE-ssp534-over-1-2-1_gr1-GMNHSH_201501-250012.nc'
atm585 = 'mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_' + \
         'ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_201501-250012.nc'

co2f_119 = netcdf(floc + atm119, 'r')
co2f_119o = co2f_119.variables
time = np.array(co2f_119o['time'][:])
Tini = datetime(2014, 12, 1)
Tend = datetime(2100, 11, 1)
Tiniinx1 = d2i(Tini, co2f_119o['time'], select='nearest')
Tendinx1 = d2i(Tend, co2f_119o['time'], select='nearest')
co2_119 = np.array(co2f_119o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_119.close()

co2f_126 = netcdf(floc + atm126, 'r')
co2f_126o = co2f_126.variables
co2_126 = np.array(co2f_126o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_126.close()

co2f_245 = netcdf(floc + atm245, 'r')
co2f_245o = co2f_245.variables
co2_245 = np.array(co2f_245o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_245.close()

co2f_370 = netcdf(floc + atm370, 'r')
co2f_370o = co2f_370.variables
co2_370 = np.array(co2f_370o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_370.close()

co2f_370_NTCF = netcdf(floc + atm370_lowNTCF, 'r')
co2f_370_NTCFo = co2f_370_NTCF.variables
co2_370_NTCF = np.array(co2f_370_NTCFo['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_370_NTCF.close()

co2f_434 = netcdf(floc + atm434, 'r')
co2f_434o = co2f_434.variables
co2_434 = np.array(co2f_434o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_434.close()

co2f_460 = netcdf(floc + atm460, 'r')
co2f_460o = co2f_460.variables
co2_460 = np.array(co2f_460o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_460.close()

co2f_534over = netcdf(floc + atm534_over, 'r')
co2f_534overo = co2f_534over.variables
co2_534over = np.array(co2f_534overo['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_534over.close()

co2f_585 = netcdf(floc + atm585, 'r')
co2f_585o = co2f_585.variables
co2_585 = np.array(co2f_585o['mole_fraction_of_carbon_dioxide_in_air'][Tiniinx1:Tendinx1, 0])
co2f_585.close()
co2_concs = pd.DataFrame({'co2_119': co2_119,
                        'co2_126': co2_126,
                        'co2_245': co2_245,
                        'co2_370': co2_370,
                        'co2_370_NTCF': co2_370_NTCF,
                        'co2_434': co2_434,
                        'co2_460': co2_460,
                        'co2_534over': co2_534over,
                        'co2_585': co2_585})
co2_concs.to_excel(floc + 'co2_atm_ssps.xlsx')
