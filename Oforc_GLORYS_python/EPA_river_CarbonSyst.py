#
######################################################################
######################################################################
#
#  Main program
#
#  Build a CROCO boundary file using GLORYS12 renanalysis data
#
######################################################################
######################################################################
#

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import num2date
from scipy.interpolate import griddata

import glob
from netCDF4 import Dataset as netcdf
from netCDF4 import date2index as d2i
from calendar import monthrange
from datetime import date, datetime
import PyCO2SYS as pyco2
import pandas as pd

import sys

# sys.path.insert(0,'/XXX/')

import croco_vgrid as vgrd
import croco_glorys as glor
from interp_Cgrid import *
from progressbar import *
from scipy.spatial import Delaunay

# import excel sheets
# point to relevant information via dataframe
# apply CO2sys to data to determine DIC
# save to excel sheet (same file?)
rivdir = '/media/dskone/CELTIC/RIVERS/'
Avoca = 'Avoca_EPA_FW.xls'
Bandon = 'Bandon_EPA_FW.xls'
Barrow = 'Barrow_EPA_FW.xls'
Blackwater = 'Blackwater_EPA_FW.xls'
Feale = 'Feale_EPA_FW.xls'
Galey = 'Galey_EPA_FW.xls'
Lee = 'Lee_EPA_FW.xls'
Maigue = 'Maigue_EPA_FW.xls'
Nore = 'Nore_EPA_FW.xls'
Shannon = 'Shannon_EPA_FW.xls'
Slaney = 'Slaney_EPA_FW.xls'
Suir = 'Suir_EPA_FW.xls'

Avo = pd.read_excel(rivdir + Avoca, sheet_name="Alk_pH_temp")
Avo_alk = np.array(Avo['Alkalinity'])/0.05
Avo_pH = np.array(Avo['pH'])
Avo_temp = np.array(Avo['Temperature'])
Avo_results = pyco2.sys(par1=Avo_alk, par1_type=1,
                        par2=Avo_pH, par2_type=3,
                        temperature=Avo_temp,
                        salinity=np.zeros_like(Avo_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Avo.insert(Avo.shape[1], "DIC", Avo_results["dic"], True)
Avo_DIC = Avo_results["dic"]
Avo.to_excel(rivdir + Avoca[:-4] + 'DIC.xls')

Ban = pd.read_excel(rivdir + Bandon, sheet_name="Alk_pH_temp")
Ban_alk = np.array(Ban['Alkalinity'])/0.05
Ban_pH = np.array(Ban['pH'])
Ban_temp = np.array(Ban['Temperature'])
Ban_results = pyco2.sys(par1=Ban_alk, par1_type=1,
                        par2=Ban_pH, par2_type=3,
                        temperature=Ban_temp,
                        salinity=np.zeros_like(Ban_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Ban.insert(Ban.shape[1], "DIC", Ban_results["dic"], True)
Ban_DIC = Ban_results["dic"]
Ban.to_excel(rivdir + Bandon[:-4] + 'DIC.xls')

Bar = pd.read_excel(rivdir + Barrow, sheet_name="Alk_pH_temp")
Bar_alk = np.array(Bar['Alkalinity'])/0.05
Bar_pH = np.array(Bar['pH'])
Bar_temp = np.array(Bar['Temperature'])
Bar_results = pyco2.sys(par1=Bar_alk, par1_type=1,
                        par2=Bar_pH, par2_type=3,
                        temperature=Bar_temp,
                        salinity=np.zeros_like(Bar_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Bar.insert(Bar.shape[1], "DIC", Bar_results["dic"], True)
Bar_DIC = Bar_results["dic"]
Bar.to_excel(rivdir + Barrow[:-4] + 'DIC.xls')

Bla = pd.read_excel(rivdir + Blackwater, sheet_name="Alk_pH_temp")
Bla_alk = np.array(Bla['Alkalinity'])/0.05
Bla_pH = np.array(Bla['pH'])
Bla_temp = np.array(Bla['Temperature'])
Bla_results = pyco2.sys(par1=Bla_alk, par1_type=1,
                        par2=Bla_pH, par2_type=3,
                        temperature=Bla_temp,
                        salinity=np.zeros_like(Bla_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Bla.insert(Bla.shape[1], "DIC", Bla_results["dic"], True)
Bla_DIC = Bla_results["dic"]
Bla.to_excel(rivdir + Blackwater[:-4] + 'DIC.xls')

Fea = pd.read_excel(rivdir + Feale, sheet_name="Alk_pH_temp")
Fea_alk = np.array(Fea['Alkalinity'])/0.05
Fea_pH = np.array(Fea['pH'])
Fea_temp = np.array(Fea['Temperature'])
Fea_results = pyco2.sys(par1=Fea_alk, par1_type=1,
                        par2=Fea_pH, par2_type=3,
                        temperature=Fea_temp,
                        salinity=np.zeros_like(Fea_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Fea.insert(Fea.shape[1], "DIC", Fea_results["dic"], True)
Fea_DIC = Fea_results["dic"]
Fea.to_excel(rivdir + Feale[:-4] + 'DIC.xls')

Gal = pd.read_excel(rivdir + Galey, sheet_name="Alk_pH_temp")
Gal_alk = np.array(Gal['Alkalinity'])/0.05
Gal_pH = np.array(Gal['pH'])
Gal_temp = np.array(Gal['Temperature'])
Gal_results = pyco2.sys(par1=Gal_alk, par1_type=1,
                        par2=Gal_pH, par2_type=3,
                        temperature=Gal_temp,
                        salinity=np.zeros_like(Gal_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Gal.insert(Gal.shape[1], "DIC", Gal_results["dic"], True)
Gal_DIC = Gal_results["dic"]
Gal.to_excel(rivdir + Galey[:-4] + 'DIC.xls')

Leef = pd.read_excel(rivdir + Lee, sheet_name="Alk_pH_temp")
Leef_alk = np.array(Leef['Alkalinity'])/0.05
Leef_pH = np.array(Leef['pH'])
Leef_temp = np.array(Leef['Temperature'])
Leef_results = pyco2.sys(par1=Leef_alk, par1_type=1,
                        par2=Leef_pH, par2_type=3,
                        temperature=Leef_temp,
                        salinity=np.zeros_like(Leef_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Leef.insert(Leef.shape[1], "DIC", Leef_results["dic"], True)
Leef_DIC = Leef_results["dic"]
Leef.to_excel(rivdir + Lee[:-4] + 'DIC.xls')

Mai = pd.read_excel(rivdir + Maigue, sheet_name="Alk_pH_temp")
Mai_alk = np.array(Mai['Alkalinity'])/0.05
Mai_pH = np.array(Mai['pH'])
Mai_temp = np.array(Mai['Temperature'])
Mai_results = pyco2.sys(par1=Mai_alk, par1_type=1,
                        par2=Mai_pH, par2_type=3,
                        temperature=Mai_temp,
                        salinity=np.zeros_like(Mai_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Mai.insert(Mai.shape[1], "DIC", Mai_results["dic"], True)
Mai_DIC = Mai_results["dic"]
Mai.to_excel(rivdir + Maigue[:-4] + 'DIC.xls')

Nor = pd.read_excel(rivdir + Nore, sheet_name="Alk_pH_temp")
Nor_alk = np.array(Nor['Alkalinity'])/0.05
Nor_pH = np.array(Nor['pH'])
Nor_temp = np.array(Nor['Temperature'])
Nor_results = pyco2.sys(par1=Nor_alk, par1_type=1,
                        par2=Nor_pH, par2_type=3,
                        temperature=Nor_temp,
                        salinity=np.zeros_like(Nor_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Nor.insert(Nor.shape[1], "DIC", Nor_results["dic"], True)
Nor_DIC = Nor_results["dic"]
Nor.to_excel(rivdir + Nore[:-4] + 'DIC.xls')

Sha = pd.read_excel(rivdir + Shannon, sheet_name="Alk_pH_temp")
Sha_alk = np.array(Sha['Alkalinity'])/0.05
Sha_pH = np.array(Sha['pH'])
Sha_temp = np.array(Sha['Temperature'])
Sha_results = pyco2.sys(par1=Sha_alk, par1_type=1,
                        par2=Sha_pH, par2_type=3,
                        temperature=Sha_temp,
                        salinity=np.zeros_like(Sha_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Sha.insert(Sha.shape[1], "DIC", Sha_results["dic"], True)
Sha_DIC = Sha_results["dic"]
Sha.to_excel(rivdir + Shannon[:-4] + 'DIC.xls')

Sla = pd.read_excel(rivdir + Slaney, sheet_name="Alk_pH_temp")
Sla_alk = np.array(Sla['Alkalinity'])/0.05
Sla_pH = np.array(Sla['pH'])
Sla_temp = np.array(Sla['Temperature'])
Sla_results = pyco2.sys(par1=Sla_alk, par1_type=1,
                        par2=Sla_pH, par2_type=3,
                        temperature=Sla_temp,
                        salinity=np.zeros_like(Sla_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Sla.insert(Sla.shape[1], "DIC", Sla_results["dic"], True)
Sla_DIC = Sla_results["dic"]
Sla.to_excel(rivdir + Slaney[:-4] + 'DIC.xls')

Sui = pd.read_excel(rivdir + Suir, sheet_name="Alk_pH_temp")
Sui_alk = np.array(Sui['Alkalinity'])/0.05
Sui_pH = np.array(Sui['pH'])
Sui_temp = np.array(Sui['Temperature'])
Sui_results = pyco2.sys(par1=Sui_alk, par1_type=1,
                        par2=Sui_pH, par2_type=3,
                        temperature=Sui_temp,
                        salinity=np.zeros_like(Sui_temp),
                        opt_pH_scale=1, opt_k_carbonic=8,
                        opt_k_bisulfate=1, opt_total_borate=1)
Sui.insert(Sui.shape[1], "DIC", Sui_results["dic"], True)
Sui_DIC = Sui_results["dic"]
Sui.to_excel(rivdir + Suir[:-4] + 'DIC.xls')

print('End')
