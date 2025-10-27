# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 16:45:00 2025

@author: u300737
"""

import os 
import glob
import sys
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import matplotlib
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

import metpy.calc as mpcalc
import metpy.constants as mpconsts
from metpy.units import units

warnings.filterwarnings("ignore")
sys.path.insert(1,os.getcwd()+"/../src/analysing/")


dates_april = ["2024-04-23","2024-04-24","2024-04-25","2024-04-26"]
dates_may   = ["2024-04-29","2024-04-30","2024-05-01","2024-05-02"] 


#import Measurement_Platforms
import Grid_data_STN
import carra_plotting as CARRA_plotting

#-----------------------------------------------------------------------------#
# Definitions
matplotlib.rcParams.update({"font.size":14})
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
STN_coords={"lat": 81.60,
            "lon": -16.65}
#-----------------------------------------------------------------------------#
# Flight definitions
flight_list={"RF12":"20240401",
             "RF13":"20240401",
             "RF14":"20240401",
             "RF1312":"20240401",
             "RF15":"20240402",
             "RF1516":"20240402",
             "RF141312":"20240401",
             "RF26":"20240410"}
flight="RF1516"
main_path=os.getcwd()+"/../"
#-----------------------------------------------------------------------------#
# Open and allocate classes
#Platforms_STN_cls=Measurement_Platforms.Measurement_Platforms_STN(
#    main_path=main_path)
#BELUGA_cls=Measurement_Platforms.BELUGA(Platforms_STN_cls)
#BELUGA_cls.get_flight_information()

#cpgn_df=BELUGA_cls.open_all_met_data(level="L2")
#if not "hor_vv" in cpgn_df.columns:
#    cpgn_df=BELUGA_cls.calc_wspeed(cpgn_df)
#    cpgn_df["TEMP"]-=273.15
#    cpgn_df[cpgn_df["TEMP"]>0]=np.nan
#%% Gridded data Reanalyses
Grid_data_cls     = Grid_data_STN.Gridded_data(main_path=main_path)
CARRA             = Grid_data_STN.CARRA(Grid_data_cls)

#%% CARRA plotting
# April data
CARRA_plotting.synoptic_carra_map(CARRA,None,main_path,dates=dates_april)
# May data
CARRA_plotting.synoptic_carra_map(CARRA,None,main_path,dates=dates_may)