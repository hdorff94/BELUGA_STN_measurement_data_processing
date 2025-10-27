# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 15:39:28 2025

@author: u300737
"""

import glob

import xarray as xr

main_data_path="C://Users//u300737//Desktop//Desktop_alter_Rechner//"+\
    "BELUGA_Leipzig_AC3"+"//Code//GIT//STN_analysis//BELUGA_data//"
ceilometer_data_path=main_data_path+"ceilometer//"


rfs_dict={"RF12":"20240401",
          "RF15":"20240402",}

# RF to analyse
flight="RF15"
flight_date=rfs_dict[flight]

ceilo_fname="A4"+flight_date[4:8]+"*.nc"
ceilo_files=glob.glob(ceilometer_data_path+ceilo_fname)

ceilo_ds=xr.open_dataset(ceilo_files[0])#xr.open_mfdataset(ceilometer_data_path+ceilo_files)