# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 10:44:19 2025

@author: u300737
"""
import glob
import os
import xarray as xr
import sys

#sys.path.insert(1, os.getcwd()+)
path="C:\\Users\\u300737\\Desktop\\Desktop_alter_Rechner\\BELUGA_Leipzig_AC3"+\
    "\\Code\\GIT\\STN_analysis\\BELUGA_data\\BELUGA_broadband_probe\\"+\
        "temporary_data\\"
version_number="v2.2"

# Open all processed BP irradiance data
file_structure="*"+version_number+".nc"
files_to_open=glob.glob(path+file_structure)
BP_ds=xr.open_mfdataset(files_to_open,combine="nested",concat_dim="time")

# Get the irradiances at the respective heights
Fnet_df=BP_ds[["z_b","F_net","segments"]].to_dataframe()
hmax_df=Fnet_df.loc[Fnet_df["segments"]=="peak"]
hmax_df=hmax_df.loc[hmax_df["z_b"]>300]
h_30_df=Fnet_df.loc[Fnet_df["z_b"].between(15,45)]
cpgn_TNIR={}
cpgn_TNIR["peak"] = hmax_df
cpgn_TNIR["h30"]  = h_30_df
# Plot the data using the beluga_plotting module
import beluga_plotting
beluga_plotting.plot_tnir_histogram(cpgn_TNIR)



