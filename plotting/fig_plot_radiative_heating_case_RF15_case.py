# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 10:22:28 2025

@author: u300737
"""
import glob
import os
import sys

import numpy  as np
import pandas as pd 
import xarray as xr

import seaborn as sns
# Information from Lonardi et al. (2024); :
# The vertical divergence of Fnet in an atmospheric layer results in local 
# temperature tendencies (cooling or warming), induced by radiation in a layer 
# between two heights zbot and ztop, quantified by radiative temperature 
# tendency rate (ζ). ρ is density of air, cp is specific heat capacity of air 
# at constant pressure. Positive values; “heating rates”, negative values;
# “cooling rates”. Magnitude of ζ depends on layer thickness over which 
# convergence or divergence of radiative fluxes occurs. 
# first: interpolation of the irradiance profiles to common vertical grid with 
# 1m resolution. Second: ζ profiles calculated using layer thickness of 10 m. 
# Such resolution allows for characterizing the temperature tendencies also 
# in the narrow cloud top region.

def vertical_uniform_regridding(profile_df,zmin=0,zmax=450):
    # Common vertical grid with 1 meter resolution from zmin to zmax
    z_1m = np.arange(zmin, zmax + 1, 1)

    # Interpolate Fnet profiles onto the common vertical grid 
    # for each time step
    def interpolate_profile(fnet_profile, z_profile, new_z):
        return np.interp(new_z, z_profile, fnet_profile)

    # Prepare an empty array to hold interpolated profiles
    # Dimensions: time x z
    #times = ds['time']
    z_profile = profile_df["z_b"]
        
    interp_fnet = interpolate_profile(profile_df["F_net"], z_profile, z_1m)
    fnet_1m     = pd.Series(data=interp_fnet,index=z_1m)
    return fnet_1m
def calc_radiative_heating_rate_zeta(f_net,layer_thickness=10):
    
    """
    Parameters
    ----------
    f_net : pandas.Series
        1 m regridded pandas Series containing the net irradiances.
    layer_thickness : float, optional
        layer thickness in m for calculating the divergence. Default is 10 m.

    Returns
    -------
    zeta : pandas.Series
        calculated heating rates using the formula of Lonardi et al. (2024)
    """
    # Constants
    cp  = 1005.0  # J/(kg·K)
    rho = 1.2922  # kg/m³
    
    common_z=f_net.index
    # 3. For each 10 m layer, compute Fnet difference
    layer_thickness = 10  # meters
    num_layers = len(common_z) - 10  # Number of 10m layers

    zeta_values     = []
    layer_midpoints  = []
    
    for i in range(num_layers):
        Fnet_start = f_net.iloc[i]
        Fnet_end = f_net.iloc[i + layer_thickness]
        dF_dz_layer = (Fnet_end - Fnet_start) / layer_thickness  
        # derivative over 10 m
        
        # Calculate ζ
        zeta_layer = (1 / (rho * cp)) * dF_dz_layer
        # Store layer midpoint height
        layer_midpoints.append(common_z[i] + layer_thickness/ 2)
        zeta_values.append((zeta_layer))
    
    # Convert to pandas  Series
    zeta_series = pd.Series(data=zeta_values, index=layer_midpoints)
    # provide zeta as rates per hour
    zeta_series*=3600
    # Resample the 1 m resolution zeta to 10 m means as actually the difference 
    # is also calculated using a layer thickness of 10 m    
    bin_width = 10
    bins = np.arange(zeta_series.index.min(), zeta_series.index.max() +\
                     bin_width, bin_width)
    # Assign each index to a bin
    bin_labels = (bins[:-1] + bin_width/2)  # Midpoints of bins for new index
    groups = pd.cut(zeta_series.index, bins=bins, labels=bin_labels)

    # Group by these bins and compute the mean
    resampled_zeta = zeta_series.groupby(groups).mean()
    #resampled_zeta  = zeta_series.rolling(10).mean()
    return resampled_zeta
flight="RF15"

data_path=path="C:\\Users\\u300737\\Desktop\\Desktop_alter_Rechner\\BELUGA_Leipzig_AC3"+\
    "\\Code\\GIT\\STN_analysis\\BELUGA_data\\BELUGA_broadband_probe\\"+\
        "temporary_data\\"
version_number="v2.5"

# Open all processed BP irradiance data
file_structure="*"+flight+"*"+version_number+".nc"
files_to_open=glob.glob(path+file_structure)
BP_ds=xr.open_dataset(files_to_open[0])
#
#%% Access first and last profile 
# (first ascent and last descent)
# Get all segments
BP_df            = BP_ds.to_dataframe()
segment_labels   = BP_df["segments"].unique()
maxima           = BP_df[BP_df["segments"]=="max"]

first_profile_start = BP_df[BP_df["segments"]=="profile 01 ascent"].index[0]
first_profile_end   = BP_df[BP_df["segments"]=="max"].index[0]

last_profile_start  = BP_df[BP_df["segments"]=="max"].index[-1]
last_profile_end    = BP_df[BP_df["segments"]=="profile 20 descent"].index[-1]

first_profile       = BP_df.loc[first_profile_start:first_profile_end].rolling("3s").mean()
last_profile        = BP_df.loc[last_profile_start:last_profile_end].rolling("3s").mean()

first_profile_1m = vertical_uniform_regridding(first_profile,zmin=0,zmax=450)
last_profile_1m  = vertical_uniform_regridding(last_profile[::-1],zmin=0,
                                               zmax=450)
first_profile_1m=first_profile_1m.rolling(10).mean()
last_profile_1m=last_profile_1m.rolling(10).mean()
# Heating rates Zeta
zeta_first = calc_radiative_heating_rate_zeta(
                first_profile_1m,layer_thickness=10)
zeta_last  = calc_radiative_heating_rate_zeta(
                last_profile_1m,layer_thickness=10)
# Quicklook
import matplotlib
import matplotlib.pyplot as plt
rad_fig=plt.figure(figsize=(16,9))
matplotlib.rcParams.update({"font.size":22})
ax1=rad_fig.add_subplot(121)
ax2=rad_fig.add_subplot(122)

ax1.plot(first_profile["F_net"],first_profile["z_b"],lw=3,color="slateblue",
         label="first profile")
ax1.plot(last_profile["F_net"],last_profile["z_b"],lw=3,color="orange",
         label="last profile")
ax1.text(0.01,0.9,"(a)",color="k",bbox=dict(
    facecolor='whitesmoke',edgecolor="black", boxstyle='round'),
    transform=ax1.transAxes)
ax1.set_xlabel("Net irradiance $F_{\mathrm{net}}$ ($\mathrm{W\,m}^{-2}$)")
ax1.set_ylabel("Barometric height (m)")
ax1.set_ylim([0,400])
ax1.set_yticks([0,100,200,300,400])
ax2.step(zeta_last,zeta_last.index,lw=3,color="orange")
ax2.step(zeta_first,zeta_first.index,lw=3,color="slateblue",zorder=2)
#Mean heating rate (withdrawn)
#ax2.step((zeta_first+zeta_last)/2,list(zeta_first.index.values),color="w", lw=6,zorder=2,where='post')#
#ax2.step((zeta_first+zeta_last)/2,list(zeta_first.index.values),color="k",
#         label="mean",lw=4,zorder=3)
ax2.axvline(x=0,ls="--",lw=3,color="grey")
ax1.axhline(y=192,lw=8,color="lightgray",zorder=0)
ax2.axhline(y=192,lw=8,color="lightgray",zorder=0)
ax1.text(-45,175,"Inversion bottom",color="dimgrey",fontsize=18)
ax2.text(0.01,0.9,"(b)",color="k",bbox=dict(
    facecolor='whitesmoke',edgecolor="black", boxstyle='round'),
    transform=ax2.transAxes)
ax2.set_xlabel("Heating rates $\zeta$ ($\mathrm{K\,h}^{-1})$")
ax2.set_ylim([0,400])
ax2.set_yticks([0,100,200,300,400])
ax1.legend(loc="upper right", fontsize=18)
for axis in ['bottom','left']:
    ax1.spines[axis].set_linewidth(2)
    ax2.spines[axis].set_linewidth(2)
ax1.yaxis.set_tick_params(width=2,length=6)
ax1.xaxis.set_tick_params(width=2,length=6)
ax2.yaxis.set_tick_params(width=2,length=6)
ax2.xaxis.set_tick_params(width=2,length=6)
ax2.set_xlim([-.6,.4])
ax1.set_xlim([-45,-20])
ax1.set_xticks([-45,-40,-35,-30,-25,-20])
ax2.set_xticks([-.6,-.3,0,.3])
ax2.legend(loc="upper right",fontsize=18)
sns.despine(offset=10)
fig_name=flight+"_Heating_rates_start_end.pdf"
plot_path=os.getcwd()+"/../plots/"
rad_fig.savefig(plot_path+fig_name,bbox_inches="tight")
print("Figure saved as:", plot_path+fig_name)
sys.exit() 
