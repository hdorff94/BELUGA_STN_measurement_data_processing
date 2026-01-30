# -*- coding: utf-8 -*-
"""
Created on Tue Aug 19 15:20:22 2025

@author: u300737
"""
import os
import glob

import metpy.calc as mpcalc
from metpy.units import units


import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

def plot_sonde_period(sonde_df, height_limit=10):
    launch_dates=np.unique(sonde_df["Theta"].index.date)
    fontsize=28
    matplotlib.rcParams.update({'font.size':fontsize})
    sonde_fig,ax=plt.subplots(3,1,sharex=True,figsize=(16,18))
    ax1=ax[0]
    ax2=ax[1]
    ax3=ax[2]
    mlw=1
    ws_max=25
    theta_max=300
    if height_limit<1:
        mlw=5
        ws_max=15
        theta_max=280
    for d, date in enumerate(launch_dates):
        sonde_day=sonde_df.loc[str(date)]
        altitude=sonde_day["altitude"]
        wspeed=np.sqrt(sonde_day["u"]**2+sonde_day["v"]**2)
        C1=ax1.scatter(np.repeat(date,altitude.shape[0]),
            altitude.values/1000,
            c=sonde_day["Theta"].values,vmin=250,vmax=theta_max,
            cmap="Spectral_r",s=500, marker="_",linewidths=mlw) # default s=100
        C2=ax2.scatter(
            np.repeat(date,altitude.shape[0]),
            altitude.values/1000,
            c=sonde_day["RH"].values,vmin=0,vmax=100,
            cmap="Blues",s=500, marker="_",linewidths=mlw)
        C3=ax3.scatter(np.repeat(date,altitude.shape[0]),
            altitude.values/1000,
            c=wspeed.values,vmin=0,vmax=ws_max,
            cmap="viridis",s=500, marker="_",linewidths=mlw)
    
    
    cax1 = sonde_fig.add_axes([0.925, 0.68, 0.01, 0.175])
    cax2 = sonde_fig.add_axes([0.925, 0.4125, 0.01, 0.175])
    cax3 = sonde_fig.add_axes([0.925, 0.15, 0.01, 0.175])
    
    cbar1=plt.colorbar(C1, cax=cax1)
    cbar1.set_label('$\Theta$ (K)', labelpad=1)
    cbar1.set_ticks([250,275,theta_max])
    if height_limit<1:
        cbar1.set_ticks([250,260,270,280])
        #cbar1.ylim([250,theta_max])
    cbar2=plt.colorbar(C2, cax=cax2)
    cbar2.set_label('RH (%)', labelpad=1)
    cbar2.set_ticks([0,50,100])
    
    cbar3=plt.colorbar(C3, cax=cax3)
    cbar3.set_label('WS ($\mathrm{m\,s}^{-1}$)', labelpad=1)
    cbar3.set_ticks([0,12.5,ws_max])
    if height_limit<1:
        cbar3.set_ticks([0,5,10,ws_max])
        #cbar3.set_ylim([0,ws_max])
    
    ax1.set_xlim([pd.to_datetime("2024-03-20"),pd.to_datetime("2024-04-12")])
    xticks=pd.to_datetime(["2024-03-21","2024-03-26","2024-04-01","2024-04-06",
                    "2024-04-11"])
    ax1.set_xticks(xticks)
    
    ax1.set_ylim([0,height_limit])
    ax2.set_ylim([0,height_limit])
    ax3.set_ylim([0,height_limit])
    ax1.set_ylabel("Height (km)")
    ax2.set_ylabel("Height (km)")
    ax3.set_ylabel("Height (km)")
    
    fig_labels=["(a)","(b)","(c)"]
    for a,ax in enumerate([ax1,ax2,ax3]):
        [ax.spines[axis].set_linewidth(2) for axis in ['bottom','left']]
        ax.spines[['right', 'top']].set_visible(False)
        ax.yaxis.set_tick_params(width=2,length=8)
        ax.xaxis.set_tick_params(width=3,length=10)
        ax.text(-0.11,0.96,fig_labels[a],color="k",fontsize=fontsize-2,
                bbox=dict(facecolor='whitesmoke',edgecolor="black",
                          boxstyle='round'),transform=ax.transAxes)

    fig_name = "fig09_radiosondes_campaign_"+str(height_limit)+"km"+".png"
    fig_path = os.getcwd()+"/../plots/campaign_overview/"
    os.makedirs(fig_path,exist_ok=True)
    sonde_fig.savefig(fig_path+fig_name, dpi=600,bbox_inches="tight")
    print("Figure saved as:", fig_path+fig_name)
    
        
    #def plot_sonde_campaign_data():
radiosonde_l2_path=os.getcwd()+"/../BELUGA_data//"+\
        "radiosondes//final_data//"

level="2"

sonde_file_style="*"+level+"*.nc"
        
file_list=glob.glob(radiosonde_l2_path+sonde_file_style)

sonde_df=xr.open_mfdataset(
    file_list,combine="nested",concat_dim="time").to_dataframe()

sonde_df["Theta"]=np.nan
# Prepare for metpy
T_profile=sonde_df["T"]+273.15
T_profile=T_profile * units.kelvin
p_hPa=sonde_df["p"]
p_hPa=p_hPa * units.hPa
Theta_profile=mpcalc.potential_temperature(
    p_hPa.values * units.hPa,T_profile.values * units.kelvin)

sonde_df["Theta"]=np.array(Theta_profile)
plot_sonde_period(sonde_df,height_limit=.5)