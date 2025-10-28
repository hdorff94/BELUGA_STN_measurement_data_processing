# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 22:50:16 2025

@author: u300737
"""

import glob
import os
import sys

import numpy  as np
import pandas as pd 
import xarray as xr

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.signal import find_peaks

def plot_profile_heights(max_points):
    matplotlib.rcParams.update({"font.size":24})
    profile_info_fig=plt.figure(figsize=(16,9))
    ax1=profile_info_fig.add_subplot(111)
    
    for idx,flight in enumerate(flights):
        print(flight)
        if flight=="RF09":
            continue
        else:
            alt_max_df=max_points.loc[max_points["flight"]==flight]
            alt_max_vals=alt_max_df["z_b"].values
            plot_idx=np.ones(alt_max_df.shape[0])*(idx+1)
            transition_type=alt_max_df["transition_type"].iloc[0]
            # check for specific transition types
            if transition_type=="Polar night to day":
                marker_color="orange"
            elif transition_type=="LLJ":
                marker_color="green"
            elif transition_type=="clear to cloudy":
                marker_color="slateblue"
            else:
                marker_color="grey"
                transition_type="none"
                
            edge_color="k"
            ax1.scatter(plot_idx,alt_max_vals,s=200,
                marker="x",color=marker_color,edgecolors="k",
                linewidths=4,label=transition_type)
            
    ax1.set_xticks(rf_indices+1,flight_labels)
    ax1.tick_params(axis='x', labelrotation=90)
    ax1.set_ylim([200,950])
    ax1.set_ylabel("Barometric height (m)")
    ax1.set_xlabel("Research flight (RF)")
    sns.despine(ax=ax1,offset=10)
    
    transition_lgd_handle=[]
    transition_lgd_handle.append(mlines.Line2D([0], [0], marker='x',
        color='orange',markeredgewidth=4,label="night$\leftrightarrow$day", 
        ls="none",markersize=20))
    transition_lgd_handle.append(mlines.Line2D([0], [0], marker='x',ls="none",
        color='slateblue',markeredgewidth=4,label="clear$\leftrightarrow$cloudy",
        markersize=20))
    transition_lgd_handle.append(mlines.Line2D([0],[0],marker="x",ls="none",
        color="green",markeredgewidth=4,label="LLJ",markersize=20))
    transition_lgd_handle.append(mlines.Line2D([0],[0],marker="x",ls="none",
        color="grey",markeredgewidth=4,label="none",markersize=20))
    [ax1.spines[axis].set_linewidth(2) for axis in ['bottom','left']]
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    ax1.legend(handles=transition_lgd_handle,loc="upper left",
               ncols=4,title="Transition types",fontsize=18)
    # Save figure
    fig_name = "fig03_scatter_profile_alts.pdf"
    fig_path = os.getcwd()+"/plots/campaign_overview/"
    os.makedirs(fig_path,exist_ok=True)
    profile_info_fig.savefig(fig_path+fig_name, dpi=300,bbox_inches="tight")
    print("Figure saved as:", fig_path+fig_name)
    return None

do_plotting=False
flights=["RF"+str(i).zfill(2) for i in np.arange(1,29)]
flight_labels=[str(i).zfill(2) for i in np.arange(1,29)]
main_data_path="C://Users//u300737//Desktop//Desktop_alter_Rechner//BELUGA_Leipzig_AC3"+\
    "//Code//GIT//STN_analysis//BELUGA_data//"
BP_data_path=main_data_path+"BELUGA_broadband_probe//"+\
        "temporary_data\\"
TMP_data_path=main_data_path+"BELUGA_TMP_met_probe//"+\
    "temporary_data\\"
version_number="v2.2"
rf_indices=np.arange(len(flights))

# Open all processed BP irradiance data
# read TMP_met first (default), and if there is no data available 
# read BP_dataset
rf_dfs=[]
for flight in flights:
    print("flight: ",flight)
    TM_file_structure="*TMP_met_"+flight+"*"+version_number+".nc"
    files_to_open=glob.glob(TMP_data_path+TM_file_structure)
    if len(files_to_open)>0:
        ds=xr.open_dataset(files_to_open[0])
        balloon_df=ds[["z_b","segments"]].to_dataframe()
        balloon_df["transition_type"]=ds.transition_type.values
        balloon_df["flight"]=ds.flight
        rf_dfs.append(balloon_df)
    else:
        BP_file_structure="*BP_"+flight+"*"+version_number+".nc"
        files_to_open=glob.glob(BP_data_path+BP_file_structure)
        if len(files_to_open)>0:
            ds=xr.open_dataset(files_to_open[0])
            balloon_df=ds[["z_b","segments"]].to_dataframe()
            balloon_df["transition_type"]=ds.transition_type.values
            balloon_df["flight"]=ds.flight
            rf_dfs.append(balloon_df)
        else:
            continue

rfs_df=pd.concat(rf_dfs)

# Get the max values of the profiles
max_points=rfs_df.loc[rfs_df["segments"]=="max"]
if do_plotting:
    plot_profile_heights(max_points)
# Get only profile periods    
profiles_df=rfs_df[rfs_df["segments"].str.contains(
    "profile",case=False,na=False)]
profiles_df["abs_height_rates"]=abs(profiles_df["z_b"].diff())

print(profiles_df["abs_height_rates"].describe())
"""
        for idx,rf in enumerate(flight_info.index):
            rf_df=BELUGA_cls.open_met_data(rf=rf)
            alt_max_vals=find_profile_peaks(rf_df,min_height=min_height)
            plot_idx=np.ones(alt_max_vals.shape[0])*(idx+1)
            transition_type=flight_info["type"].loc[rf]
            # check for specific transition types
            if transition_type=="Polar night to day":
                marker_color="orange"
            elif transition_type=="LLJ":
                marker_color="green"
            elif transition_type=="clear to cloudy":
                marker_color="slateblue"
            else:
                marker_color="grey"
                transition_type="none"
                
            edge_color="k"
            ax1.scatter(plot_idx,alt_max_vals,s=200,marker="x",color=marker_color,
                       edgecolors="k",linewidths=4,label=transition_type)
            
        ax1.set_xticks(rf_indices+1,flight_info.index)
        ax1.tick_params(axis='x', labelrotation=90)
        ax1.set_ylim([200,950])
        ax1.set_ylabel("Height AGL (m)")
        ax1.set_xlabel("Research flight")
        #ax1.legend(loc="upper left",ncol=4)
        sns.despine(ax=ax1,offset=10)
        transition_lgd_handle=[]
        
        transition_lgd_handle.append(mlines.Line2D([0], [0], marker='x',
                color='orange',markeredgewidth=4,label="night$\leftrightarrow$day", ls="none",
                markersize=20))
        transition_lgd_handle.append(mlines.Line2D([0], [0], marker='x',ls="none",
                color='slateblue',markeredgewidth=4,label="Clear$\leftrightarrow$cloudy",
                markersize=20))
        transition_lgd_handle.append(mlines.Line2D([0],[0],marker="x",ls="none",
                color="green",markeredgewidth=4,label="LLJ",markersize=20))
        transition_lgd_handle.append(mlines.Line2D([0],[0],marker="x",ls="none",
                color="grey",markeredgewidth=4,label="none",markersize=20))
        
        [ax1.spines[axis].set_linewidth(2) for axis in ['bottom','left']]
            
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)
        ax1.legend(handles=transition_lgd_handle,loc="upper left",ncols=4,
                   title="Transition types",fontsize=18)
        # Save figure
        fig_name = "scatter_profile_alts.pdf"
        fig_path = os.getcwd()+"/plots/campaign_overview/"
        os.makedirs(fig_path,exist_ok=True)
        profile_info_fig.savefig(fig_path+fig_name, dpi=300,bbox_inches="tight")
        print("Figure saved as:", fig_path+fig_name)
    """