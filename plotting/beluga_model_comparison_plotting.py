# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 15:11:52 2025

@author: u300737
"""
import os
import pandas as pd
import numpy as np 
import xarray as xr

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.dates as mdates
import matplotlib.ticker as mticker

from matplotlib.gridspec import GridSpec

import metpy.calc as mpcalc
import metpy.constants as mpconsts
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

import seaborn as sns

from metpy.units import units

def verticalprofile_ERA5_for_specific_flight(STN_era5,df,rf="RF12"):
    flight_date            = str(df.index[0].date())
    era5_theta_df          = STN_era5["theta"]#.to_dataframe()
    era5_wspeed_df         = STN_era5["wspeed"]#.to_dataframe()
    era5_rh_df             = STN_era5["r"]#.to_dataframe()
    
    quick_scatter=plt.figure(figsize=(16,12))
    matplotlib.rcParams.update({"font.size":24})
    ax1=quick_scatter.add_subplot(131)
    ax2=quick_scatter.add_subplot(132)
    ax3=quick_scatter.add_subplot(133)
    ax1.scatter(df["Theta"],df["PRES"],color="salmon",s=1)
    pressure=np.array(era5_theta_df.columns)
    pressure=pressure[np.newaxis,:]
    pressure_array=np.repeat(pressure,era5_theta_df.shape[0],axis=0)
    
    ax1.scatter(era5_theta_df.values,
                pressure_array,
                color="w",s=15,edgecolor="darkred",lw=2)
    ax1.set_ylim([920,1030])
    ax1.set_xlim([int(df["Theta"].dropna().min())-2,
                  int(df["Theta"].loc[lambda v: v<np.Inf].max())+2])
    ax1.set_ylabel("Pressure (hPa)")
    ax1.set_xlabel("$\Theta$ (K)")
    #-------------------------------------------------------------------------#
    # Wind speed
    ax2.scatter(df["hor_vv"],df["PRES"],color="lightgreen",s=1)
    ax2.scatter(era5_wspeed_df.values,pressure_array,
                color="w",s=15,edgecolor="darkgreen",lw=2)
    ax2.set_ylim([920,1030])
    ax2.set_xlim([0,6])
    if df["hor_vv"].max()>6:
        ax2.set_xlim([int(df["hor_vv"].dropna().min())-2,
                  int(df["hor_vv"].loc[lambda v: v<np.Inf].max())+2])
    ax2.set_xlabel("Wind speed (m/s)")
    #-------------------------------------------------------------------------#
    # Relative humidity
    ax3.scatter(df["RH"],df["PRES"],color="lightblue",s=1)
    ax3.scatter(era5_rh_df.values,pressure_array,
                color="w",s=15,edgecolor="darkblue",lw=2)
    ax3.set_ylim([920,1030])
    ax3.set_xlim([0,100])
    ax3.set_xlabel("Relative Humidity (%)")
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
        ax2.spines[axis].set_linewidth(2)
        ax3.spines[axis].set_linewidth(2)
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)
        ax2.yaxis.set_tick_params(width=2,length=6)
        ax2.xaxis.set_tick_params(width=2,length=6)
        ax3.yaxis.set_tick_params(width=2,length=6)
        ax3.xaxis.set_tick_params(width=2,length=6)
    
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
    sns.despine(offset=10)
    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_name="Vertical_profiles_BELUGA_ERA5_"+rf+"_"+flight_date+".png"
    quick_scatter.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)

def verticalprofile_CARRA_for_specific_flight(STN_carra,data_df,rf="RF12",
                                              plot_flight=False,
                                              moisture="rh", icon_dict={}):
    
    matplotlib.rcParams.update({"font.size":24})
    
    from matplotlib.lines import Line2D
    # Correct offset in different heights. BELUGA measures height above msl, 
    # whereas CARRA height is provided in height above ground level
    df=data_df.copy()
    ground_alt=30
    df["ALT"]-=ground_alt
    df.name=rf
    ylim=700
    if df["ALT"].max()<550:
        ylim=550
    yticks=np.arange(0,ylim+50,100)
    ##### Reassign data
    flight_date             = str(df.index[0].date())
    if not STN_carra=={}:
        carra_theta_df          = STN_carra["theta"]#.to_dataframe()
        carra_wspeed_df         = STN_carra["wspeed"]#.to_dataframe()
        carra_moisture_df       = STN_carra[moisture]#.to_dataframe()
        alt=np.array(carra_theta_df.columns)
        alt=alt[np.newaxis,:]
        alt_array=np.repeat(alt,carra_theta_df.shape[0],axis=0)
        
    if moisture=="rh":
        moisture_label = "RH"
        moisture_unit  = " (%)"
    elif moisture=="q":
        moisture_label = "q"
        moisture_unit  = " ($\mathrm{g\,kg}^{-1})$"
    
    quick_scatter=plt.figure(figsize=(16,12))
    if not plot_flight:
        gs = GridSpec(1, 3, figure=quick_scatter)    
        ax1=quick_scatter.add_subplot(gs[0])
        ax2=quick_scatter.add_subplot(gs[1])
        ax3=quick_scatter.add_subplot(gs[2])
        sub_fig_labels=["(a)","(b)","(c)"]
    else:
        import beluga_plotting
        gs = GridSpec(2, 3, figure=quick_scatter,
                      height_ratios=[0.35,1],wspace=.2)    
        ax0=quick_scatter.add_subplot(gs[0,:])
        beluga_plotting.plot_RF_height(df,ax_obj=ax0)      
        ax1=quick_scatter.add_subplot(gs[1,0])
        ax2=quick_scatter.add_subplot(gs[1,1])
        ax3=quick_scatter.add_subplot(gs[1,2])
        sub_fig_labels=["(b)","(c)","(d)"]
    # Create a proxy artist for the legend with a larger marker size
    legend_marker_size = 20  # Larger size for the legend
    beluga_lgd_handle = Line2D([0], [0], marker='o', color='salmon',
            markersize=legend_marker_size, label='BELUGA')
    carra_lgd_handle = Line2D([0], [0], marker='o', color='white',
            markersize=legend_marker_size, markeredgecolor="darkred",
            markeredgewidth=3,label='CARRA')
    icon_lgd_handle = Line2D([0], [0], marker='v', color='grey',
            markersize=legend_marker_size+5, markeredgecolor="k",
            markeredgewidth=3,label='ICON')
    
    # Create legend objects
    
    # Scatter plots
    ax1.scatter(df["Theta"],df["ALT"],color="salmon",s=1) # BELUGA
    if not STN_carra=={}:
        ax1.scatter(carra_theta_df.values,alt_array,color="w",s=50,
                edgecolor="darkred",lw=2) # CARRA
    
    
    ax1.set_yticks(yticks)
    ax1.set_xlim([int(df["Theta"].dropna().min())-2,
                  int(df["Theta"].loc[lambda v: v<np.Inf].max())+2])
    ax1.set_ylabel("Height AGL (m)")
    ax1.set_xlabel("$\Theta$ (K)")
    ax1.legend(handles=[beluga_lgd_handle],
           loc="center left",fontsize=16)

    if not STN_carra=={}:
        ax1.legend(handles=[beluga_lgd_handle,carra_lgd_handle],
               loc="center left",fontsize=16)
    if not len(icon_dict)==0:
        ax1.legend(handles=[beluga_lgd_handle,carra_lgd_handle,icon_lgd_handle],
                   loc="center left",fontsize=16)
            
    #-------------------------------------------------------------------------#
    # Wind speed
    ax2.scatter(df["hor_vv"],df["ALT"],color="lightgreen",s=1)
    if not STN_carra=={}:
        ax2.scatter(carra_wspeed_df.values,alt_array,
                color="w",s=50,edgecolor="darkgreen",lw=2)
    
    ax2.set_yticks(yticks)
    ax2.set_xlim([int(df["hor_vv"].dropna().min())-2,
                  int(df["hor_vv"].loc[lambda v: v<np.Inf].max())+2])
    if df["hor_vv"].max()<7.5:
        ax2.set_xlim([0,7.5])
    else:
        ax2.set_xlim([0,10])
    ax2.set_xlabel("Wind speed ($\mathrm{m\,s}^{-1}$)")
    #-------------------------------------------------------------------------#
    # Relative humidity
    ax3.scatter(df[moisture_label],df["ALT"],color="lightblue",s=1)
    if not STN_carra=={}:
        ax3.scatter(carra_moisture_df.values,alt_array,
                color="w",s=50,edgecolor="darkblue",lw=2)
    
    if not len(icon_dict)==0:
        icon_theta_df          = icon_dict["theta"]#.to_dataframe()
        icon_wspeed_df         = icon_dict["wspeed"]#.to_dataframe()
        icon_moisture_df       = icon_dict[moisture]#.to_dataframe()
        if moisture=="rh":
            icon_moisture_df*=100
        icon_alt=np.array(icon_theta_df.columns)
        icon_alt=icon_alt[np.newaxis,:]
        
        icon_alt_array=np.repeat(icon_alt,icon_theta_df.shape[0],axis=0)
        
        ax1.scatter(icon_theta_df.values,icon_alt_array,#-ground_alt,
                    marker="v",color="dimgrey",s=55,
                    edgecolor="k",lw=2) # ICON
        ax2.scatter(icon_wspeed_df.values,icon_alt_array,#-ground_alt,
                    marker="v",color="dimgrey",s=55,
                    edgecolor="k",lw=2)
        ax3.scatter(icon_moisture_df.values,icon_alt_array,#-ground_alt,
                    marker="v",color="dimgrey",s=55,
                    edgecolor="k",lw=2)
        
    ax1.set_ylim([0,ylim])
    ax2.set_ylim([0,ylim])
    ax3.set_ylim([0,ylim])
    ax3.set_yticks(yticks)
    if moisture=="rh":
        ax3.set_xlim([0,100])
    elif moisture=="q":
        ax3.set_xlim([0,2])
    ax3.set_xlabel(moisture_label+moisture_unit)
        
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
        ax2.spines[axis].set_linewidth(2)
        ax3.spines[axis].set_linewidth(2)
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)
        ax2.yaxis.set_tick_params(width=2,length=6)
        ax2.xaxis.set_tick_params(width=2,length=6)
        ax3.yaxis.set_tick_params(width=2,length=6)
        ax3.xaxis.set_tick_params(width=2,length=6)
    
    ax1.text(0.03,0.95,sub_fig_labels[0],transform=ax1.transAxes,fontsize=20)
    ax2.text(0.03,0.95,sub_fig_labels[1],transform=ax2.transAxes,fontsize=20)
    ax3.text(0.03,0.95,sub_fig_labels[2],transform=ax3.transAxes,fontsize=20)
    
    plt.subplots_adjust(wspace=0.5,hspace=0.3)
    sns.despine(offset=10)
    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    major_name="Vertical_profiles_BELUGA_"
    if not STN_carra=={}:
        major_name+="CARRA_"
    if not len(icon_dict)==0:
        major_name+="ICON_"
    fig_name=major_name+moisture_label+"_"+\
        rf+"_"+flight_date+".png"
    quick_scatter.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)
    
def contour_profile_CARRA_ICON_BELUGA(beluga_dataframe,rf,carra_dict,icon_dict,
                                      var_to_plot="Theta"):
    ms=12 # marker size
    x_carra,y_carra = np.meshgrid(carra_dict[var_to_plot.lower()].index,
                                  carra_dict[var_to_plot.lower()].columns)
    if icon_dict:
        x_icon,y_icon = np.meshgrid(icon_dict[var_to_plot.lower()].index,
                                  icon_dict[var_to_plot.lower()].columns)  
    
    colormap={"theta":"Spectral_r",
              "theta_e":"Spectral_r",
              "q":"Blues",
              "rh":"Greys",
              "wspeed":"viridis"} 
    
    units   = {"theta":" (K)","theta_e":" (K)",
              "q":" ($\mathrm{g\,kg}^{-1}$)",
              "rh":"%","wspeed":" ($\mathrm{m\,s}^{-1}$)"}
    
    c_ticks= {"theta":[255,260,265,270],
              "theta_e":[260,270,280,290,300],
              "q":[0,0.5,1,1.5,2],
              "rh":[20,60,100],
              "wspeed":[0,2,4,6,8]}
    
    beluga_df=beluga_dataframe.copy()
    beluga_df["ALT"]-=31
    if rf=="RF1516":
        # interpolate theta
        #beluga_df=
        beluga_df["Theta"]=beluga_df["Theta"].interpolate("linear")
        beluga_df["ALT"]=beluga_df["ALT"].interpolate("linear")
        beluga_df["ALT"].loc["2024-04-02 10:11:18"]=np.nan
    if var_to_plot=="wspeed":
        beluga_df["Wspeed"]=beluga_df["hor_vv"].copy()
    ### Plotting
    matplotlib.rcParams.update({"font.size":24})
    profile_comparison_fig=plt.figure(figsize=(16,12))
    # check if ICON dictionary is empty or not
    if not icon_dict:
        ax1=profile_comparison_fig.add_subplot(111)
    else:
        ax1=profile_comparison_fig.add_subplot(211)
    C1=ax1.contourf(x_carra,y_carra,carra_dict[var_to_plot.lower()].T.values,
                    levels=np.linspace(c_ticks[var_to_plot.lower()][0],
                                       c_ticks[var_to_plot.lower()][-1],31),
                    cmap=colormap[var_to_plot.lower()],
                    extend="both")
    
    cax1 = profile_comparison_fig.add_axes([.95, 0.25, 0.01, 0.5])
    cbar=plt.colorbar(C1, cax=cax1)
    cbar.set_label(var_to_plot+units[var_to_plot.lower()], labelpad=1)
    cbar.set_ticks(c_ticks[var_to_plot.lower()])
    
    ax1.set_ylabel("Height AGL (m)")
    ax1.plot(beluga_df.index,beluga_df["ALT"],color="whitesmoke",
             lw=7,zorder=2)
    
    try:
        ax1.scatter(beluga_df.index,beluga_df["ALT"],
                c=beluga_df[var_to_plot.capitalize()],s=ms,
                vmin=c_ticks[var_to_plot.lower()][0],
                vmax=c_ticks[var_to_plot.lower()][-1],
                cmap=colormap[var_to_plot.lower()],zorder=3)
    except:
        try:
            ax1.scatter(beluga_df.index,beluga_df["ALT"],
                c=beluga_df[var_to_plot.upper()],s=ms,
                vmin=c_ticks[var_to_plot.lower()][0],
                vmax=c_ticks[var_to_plot.lower()][-1],
                cmap=colormap[var_to_plot.lower()],zorder=3)
        except:
            ax1.scatter(beluga_df.index,beluga_df["ALT"],
                c=beluga_df[var_to_plot.lower()],s=ms,
                vmin=c_ticks[var_to_plot.lower()][0],
                vmax=c_ticks[var_to_plot.lower()][-1],
                cmap=colormap[var_to_plot.lower()],zorder=3)
    if not icon_dict:
        ylim=550
    else:
        ylim=800
    ax1.set_ylim([0,ylim])
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax1.text(.05,.95,"(a) CARRA",fontsize=22,transform=ax1.transAxes)
    if icon_dict:
        ax2=profile_comparison_fig.add_subplot(212)
        C2=ax2.contourf(x_icon,y_icon,icon_dict[var_to_plot.lower()].T.values,
                    levels=np.linspace(c_ticks[var_to_plot.lower()][0],
                                       c_ticks[var_to_plot.lower()][-1],31),
                    cmap=colormap[var_to_plot.lower()],extend="both")
    
        ax2.plot(beluga_df.index,beluga_df["ALT"],color="whitesmoke",
             lw=7,zorder=2)
        try:
            ax2.scatter(beluga_df.index,beluga_df["ALT"],
                c=beluga_df[var_to_plot.capitalize()],s=ms,
                vmin=c_ticks[var_to_plot.lower()][0],
                vmax=c_ticks[var_to_plot.lower()][-1],
                cmap=colormap[var_to_plot],zorder=3)
        except:
            try:
                ax2.scatter(beluga_df.index,beluga_df["ALT"],
                            c=beluga_df[var_to_plot.upper()],s=ms,
                            vmin=c_ticks[var_to_plot.lower()][0],
                            vmax=c_ticks[var_to_plot.lower()][-1],
                            cmap=colormap[var_to_plot.lower()],zorder=3)
            except:
                ax2.scatter(beluga_df.index,beluga_df["ALT"],
                            c=beluga_df[var_to_plot.lower()],s=ms,
                            vmin=c_ticks[var_to_plot.lower()][0],
                            vmax=c_ticks[var_to_plot.lower()][-1],
                            cmap=colormap[var_to_plot.lower()],zorder=3)
    
        ax2.set_ylim([0,ylim])
        ax2.set_ylabel("Height AGL (m)")
    
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        for axis in ['bottom','left']:
            ax2.spines[axis].set_linewidth(2)
        ax2.yaxis.set_tick_params(width=2,length=6)
        ax2.xaxis.set_tick_params(width=2,length=6)
    
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.text(.05,.9,"(b) ICON-2km",fontsize=22,transform=ax2.transAxes)
        ax2.set_xlabel("2024-04-01 (HH:MM)")
    else:
        ax1.set_xlabel(str(beluga_df.index[0].date())+" (HH:MM)")
    fig_name=rf+"_transition_event"+var_to_plot+".png"
    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_path=plot_path+fig_name
    profile_comparison_fig.savefig(fig_path,dpi=300,
                                   bbox_inches="tight")
    print("Figure saved as:", fig_path)
    
    