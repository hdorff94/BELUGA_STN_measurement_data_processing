import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

import seaborn as sns

#------------------------------------------------------------------------------------#
# Turbulent meteorological probe (meteorological component)
#------------------------------------------------------------------------------------#
def TMP_met_plot_series_all_vars(TM_ds,rf):
    """
    This routines plots time series of major meteorological quantities 
    measured by TMP_met. 
    
    Input:
        Parameters
        ----------
        TM_ds : xr.Dataset
            dataset of TMP_met component.
        rf : str
            research flight.

        Returns
        -------
        None.
    """
    
    TM_df=TM_ds[["z_b","T","rh","vv","dd"]].to_dataframe()
    flight_time=TM_df.index
    matplotlib.rcParams.update({"font.size":16})
    # Set the style and colormap
    sns.set_style('whitegrid')
    cmap = sns.color_palette("icefire", as_cmap=True)
    #norm = plt.Normalize(flight_time.min(), flight_time.max())
    # Create figure and axes objects
    fig, axes = plt.subplots(5, 1, figsize=(8, 14), sharex=True)

    # Plot z_b (barometric height)
    axes[0].scatter(flight_time,TM_df['z_b'],s=3, c=flight_time,cmap=cmap)
    axes[0].set_ylabel('Barometric Height (m)')
    # Plot temperature
    sc = axes[1].scatter(flight_time,TM_df['T'],s=3, c=flight_time, cmap='icefire')
    axes[1].set_ylabel('Temperature (K)')
    # Plot relative humidity
    sc = axes[2].scatter(flight_time,TM_df['rh'],s=3, c=flight_time, cmap='icefire')
    axes[2].set_ylabel('Relative Humidity (%)')
    axes[2].set_xlabel('')
    axes[2].set_title('Relative Humidity over Flight Time')
    # Wind speed
    sc = axes[3].scatter(flight_time,TM_df['vv'],s=3, c=flight_time, cmap='icefire')
    axes[3].set_ylabel('Wind Speed (m/s)')

    # Wind direction
    sc = axes[4].scatter(flight_time,TM_df['dd'],s=3, c=flight_time, cmap='icefire')
    axes[4].set_ylabel('Wind Direction (degrees)')
    axes[4].set_xlabel('Flight Time')
    axes[4].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        
    plt.tight_layout()
    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok="True")
    figname="TMP_met_overview_"+rf+".png"
    fig.savefig(plot_path+figname,dpi=500,bbox_inches="tight")
    print("Fig. saved as:", plot_path+figname)
    return None

def TMP_met_plot_single_var_time_series(TM_ds,var,rf):
    """
    This routines plots time series of a specified meteorological variable 
    measured by TMP_met. 
    
    Input:
        Parameters
        ----------
        TM_ds : xr.Dataset
            dataset of TMP_met component.
        rf    : str
               research flight.
        var   : str
            variable to plot. Available vars: 
            "z_b","T","sonic_T","rh","vv","dd" 
        
        Returns
        -------
        None.
    """

    TM_df=TM_ds[["z_b","T","sonic_T","rh","vv","dd"]].to_dataframe()
    var_name={"z_b":"Barometric altitude (m)",
              "T":"Air temperature (K)",
              "sonic_T":"Sonic temperature",
              "rh": "Relative Humidity (%)",
              "vv":"Wind speed ($\mathrm{m\,s}^{-1}$)",
              "dd":"Wind direction (°)"}
    
    flight_time=TM_df.index
    matplotlib.rcParams.update({"font.size":16})
    # Set the style and colormap
    sns.set_style('whitegrid')
    cmap = sns.color_palette("icefire", as_cmap=True)
    # Create figure and axes objects
    fig, ax = plt.subplots(2, 1, figsize=(9, 8), height_ratios=[0.5,1],sharex=True)
    # Plot z_b (barometric height)
    ax[0].scatter(flight_time,TM_df['z_b'],s=3, c=flight_time,cmap=cmap)
    ax[0].set_ylabel('Barometric Height (m)')
    # Plot temperature
    sc = ax[1].scatter(flight_time,TM_df[var],s=2,marker="o",
                    c=flight_time, cmap='icefire',label="Sonde")
    if var=="T":
        ax[1].plot(flight_time,TM_df["sonic_T"],
            lw=1,color="grey",label="Sonic T")
    ax[1].set_ylabel(var_name[var])
    ax[1].set_xlabel('Flight Time')
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax[1].legend(fontsize=12)
    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok="True")
    figname="TMP_met_"+var+"_"+rf+".png"
    fig.savefig(plot_path+figname,dpi=500,bbox_inches="tight")
    return None

def TMP_met_plot_single_var_profile(TM_ds,var,rf):
    """
    This routines plots profiles (f(z_b)) of a specified meteorological variable 
    measured by TMP_met. 
    
    Input:
        Parameters
        ----------
        TM_ds : xr.Dataset
            dataset of TMP_met component.
        rf    : str
               research flight.
        var   : str
            variable to plot. Available vars: 
            "z_b","T","sonic_T","rh","vv","dd" 
        
        Returns
        -------
        None.
    """

    TM_df=TM_ds[["z_b","T","rh","vv","dd"]].to_dataframe()
    var_name={"z_b":"Barometric altitude (m)",
              "T":"Air temperature (K)",
              "rh": "Relative Humidity (%)",
              "vv":"Wind speed ($\mathrm{m\,s}^{-1}$)",
              "dd":"Wind direction (°)"}
    flight_time=TM_df.index
    matplotlib.rcParams.update({"font.size":16})
    # Set the style and colormap
    sns.set_style('whitegrid')
    cmap = sns.color_palette("icefire", as_cmap=True)
    # Create figure and axes objects
    fig, ax = plt.subplots(2, 1, figsize=(7, 8),
                           height_ratios=[0.3,1.0])
    # Plot z_b (barometric height)
    ax[0].scatter(flight_time,TM_df['z_b'],s=3, c=flight_time,cmap=cmap)
    ax[0].set_ylabel(' Height (m)')
    ax[0].set_ylim([0,900])
    ax[1].scatter(TM_df[var],
                  TM_df['z_b'],s=3, c=flight_time,cmap=cmap)
    ax[1].set_ylabel('Barometric Height (m)')
    ax[1].set_xlabel(var_name[var])
    ax[1].set_ylim([0,900])
    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok="True")
    figname="TMP_met_profile_"+var+"_"+rf+".png"
    fig.savefig(plot_path+figname,dpi=500,bbox_inches="tight")
    return None  
#------------------------------------------------------------------------------------#
# Turbulent meteorological probe (turbulent component)
#------------------------------------------------------------------------------------#
def plot_TMP_turb_windspeed(TU_ds,rf):
    """
    This routines plots time series of high resolution 
    wind speed as measured by TMP_turb. 
    
    Input:
        Parameters
        ----------
        TU_ds : xr.Dataset
            dataset of TMP_turb component.
        rf    : str
               research flight.
        
        Returns
        -------
        None.
    """
    ws_series=TU_ds["U"].to_series()
    date=ws_series.index.date[0]
    import seaborn as sns
    tu_fig=plt.figure(figsize=(8,4))
    matplotlib.rcParams.update({"font.size":12})
    resampled_u=pd.DataFrame()
    resampled_u["mean"] = ws_series.resample("10s").mean()
    resampled_u["std"]  = ws_series.resample("10s").std()
    
    ax1=tu_fig.add_subplot(111)
    ax1.plot(ws_series.index,ws_series.values,lw=.5,color="lightgrey",zorder=0)
    ax1.fill_between(resampled_u.index,resampled_u["mean"]-resampled_u["std"],
                     resampled_u["mean"]+resampled_u["std"],color="lightgreen",label="+-std",zorder=1)
    ax1.plot(resampled_u.index,resampled_u["mean"],label="10 s mean wind speed",color="darkgreen",lw=1,zorder=2)
    ax1.set_ylabel("Wind speed in $\mathrm{m\,s}^{-1}$")
    ax1.set_xlabel(str(date)+': Flight Time (UTC)')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax1.legend()
    ax1.set_ylim([0,8])
    print(ws_series.max())
    if float(ws_series.max())>8:
        ax1.set_ylim([0,10])
    if float(ws_series.max())>10:
        ax1.set_ylim([0,12])
    sns.despine(offset=10)
    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok="True")
    figname="TMP_turb_series_windspeed_"+rf+".png"
    tu_fig.savefig(plot_path+figname,dpi=500,bbox_inches="tight")
    
    return None

#------------------------------------------------------------------------------------#
# Broadband radiation probe
#------------------------------------------------------------------------------------#
def plot_BP_up_and_down_ward_terr_radiation(BP_ds,rf):
    """
    This function creates two subplots showing 
    the vertical profiles of F_up and F_down
    """
    matplotlib.rcParams.update({"font.size":18})
    df=BP_ds[["z_b","F_up","F_down","F_net"]].to_dataframe()
    df["ALT"]=df["z_b"]
    flight_time=df.index
    cmap="icefire"
    max_height=df["z_b"].max()    
    quick_scatter=plt.figure(figsize=(11,11))
    gs = GridSpec(2, 2, figure=quick_scatter,
        height_ratios=[0.35,1],wspace=.2)    
    ax0=quick_scatter.add_subplot(gs[0,:])
    # Plot z_b (barometric height)
    ax0.scatter(flight_time,df['z_b'],s=3, c=flight_time,cmap=cmap)
    ax0.set_ylabel('Barometric Height (m)')
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
         
    ax1=quick_scatter.add_subplot(gs[1,0])
    ax2=quick_scatter.add_subplot(gs[1,1])
    sub_fig_labels=["(b)","(c)"]
    # Create a proxy artist for the legend with a larger marker size
    # Scatter plots
    sns_cmap=sns.color_palette("icefire", as_cmap=True)
    ## Normalize time index for color mapping
    time=df.index.values
    norm_main = mcolors.Normalize(vmin=time.min(), vmax=time.max())
    # Plot 3: Balloon    
    ax1.scatter(df["F_up"],df["z_b"],c=df.index,cmap=sns_cmap,
                   norm=norm_main,s=3,label=rf) # BELUGA
    ax1.set_ylabel("Barometric height AGL (m)")
    ax2.scatter(df["F_down"],df["z_b"],c=df.index,
                cmap=sns_cmap,norm=norm_main,s=3,label=rf)
    ax2.set_yticks([0,100,200,300,400,500,600,700])
    #-------------------------------------------------------------------------#
    if max_height<500:
       ax1.set_ylim([0,500])
       ax2.set_ylim([0,500])
       ax1.set_yticks([0,250,500])
       ax2.set_yticks([250,500])
    elif 500<max_height<700:
        ax1.set_ylim([0,700])
        ax2.set_ylim([0,700])
        ax1.set_yticks([0,100,200,300,400,500,600,700])
        ax2.set_yticks([0,100,200,300,400,500,600,700])
    elif 700<max_height<800:
        ax1.set_ylim([0,800])
        ax2.set_ylim([0,800])
        ax1.set_yticks([0,200,400,600,800])
        ax2.set_yticks([0,200,400,600,800])
    else:
        ax1.set_ylim([0,1000])
        ax2.set_ylim([0,1000])
        ax1.set_yticks([0,200,400,600,800,1000])
        ax2.set_yticks([0,200,400,600,800,1000])
    
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
        ax2.spines[axis].set_linewidth(2)
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)
        ax2.yaxis.set_tick_params(width=2,length=6)
        ax2.xaxis.set_tick_params(width=2,length=6)
        
    ax1.text(0.03,0.95,sub_fig_labels[0],transform=ax1.transAxes,
             fontsize=20)
    ax2.text(0.03,0.95,sub_fig_labels[1],transform=ax2.transAxes,
             fontsize=20)
    ax1.set_xlabel('Upward thermal near-infrared \nirradiance ($\mathrm{W\,m}^{-2})$')
    ax2.set_xlabel('Downward thermal near-infrared \nirradiance \nRadiation ($\mathrm{W\,m}^{-2})$')
    f_up_p5,f_up_p95=df["F_up"].quantile([0.05,0.95])
    f_down_p5,f_down_p95=df["F_down"].quantile([0.05,0.95])
    
    ax1.set_xlim([f_up_p5//10*10-10,f_up_p95//10*10+10])
    ax2.set_xlim([f_down_p5//10*10-10,f_down_p95//10*10+10])
    
    plt.subplots_adjust(wspace=0.5,hspace=0.3)
    sns.despine(offset=10)
    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_name="Vertical_profiles_F_up_F_down_"+rf+".png"
    fname=plot_path+fig_name
    quick_scatter.savefig(fname,dpi=300,bbox_inches="tight")
    return None

def plot_BP_net_terr_radiation(BP_ds,rf):
    """
    this routines plot the net thermal near-infrared irradiance (F_down-F_up) 
    as a function of height.
    """
    matplotlib.rcParams.update({"font.size":16})
    df=BP_ds[["z_b","F_up","F_down","F_net"]].to_dataframe()
    df["ALT"]=df["z_b"]
    flight_time=df.index
    cmap="icefire"
    max_height=df["z_b"].max()    
    fig, ax = plt.subplots(2, 1, figsize=(7, 8),
                           height_ratios=[0.3,1.0])
    # Plot z_b (barometric height)
    ax[0].scatter(flight_time,df['z_b'],s=3, c=flight_time,cmap=cmap)
    ax[0].set_ylabel('Height (m)')
    ax[0].set_ylim([0,900])
    ax[0].xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M"))
    
    # Create a proxy artist for the legend with a larger marker size
    # Scatter plots
    sns_cmap=sns.color_palette("icefire", as_cmap=True)
    ## Normalize time index for color mapping
    time=df.index.values
    norm_main = mcolors.Normalize(vmin=time.min(), vmax=time.max())
    # Plot 3: Balloon    
    ax[1].scatter(df["F_net"],df["z_b"],c=df.index,
        cmap=sns_cmap,
        norm=norm_main,s=3,label=rf) 
    ax[1].set_ylabel("Barometric height (m)")
    #-------------------------------------------------------------------------#
    if max_height<500:
       ax[1].set_ylim([0,500])
       ax[1].set_yticks([0,250,500])   
    else:
        ax[1].set_ylim([0,900])
        ax[1].set_yticks([0,300,600,900])
    
        
    ax[1].set_xlabel('Net thermal near-infrared \nirradiance ($\mathrm{W\,m}^{-2})$')
    f_net_p5,f_net_p95=df["F_net"].quantile([0.05,0.95])
    ax[1].set_xlim([f_net_p5//10*10,f_net_p95//10*10+10])
    plt.subplots_adjust(wspace=0.5,hspace=0.3)
    sns.despine(offset=10)
    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    fig_name="Vertical_profiles_F_net_"+rf+".png"
    fname=plot_path+fig_name
    fig.savefig(fname,dpi=300,bbox_inches="tight")

#------------------------------------------------------------------------------------#
# Radiosonde
#------------------------------------------------------------------------------------#
def plot_radiosonde(ds,rf,date):
    """
    this routine plots a quicklook showing the vertical profiles of the 
    radiosonde measurements for single launches
    """
    matplotlib.rcParams.update({"font.size":12})
    fig, ax = plt.subplots(1, 8, figsize=(20, 9), sharey=True)
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'k', 'k', 'k']
    ds.plot.scatter(y='altitude', x='RH', ax=ax[0], edgecolor='k',
                    s=5, lw=0, color=colors[0])
    ds.plot.scatter(y='altitude', x='T', ax=ax[1], edgecolor='k',
                    s=5, lw=0, color=colors[1])
    ds.plot.scatter(y='altitude', x='Tu_Cap', ax=ax[2], edgecolor='k',
                    s=5, lw=0, color=colors[2])
    ds.plot.scatter(y='altitude', x='horizontal_speed', 
                    ax=ax[3], edgecolor='k', s=5, lw=0, color=colors[3])
    ds.plot.scatter(y='altitude', x='vertical_speed', ax=ax[4], 
                    edgecolor='k', s=5, lw=0, color=colors[4])
    ds.plot.scatter(y='altitude', x='distance', ax=ax[5],
                    edgecolor='k', s=5, lw=0, color=colors[5])
    ds.plot.scatter(y='altitude', x='longitude', ax=ax[6],
                    edgecolor='k', s=5, lw=0, color=colors[6])
    ds.plot.scatter(y='altitude', x='latitude', ax=ax[7],
                    edgecolor='k', s=5, lw=0, color=colors[7])

    for i, var in enumerate(['RH', 'T', 'Tu_Cap', 'horizontal_speed',
                             'vertical_speed', 'distance', 'longitude',
                             'latitude']):
        a = ax[i]
        a.text(0.05, 1, var, ha='left', va='center',
               transform=a.transAxes, color=colors[i],
               fontsize=10, weight='bold',)

    for a in ax:
        a.set_ylabel('')
        a.grid(True)
        a.spines[['top', 'right']].set_visible(False)
    ax[0].set_ylabel("Height (m)")
    date_str = ds['time'].dt.strftime('%Y-%m-%d').values[0]
    plt.savefig(os.getcwd()+"/../plots/"+rf+'//radiosonde_{date_str}.png',
                dpi=300)

def plot_flight_segmentation(df_initial,rf,probe="TMP_met"):
    """
    This routines plots in a colour-coded timeseries of BELUGA height 
    showing how the flight segments are categorised.
    
    """
    df=df_initial.copy()
    df["z_b"][df["z_b"]==np.Inf]=np.nan
    df["segments"]=df["flight_segments"].values
    # Do not separate between profile numbers anymore
    df["segments"][df_initial["flight_segments"]>=100]  = 'ascent'
    df["segments"][df_initial["flight_segments"]<=-100] = 'descent'
    df["segments"][df_initial["flight_segments"]==0]   = 'near ground'
    df["segments"][df_initial["flight_segments"]==1]   = 'const altitude'
    df["segments"][df_initial["flight_segments"]==2]   = 'peak'
    df["segments"][df_initial["flight_segments"]==3]   = 'max'
    #---------------------------------------------------------------------#
    matplotlib.rcParams.update({"font.size":16})
    #% BELUGA Plotting
    fig=plt.figure(figsize=(10,6))
    
    ax1=fig.add_subplot(111)
    ax1.plot(df["z_b"],lw=4,color="w")
    ax1.plot(df["z_b"],lw=2,color="grey")
    # Colour-code flight segmentation
    segmentation_classes=["ascent","descent","peak",
                          "near ground", "const altitude"]
    segmentation_legend=["ascent","descent","peak",
                          "near ground", "constant altitude"]
    segmentation_cls_colors=["blue","orange","red",
                             "sienna","green"]
    
    for s, seg in enumerate(segmentation_classes):
        seg_df=df["z_b"][df["segments"]==seg]
        ax1.scatter(seg_df.index,seg_df,marker="s",s=40,
                    color=segmentation_cls_colors[s],
                    zorder=2,label=segmentation_legend[s])
    
    # Add extra marker for maximum heights    
    max_max=df["z_b"].dropna().max()
    ax1.scatter(df["z_b"].idxmax(),max_max,
                marker="o",s=200,color="darkred",
                edgecolor="k",zorder=4)
    ax1.axhline(y=200, color="darkgrey", lw=3,ls="--")
    if max_max>600:
        ylim_max=900
    else:
        ylim_max=600
    
    ax1.text(df["z_b"].idxmax(),max_max+20,
             str(int(max_max))+" m")
    ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
    ax1.set_ylabel("Barometric height (m)")
    ax1.set_ylim([0,ylim_max])
    ax1.legend(fontsize=12)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    
    sns.despine(offset=10)
    
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    plot_path=os.getcwd()+"/../plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_name="BELUGA_"+probe+"_flight_segmentation_"+rf
    file_end=".png"
    fig_name+=file_end
    fname=plot_path+fig_name
    
    fig.savefig(fname, bbox_inches="tight")
    print("Figure saved as:", fname)

def show_data_quality(ds,probe_name):
    """
    This routine creates a bar plot showing the data quality of the instrument probe
    for the given research flight. It refers to the "data_quality" DataArray which is
    assigned during post-processing.
    Input:
        Parameters
        ----------
        ds : xr.Dataset
            dataset of instrument probe for research flight.
        probe_name : str
            name of probe given as plot title.

        Returns
        -------
        None.
    """
    import matplotlib.pyplot as plt
    quality_series=ds["quality_flag"].to_series()
    # total number of entries
    total = quality_series.shape[0]

    # Count how many entries fall into each category
    good_prop = (quality_series == 'good').sum() / total
    ok_prop   = (quality_series.str.startswith('ok_')).sum() / total
    bad_prop  = (quality_series.str.startswith('bad_')).sum() / total

    props = [good_prop, ok_prop, bad_prop]
    labels = ['Good', 'Ok', 'Bad']
    colors = ['darkgreen', 'orange', 'darkred']

    # Plot
    fig, ax = plt.subplots(figsize=(4, 1))
    #ax.bar(0, 1, bottom=0, width=0.5, color='white')  # optional: background
    ax.barh(0, props[0], height=0.5, color=colors[0],
            label=labels[0]+":"+str(int(round(good_prop*100)))+"%")
    ax.barh(0, props[1], left=props[0], height=0.5, color=colors[1],
            label=labels[1]+":"+str(int(round(ok_prop*100)))+"%")
    ax.barh(0, props[2], left=props[0]+props[1], height=0.5, color=colors[2],
            label=labels[2]+":"+str(int(round(bad_prop*100)))+"%")

    # Remove axes for a cleaner look
    ax.set_xticks([])
    ax.set_yticks([])

    # Optional: add labels or legend
    plt.legend(loc='upper left')

    plt.title(probe_name+" : Data Quality, relative frequency")

def plot_series_quality_flag(ds,probe="TMP_met"):  
    """
    This routine shows for the given flight altitude time series 
    the assigned quality flag.
    """
    quality_flag_plot=plt.figure(figsize=(12,9))
    ax1=quality_flag_plot.add_subplot(111)
    df=ds.to_dataframe()
    if probe=="BP":
        #bad_data  = pd.concat([df["z_b"].loc[df["quality_flag"]==3],
        #                      df["z_b"].loc[df["quality_flag"]==4]])
        bad_data  = df["z_b"].loc[df["quality_flag"]>2.5]
        ok_data   = df["z_b"].loc[df["quality_flag"].between(0.5,2.5)]
        good_data = df["z_b"].loc[df["quality_flag"]==0]
    if probe=="TMP_met":
        bad_data  = df["z_b"].loc[df["quality_flag"]>2.5]
        ok_data   = df["z_b"].loc[df["quality_flag"].between(0.5,2.5)]
        good_data = df["z_b"].loc[df["quality_flag"]==0]
    print("bad data:",bad_data.shape[0]," of ", 
          ds["z_b"].shape[0], " time steps")
    print("ok data:",ok_data.shape[0]," of ",
          ds["z_b"].shape[0], " time steps")
    print("good data:",good_data.shape[0]," of ",
          ds["z_b"].shape[0]," time steps")
    ax1.scatter(good_data.index,good_data.values,color="green",zorder=0)
    ax1.scatter(ok_data.index,ok_data.values,color="yellow",zorder=1)
    ax1.scatter(bad_data.index,bad_data.values,color="red",zorder=2)
    ax1.set_ylim([0,900])
    ax1.set_ylabel("Height (m)")
    ax1.set_xlabel(str(pd.DatetimeIndex(ds.time[0]).date)+" Time (UTC)")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    sns.despine(offset=10)
    #plt.savefig(