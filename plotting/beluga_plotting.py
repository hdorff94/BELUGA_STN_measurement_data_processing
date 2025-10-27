# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 20:29:24 2025

@author: u300737
"""
import os
import matplotlib
import matplotlib.pyplot as plt

import matplotlib.dates  as mdates
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors

import metpy.calc as mpcalc
import metpy.constants as mpconsts

import numpy as np
import pandas as pd
import seaborn as sns
import sys
import Performance
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
#import Performance
#performance=Performance.performance()

def quicklook_beluga_temp(df):
    """
    

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    df[["T_SENS_OLA_filtered","T_OLA_filtered","T_AMB_SENS_OLA_filtered",
        "T_filtered","TEMP"]].plot(ylim=[-25,5],ylabel="Temperature")
    return None

def plot_beluga_temp(df,rf,sensor1="TEMP",sensor2="T_SENS_OLA_filtered"):
    """
    This function creates a timeseries of the BELUGA temperature 
    measurements for an entire flight.  In contrast to the 
    quicklook routine only temperatures from two sensors are compared
    and more effort is put on plotting style.

    Parameters
    ----------
    df : pd.DataFrame
        BELUGA data containing temperature values to be plotted.
    rf : str
        BELUGA research flight name (RFxy)
    Returns
    -------
    None.

    """
    if sensor1=="TEMP":
        sensor1_name="Sonde_TEMP"
    else:
        # to be diversified in the future
        sensor1_name=sensor1
    #% BELUGA Plotting
    fig=plt.figure(figsize=(12,8))
    ax1=fig.add_subplot(111)
    ax1.plot(df[sensor1],lw=4,color="w")
    ax1.plot(df[sensor1],label=sensor1_name,lw=2,color="darkred")

    ax1.plot(df[sensor2],lw=4,color="w")
    ax1.plot(df[sensor2],label=sensor2,lw=2,color="orange")
    ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
    ax1.set_ylim([-25,0])
    ax1.set_ylabel("Temperature (°C)")
    ax12=ax1.twinx()
    ax1.legend(ncol=2)
    ax12.set_ylim([0,700])
    ax12.set_ylabel("Height (m)")
    ax12.plot(df["ALT"],color="grey",ls="--")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    sns.despine(offset=10)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_name="BELUGA_Temperature_timeseries_"+sensor1+"_"+sensor2
    file_end=".pdf"
    fig_name+=file_end
    fig.savefig(plot_path+fig_name, bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)

def plot_basic_flight_infos(df,max_vals):
    #% BELUGA Plotting
    fig=plt.figure(figsize=(12,8))
    rf=df.name
    ax1=fig.add_subplot(111)
    ax1.plot(df["ALT"],lw=4,color="w")
    ax1.plot(df["ALT"],lw=2,color="blue")
    ax1.scatter(max_vals.index,max_vals, marker="o", 
                s=100,color="red",edgecolor="k",zorder=3)
    max_max=df["ALT"].max()
    ax1.scatter(df["ALT"].idxmax(),max_max,marker="o",s=200,color="darkred",
                edgecolor="k",zorder=4)
    ax1.axhline(y=200, color="darkgrey", lw=3,ls="--")
    ax1.text(df["ALT"].idxmax(),max_max+20,
             str(int(max_vals.max()))+" m")
    ax1.text(0.1,0.9,"Number of profiles: "+str(2*len(max_vals)),
             transform=ax1.transAxes)
    ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
    ax1.set_ylim([0,900])
    ax1.set_ylabel("Height AGL (m)")
    ax1.set_title(rf)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    sns.despine(offset=10)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_name="BELUGA_Profiles_timeseries_"+rf
    file_end=".png"
    fig_name+=file_end
    fig.savefig(plot_path+fig_name,dpi=300, bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)

def find_profile_peaks(rf_df,min_height=200):
    from scipy.signal import find_peaks
    min_height=min_height
    alt_max_idx,_ = find_peaks(rf_df["ALT"],height=min_height,
        distance=500,width=2)
    alt_max_vals=rf_df["ALT"].iloc[alt_max_idx]
    
    return alt_max_vals 

def segment_flight_sections(rf_df,rate_threshold=0.7):
    
    # Classes of segmentation:
        # ascent, descent, peak, near-ground
    print("Perform flight segmentation:",rf_df.name)
    performance=Performance.performance()
    segmentation_classes=["ascent","descent","peak","max"
                          "near_ground", "constant_alt"]
    # find profile peaks as first separation
    profile_peaks=find_profile_peaks(rf_df)
    # get gradients of altitude
    alt_grad=rf_df["ALT"].diff()
    alt_grad_30s=alt_grad.rolling("30s").mean()
    
    rf_df["segments"]="none"
    alt_grad_threshold=rate_threshold
    asc_index=alt_grad_30s[alt_grad_30s>alt_grad_threshold].index
    dsc_index=alt_grad_30s[alt_grad_30s<-1*alt_grad_threshold].index
    con_index=alt_grad_30s[alt_grad_30s.between(-1*alt_grad_threshold,
                                                alt_grad_threshold)].index
    gnd_index=con_index.intersection(rf_df[rf_df["ALT"]<50].index)
    
    # Assign segment flags
    rf_df["segments"].loc[asc_index]="ascent"
    rf_df["segments"].loc[dsc_index]="descent"
    rf_df["segments"].loc[con_index]="const altitude"
    rf_df["segments"].loc[gnd_index]="near ground"
    
    # Update flight segmentation with added profile information
    # Merge flag ascent/descent as boolean
    added_rf_df=pd.DataFrame(data=rf_df[["ALT","segments"]],
        columns=["ALT","segments"],index=rf_df.index)
    added_rf_df['is_ascent_or_descent'] = added_rf_df['segments'].isin(
                        ['ascent', 'descent'])
    added_rf_df["DATETIME"]=added_rf_df.index
    # Assign a group id to each connected ascent/descent segment
    added_rf_df['profile_group'] = (~added_rf_df['is_ascent_or_descent']).cumsum()
    # Calculate height change over each group
    def compute_height_change(group):
        if group['is_ascent_or_descent'].any():
            height_diff = group['ALT'].iloc[-1] - group['ALT'].iloc[0]
            segment_length = len(group)  
            return pd.Series({
            'height_change': height_diff,
            'segment_length': segment_length
            })
        else:
            return pd.Series({'height_change': 0, 'segment_length': len(group)})

    segment_stats = added_rf_df.groupby('profile_group').apply(
                    compute_height_change)

    # Step 4: Merge segment stats back into the original DataFrame
    added_rf_df = added_rf_df.merge(segment_stats, on='profile_group')

    # Step 5: Flag segments longer than 100 m in height change as 'profile'
    def assign_profile(row):
        if row['is_ascent_or_descent']:
            if row['height_change'] > 100: 
                return 'profile_ascent'
            elif row['height_change'] < -100:
                return 'profile_descent'
            else:
                return row['segments']
        else:
            return row['segments']

    added_rf_df['segment_category'] = added_rf_df.apply(assign_profile, axis=1)

    # Profile numbers to each labeled 'profile' segment
    # First, create a boolean mask for 'profile' segments
    added_rf_df['is_profile'] = \
        added_rf_df['segment_category'].str.startswith('profile')
    
    # Then, assign a cumulative count to each profile
    added_rf_df['profile_number'] = (added_rf_df['is_profile'] & \
                (~added_rf_df['is_profile'].shift().fillna(False))).cumsum()
    
    # For all 'profile' segments, 'profile_number' indicates their id
    # For non-profile segments, 'profile_number' will be NaN or 0; set as needed
    added_rf_df.loc[~added_rf_df['is_profile'], 'profile_number'] = pd.NA
    splitted=added_rf_df["segment_category"].str.split("_")
    # Make column naming as e.g. "profile 0i ascent" 
    for i in range(len(splitted)):
        performance.updt(len(splitted),i)
        if len(splitted.iloc[i]) == 2:
            added_rf_df.at[i, "segment_category"] = splitted.iloc[i][0]+\
                " " +str(int(added_rf_df["profile_number"].iloc[i])).zfill(2) +\
                    " "+splitted.iloc[i][1]
    added_rf_df=added_rf_df.dropna(subset="profile_number")
    added_rf_df.index=added_rf_df["DATETIME"]            
    rf_df["segments"].loc[added_rf_df.index]=added_rf_df["segment_category"]
    #---------------------------------------------------------#
    for peak in profile_peaks.index:
        peak_idx=rf_df.index.get_loc(peak)
        rf_df["segments"].iloc[peak_idx-30:peak_idx+30]="peak"
    rf_df["segments"].loc[profile_peaks.index]="max"
    
    return rf_df,profile_peaks

def plot_RF_height(df,alt_var="ALT",ax_obj=""):
    rf=df.name
    sns_cmap=sns.color_palette("icefire", as_cmap=True)
    ## Normalize time index for color mapping
    time=df.index.values
    norm_main = mcolors.Normalize(vmin=time.min(), vmax=time.max())
    if ax_obj=="":
        fig=plt.figure(figsize=(9,6))
        ax_obj=fig.add_subplot(111)
    else:
        pass
    # Plot 3: Balloon    
    ax_obj.scatter(df.index,df[alt_var],c=df.index,cmap=sns_cmap,
                   norm=norm_main,s=6,label=rf)
    ax_obj.set_xlabel(str(df.index.date[0])+" Time (UTC)",fontsize=18)
    if alt_var=="ALT":
        ax_obj.set_ylabel("Height AGL (m)", fontsize=18)
    elif alt_var=="z_b":
        ax_obj.set_ylabel("Barom. height (m)", fontsize=18)
    else:
        ax_obj.set_ylabel("Height (m)", fontsize=18)
    
    ax_obj.tick_params(axis='both', which='major', labelsize=16)
    ylim_max=700
    ax_obj.set_yticks([0,350,ylim_max])
    
    if df[alt_var].max()>700:
        ylim_max=900
        ax_obj.set_yticks([0,300,600,ylim_max])
    ax_obj.set_ylim([0,ylim_max])
    ax_obj.text(0.1,0.9,rf,transform=ax_obj.transAxes,fontsize=20)
    #ax_obj.legend(loc="center left",fontsize=18)
    sns.despine(offset=10)
    for axis in ['bottom','left']:
        ax_obj.spines[axis].set_linewidth(2)
    ax_obj.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    
    ax_obj.yaxis.set_tick_params(width=2,length=6)
    ax_obj.xaxis.set_tick_params(width=2,length=6)

    ax_obj.spines['top'].set_visible(False)
    ax_obj.spines['right'].set_visible(False)
    
    ax_obj.text(0.01,0.9,"(a)",fontsize=20,transform=ax_obj.transAxes)
    
def plot_flight_segmentation(df,profile_peaks,add_peak_infos=False):
    #% BELUGA Plotting
    fig=plt.figure(figsize=(9,6))
    rf=df.name
    ax1=fig.add_subplot(111)
    ax1.plot(df["ALT"],lw=4,color="w")
    ax1.plot(df["ALT"],lw=2,color="grey")
    
    # Colour-code flight segmentation
    segmentation_classes=["ascent","descent","peak",
                          "near_ground", "const_alt"]
    
    segmentation_cls_colors=["blue","orange","red",
                             "sienna","green"]
    
    for seg,clr in zip(segmentation_classes,segmentation_cls_colors):
        seg_df=df["ALT"][df["segments"]==seg]
        ax1.scatter(seg_df.index,seg_df,marker="s",s=40,
                    color=clr,zorder=2,label=seg)
    # Add extra marker for maximum heights    
    ax1.scatter(profile_peaks.index,profile_peaks, marker="o", 
                s=100,color="red",edgecolor="k",zorder=3)
    max_max=df["ALT"].max()
    ax1.scatter(df["ALT"].idxmax(),max_max,
                marker="o",s=200,color="darkred",
                edgecolor="k",zorder=4)
    ax1.axhline(y=200, color="darkgrey", lw=3,ls="--")
    # Add profile peak infos
    
    if add_peak_infos:
        ax1.text(df["ALT"].idxmax(),max_max+20,
             str(int(max_max)+" m"))
        ax1.text(0.1,0.9,"Number of profiles: "+str(2*len(profile_peaks)),
             transform=ax1.transAxes)
    
    ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
    ax1.set_ylim([0,900])
    ax1.set_ylabel("Height AGL (m)")
    ax1.set_title(rf)
    ax1.legend()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    sns.despine(offset=10)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig_name="BELUGA_flight_segmentation_"+rf
    file_end=".png"
    fig_name+=file_end
    fig.savefig(plot_path+fig_name,dpi=300, bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)

def find_max_rf_height(rf_df):
    return rf_df["ALT"].max()
    
def basic_flight_information(cpgn_df,BELUGA_cls,min_height=200,
                             add_profile_infos_in_file=False):
    
    flight_infos=BELUGA_cls.flight_infos
    if add_profile_infos_in_file:
        flight_infos["No_of_profiles"] = np.nan
        flight_infos["Max_height"]     = np.nan
        
    for rf in flight_infos.index:
        print("Plot ", rf)
        rf_df=BELUGA_cls.open_met_data(rf=rf)
        rf_df.name=rf
        # Get profile maxima from research flight
        alt_max_vals = find_profile_peaks(rf_df,min_height=min_height)
        # Plotting
        plot_basic_flight_infos(rf_df, alt_max_vals)
        
        if add_profile_infos_in_file:    
            flight_infos["No_of_profiles"].loc[rf] = 2*len(alt_max_vals)
            flight_infos["Max_height"].loc[rf]     = find_max_rf_height(rf_df)
    if add_profile_infos_in_file:
        flight_info_fname=BELUGA_cls.main_path+"flight_schedule.csv"
        flight_infos.to_csv(flight_info_fname)
        print("flight info file saved as:", flight_info_fname)
    return alt_max_vals

def advanced_flight_information(cpgn_df,BELUGA_cls,min_height=200,
                             add_profile_infos_in_file=False,
                             plot_segments=True):
    
    merged_rfs=[]
    flight_infos=BELUGA_cls.flight_infos
    if add_profile_infos_in_file:
        flight_infos["No_of_profiles"] = np.nan
        flight_infos["Max_height"]     = np.nan
        
    for rf in flight_infos.index:
        rf_df=BELUGA_cls.open_met_data(rf=rf)
        # Drop interrupted indexes (NAT) and duplicates
        rf_df = rf_df.loc[rf_df.index.dropna()]
        rf_df = rf_df[~rf_df.index.duplicated()]
        rf_df.name=rf
        # Get profile maxima from research flight
        rf_df,profile_peaks = segment_flight_sections(rf_df,rate_threshold=0.1)
        # Plot flight segmentation
        if plot_segments:
            plot_flight_segmentation(rf_df, profile_peaks)
        if add_profile_infos_in_file:    
            flight_infos["No_of_profiles"].loc[rf] = 2*len(profile_peaks)
            flight_infos["Max_height"].loc[rf]     = find_max_rf_height(rf_df)
        BELUGA_cls.met_df=rf_df
        BELUGA_cls.save_met_data()
        merged_rfs.append(rf_df)
    #-------------------------------------------------------------------------#
    # Prepare output
    
    # merge Rfs in single DataFrame
    merged_rfs_df=pd.concat(merged_rfs)    
    
    #  
    if add_profile_infos_in_file:
        flight_info_fname=BELUGA_cls.main_path+"flight_schedule.csv"
        flight_infos.to_csv(flight_info_fname)
        print("flight info file saved as:", flight_info_fname)
    
    return flight_infos,merged_rfs_df

def plot_flight_segments_histogram(df):
    # Colour-code flight segmentation
    #Reference
    #segmentation_classes=["ascent","descent","peak","none"
    #                      "near_ground", "const_alt"]
    # 
    #segmentation_cls_colors=["blue","orange","red","grey"
    #                         "sienna","green"]
    seg_colors=["orange","blue","sienna","red","green","grey","darkred"]
    hist_fig=plt.figure(figsize=(7,4))
    # Use the figure to create axes
    ax_new = hist_fig.add_subplot(111)
    segment_flags=df["segments"]
    # Count the occurrences of each flag
    counts = segment_flags.value_counts()

    # Calculate the percentage
    percentages = (counts / len(segment_flags)) * 100

    # Plot the histogram
    percentages.plot(kind='bar', color= seg_colors,ax=ax_new)

    # Add labels and title
    ax_new.set_ylabel('Relative flight duration (%)')
    ax_new.set_xlabel('Flight segment flag')
    # Rotate x-axis labels by 45 degrees
    sns.despine(offset=10,ax=ax_new)
    ax_new.set_xticklabels(percentages.index, rotation=45)
    # Save figure
    fig_name = "hist_segments.pdf"
    fig_path = os.getcwd()+"/plots/"
    hist_fig.savefig(fig_path+fig_name, dpi=300,bbox_inches="tight")
    print("Figure saved as:", fig_path+fig_name)      
    return percentages

def plot_profile_infos_RFs(BELUGA_cls,min_height=200):
    import matplotlib.lines as mlines
    from scipy.signal import find_peaks
    
    matplotlib.rcParams.update({"font.size":24})
    
    flight_info=BELUGA_cls.flight_infos
    
    profile_info_fig=plt.figure(figsize=(16,9))
    ax1=profile_info_fig.add_subplot(111)
    rf_indices=np.arange(flight_info.shape[0])
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

def plot_tnir_histogram(cpgn_tnir):
    
    """
    Parameters
    ----------
    cpgn_tnir_df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """ 
    matplotlib.rcParams.update({"font.size":28})
    if not isinstance(cpgn_tnir, dict):
        import Measurement_Platforms
        BELUGA_cls=Measurement_Platforms.BELUGA
        timesteps=3600*4*25 #assuming 25 RFs, recording 4 hours with 1Hz res     
        
        cpgn_tnir_df=BELUGA_cls.create_artificial_tnirs(timesteps)
        h30_df=cpgn_tnir_df["30m"]
        hmax_var="h$_{\mathrm{max}}$"
        hmax_df=cpgn_tnir_df[hmax_var].dropna()
    else:
        h30_df  = cpgn_tnir["h30"]["F_net"]
        hmax_df = cpgn_tnir["peak"]["F_net"]
        heights = cpgn_tnir["peak"]["z_b"]
        print(heights.describe())
    
    # Plotting the histograms
    tnir_fig=plt.figure(figsize=(16,9))
    
    # Calculate histogram for 30 m height
    count_30m, bins_30m = np.histogram(h30_df, bins=np.linspace(-100,40,30),
                                           density=True)
    # Calculate histogram for hmax
    count_hmax, bins_hmax = np.histogram(hmax_df, bins=np.linspace(-100,40,30),
                                         density=True)
    # Calculate the bin centers for plotting
    bin_centers_30m = 0.5 * (bins_30m[1:] + bins_30m[:-1])
    bin_centers_hmax = 0.5 * (bins_hmax[1:] + bins_hmax[:-1])
    
    # Normalize the counts to percentage
    count_30m_percent = count_30m * np.diff(bins_30m) * 100
    count_hmax_percent = count_hmax * np.diff(bins_hmax) * 100
    
    # Plotting the histograms
    ax1=tnir_fig.add_subplot(111)
    
    ax1.bar(bin_centers_30m, count_30m_percent, 
            width=(bins_30m[1] - bins_30m[0]), color='slateblue', alpha=0.6,
            label='30m', edgecolor='black')
    
    ax1.bar(bin_centers_hmax, count_hmax_percent, 
            width=(bins_hmax[1] - bins_hmax[0]), color='orange', alpha=0.6,
            label="profile peaks", edgecolor='black')
    
    # Add labels and a title to the plot
    ax1.set_xlabel('Net irradiance $F_{\mathrm{net}}$ (W m$^{-2}$)')
    ax1.set_ylabel('Relative Occurence (%)')
    ax1.axhline(0, color='black',linewidth=0.8, ls='--')
    ax1.legend()
    
    # Show grid for better readability
    #ax1.grid(axis='y', linestyle='--', alpha=0.7)
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Show the plot
    ax1.set_xlim(-100, 50)
    #ax1.set_ylim(0, max(count_30m_percent.max(), 
    #                    count_hmax_percent.max()) + 5)
    ax1.set_ylim(0,25)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    sns.despine(ax=ax1,offset=10)
        
    fig_name="BELUGA_STN_TNIR_net_histogram_campaign"
    file_end=".png"
    fig_name+=file_end
    plot_path=os.getcwd()+"/../plots/campaign_overview/"
    os.makedirs(plot_path,exist_ok=True)
    
    tnir_fig.savefig(plot_path+fig_name,dpi=600,bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)
    return None

def plot_met_var_statistics(cpgn_df,as_scatter=True):
    matplotlib.rcParams.update({"font.size":24})
    stat_fig, (ax1,ax2,ax3,ax4) = plt.subplots(1, 4, figsize=(18,12),
                            gridspec_kw={'width_ratios': [1,2,2,2]},sharey=True)
    # Change height above sea level to height above ground
    #cpgn_df["ALT"]+=-31
    # Some u,v values cause very high wind speeds
    cpgn_df[cpgn_df["hor_vv"]>12]=np.nan
    cpgn_df[cpgn_df["ALT"]<0]=0.5
    altitudes=cpgn_df["ALT"].dropna()
    #
    #alt_hist= np.histogram(cpgn_df["ALT"].dropna()-31,range=(0,800),bins=16)

    # Define the bins and labels
    bins = np.arange(0, 850, 50)  # Bins from 0 to 800, exclusive of 800
    labels = [f"{i}-{i+50}" for i in range(0, 800, 50)]
    bin_midpoints = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]

    # Cut the altitudes into bins and count frequencies
    altitude_counts = pd.cut(altitudes, bins=bins, right=False)
    frequency = altitude_counts.value_counts().sort_index()

    # Create a horizontal bar plot
    #plt.figure(figsize=(10, 6))
    ax1.barh(bin_midpoints, frequency/3600, height=48, color='skyblue')
    #ax1.barh(frequency.values, width=50)
    ax1.set_yticks(np.linspace(0,800,5))
    ax1.set_xticks([0,5],["0","5"])
    ax1.text(1.25,20,str(round(frequency[0]/3600,1))+" h",color="white",
             fontweight="bold",zorder=3)
    ax1.set_xlim([0,7])
    ax1.set_ylim([0,800])
    ax1.set_xlabel("Flight hours (h)")
    ax1.set_ylabel("Height AGL (m)")
    ax1.text(0.1,0.97,"(a)",color="k",bbox=dict(facecolor='whitesmoke',
        edgecolor="black", boxstyle='round'),fontsize=22, 
             transform=ax1.transAxes)
    # Move y-ticks to the right spine
    plt_cpgn_df=cpgn_df[["TEMP","RH","ALT","hor_vv"]]
    plt_cpgn_df=plt_cpgn_df.dropna(how="any")
    if as_scatter:
        ax2.scatter(plt_cpgn_df["TEMP"],plt_cpgn_df["ALT"],color="salmon",s=2)
        ax3.scatter(plt_cpgn_df["RH"],plt_cpgn_df["ALT"],color="deepskyblue",s=2)
        ax4.scatter(plt_cpgn_df["hor_vv"],plt_cpgn_df["ALT"],color="teal",s=2)
    else:
        # Create 2D histograms for Temperature, RH, and Horizontal Wind Velocity
        hist_t  = ax2.hist2d(plt_cpgn_df["TEMP"], plt_cpgn_df["ALT"],
                             bins=(16, 16), cmap='Reds',density=True)
        hist_rh = ax3.hist2d(plt_cpgn_df["RH"], plt_cpgn_df["ALT"],
                             bins=(16, 16), cmap='Blues',density=True)
        hist_ws = ax4.hist2d(plt_cpgn_df["hor_vv"], plt_cpgn_df["ALT"], 
                             bins=(16, 16), cmap='Greens',density=True)
        cax_t = ax2.inset_axes([0.1, 0.9, 0.85, 0.03])
        cbar_t=stat_fig.colorbar(hist_t[3], ticks=[0,0.0002],
                                 cax=cax_t, orientation='horizontal')
        cbar_t.ax.set_xticklabels(['0', '2e-4'])
        # Add colorbar for Temperature
        cax_rh = ax3.inset_axes([0.1, 0.9, 0.85, 0.03])
        cbar_rh=plt.colorbar(hist_rh[3], cax=cax_rh, ticks=[0,.0001],
                             orientation='horizontal')  
        cbar_rh.ax.set_xticklabels(['0', '1e-4'])
        cbar_rh.set_label('Density', rotation=0, labelpad=10, ha='center')  # Set a title
        cbar_rh.ax.xaxis.set_label_position('top')
        # Add colorbar for Temperature
        cax_ws = ax4.inset_axes([0.1, 0.9, 0.85, 0.03])
        cbar_ws=plt.colorbar(hist_ws[3], cax=cax_ws, ticks=[0,.001],
                             orientation='horizontal')  
        cbar_ws.ax.set_xticklabels(['0', '1e-3'])
         
        #cbar_ws.set_label('Density', rotation=0, labelpad=25, ha='center')  # Set a title
        
        # Add colorbar for Temperature
        # Adjust the size of the colorbars if necessary
        
    #ax1.barh(np.linspace(0,800,17),cpgn_df["ALT"].dropna()-31, color="grey")
    ax2.set_xlim([-30,0])
    ax3.set_xlim([0,100])
    ax4.set_xlim([0,8])
    ax1.set_ylim([0,800])
    ax2.set_ylim([0,800])
    ax3.set_ylim([0,800])
    ax4.set_ylim([0,800])

    ax2.set_xlabel("Temperature (°C)")
    ax3.set_xlabel("RH (%)")
    ax4.set_xlabel("Wind speed ($\mathrm{m\,s}^{-1})$")
    ax4.set_xticks([0,4,8])
    ax4.yaxis.tick_right()  # Move the ticks to the right
    ax4.yaxis.set_label_position("right")  # Move the y-axis label to the right as well
    ax4.set_yticks([0,200,400,600,800])
    ax4.set_yticklabels(["0","200","400","600","800"])
    for axis in ['bottom','left']:
        ax2.spines[axis].set_linewidth(2)
        ax3.spines[axis].set_linewidth(2)
    for axis in ['bottom','right']:
        ax1.spines[axis].set_linewidth(2)
        ax4.spines[axis].set_linewidth(2)

    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    ax2.yaxis.set_tick_params(width=2,length=6)
    ax2.xaxis.set_tick_params(width=2,length=6)
    ax3.yaxis.set_tick_params(width=2,length=6)
    ax3.xaxis.set_tick_params(width=2,length=6)
    ax4.yaxis.set_tick_params(width=2,length=6)
    ax4.xaxis.set_tick_params(width=2,length=6)

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    ax4.spines['top'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    
    ax2.text(0.05,0.97,"(b)",color="k",bbox=dict(
        facecolor='whitesmoke',edgecolor="black", boxstyle='round'),
        fontsize=22,transform=ax2.transAxes)
    ax3.text(0.05,0.97,"(c)",color="k",bbox=dict(facecolor='whitesmoke',
        edgecolor="black", boxstyle='round'),
        fontsize=22,transform=ax3.transAxes)
    ax4.text(0.05,0.97,"(d)",color="k",bbox=dict(facecolor='whitesmoke',
        edgecolor="black", boxstyle='round'),
        fontsize=22,transform=ax4.transAxes)


    ax4.yaxis.tick_right()  # Move the ticks to the right
    ax4.yaxis.set_label_position("right")  # Move the y-axis label to the right as well
    plt.subplots_adjust(wspace=0.5)
    fig_name="fig05_BELUGA_STN_campaign_stats_"
    if as_scatter:
        fig_name+="_scatter"
    else:
        fig_name+="_hist"
    file_end=".png"
    fig_name+=file_end
    plot_path=os.getcwd()+"/plots/campaign_overview/"
    os.makedirs(plot_path,exist_ok=True)
    stat_fig.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
    print("Figure saved as: ",plot_path+fig_name)