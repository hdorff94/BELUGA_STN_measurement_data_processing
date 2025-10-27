# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 23:12:39 2025

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
import metpy.calc as mpcalc
import metpy.constants as mpconsts
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

import seaborn as sns

from metpy.units import units

class StadiaStamen(cimgt.Stamen):
    def _image_url(self, tile):
         x,y,z = tile
         url =  f"https://tiles.stadiamaps.com/tiles/stamen_terrain_background/{z}/{x}/{y}.jpg?api_key=0963bb5f-6e8c-4978-9af0-4cd3a2627df9"
         return url

def plot_campaign_overview(carra_STN_dict,BELUGA_cls,moisture="RH",
                           include_theta_contours=False,theta_var="theta"):
    """
    This routine plots the CARRA data interpolated onto the STN location over 
    the period of the STN BELUGA measurement campaign by timeseries for 
    different meteorological quantities in a multipanel plot. 
    Illustrated are potential temperature for two levels as a proxy for 
    stability and near-surface wind speeds (a). Furthermore the RF categorised
    by defined ABL transition events are specified as vertical lines. 

    Parameters
    ----------
    carra_STN_dict : dict
        CARRA vertical data interpolated to STN location.
    BELUGA_cls : class
        BELUGA data class object for executing all functions of BELUGA.
    moisture : str, optional
        Moisture label to be plotted, this can be either relative humidity (RH)
        or specific humidity (q). The default is "RH".
    include_theta_contours : boolean, optional
        Switch for either line plot of selected theta heights (see manuscript)
        or contourplots of temperature variable as for RH and wind speed below.
        Default is False.
    
    theta_var: str, optional
        Variable name of temperature quantity to be selected if 
        include_theta_contours is True. Default is "theta".
    Returns
    -------
    None.

    """
    # Campaign overview 
    matplotlib.rcParams.update({"font.size":22})
    cmpgn_ts_plot=plt.figure(figsize=(16,12))
    ax1=cmpgn_ts_plot.add_subplot(311)
    ax1.spines['top'].set_visible(False)
    if not include_theta_contours:
        ax1.plot(carra_STN_dict["theta"][50], color="w",
             lw=3)
        ax1.plot(carra_STN_dict["theta"][500], color="w",
             ls="-",lw=3)

        ax1.plot(carra_STN_dict["theta"][50], color="darkred",
             label="$\Theta_{\mathrm{50\,m}}$",zorder=2)
        ax1.plot(carra_STN_dict["theta"][500], color="salmon",
             ls="-",label="$\Theta_{\mathrm{500\,m}}$",zorder=2)

        ax12=ax1.twinx()
        ax12.plot(carra_STN_dict["wspeed"][50],color="w",
                  ls="-",lw=3)
        ax12.plot(carra_STN_dict["wspeed"][50],color="grey",
                  ls="--",label="$WS_{50\,m}$",zorder=2)
        ax12.set_ylabel("Wind speed ($\mathrm{m\,s}^{-1}$)")
        ax12.set_ylim([0,10])
        for f in BELUGA_cls.flight_infos.index:
            ax1.axvspan(pd.Timestamp(
                BELUGA_cls.flight_infos["Start Time"].loc[f]),
                pd.Timestamp(BELUGA_cls.flight_infos["End Time"].loc[f]),
                facecolor="grey",alpha=.5)
            
        for l in range(BELUGA_cls.LLJ_flights.shape[0]):
            label=""
            if l==BELUGA_cls.LLJ_flights.shape[0]-1:
                label="LLJ"
            index_val=BELUGA_cls.LLJ_flights.index[l]
            ax1.axvspan(
                pd.Timestamp(BELUGA_cls.LLJ_flights["Start Time"].loc[index_val]),
                pd.Timestamp(BELUGA_cls.LLJ_flights["End Time"].loc[index_val]),
                facecolor='mediumseagreen',alpha=.5,label=label)

        for c in range(BELUGA_cls.cloud_flights.shape[0]):
            label=""
            if c==BELUGA_cls.cloud_flights.shape[0]-1:
                label="clear$\leftrightarrow$cloudy"
            index_val=BELUGA_cls.cloud_flights.index[c]
            ax1.axvspan(
                pd.Timestamp(BELUGA_cls.cloud_flights["Start Time"].loc[index_val]),
                pd.Timestamp(BELUGA_cls.cloud_flights["End Time"].loc[index_val]),
                alpha=0.5, facecolor='steelblue',
                label=label)    
        
        for n in range(BELUGA_cls.night_flights.shape[0]):
            label=""
            if n==BELUGA_cls.night_flights.shape[0]-1:
                label="day$\leftrightarrow$night"
            index_val=BELUGA_cls.night_flights.index[n]
            
            ax1.axvspan(pd.Timestamp(
                BELUGA_cls.night_flights["Start Time"].loc[index_val]),
                pd.Timestamp(BELUGA_cls.night_flights["End Time"].loc[index_val]),
                alpha=0.5, facecolor='orange',label=label)    

        ax1.spines['top'].set_visible(False)

        ax1.legend(ncol=6,fontsize=18,bbox_to_anchor=(0.04,1.02),
                   bbox_transform=ax1.transAxes)
        ax12.legend(loc="upper right",fontsize=16.5)
        ax12.spines['top'].set_visible(False)

        
        ax1.set_ylabel("$\Theta$ (K)")
        ax1.set_ylim([240,280])
    else:
        x,y=np.meshgrid(carra_STN_dict[theta_var].index,
                        carra_STN_dict[theta_var].columns)

        C0=ax1.contourf(x,y,carra_STN_dict[theta_var].T.values,
                        levels=np.linspace(240,270,31),cmap="Spectral_r",
                        extend="both")
        cax0 = cmpgn_ts_plot.add_axes([0.95, 0.65, 0.01, 0.2])
        cbar=plt.colorbar(C0, cax=cax0)
        cbar.set_label(theta_var+'(K)', labelpad=1)
        cbar.set_ticks([240,250,260,270])
        ax1.set_ylabel("Height AGL (m)")
        ax1.set_ylim([0,500])
        
    
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    
    # x-axis
    ax1.set_xticks(["2024-03-24","2024-03-26","2024-03-28","2024-03-30",
                "2024-04-01","2024-04-03","2024-04-05","2024-04-07",
                "2024-04-09","2024-04-11"])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
   
    ax1.set_xlim([pd.Timestamp("2024-03-24"),pd.Timestamp("2024-04-11")])
    #-----------------------------------------------------------------------------#
    ax2=cmpgn_ts_plot.add_subplot(312)
    x,y=np.meshgrid(carra_STN_dict["wspeed"].index,
                    carra_STN_dict["wspeed"].columns)

    C1=ax2.contourf(x,y,carra_STN_dict["wspeed"].T.values,vmin=0,vmax=15,
                    extend="max")
    cax1 = cmpgn_ts_plot.add_axes([0.95, 0.4, 0.01, 0.2])
    cbar=plt.colorbar(C1, cax=cax1)
    cbar.set_label('Wind speed ($\mathrm{m\,s}^{-1}$)', labelpad=1)
    cbar.set_ticks([0,5,10,15])
    ax2.set_ylabel("Height AGL (m)")
    ax2.set_ylim([0,500])
    ax2.set_xticks(["2024-03-24","2024-03-26","2024-03-28","2024-03-30",
                    "2024-04-01","2024-04-03","2024-04-05","2024-04-07",
                    "2024-04-09","2024-04-11"])
    
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax2.set_xlim([pd.Timestamp("2024-03-24"),pd.Timestamp("2024-04-11")])

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    for axis in ['bottom','left']:
        ax2.spines[axis].set_linewidth(2)
    ax2.yaxis.set_tick_params(width=2,length=6)
    ax2.xaxis.set_tick_params(width=2,length=6)

    #-------------------------------------------------------------------------#
    if moisture=="rh":
        moisture_label = "RH"
        moisture_unit  = " (%)"
        moisture_cmap  = "Blues"
    elif moisture=="q":
        moisture_label = "q"
        moisture_unit  = " ($\mathrm{g\,kg}^{-1}$"
        moisture_cmap  = "Blues_r"
    ax3=cmpgn_ts_plot.add_subplot(313)
    C2=ax3.contourf(x,y,carra_STN_dict[moisture].T.values,cmap=moisture_cmap,
                    vmin=0,vmax=100,extend="max")
    cax2 = cmpgn_ts_plot.add_axes([0.95, 0.13, 0.01, 0.2])

    cbar2=plt.colorbar(C2, cax=cax2)
    cbar2.set_label(moisture_label+moisture_unit)
    cbar2.set_ticks([20,40,60,80,100])
    
    ax3.set_ylabel("Height AGL (m)")
    ax3.set_xticks(["2024-03-24","2024-03-26","2024-03-28","2024-03-30",
                    "2024-04-01","2024-04-03","2024-04-05","2024-04-07",
                    "2024-04-09","2024-04-11"])
    
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_ylim([0,500])
    ax3.set_xlim([pd.Timestamp("2024-03-24"),pd.Timestamp("2024-04-11")])
    for axis in ['bottom','left']:
        ax3.spines[axis].set_linewidth(2)
    ax3.yaxis.set_tick_params(width=2,length=6)
    ax3.xaxis.set_tick_params(width=2,length=6)
    ax3.set_xlabel("Date in 2024 (MM-DD)")
    # Make vertical ticks fixed
    ax2.set_yticks([0,250,500])
    ax3.set_yticks([0,250,500])
    ax1.text(0.01,0.9,"(a)",color="k",transform=ax1.transAxes,
             fontsize=16.5,bbox=dict(facecolor='whitesmoke',
                 edgecolor="black", boxstyle='round'))
    ax2.text(0.01,0.9,"(b)",color="k",transform=ax2.transAxes,
             fontsize=16.5,bbox=dict(facecolor='whitesmoke',
                       edgecolor="black", boxstyle='round'))
    ax3.text(0.01,0.9,"(c)",color="k",transform=ax3.transAxes,
             fontsize=16.5,bbox=dict(facecolor='whitesmoke',
                 edgecolor="black", boxstyle='round'))
    
    plt.subplots_adjust(wspace=0.4,hspace=0.4)

    fig_name="CARRA_STN_campaign_series.png"
    if include_theta_contours:
        fig_name=theta_var+"_"+fig_name
    plot_path=os.getcwd()+"/plots/campaign_overview/"
    os.makedirs(plot_path,exist_ok=True)

    cmpgn_ts_plot.savefig(plot_path+fig_name, dpi=600, bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)

def synoptic_carra_map(CARRA_cls,carra_ds,main_path,
        dates=['2024-03-23', '2024-03-30', '2024-04-05', '2024-04-10'],
        theta_height=850,wind_height=1000):
    """
    This routine plots maps of Theta_e at 850 hPa for four specific days during
    the STN BELUGA measurement campaign in Spring 2024 as 2x2 multipanels,
    referring to the CARRA west domain data. Isolines of 500 hPa geopotential
    meters are superimposed. 

    Parameters
    ----------
    CARRA_cls : class
        class object of the CARRA reanalysis.
    carra_ds : xr.Dataset
        CARRA dataset for the campaign period.
    main_path : str
        main path to navigate towards plot data paths.
    theta_height : int, optional
        Pressure height to consider Theta_e from. The default is 850.
    wind_height : int, optional
        Pressure height to consider wind speed. The default is 1000.
        Note, currently wind data is not shown and neglected.

    Returns
    -------
    None.

    """
    
    # Define the list of dates
    
    if carra_ds==None:
        # Just get the grid region
        from scipy.spatial import ConvexHull
        #---------------------------------------------------------------------#
        #Collocation to Station Noord coordinates
        STN_coords={"lat": 81.60,
                "lon": -16.65}
    
        lon_lim=[-22,-14]
        lat_lim=[81,82.]
        # Open CARRA grid
        carra_grid_df=CARRA_cls.get_carra_grid()
        grid_points             = carra_grid_df.values
        # Step 1: Calculate the convex hull
        hull = ConvexHull(grid_points)

        # Step 2: Get the vertices of the convex hull for plotting
        hull_points = grid_points[hull.vertices]
        hull_edges=pd.DataFrame(data=np.nan,index=range(5),
                                columns=["lat","lon"])
        hull_edges["lat"].iloc[0:4]=hull_points[0:4, 0]
        hull_edges["lon"].iloc[0:4]=hull_points[0:4,1]
        hull_edges["lat"][4]=hull_edges["lat"].iloc[0]
        hull_edges["lon"][4]=hull_edges["lon"].iloc[0]
    #-------------------------------------------------------------------------#
    #Collocation to Station Noord coordinates
    STN_coords={"lat": 81.60,
                "lon": -16.65}
    x1=STN_coords["lon"]
    y1=STN_coords["lat"]
    STN_color="mediumorchid"
    # Lat/Lon
    lon_lim=[-20,-14]
    lat_lim=[81,82]
    central_lon=-30
    theta_levels=np.linspace(250,290,21)
    
    big_extent  =[-57,-7,50,89]
    stamen_terrain = StadiaStamen('terrain-background')
        
    # Prepare the figure and axes
    fig, axs = plt.subplots(2, 2, figsize=(16, 12),
                subplot_kw={"projection":\
                    ccrs.NorthPolarStereo(central_longitude=central_lon)})
        
    # Set common features for the plots
    panel_labels=["(a)","(b)","(c)","(d)"]
    
    for i, ax in enumerate(axs.flatten()):
        
        
        # Open CARRA data
        if dates==['2024-03-23', '2024-03-30', '2024-04-05', '2024-04-10']:
            is_additional=False
            # Define model data month
            if dates[i][5:7]=="03":
                carra_fname="CARRA_West_p_levels_00_12_March_2024.nc"
            else:
                carra_fname="CARRA_West_p_levels_00_12_April_2024.nc"
            
            carra_file=main_path+carra_fname
            carra=xr.open_dataset(carra_file)
            carra_ds=carra.sel({"valid_time":dates[i]+" 12:00"})
            carra_ds["pres"]=carra_ds["pressure_level"].copy()
            carra_ds=CARRA_cls.calculate_carra_theta_e(var_to_use=carra_ds)
            carra_ds=CARRA_cls.calculate_carra_wspeed_and_wdir(var_to_use=carra_ds)
            carra_ds["longitude"]=carra_ds["longitude"].where(
                carra_ds["longitude"]<180,carra_ds["longitude"]-360)
        else:
            is_additional=True
            theta_levels=np.linspace(260,300,21)
            
            import glob
            CARRA_path=os.getcwd()+"/../../CARRA_STN/"
            CARRA_fname="CARRA_west_STN*"
            CARRA_files=glob.glob(CARRA_path+CARRA_fname)
            carra_ds=xr.open_mfdataset(CARRA_files,combine="nested",
                           concat_dim="valid_time")
            carra_ds=carra_ds.sel({"valid_time":dates[i]+" 12:00"})
            
            carra_ds["pres"]=carra_ds["pressure_level"].copy()
            carra_ds=CARRA_cls.calculate_carra_theta_e(var_to_use=carra_ds)
            carra_ds=CARRA_cls.calculate_carra_wspeed_and_wdir(var_to_use=carra_ds)
            carra_ds["longitude"]=carra_ds["longitude"].where(
                carra_ds["longitude"]<180,carra_ds["longitude"]-360)

        theta_e  = carra_ds["theta_e"].sel({"pressure_level":theta_height})
        geopot = carra_ds["z"].sel({"pressure_level":500})
        
        # Plot map
        ax.coastlines(resolution="50m")
        ax.add_feature(cartopy.feature.BORDERS)
        ax.set_extent(big_extent, crs=ccrs.Geodetic())
        ax.add_image(stamen_terrain, 2) #5)
        ## Iterate over the dates to create the subpanels
        # Extract data and Plotting theta
        lat_1d=np.array(carra_ds.latitude).flatten()
        lon_1d=np.array(carra_ds.longitude).flatten()
        
        #Plot Theta e
        theta_1d=np.array(theta_e.values).flatten()
        C1 = ax.scatter(lon_1d, lat_1d, c=theta_1d,vmin=theta_levels[0],
            vmax=theta_levels[-1],cmap="Spectral_r",s=2,
            transform=ccrs.PlateCarree())
        
        
        if i == 0:
            ax.text(x1-25, y1 - 2, "STN",
            fontsize=14,transform=ccrs.PlateCarree(),color=STN_color,
            ha="center",bbox=dict(facecolor='whitesmoke',edgecolor="black",
                                  boxstyle='round',alpha=0.6))
        
        ax.plot(x1, y1, '*',color=STN_color, markersize=20,markeredgecolor="k",
             transform=ccrs.PlateCarree(),zorder=20)
        
        ax.text(0.05,0.91,panel_labels[i],color="k",fontsize=20,
                     transform=ax.transAxes,zorder=30,
                     bbox=dict(facecolor='whitesmoke',edgecolor="black",
                                           boxstyle='round',alpha=.6))
        
        ax.text(0.25,0.04,dates[i],fontsize=20,transform=ax.transAxes,
                zorder=30,bbox=dict(facecolor="whitesmoke",edgecolor="black",
                                    boxstyle="round"))
        
        if not carra_ds:
            ax.plot([hull_edges["lon"][0],hull_edges["lon"][1]],
                      [hull_edges["lat"][0],hull_edges["lat"][1]],
                     color='coral', marker="s",markersize=10, 
                     markeredgecolor="k",linestyle="--", lw=4,
                     transform=ccrs.Geodetic())
            ax.plot([hull_edges["lon"][1],hull_edges["lon"][2]],
                      [hull_edges["lat"][1],hull_edges["lat"][2]],
                     color='coral', marker="s",markeredgecolor="k",
                     markersize=10, linestyle="--", lw=4,
                     transform=ccrs.Geodetic())
            ax.plot([hull_edges["lon"][1],hull_edges["lon"][2]],
                      [hull_edges["lat"][1],hull_edges["lat"][2]],
                     color='coral', marker="s",markersize=10,
                     markeredgecolor="k",mew=2,linestyle="--",
                     lw=4,transform=ccrs.Geodetic())
            ax.plot([hull_edges["lon"][2],hull_edges["lon"][3]],
                      [hull_edges["lat"][2],hull_edges["lat"][3]],
                     color='coral', marker="s",markersize=10, 
                     markeredgecolor="k",mew=2,linestyle="--",
                     lw=4,transform=ccrs.Geodetic())
            ax.plot([hull_edges["lon"][3],hull_edges["lon"][4]],
                      [hull_edges["lat"][3],hull_edges["lat"][4]],
                     color='coral', marker="s",markersize=10,
                     markeredgecolor="k",mew=2,linestyle="--",
                     lw=4,transform=ccrs.Geodetic())
        
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                              x_inline=False, y_inline=False)
        
        if i==0:    
            gl.bottom_labels = False
            gl.top_labels    = True
            gl.right_labels  = False
        if i==1:
            gl.bottom_labels = False
            gl.top_labels    = True
            gl.right_labels  = True
            gl.left_labels   = False
        if i==2:
            gl.bottom_labels = True
            gl.top_labels    = False
            gl.right_labels  = False
            gl.left_labels   = True            
        if i==3:
            gl.bottom_labels = True
            gl.top_labels    = False
            gl.right_labels  = True
            gl.left_labels   = False
        
        gl.xlabel_style = {'size': 18}
        gl.ylabel_style = {'size': 18}
        gl.xlocator= mticker.FixedLocator([-100,-90,-80,-70,-60,-50,
                                           -40,-30,-20,-10,0,10,20,30])
        gl.ylocator = mticker.FixedLocator([55, 60,65,70,75,80])
    
        # Plotting isolines Geopot gpdm
        z_levels=np.linspace(4800, 5600, 17)
        lon=carra_ds["longitude"].values
        lat=carra_ds["latitude"].values
        z=geopot.values/10
        Cs=ax.contour(lon, lat,z,levels=z_levels,linestyles="--",
            colors='k', linewidths=1.5,transform=ccrs.PlateCarree())
        ax.clabel(Cs,z_levels[::2], fmt='%04d', fontsize=12,
                   inline=True,colors="k")
    
        # To be included for wind field illustration by quiver plots
        #    lon_step=2
        #    lat_step=1
        #    quiver_lon=np.array(era5_u["longitude"][::lon_step])
        #    quiver_lat=np.array(era5_u["latitude"][:])
        #    u=np.array(era5_u[12,:,::lon_step])
        #    v=np.array(era5_v[12,:,::lon_step])
            
        #    quiver=ax.quiver(quiver_lon,quiver_lat,
        #                      u,v,color="mintcream",
        #                      edgecolor="k",lw=1,
        #                      scale=8,scale_units="inches",
        #                      pivot="mid",width=0.02,zorder=20,
        #                      transform=ccrs.PlateCarree())
        # Add a color bar for theta
        
    cbar_ax = fig.add_axes([0.4, 0.02, 0.25, 0.02])
    cbar=fig.colorbar(C1, cax=cbar_ax,orientation="horizontal",
                      extend="both")
    cbar_ticks=[250,260,270,280,290]
    if is_additional:
        cbar_ticks=[260,270,280,290,300]
    cbar.set_ticks(cbar_ticks)
    cbar.ax.tick_params(labelsize=18)
    cbar_ax.set_xlabel("$\Theta_{e,\mathrm{"+str(theta_height)+\
                       "hPa}}$ (K)",fontsize=18)
    
    plt.subplots_adjust(wspace=-0.45,hspace=.04)
    #######################################################################
    if dates==['2024-03-23', '2024-03-30', '2024-04-05', '2024-04-10']:
        fig_name="CARRA_theta_"+str(theta_height)+"hPa_region_wind_"+\
                    str(wind_height)+"hPa_campaign_overview"
        plot_path=os.getcwd()+"/plots/campaign_overview/"
        
    elif dates==["2024-04-23","2024-04-24","2024-04-25","2024-04-26"]:
        fig_name="CARRA_theta_"+str(theta_height)+"hPa_region_wind_"+\
                    "hPa_end_April_overview"
        plot_path=os.getcwd()+"/../../CARRA_STN/"
    elif dates==["2024-04-29","2024-04-30","2024-05-01","2024-05-02"]:
        fig_name="CARRA_theta_"+str(theta_height)+"hPa_region_wind_"+\
                    "hPa_begin_May_overview"
        plot_path=os.getcwd()+"/../../CARRA_STN/"
    file_end=".png"
    fig_name+=file_end
    fig.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)
    
