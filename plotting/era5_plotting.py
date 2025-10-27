# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 09:02:35 2025

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

from metpy.units import units

class StadiaStamen(cimgt.Stamen):
    def _image_url(self, tile):
         x,y,z = tile
         url =  f"https://tiles.stadiamaps.com/tiles/stamen_terrain_background/{z}/{x}/{y}.jpg?api_key=0963bb5f-6e8c-4978-9af0-4cd3a2627df9"
        #https://tiles.stadiamaps.com/tiles/stamen_terrain_background/{z}/{x}/{y}.png?api_key={API_KEY}"
         return url

def map_STN(era_grid, sea_ice_ds, carra_seaice, main_path,
            with_miniplot=False,start_hour=12):
    from scipy.spatial import ConvexHull
    
    fig_name="fig01_STN_map_region"
    #-------------------------------------------------------------------------#
    #Collocation to Station Noord coordinates
    STN_coords={"lat": 81.60,
                "lon": -16.65}
    
    lon_lim=[-22,-14]
    lat_lim=[81,82.]
    central_lon=-30
    
    #location_plot=plt.figure(figsize=(10,10))
    #era_grid=np.meshgrid(era_ds.longitude,era_ds.latitude)
    # Open CARRA grid
    carra_grid              = xr.open_dataset(main_path+"CARRA_grid.nc")
    carra_lat               = carra_grid.latitude.values.flatten()
    carra_lon               = carra_grid.longitude.values.flatten()
    carra_grid_df           = pd.DataFrame(data=np.nan, index=range(len(carra_lat)),
                                       columns=["lat","lon"])
    carra_grid_df["lat"]    = carra_lat
    carra_grid_df["lon"]    = carra_lon
    carra_grid_df["lon"][carra_grid_df["lon"]>180]-=360
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
    
    # Plotting
    domain_fig,ax1=plt.subplots(1,1,figsize=(10,10),subplot_kw={"projection":\
        ccrs.NorthPolarStereo(central_longitude=central_lon)})
    
    ax1.plot([hull_edges["lon"][0],hull_edges["lon"][1]],
              [hull_edges["lat"][0],hull_edges["lat"][1]],
             color='coral', marker="s",markersize=10, markeredgecolor="k",
             linestyle="--", lw=4,
             transform=ccrs.Geodetic(),label="CARRA West \ndomain")
    ax1.plot([hull_edges["lon"][1],hull_edges["lon"][2]],
              [hull_edges["lat"][1],hull_edges["lat"][2]],
             color='coral', marker="s",markeredgecolor="k",
             markersize=10, linestyle="--", lw=4,
             transform=ccrs.Geodetic())
    ax1.plot([hull_edges["lon"][1],hull_edges["lon"][2]],
              [hull_edges["lat"][1],hull_edges["lat"][2]],
             color='coral', marker="s",markersize=10,markeredgecolor="k",mew=2,
             linestyle="--", lw=4,transform=ccrs.Geodetic())
    ax1.plot([hull_edges["lon"][2],hull_edges["lon"][3]],
              [hull_edges["lat"][2],hull_edges["lat"][3]],
             color='coral', marker="s",markersize=10, markeredgecolor="k",mew=2,
             linestyle="--", lw=4,transform=ccrs.Geodetic())
    ax1.plot([hull_edges["lon"][3],hull_edges["lon"][4]],
              [hull_edges["lat"][3],hull_edges["lat"][4]],
             color='coral', marker="s",markersize=10,markeredgecolor="k",mew=2,
             linestyle="--", lw=4,transform=ccrs.Geodetic())
    
    small_extent=[lon_lim[0],lon_lim[1],lat_lim[0],lat_lim[1]]#[-40,30,55,90]
    big_extent  =[-57,-7,50,89]
    ax1.set_extent(big_extent, crs=ccrs.Geodetic())
    stamen_terrain = StadiaStamen('terrain-background')
    ax1.add_image(stamen_terrain, 5)
    
    ax1.coastlines(resolution="50m")
    ax1.add_feature(cartopy.feature.BORDERS)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                          x_inline=False, y_inline=False)
    #-------------------------------------------------------------------------#
    # Add locations as text in plots
    text_color=(0.031, 0.188, 0.419)
    
    x1=STN_coords["lon"]
    y1=STN_coords["lat"]
    STN_color="mediumorchid"
    ax1.text(x1-35, y1-1, "Station \nNord (STN)",
             fontsize=16,transform=ccrs.PlateCarree(),
             color=STN_color,ha="center",bbox=dict(
                 facecolor='whitesmoke',edgecolor="black", boxstyle='round'))
    ax1.plot(x1, y1, 'X',color=STN_color, markersize=20,markeredgecolor="w",
         mew=3,transform=ccrs.PlateCarree(),zorder=20)
    
    #-------------------------------------------------------------------------#
    # Sea ice regions
    orig_map = plt.cm.get_cmap('Blues') # getting the original colormap using cm.get_cmap() function
    reversed_map = orig_map.reversed()
    sea_text_size=16
    # # open old sea ice grid
    old_sea_ice_path="C:\\Users\\u300737\Desktop\\Desktop_alter_Rechner//"+\
        "PhD_UHH_WIMI\Work\GIT_Repository\hamp_processing_py//"+\
            "hamp_processing_python\Flight_Data\HALO_AC3\sea_ice\\"
    exemplaric_file="asi-AMSR2-n6250-20220316-v5.4.nc"
    grid_file=xr.open_dataset(old_sea_ice_path+exemplaric_file)
    #
    C1=ax1.pcolormesh(grid_file.lon,grid_file.lat,
           np.array(sea_ice_ds["z"][:]), 
           transform=ccrs.PlateCarree(),
               cmap=reversed_map)    
    
    gl.top_labels   = False
    gl.right_labels = False
    gl.left_labels  = True
    gl.xlabel_style = {'size': 18}
    gl.ylabel_style = {'size': 18}
    gl.xlocator= mticker.FixedLocator([-100,-90,-80,-70,-60,-50,
                                       -40,-30,-20,-10,0,10,20,30])
    gl.ylocator = mticker.FixedLocator([55, 60,65,70,75,80])
    
    
    y_box=[lat_lim[0],lat_lim[1],lat_lim[1],lat_lim[0],lat_lim[0]]
    x_box=[lon_lim[0],lon_lim[0],lon_lim[1],lon_lim[1],lon_lim[0]]
    
    
    ax1.plot([x_box[2],35],[y_box[2],68],lw=2,ls="-",color="k",
             transform=ccrs.Geodetic())
    ax1.plot([x_box[0],0],[y_box[0],55],lw=2,ls="-",color="k",
             transform=ccrs.PlateCarree())
    ax1.plot(x_box,y_box,color="k",ls="--",lw=2,marker="s",
             transform=ccrs.Geodetic(),zorder=4)
    ax1.legend(loc="right")
    #######################################################################
    mini_plot_axes=[.86, .25, .5, .5]
    if not with_miniplot:
        ax2= plt.axes(mini_plot_axes,alpha=0)
        ax2.patch.set_facecolor('none')
        ax2.spines[['right', 'bottom','left','top']].set_visible(False)
        ax2.axes.xaxis.set_ticklabels([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.axes.yaxis.set_ticklabels([])
        
    else:
        ax2= plt.axes(mini_plot_axes,
                                    projection=ccrs.NorthPolarStereo(
                                        central_longitude=central_lon)) 
        
        #######################################################################
        gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                              x_inline=False, y_inline=False,zorder=30)
        
        gl2.bottom_labels = True
        gl2.top_labels    = False
        gl2.left_labels   = False
        gl2.right_labels  = True
        gl2.xlabel_style = {'size': 18}
        gl2.ylabel_style = {'size': 18}
        
        #gl2.ylocator = mticker.FixedLocator([81,81.25,81.5,81.75])
        #gl2.xlocator= mticker.FixedLocator([-18,-16,-14])
        ax2.coastlines(resolution="10m")
        ax2.add_feature(cartopy.feature.BORDERS)
        #Add your line modifications here
        ax1.spines['geo'].set_linewidth(2)
        ax2.spines['geo'].set_linewidth(2)
        ax2.spines['geo'].set_linestyle('--')
        
        ax1.text(0.02,0.95,"(a)",color="k",fontsize=22,
                 transform=ax1.transAxes,bbox={"boxstyle":"round", 
                    "facecolor":"white","alpha":.9}, zorder=30)
        ax2.text(0.05,0.90,"(b)",color="k",fontsize=22,bbox={"boxstyle":"round", 
           "facecolor":"white","alpha":.9},
                 transform=ax2.transAxes,zorder=30)
        #70N - 85N 20W - 30E
        ax2.set_extent([small_extent[0]+1,small_extent[1],
                        small_extent[2],small_extent[3]],
                       crs=ccrs.Geodetic())
        ax2.add_image(stamen_terrain, 6)
        # AMSR-2 sea-ice is too coarse for this illustration
        #C2=ax2.pcolormesh(grid_file.lon,grid_file.lat,
        #       np.array(sea_ice_ds["z"][:]), 
        #       transform=ccrs.PlateCarree(),
        #           cmap=reversed_map)    
        # Use CARRA instead
        carra_masked=carra_seaice.copy()#where(carra_seaice != 0)
        C2=ax2.pcolormesh(carra_seaice.longitude,carra_seaice.latitude,
               np.array(carra_masked.values[:]*100), 
               transform=ccrs.PlateCarree(),
                   cmap=reversed_map)    
        
        cbar_ax = domain_fig.add_axes([0.96, 0.78, 0.3, 0.04])
        cbar    = domain_fig.colorbar(C2, cax=cbar_ax,orientation="horizontal")
        cbar.set_ticks([0,50,100])
        cbar.ax.set_title("AMSR-2 mean \nsea-ice fraction / %",fontsize=16)
        cbar.ax.tick_params(labelsize=16)

        ax2.scatter(STN_coords["lon"],STN_coords["lat"], marker="X",s=700,
                    color=STN_color,edgecolor="k", label="STN",
                    transform=ccrs.PlateCarree())
        ax2.scatter(grid_points[:, 1], grid_points[:, 0],color='orange', s=2,
                    lw=0.5,label='CARRA grid',
                    transform=ccrs.PlateCarree())
        
        ax2.scatter(era_grid[0],era_grid[1],s=10,color="black",
                    label="ERA5 grid",transform=ccrs.PlateCarree())
        
        ax2.gridlines()
        ax2.legend(loc="lower left")
        fig_name+="_with_subregion"
    file_end=".png"
    fig_name+=file_end
    plot_path=os.getcwd()+"/plots/campaign_overview/"
    os.makedirs(plot_path,exist_ok=True)
    domain_fig.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)
    
    
    return None

def plot_era5_theta_in_region(era_ds,rf, with_miniplot=True,start_hour=12,
                              theta_height=1000,wind_height=1000):
    
    # Theta 850 hPa
    theta_levels=np.linspace(240,270,31)
    z_levels=np.linspace(5000,5600,13)
    era5_theta=era_ds["theta"].sel({"pressure_level":theta_height})
    # Geopotential 500 hPa
    era5_z_500=era_ds["z"].sel({"pressure_level":500})/9.81
    era5_u=era_ds["u"].sel({"pressure_level":wind_height})
    era5_v=era_ds["v"].sel({"pressure_level":wind_height})
    #-------------------------------------------------------------------------#
    #Collocation to Station Noord coordinates
    STN_coords={"lat": 81.60,
                "lon": -16.65}
    lon_lim=[-20,-14]
    lat_lim=[81,82]
    central_lon=-17
    
    location_plot=plt.figure(figsize=(10,10))
    era_grid=np.meshgrid(era_ds.longitude,era_ds.latitude)
    fig, ax1=plt.subplots(1,1,figsize=(10,10),subplot_kw={"projection":
        ccrs.NorthPolarStereo(central_longitude=central_lon)})
    small_extent=[lon_lim[0],lon_lim[1],lat_lim[0],lat_lim[1]]#[-40,30,55,90]
    big_extent  =[-30,10,75,85]
    ax1.set_extent(big_extent, crs=ccrs.Geodetic())
    stamen_terrain = StadiaStamen('terrain-background')
    ax1.add_image(stamen_terrain, 5)
    
    C1=ax1.contourf(era5_theta.longitude,era5_theta.latitude,
                 era5_theta.values[12,:,:],levels=theta_levels,
                 cmap="Spectral_r",transform=ccrs.PlateCarree(),alpha=0.8)
    
    ax1.coastlines(resolution="50m")
    ax1.add_feature(cartopy.feature.BORDERS)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                          x_inline=False, y_inline=False)
    #-------------------------------------------------------------------------#
    # Add locations as text in plots
    text_color=(0.031, 0.188, 0.419)
    
    x1=STN_coords["lon"]
    y1=STN_coords["lat"]
    STN_color="mediumorchid"
    ax1.text(x1-12, y1 + 0.3, "Station \nNoord (STN)",
             fontsize=16,transform=ccrs.PlateCarree(),
             color=STN_color,ha="center",bbox=dict(
                 facecolor='whitesmoke',edgecolor="black", boxstyle='round'))
    ax1.plot(x1, y1, '*',color=STN_color, markersize=20,markeredgecolor="k",
         transform=ccrs.PlateCarree(),zorder=20)
    
    # Sea regions
    sea_text_size=16
    #C1=ax1.pcolormesh(sea_ice.lon,sea_ice.lat,
    #       np.array(sea_ice[:]), 
    #       transform=ccrs.PlateCarree(),
    #           cmap=reversed_map)
    
    gl.bottom_labels = True
    gl.top_labels   = True
    gl.right_labels = False
    gl.xlabel_style = {'size': 18}
    gl.ylabel_style = {'size': 18}
    gl.xlocator= mticker.FixedLocator([-40,-20,0,20,40])
    gl.ylocator = mticker.FixedLocator([76,78, 80,82,84])
    
    
    y_box=[lat_lim[0],lat_lim[1],lat_lim[1],lat_lim[0],lat_lim[0]]
    x_box=[lon_lim[0],lon_lim[0],lon_lim[1],lon_lim[1],lon_lim[0]]
    ax1.plot([16,40],[74.0,71],lw=1,ls="-",color="k",
             transform=ccrs.PlateCarree())        
    
    ax1.plot([x_box[2],19],[y_box[2],76],lw=2,ls="-",color="k",
             transform=ccrs.PlateCarree())
    ax1.plot([x_box[0],-8],[y_box[0],74],lw=2,ls="-",color="k",
             transform=ccrs.PlateCarree())
    ax1.plot(x_box,y_box,color="k",ls="--",lw=2,marker="s",
             transform=ccrs.PlateCarree(),zorder=4)
    ax1.contour(era5_z_500.longitude,era5_z_500.latitude,
                era5_z_500[12,:,:].values,
                levels=z_levels,linewidths=5,colors="white",zorder=10,
                transform=ccrs.PlateCarree())
    
    Cs=ax1.contour(era5_z_500.longitude,era5_z_500.latitude,
                era5_z_500[12,:,:].values,
                levels=z_levels,linewidths=3,cmap="Greys",
                zorder=12,transform=ccrs.PlateCarree(),label="Geopotential")
    ax1.clabel(Cs, z_levels, fmt='%04d m', fontsize=20, inline=True,colors="k")
    #ax1.legend(loc="upper right",fontsize=20)
    cbar_ax = fig.add_axes([0.2, -0.05, 0.25, 0.04])
    cbar=fig.colorbar(C1, cax=cbar_ax,orientation="horizontal")
    cbar.set_ticks([240,250,260,270])
    cbar.ax.set_title("ERA5 $\Theta_{\mathrm{"+str(theta_height)+"hPa}}$ (K)",
                      fontsize=18)
    #######################################################################
    fig_name="ERA5_theta_"+str(theta_height)+"hPa_region_wind_"+\
                    str(wind_height)+"hPa"
    if not with_miniplot:
        ax2= plt.axes([.55, -.15, .39, .45],alpha=0)
        ax2.patch.set_facecolor('none')
        ax2.spines[['right', 'bottom','left','top']].set_visible(False)
        ax2.axes.xaxis.set_ticklabels([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.axes.yaxis.set_ticklabels([])
        
    else:
        ax2= plt.axes([.55, -.15, .39, .45],
                                    projection=ccrs.NorthPolarStereo(
                                        central_longitude=central_lon)) 
        
        #######################################################################
        #
        gl2 = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                              x_inline=False, y_inline=False,zorder=30)
        
        gl2.bottom_labels = True
        gl2.top_labels    = False
        gl2.left_labels   = False
        gl2.right_labels  = True
        gl2.xlabel_style = {'size': 18}
        gl2.ylabel_style = {'size': 18}
        
        gl2.ylocator = mticker.FixedLocator([81,81.25,81.5,81.75])
        gl2.xlocator= mticker.FixedLocator([-18,-16,-14])
        
        ax2.coastlines(resolution="10m")
        ax2.add_feature(cartopy.feature.BORDERS)
        #Add your line modifications here
        ax1.spines['geo'].set_linewidth(2)
        ax2.spines['geo'].set_linewidth(2)
        ax2.spines['geo'].set_linestyle('--')
        
        ax1.text(0.02,0.95,"(a)",color="k",fontsize=22,
                 transform=ax1.transAxes,zorder=30)
        ax2.text(0.05,0.90,"(b)",color="k",fontsize=22,
                 transform=ax2.transAxes,zorder=30)
        #70N - 85N 20W - 30E
        ax2.set_extent([small_extent[0]+1,small_extent[1],
                        small_extent[2],small_extent[3]],
                       crs=ccrs.Geodetic())
        ax2.add_image(stamen_terrain, 6)
        
        #ax2.legend(loc="lower left",fontsize=18,ncol=1)
        ax2.scatter(STN_coords["lon"],STN_coords["lat"], marker="*",s=1000,
                    color=STN_color,edgecolor="k", label="STN",
                    transform=ccrs.PlateCarree())
        ax2.scatter(era_grid[0],era_grid[1],c=era5_theta[12,:,:].values,
                    s=200,cmap=matplotlib.cm.Spectral_r,
                    vmin=theta_levels[0],vmax=theta_levels[-1],
                    edgecolor="k",label="ERA5 grid",
                    transform=ccrs.PlateCarree())
        lon_step=2
        lat_step=1
        quiver_lon=np.array(era5_u["longitude"][::lon_step])
        quiver_lat=np.array(era5_u["latitude"][:])
        u=np.array(era5_u[12,:,::lon_step])
        v=np.array(era5_v[12,:,::lon_step])
        
        quiver=ax2.quiver(quiver_lon,quiver_lat,
                          u,v,color="mintcream",
                          edgecolor="k",lw=1,
                          scale=8,scale_units="inches",
                          pivot="mid",width=0.02,zorder=20,
                          transform=ccrs.PlateCarree())
        
        ax2.gridlines()
        ax2.legend()
        fig_name+="_with_subregion"
    file_end=".png"
    fig_name+=file_end
    plot_path=os.getcwd()+"/plots/"+rf+"/"
    os.makedirs(plot_path,exist_ok=True)
    fig.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)
    return None

def plot_campaign_overview(era_STN_dict,BELUGA_cls):
    # Campaign overview 
    
    cmpgn_ts_plot=plt.figure(figsize=(16,12))
    ax1=cmpgn_ts_plot.add_subplot(311)
    ax1.spines['top'].set_visible(False)
    ax1.plot(era_STN_dict["theta"][1000], color="darkred",
             label="$\Theta_{\mathrm{1000\,hPa}}$")
    ax1.plot(era_STN_dict["theta"][925], color="salmon",
             ls="-",label="$\Theta_{\mathrm{925\,hPa}}$")

    ax1.set_ylim([240,280])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax1.set_ylabel("$\Theta$ (K)")
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(2)
    ax1.yaxis.set_tick_params(width=2,length=6)
    ax1.xaxis.set_tick_params(width=2,length=6)
    ax1.set_xlim([pd.Timestamp("2024-03-16"),pd.Timestamp("2024-04-16")])
    ax12=ax1.twinx()
    ax12.plot(era_STN_dict["wspeed"][1000],color="grey",
              ls="--",label="$WS_{1000\,hPa}$")
    ax12.set_ylabel("Wind speed (m/s)")
    ax12.set_ylim([0,10])
    for f in BELUGA_cls.flight_infos.index:
        ax1.axvline(pd.Timestamp(BELUGA_cls.flight_infos["Start Time"].loc[f]),
                    color="grey",ls=":")
    for l in range(BELUGA_cls.LLJ_flights.shape[0]):
        label=""
        if l==BELUGA_cls.LLJ_flights.shape[0]-1:
            label="LLJ RF"
        index_val=BELUGA_cls.LLJ_flights.index[l]
        ax1.axvspan(pd.Timestamp(BELUGA_cls.LLJ_flights["Start Time"].loc[index_val]),
                   pd.Timestamp(BELUGA_cls.LLJ_flights["End Time"].loc[index_val]),
                   facecolor='lightsteelblue',alpha=.5,label=label)

    for c in range(BELUGA_cls.cloud_flights.shape[0]):
        label=""
        if c==BELUGA_cls.cloud_flights.shape[0]-1:
            label="clear to cloudy RF"
        index_val=BELUGA_cls.cloud_flights.index[c]
        ax1.axvspan(pd.Timestamp(BELUGA_cls.cloud_flights["Start Time"].loc[index_val]),
                   pd.Timestamp(BELUGA_cls.cloud_flights["End Time"].loc[index_val]),
                   alpha=0.5, facecolor='grey',
                   label=label)    
    for n in range(BELUGA_cls.night_flights.shape[0]):
        label=""
        if n==BELUGA_cls.night_flights.shape[0]-1:
            label="Night to day RF"
        index_val=BELUGA_cls.night_flights.index[n]
        ax1.axvspan(pd.Timestamp(BELUGA_cls.night_flights["Start Time"].loc[index_val]),
                   pd.Timestamp(BELUGA_cls.night_flights["End Time"].loc[index_val]),
                   alpha=0.5, facecolor='orange',label=label)    

    ax1.spines['top'].set_visible(False)

    ax1.legend(ncol=5)
    ax12.legend(loc="upper right")
    ax12.spines['top'].set_visible(False)

    #-----------------------------------------------------------------------------#
    ax2=cmpgn_ts_plot.add_subplot(312)
    x,y=np.meshgrid(era_STN_dict["wspeed"].index,era_STN_dict["wspeed"].columns)

    C1=ax2.contourf(x,y,era_STN_dict["wspeed"].T.values)
    cax1 = cmpgn_ts_plot.add_axes([0.95, 0.4, 0.01, 0.2])
    cbar=plt.colorbar(C1, cax=cax1)
    cbar.set_label('Wind speed (m/s)', labelpad=1)
    ax2.set_ylabel("Pressure (hPa)")
    ax2.set_ylim([1000,700])
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax2.set_xlim([pd.Timestamp("2024-03-16"),pd.Timestamp("2024-04-16")])

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    for axis in ['bottom','left']:
        ax2.spines[axis].set_linewidth(2)
    ax2.yaxis.set_tick_params(width=2,length=6)
    ax2.xaxis.set_tick_params(width=2,length=6)

    #-----------------------------------------------------------------------------#
    ax3=cmpgn_ts_plot.add_subplot(313)
    C2=ax3.contourf(x,y,era_STN_dict["r"].T.values,cmap="Blues")
    cax2 = cmpgn_ts_plot.add_axes([0.95, 0.1, 0.01, 0.2])

    cbar2=plt.colorbar(C2, cax=cax2)
    cbar2.set_label('RH (%)')

    ax3.set_ylabel("Pressure (hPa)")
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))

    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_ylim([700,1000])
    ax3.set_xlim([pd.Timestamp("2024-03-16"),pd.Timestamp("2024-04-16")])

    for axis in ['bottom','left']:
        ax3.spines[axis].set_linewidth(2)
    ax3.yaxis.set_tick_params(width=2,length=6)
    ax3.xaxis.set_tick_params(width=2,length=6)
    ax3.invert_yaxis()
    ax3.set_xlabel("Date in 2024 (MM-DD)")
    plt.subplots_adjust(wspace=0.4,hspace=0.4)

    fig_name="ERA5_STN_campaign_series.png"
    plot_path=os.getcwd()+"/plots/campaign_overview/"
    os.makedirs(plot_path,exist_ok=True)

    cmpgn_ts_plot.savefig(plot_path+fig_name, dpi=600, bbox_inches="tight")
    print("Figure saved as:", plot_path+fig_name)
