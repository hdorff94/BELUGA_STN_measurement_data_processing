# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 17:43:30 2025

@author: u300737
"""
import glob
import os
import sys
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(1,os.getcwd()+"/plotting/")
sys.path.insert(2,os.getcwd()+"/src/analysing/")

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import matplotlib
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt

import metpy.calc as mpcalc
import metpy.constants as mpconsts
from metpy.units import units


import Measurement_Platforms
import Grid_data_STN

#-----------------------------------------------------------------------------#
# Definitions
matplotlib.rcParams.update({"font.size":14})
#-----------------------------------------------------------------------------#
# Switches for data sources to plot
do_introductory_plot = False # Simple map of Station Noord location
do_campaign_overview = False # Campaign overview

do_beluga_plots      = False
do_era5_plots        = False #ERA5
compare_era5_beluga  = False
do_carra_plots       = True  #CARRA
compare_carra_beluga = True
compare_icon_beluga  = False
#-----------------------------------------------------------------------------#
STN_coords={"lat": 81.60,
            "lon": -16.65}
#-----------------------------------------------------------------------------#
# Flight definitions
flight_list={"RF12":"20240401",
             "RF13":"20240401",
             "RF14":"20240401",
             "RF1312":"20240401",
             "RF15":"20240402",
             "RF1516":"20240402",
             "RF141312":"20240401",
             "RF22":"20240407",
             "RF26":"20240410"}
flight="RF22"
main_path=os.getcwd()+"/BELUGA_data/"
#-----------------------------------------------------------------------------#
# Open and allocate classes
Platforms_STN_cls=Measurement_Platforms.Measurement_Platforms_STN(
    main_path=main_path)
BELUGA_cls=Measurement_Platforms.BELUGA(Platforms_STN_cls)
BELUGA_cls.get_flight_information()

#print(BELUGA_cls.flight_infos)
if not do_campaign_overview:
    # Open Beluga data
    if not flight=="RF141312" and not flight=="RF1516":
        df=BELUGA_cls.open_met_data(rf=flight)
    elif flight=="RF141312":
        df_rf12=BELUGA_cls.open_met_data(rf="RF12")
        df_rf13=BELUGA_cls.open_met_data(rf="RF13")
        df_rf14=BELUGA_cls.open_met_data(rf="RF14")
        df=pd.concat([df_rf12,df_rf13,df_rf14])
        
        df["TEMP"][df["TEMP"]>0]=np.nan
    
        # Calculate Auxiliary meteorological quantities
        df = BELUGA_cls.calc_wspeed(df)
        df = BELUGA_cls.calc_theta(df)
        df = BELUGA_cls.calc_theta_e(df)
        df = BELUGA_cls.calc_q_from_rh(df)
    elif flight=="RF1516":
        df_rf15=BELUGA_cls.open_met_data(rf="RF15")
        df_rf16=BELUGA_cls.open_met_data(rf="RF16")
        df=pd.concat([df_rf15,df_rf16])
        
        df["TEMP"][df["TEMP"]>0]=np.nan
    
    # Calculate Auxiliary meteorological quantities
    df = BELUGA_cls.calc_wspeed(df)
    df = BELUGA_cls.calc_theta(df)
    df = BELUGA_cls.calc_theta_e(df)
    df = BELUGA_cls.calc_q_from_rh(df)
        
        # single specific day
if not flight=="RF1312" and not flight=="RF141312" and not flight=="RF1516":
    start_hour = BELUGA_cls.flight_infos["Start Time"][flight].hour
    end_hour   = BELUGA_cls.flight_infos["End Time"][flight].hour+2 # for correct indexing
elif flight=="RF1312":
    start_hour                = BELUGA_cls.flight_infos.loc["RF12"][\
                                        "Start Time"].hour
    end_hour                  = BELUGA_cls.flight_infos.loc["RF13"][\
                                        "End Time"].hour+2
elif flight=="RF141312":
    start_hour                = BELUGA_cls.flight_infos.loc["RF12"][\
                                        "Start Time"].hour
    end_hour                  = BELUGA_cls.flight_infos.loc["RF14"][\
                                        "End Time"].hour+2
elif flight=="RF1516":
    start_hour                = BELUGA_cls.flight_infos.loc["RF15"][\
                                        "Start Time"].hour
    end_hour                  = BELUGA_cls.flight_infos.loc["RF16"][\
                                        "End Time"].hour+1
#%% All campaign measurements
if do_campaign_overview:
    cpgn_df=BELUGA_cls.open_all_met_data(level="L2")
    if not "hor_vv" in cpgn_df.columns:
        cpgn_df=BELUGA_cls.calc_wspeed(cpgn_df)
    cpgn_df["TEMP"]-=273.15
    cpgn_df[cpgn_df["TEMP"]>0]=np.nan
#%% BELUGA Plotting
else:
    if do_beluga_plots:
        import beluga_plotting
        print("Plot Beluga met data.")
        beluga_plotting.quicklook_beluga_temp(df)
        beluga_plotting.plot_beluga_temp(df,flight)
    else:
        pass
#%% Gridded data Reanalyses
Grid_data_cls     = Grid_data_STN.Gridded_data(main_path=os.getcwd())
ERA5              = Grid_data_STN.ERA5(Grid_data_cls)
CARRA             = Grid_data_STN.CARRA(Grid_data_cls)
ICON              = Grid_data_STN.ICON(Grid_data_cls)

#%% Introductory plots
if do_introductory_plot:
    import era5_plotting
    # AMSR-2 Sea ice
    SEAICE          = Grid_data_STN.AMSR2_seaice(Grid_data_cls)
    SEAICE.open_sea_ice()
    sea_ice_ds      = SEAICE.sea_ice
    CARRA_seaice_ds = CARRA.open_sea_ice_2024()
    carra_seaice    = CARRA_seaice_ds["siconc"]
    ERA5.open_ERA5_single_day(BELUGA_cls,flight)
    era_ds=ERA5.era5_day_ds
    era_grid=np.meshgrid(era_ds.longitude,era_ds.latitude)
    era5_plotting.map_STN(era_grid,sea_ice_ds,carra_seaice,
                          main_path,with_miniplot=True)
#%% ERA5 Plotting
if do_era5_plots:
    import era5_plotting
    ERA5.open_ERA5_single_day(BELUGA_cls,flight)
    era_ds=ERA5.era5_day_ds
    era5_plotting.plot_era5_theta_in_region(era_ds,flight,
                    start_hour=start_hour,with_miniplot=False)
    era5_plotting.plot_era5_theta_in_region(era_ds,flight,
                    start_hour=start_hour,with_miniplot=True)
    if do_campaign_overview:
        ERA5.open_era5_plevel_period_data(take_as_major_era5_file=False)
        STN_ERA5_dict={}
        STN_ERA5_dict["theta"] = pd.DataFrame(
            data=ERA5.STN_era5_spring["theta"].values,
            index=pd.DatetimeIndex(np.array(
                ERA5.STN_era5_spring.valid_time[:])),
            columns=ERA5.STN_era5_spring.pressure_level)                                                          

        STN_ERA5_dict["wspeed"]  = pd.DataFrame(
            data=ERA5.STN_era5_spring["wspeed"].values,
            index=pd.DatetimeIndex(np.array(
                ERA5.STN_era5_spring.valid_time[:])),
            columns=ERA5.STN_era5_spring.pressure_level)                                                          

        STN_ERA5_dict["r"]   = pd.DataFrame(
            data=ERA5.STN_era5_spring["r"].values,
            index=pd.DatetimeIndex(np.array(
                ERA5.STN_era5_spring.valid_time[:])),
            columns=ERA5.STN_era5_spring.pressure_level)             
        
###############################################################################
if compare_era5_beluga:
    ERA5.open_ERA5_single_day(BELUGA_cls,flight)
    ERA5.era5_ds=ERA5.era5_day_ds
    ERA5.cut_era5_to_beluga_rf(start_hour,end_hour)
    
    STN_ERA5_rf_dict={}
    STN_ERA5_rf_dict["theta"] = pd.DataFrame(
        data=ERA5.STN_era5_flight["theta"].values,
        index=pd.DatetimeIndex(np.array(ERA5.STN_era5_flight.valid_time[:])),
        columns=ERA5.STN_era5_flight.pressure_level)                                                          

    STN_ERA5_rf_dict["wspeed"]  = pd.DataFrame(
        data=ERA5.STN_era5_flight["wspeed"].values,
        index=pd.DatetimeIndex(np.array(ERA5.STN_era5_flight.valid_time[:])),
        columns=ERA5.STN_era5_flight.pressure_level)                                                          

    STN_ERA5_rf_dict["r"]   = pd.DataFrame(
        data=ERA5.STN_era5_flight["r"].values,
        index=pd.DatetimeIndex(np.array(ERA5.STN_era5_flight.valid_time[:])),
        columns=ERA5.STN_era5_flight.pressure_level)                                                          
    import beluga_model_comparison_plotting as BELUGA_ERA5_comp
    BELUGA_ERA5_comp.verticalprofile_ERA5_for_specific_flight(
        STN_ERA5_rf_dict,df,rf=flight)

#%% CARRA plotting
if do_carra_plots:
    # Download CARRA data for entire campaign period as given in data paper
    moisture_var="rh"
    CARRA.open_STN_carra_period_heights()
    carra_ds=CARRA.carra_STN_ds
    carra_ds=carra_ds.sortby("valid_time")
    #CARRA.cut_carra_to_beluga_rf(start_hour,end_hour)
    if not flight=="RF1312" and not flight=="RF141312" and not flight=="RF1516":
        times=BELUGA_cls.flight_infos.loc[flight]
    
        carra_start_time                = times["Start Time"].round("60min")-\
                                        pd.Timedelta("60min")
        carra_end_time                  = times["End Time"].round("60min")+\
                                        pd.Timedelta("60min")
    elif flight=="RF1312":
        carra_start_time                = BELUGA_cls.flight_infos.loc["RF12"][\
                                            "Start Time"].round("60min")-\
                                            pd.Timedelta("60min")
        carra_end_time                  = BELUGA_cls.flight_infos.loc["RF13"][\
                                            "End Time"].round("60min")+\
                                            pd.Timedelta("60min")
    elif flight=="RF141312":
        carra_start_time                = BELUGA_cls.flight_infos.loc["RF12"][\
                                            "Start Time"].round("60min")-\
                                            pd.Timedelta("60min")
                                            
        carra_end_time                  = BELUGA_cls.flight_infos.loc["RF14"][\
                                            "End Time"].round("60min")+\
                                            pd.Timedelta("60min")
    elif flight=="RF1516":
        carra_start_time                = BELUGA_cls.flight_infos.loc["RF15"][\
                                            "Start Time"].round("60min")-\
                                            pd.Timedelta("60min")
                                            
        carra_end_time                  = BELUGA_cls.flight_infos.loc["RF16"][\
                                            "End Time"].round("60min")+\
                                            pd.Timedelta("60min")
                                            
    cutted_carra_rf_ds              = carra_ds.sel(valid_time=\
                                        slice(carra_start_time,carra_end_time))
    cutted_carra_rf_ds              = CARRA.calculate_carra_theta_e(
                                        var_to_use=cutted_carra_rf_ds)
    cutted_carra_rf_ds              = CARRA.calculate_carra_q_from_rh(
                                        var_to_use=cutted_carra_rf_ds)
    # Dictionary of relevant STN collocated data variables
    STN_CARRA_rf_dict               = {}
    STN_CARRA_rf_dict["theta"]      = pd.DataFrame(
                        data=np.array(cutted_carra_rf_ds["theta"].values[:]),
                        index=pd.DatetimeIndex(cutted_carra_rf_ds.valid_time[:]),
                        columns=np.array(cutted_carra_rf_ds.heightAboveGround[:]))
    
    STN_CARRA_rf_dict["theta_e"]    = pd.DataFrame(
                        data=np.array(cutted_carra_rf_ds["theta_e"].values[:]),
                        index=pd.DatetimeIndex(cutted_carra_rf_ds.valid_time[:]),
                        columns=np.array(cutted_carra_rf_ds.heightAboveGround[:]))
    
    STN_CARRA_rf_dict["rh"]         = pd.DataFrame(
                        data=np.array(cutted_carra_rf_ds["r"].values[:]),
                        index=pd.DatetimeIndex(cutted_carra_rf_ds.valid_time[:]),
                        columns=np.array(cutted_carra_rf_ds.heightAboveGround[:]))
    STN_CARRA_rf_dict["wspeed"]     = pd.DataFrame(
                        data=np.array(cutted_carra_rf_ds["ws"].values[:]),
                        index=pd.DatetimeIndex(cutted_carra_rf_ds.valid_time[:]),
                        columns=np.array(cutted_carra_rf_ds.heightAboveGround[:]))
   
    STN_CARRA_rf_dict["q"]          = pd.DataFrame(
                        data=np.array(
                            cutted_carra_rf_ds["specific_humidity"].values[:])*1000,
                        index=pd.DatetimeIndex(cutted_carra_rf_ds.valid_time[:]),
                        columns=np.array(cutted_carra_rf_ds.heightAboveGround[:]))
    STN_CARRA_rf_dict["pressure"]   = pd.DataFrame(
                        data=np.array(cutted_carra_rf_ds["pres"])) 
    
    if not do_campaign_overview:
        # Along single flight
        import beluga_model_comparison_plotting as BELUGA_CARRA_comp
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            {},df,rf=flight,moisture="q",icon_dict={},
            plot_flight=True)
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            STN_CARRA_rf_dict,df,rf=flight,moisture="q",icon_dict={},
            plot_flight=True)
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            STN_CARRA_rf_dict,df,rf=flight,moisture=moisture_var,icon_dict={},
            plot_flight=True)
    
    if compare_icon_beluga:
        if not flight=="RF1312"and not flight=="RF141312":
            times=BELUGA_cls.flight_infos.loc[flight]
        
            icon_start_time     = times["Start Time"].round("60min")-\
                                            pd.Timedelta("60min")
            icon_end_time       = times["End Time"].round("60min")+\
                                            pd.Timedelta("60min")
        elif  flight=="RF1312":
            icon_start_time     = BELUGA_cls.flight_infos.loc["RF12"][\
                                                "Start Time"].round("60min")-\
                                                pd.Timedelta("60min")
            icon_end_time       = BELUGA_cls.flight_infos.loc["RF13"][\
                                                "End Time"].round("60min")+\
                                                pd.Timedelta("60min")
        elif flight=="RF141312":
            icon_start_time     = BELUGA_cls.flight_infos.loc["RF12"][\
                                                "Start Time"].round("60min")-\
                                                pd.Timedelta("60min")
                                                
            icon_end_time       = BELUGA_cls.flight_infos.loc["RF14"][\
                                                "End Time"].round("60min")+\
                                                pd.Timedelta("60min")
            
        
        ICON.open_icon_STN_file(BELUGA_cls,flight)
        date=str(icon_start_time.date())
        ICON.adapt_icon_time_index(date)
        icon_ds=ICON.STN_ds
        
        cutted_icon_rf_ds              = icon_ds.sel(time=\
                                            slice(icon_start_time,
                                                  icon_end_time))
        cutted_icon_rf_ds              = ICON.calculate_icon_rh_from_q(
                                            var_to_use=cutted_icon_rf_ds)
        cutted_icon_rf_ds              = ICON.calculate_icon_theta_e(
                                            var_to_use=cutted_icon_rf_ds)
        cutted_icon_rf_ds              = ICON.calculate_icon_wspeed_and_wdir(
                                            var_to_use=cutted_icon_rf_ds)
        icon_ds=cutted_icon_rf_ds
        height_offset=31
        STN_ICON_rf_dict               = {}
        STN_ICON_rf_dict["theta"]      = pd.DataFrame(
                            data=np.array(icon_ds["theta"].values[:]),
                            index=pd.DatetimeIndex(icon_ds.time[:]),
                            columns=np.array(icon_ds.z_mc[:]-height_offset))
        
        STN_ICON_rf_dict["theta_e"]    = pd.DataFrame(
                            data=np.array(icon_ds["theta_e"].values[:]),
                            index=pd.DatetimeIndex(icon_ds.time[:]),
                            columns=np.array(icon_ds.z_mc[:]-height_offset))
        
        STN_ICON_rf_dict["rh"]          = pd.DataFrame(
                            data=np.array(icon_ds["rh"].values[:]),
                            index=pd.DatetimeIndex(icon_ds.time[:]),
                            columns=np.array(icon_ds.z_mc[:]-height_offset))
        
        STN_ICON_rf_dict["q"]          = pd.DataFrame(
                            data=np.array(icon_ds["qv"].values[:])*1000,
                            index=pd.DatetimeIndex(icon_ds.time[:]),
                            columns=np.array(icon_ds.z_mc[:]-height_offset))
        
        STN_ICON_rf_dict["wspeed"]     = pd.DataFrame(
                            data=np.array(icon_ds["wspeed"].values[:]),
                            index=pd.DatetimeIndex(icon_ds.time[:]),
                            columns=np.array(icon_ds.z_mc[:]-height_offset))
        
        STN_ICON_rf_dict["pressure"]   = pd.DataFrame(
                            data=np.array(icon_ds["pres"]))
    
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            STN_CARRA_rf_dict,df,rf=flight,moisture=moisture_var,
            icon_dict={},plot_flight=True)
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            STN_CARRA_rf_dict,df,rf=flight,moisture=moisture_var,
            icon_dict=STN_ICON_rf_dict,plot_flight=True)
        
        # additional moisture var
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            STN_CARRA_rf_dict,df,rf=flight,moisture="q",
            icon_dict={},plot_flight=True)
        BELUGA_CARRA_comp.verticalprofile_CARRA_for_specific_flight(
            STN_CARRA_rf_dict,df,rf=flight,moisture="q",
            icon_dict=STN_ICON_rf_dict,plot_flight=True)
    else:
        STN_ICON_rf_dict={}
    if flight=="RF1312" or flight=="RF141312" or flight=="RF1516":
        BELUGA_CARRA_comp.contour_profile_CARRA_ICON_BELUGA(
                df,flight,STN_CARRA_rf_dict,STN_ICON_rf_dict,var_to_plot="theta")
        BELUGA_CARRA_comp.contour_profile_CARRA_ICON_BELUGA(
                df,flight,STN_CARRA_rf_dict,STN_ICON_rf_dict,var_to_plot="theta_e")
        BELUGA_CARRA_comp.contour_profile_CARRA_ICON_BELUGA(
                df,flight,STN_CARRA_rf_dict,STN_ICON_rf_dict,var_to_plot="RH")
        BELUGA_CARRA_comp.contour_profile_CARRA_ICON_BELUGA(
                df,flight,STN_CARRA_rf_dict,STN_ICON_rf_dict,var_to_plot="wspeed")
        BELUGA_CARRA_comp.contour_profile_CARRA_ICON_BELUGA(
                df,flight,STN_CARRA_rf_dict,STN_ICON_rf_dict,var_to_plot="q")
        
    #entire campaign times
    import carra_plotting as CARRA_plotting
    print("Calculate entire Theta values. Takes a while")
    carra_ds=CARRA.calculate_carra_theta_e(var_to_use=carra_ds)
    carra_ds=CARRA.calculate_carra_q_from_rh(var_to_use=carra_ds)
    STN_CARRA_period_dict               = {}
    STN_CARRA_period_dict["theta"]      = pd.DataFrame(
                        data=np.array(carra_ds["theta"].values[:]),
                        index=pd.DatetimeIndex(carra_ds.valid_time[:]),
                        columns=np.array(carra_ds.heightAboveGround[:]))
    
    STN_CARRA_period_dict["theta_e"]    = pd.DataFrame(
                        data=np.array(carra_ds["theta_e"].values[:]),
                        index=pd.DatetimeIndex(carra_ds.valid_time[:]),
                        columns=np.array(carra_ds.heightAboveGround[:]))
    
    STN_CARRA_period_dict["rh"]          = pd.DataFrame(
                        data=np.array(carra_ds["r"].values[:]),
                        index=pd.DatetimeIndex(carra_ds.valid_time[:]),
                        columns=np.array(carra_ds.heightAboveGround[:]))
    
    STN_CARRA_period_dict["q"]          = pd.DataFrame(
                        data=np.array(carra_ds["specific_humidity"].values[:]),
                        index=pd.DatetimeIndex(carra_ds.valid_time[:]),
                        columns=np.array(carra_ds.heightAboveGround[:]))
    
    STN_CARRA_period_dict["wspeed"]     = pd.DataFrame(
                        data=np.array(carra_ds["ws"].values[:]),
                        index=pd.DatetimeIndex(carra_ds.valid_time[:]),
                        columns=np.array(carra_ds.heightAboveGround[:]))
    
    STN_CARRA_period_dict["pressure"]   = pd.DataFrame(
                        data=np.array(carra_ds["pres"])) 

        
if do_campaign_overview:
    if do_beluga_plots:
        import beluga_plotting
        #beluga_plotting.basic_flight_information(cpgn_df,BELUGA_cls,
        #                            add_profile_infos_in_file=True)
        #beluga_plotting.plot_profile_infos_RFs(BELUGA_cls)
        #tnir_df=pd.DataFrame()
        #beluga_plotting.plot_tnir_histogram(tnir_df,BELUGA_cls)
        beluga_plotting.plot_met_var_statistics(cpgn_df, as_scatter=False)
    if do_era5_plots:
        era5_plotting.plot_campaign_overview(STN_ERA5_dict,BELUGA_cls)
    if do_carra_plots:
        CARRA_plotting.synoptic_carra_map(CARRA,None,main_path)
        #CARRA_plotting.plot_campaign_overview(STN_CARRA_period_dict,
        #        BELUGA_cls,moisture="RH")
        #CARRA_plotting.plot_campaign_overview(STN_CARRA_period_dict,
        #        BELUGA_cls,moisture="RH",include_theta_contours=True)
        #CARRA_plotting.plot_campaign_overview(STN_CARRA_period_dict,
        #        BELUGA_cls,moisture="rh",theta_var="theta_e",
        #        include_theta_contours=True)

