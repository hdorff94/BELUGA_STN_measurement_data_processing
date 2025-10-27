# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:00:46 2025

@author: u300737
"""
import glob
import os

import numpy as np
import pandas as pd
import xarray as xr

import metpy

from metpy.units import units
import metpy.calc as mpcalc

class Gridded_data():
    def __init__(self,main_path=os.getcwd()+"/../"):
        self.station_name       = "Station_Noord"
        self.coordinates        = STN_coords={"lat": 81.60,"lon": -16.65}
        self.main_path          = main_path
         
        print("Main path is:",main_path)
        
class AMSR2_seaice(Gridded_data):
    
    def __init__(self,Gridded_data):
        self.name           = "AMSR2 sea ice"
        self.station_name   = Gridded_data.station_name        
        self.coordinates    = Gridded_data.coordinates
        self.main_path      = Gridded_data.main_path
        self.sea_ice_path   = self.main_path+"/BELUGA_data/"
        self.nearest_lat    = 81.5
        self.nearest_lon    = -16.75
    
    def open_sea_ice(self,entire_campaign=True):
        import glob
        sea_ice_file_list=glob.glob(self.sea_ice_path+"*AMSR2*.nc")
        file_date_str_list=[file[-16:-8] for file in sea_ice_file_list]
        sea_ice_ds=xr.open_mfdataset(sea_ice_file_list,combine="nested",
                           concat_dim='time',
                           preprocess=self.add_time_dim)
        
        file_dates=pd.DatetimeIndex(file_date_str_list)
        sea_ice_ds["time"]=xr.DataArray(file_dates,dims="time")
        self.sea_ice=sea_ice_ds.mean(dim="time").compute()
    

    def add_time_dim(self,xda):
        from datetime import datetime
        xda = xda.expand_dims(time = [datetime.now()])
        return xda    

class ICON(Gridded_data):
    def __init__(self,Gridded_data):
        self.name               = "ERA5"
        self.station_name       = Gridded_data.station_name        
        self.coordinates        = Gridded_data.coordinates
        self.main_path          = Gridded_data.main_path+"/"
        self.data_path           = self.main_path+"/BELUGA_data/"
        self.nearest_lat        = 81.5
        self.nearest_lon        = -16.75
    
    def open_icon_STN_file(self,BELUGA_cls,flight):
        if not flight=="RF1312"and not flight=="RF141312":
            flight_date=str(
                BELUGA_cls.flight_infos["Start Time"].loc[flight].date())
        else:
            flight_date="2024-04-01"
        icon_file="ICON_STN_"+flight_date+".nc"
        self.STN_ds=xr.open_dataset(self.data_path+icon_file)
    
    def adapt_icon_time_index(self,date):
        """ 
        The ICON time coordinates are given in ordinal format, 
        this function changes the index into the typical format 
        %YYYY-%MM-%DD %HH:%MM (so rounded to minutes). This is required 
        for the precedent interpolation in time to minutely frequencies.
        
        date : str
            str with 8 letters specifying the date, format: YYYYMMDD
        
        Returns
        -------
        da : xr.DataArray
            ICON-Variable as DataArray (xarray) with updated Time Coordinates.
        """
        #Create index, fixed for 10-minute data
        if "-" in date:
            start_date=date[0:4]+date[5:7]+date[8:10]
        else:
            start_date=date
        new_time_index=pd.to_datetime(abs(int(start_date)-\
                        np.array(self.STN_ds.time)),unit="d",origin=start_date)
        new_time_index=new_time_index.round("min")
        
        #Assign this index to the DataArray
        self.STN_ds=self.STN_ds.assign_coords(time=new_time_index)
        print("Index (Time Coordinate) of DataArray changed")
        
    def calculate_icon_theta_e(self,var_to_use="self"):
        """
        Updated ICON dataset with new variables Theta and Theta_e.

        Raises
        ------
        Exception
            raises exception if variables needed to calculate Theta_e,
            such as T or RH, are missing.

        """
        if isinstance(var_to_use, str):
            if var_to_use=="self":
                icon_var=self.STN_ds
        else:
            icon_var=var_to_use
        
        if not [i in list(icon_var.keys()) for i in ["r","t"]]:
            raise Exception("Some variables are not included in the dataset.",
                            "Re-download the correct Reanalysis data")
        print("Start calculating the different temperatures")
        T_profile   = icon_var["temp"]
        RH_profile  = icon_var["rh"]
        p_hPa       = icon_var["pres"]
        if float(p_hPa.max())>10000:
            p_hPa       /= 100
        
        p_hPa=p_hPa * units.hPa
        
        print("Calculate Potential Temperature")
        Theta_profile=mpcalc.potential_temperature(p_hPa, T_profile)
        icon_var["theta"]=xr.DataArray(np.array(Theta_profile),
            coords=icon_var["temp"].coords,attrs={'units': "K",
                'long_name':"Potential Temperature"})
        
        print("Calculate Dewpoint")
        Td_profile=mpcalc.dewpoint_from_relative_humidity(T_profile,RH_profile)
        
        print("Calculate Equivalent Potential Temperature -- takes a while")
        Te_profile=mpcalc.equivalent_potential_temperature(p_hPa,
                                            T_profile,Td_profile)
        icon_var["theta_e"]=xr.DataArray(
            np.array(Te_profile),coords=icon_var["temp"].coords,
            attrs={'units': "K",
                  'long_name':"Equivalent Potential Temperature"})
        
        print("Theta_e calculated")
        return icon_var   
        
    def calculate_icon_rh_from_q(self,var_to_use="self"):
        """
        
        
        Returns
        -------
        None.
            
        """
        print("Calculate RH from q")
        #Using metpy functions to calculate specific humidity from RH
        if isinstance(var_to_use, str):
            if var_to_use=="self":
                icon_var=self.STN_ds
        else:
            icon_var=var_to_use
        
        p_hPa       = icon_var["pres"]/100
        p_hPa=p_hPa * units.hPa
        
    
        q=icon_var["qv"].data
        temperature=icon_var["temp"].data * units.K
        
        rh_data=mpcalc.relative_humidity_from_specific_humidity(
            p_hPa, temperature, q)
        relative_humidity=xr.DataArray(
            np.array(rh_data),
                dims=["time","height"])
        icon_var=icon_var.assign({"rh":relative_humidity})
        return icon_var
    def calculate_icon_q_from_rh(self,var_to_use="self"):
        """
        
        
        Returns
        -------
        None.
            
        """
        print("Calculate q from RH")
        #Using metpy functions to calculate specific humidity from RH
        if isinstance(var_to_use, str):
            if var_to_use=="self":
                icon_var=self.STN_ds
        else:
            icon_var=var_to_use
        
        p_hPa       = icon_var["pres"]/100
        p_hPa=p_hPa * units.hPa
        
    
        rh=icon_var["r"].data/100
        temperature=icon_var["t"].data * units.K
        
        mixing_ratio=mpcalc.mixing_ratio_from_relative_humidity(
            p_hPa,temperature,rh)
        q=mpcalc.specific_humidity_from_mixing_ratio(mixing_ratio)
        specific_humidity=xr.DataArray(
            np.array(q),
                dims=["valid_time","heightAboveGround"])
    
        icon_var=icon_var.assign({"specific_humidity":specific_humidity})
        return icon_var
    
    def calculate_icon_wspeed_and_wdir(self,var_to_use="self"):
        """
        Calculates wspeed and wdir from u and v component.

        """
        #Wind speed
        if isinstance(var_to_use,str):
            if var_to_use=="self":
                icon_var=self.STN_ds
        else:
            icon_var=var_to_use
            icon_var["wspeed"]=np.sqrt(icon_var["u"]**2+\
                                       icon_var["v"]**2)
        print("Wspeed calculated")
        # Wind direction
        print("wind direction not yet calculated. Tbd.")
        return icon_var

class ERA5(Gridded_data):
    def __init__(self,Gridded_data):
        self.name               = "ERA5"
        self.station_name       = Gridded_data.station_name        
        self.coordinates        = Gridded_data.coordinates
        self.main_path          = Gridded_data.main_path
        self.era_path           = self.main_path+"/BELUGA_data/"
        self.nearest_lat        = 81.5
        self.nearest_lon        = -16.75

    def get_info(self):
        print("This is the class for ",self.station_name,
              "gridded data from simulations (models, reanalyses).",
              " It uses the representation of ",self.name)
    def open_era5_plevel_period_data(self,take_as_major_era5_file=True,
                                     calc_for_all=False):
        era5_march=xr.open_dataset(self.era_path+\
                "ERA5_Vertical_Pressure_Levels_March_2024.nc")
        era5_april=xr.open_dataset(self.era_path+\
                "ERA5_Vertical_Pressure_Levels_April_2024.nc")
        print("Merge ERA5 files")
        era5_spring=xr.concat([era5_march,era5_april],dim="valid_time")
        del era5_march, era5_april
        if calc_for_all:    
            era5_spring=self.calculate_era5_theta_e(era5_spring)
            era5_spring=self.calculate_era5_wspeed_and_wdir(era5_spring)
        else:
            self.STN_era5_spring=era5_spring.sel({"longitude":self.nearest_lon,
                                                  "latitude":self.nearest_lat})
            self.STN_era5_spring=\
                self.calculate_era5_theta_e(self.STN_era5_spring)
            self.STN_era5_spring=\
                self.calculate_era5_wspeed_and_wdir(self.STN_era5_spring)
        
        if take_as_major_era5_file:
            self.era5_ds=self.STN_era5_spring.copy()
        
    def open_ERA5_single_day(self,BELUGA_cls,flight):
        
        flight_date=str(BELUGA_cls.flight_dates[flight])
        date_dashed=flight_date[0:4]+"-"+flight_date[4:6]+\
            "-"+flight_date[6:8]
        era_file=self.era_path+"ERA5*"+date_dashed+"*"

        file_list=glob.glob(era_file)
        if not len(file_list)==1:
            raise FileExistsError("Either no or several files label with",
                                      flight, "exist!")
        else:
            print(file_list[0])
            era_ds=xr.open_dataset(file_list[0])
            # Cute ERA5 to desired date
            era_rf_ds=era_ds.sel({"valid_time":date_dashed})
        era_ds=self.calculate_era5_theta_e(era_ds)
        era_ds=self.calculate_era5_wspeed_and_wdir(era_ds)
        self.era5_day_ds=era_ds
    
    def get_STN_nearest_ERA5_cells(self,ds):
        era5_STN=ds.sel({"longitude": self.nearest_lon,
                             "latitude" : self.nearest_lat})
        return era5_STN
    
    def cut_era5_to_beluga_rf(self,start_hour,end_hour):
        cutted_era_ds   = self.era5_ds.isel(valid_time=\
                                            slice(start_hour, end_hour))
        STN_era5_flight = self.get_STN_nearest_ERA5_cells(cutted_era_ds)
        STN_era5_flight = self.calculate_era5_theta_e(STN_era5_flight)
        STN_era5_flight = self.calculate_era5_wspeed_and_wdir(STN_era5_flight)
        self.STN_era5_flight=STN_era5_flight
        
    def calculate_era5_theta_e(self,var_to_use="self"):
        """
        Updated ERA5 dataset with new variables Theta and Theta_e.

        Raises
        ------
        Exception
            raises exception if variables needed to calculate Theta_e,
            such as T or RH, are missing.

        """
        if isinstance(var_to_use, str):
            if var_to_use=="self":
                era5_var=self.era5_ds
        else:
            era5_var=var_to_use
        
        if not [i in list(era5_var.keys()) for i in ["r","t"]]:
            raise Exception("Some variables are not included in the dataset.",
                            "Re-download the correct Reanalysis data")
        print("Start calculating the different temperatures")
        T_profile=era5_var["t"]#[:,:,50,50]
        RH_profile=era5_var["r"]#[:,:,50,50]
        
        #Replicate p-levels to ndarray
        if not hasattr(era5_var, "name"):    
            p_hPa=np.tile(era5_var["pressure_level"],
                          (era5_var["t"].shape[0],1))
        else:
            p_hPa=np.tile(float(str(era5_var["name"].values)[:-3]),
                          (era5_var["t"].shape[0],
                           era5_var["t"].shape[1]))
        
        #if not len(p_hPa.shape)==3:
        #    p_hPa=np.expand_dims(p_hPa,axis=2)
        
        #p_hPa=np.repeat(p_hPa,ds["t"].shape[1],axis=1)
        #try:
        #    p_hPa=np.repeat(p_hPa, ds["t"].shape[2],axis=2)
        #except:
        #    pass
        if len(era5_var["t"].shape)==4:
            p_hPa=np.expand_dims(p_hPa,axis=2)
            p_hPa=np.repeat(p_hPa, era5_var["t"].shape[2],axis=2)
            p_hPa=np.expand_dims(p_hPa,axis=3)
            p_hPa=np.repeat(p_hPa, era5_var["t"].shape[3],axis=3)
        p_hPa=p_hPa * units.hPa
        
        print("Calculate Potential Temperature")
        Theta_profile=mpcalc.potential_temperature(p_hPa, T_profile)
        era5_var["theta"]=xr.DataArray(np.array(Theta_profile),
            coords=era5_var["t"].coords,attrs={'units': "K",
                'long_name':"Potential Temperature"})
        
        print("Calculate Dewpoint")
        Td_profile=mpcalc.dewpoint_from_relative_humidity(T_profile,
                                                          RH_profile)
        
        print("Calculate Equivalent Potential Temperature -- takes a while")
        Te_profile=mpcalc.equivalent_potential_temperature(p_hPa,
                                                           T_profile,
                                                           Td_profile)
        era5_var["theta_e"]=xr.DataArray(
            np.array(Te_profile),coords=era5_var["t"].coords,
            attrs={'units': "K",
                  'long_name':"Equivalent Potential Temperature"})
        
        print("Theta_e calculated")
        return era5_var   
    
    def calculate_era5_wspeed_and_wdir(self,var_to_use="self"):
        """
        Calculates wspeed and wdir from u and v component.

        """
        #Wind speed
        if isinstance(var_to_use,str):
            if var_to_use=="self":
                era5_var=self.era5_ds
        else:
            era5_var=var_to_use
            era5_var["wspeed"]=np.sqrt(era5_var["u"]**2+\
                                       era5_var["v"]**2)
        print("Wspeed calculated")
        # Wind direction
        print("wind direction not yet calculated. Tbd.")
        return era5_var
        

class CARRA(Gridded_data,):
    def __init__(self,Gridded_data):        
        self.name               = "CARRA"
        self.station_name       = Gridded_data.station_name        
        self.coordinates        = Gridded_data.coordinates
        self.main_path          = Gridded_data.main_path
        self.carra_path=self.main_path+"/BELUGA_data/"
        
    def get_info(self):
        print("This is the class for ",self.station_name,
              "gridded data from simulations (models, reanalyses).",
              " It uses the representation of ",self.name)
    def get_carra_grid(self):
        carra_grid              = xr.open_dataset(self.carra_path+\
                                                  "CARRA_grid.nc")
        carra_lat               = carra_grid.latitude.values.flatten()
        carra_lon               = carra_grid.longitude.values.flatten()
        carra_grid_df           = pd.DataFrame(data=np.nan, 
                            index=range(len(carra_lat)),columns=["lat","lon"])
        carra_grid_df["lat"]    = carra_lat
        carra_grid_df["lon"]    = carra_lon
        carra_grid_df["lon"][carra_grid_df["lon"]>180]-=360
        return carra_grid_df
        
    def get_carra_grid_convex_hull(self,grid_ds):
        self.CARRA_grid = grid_ds
        
    def open_STN_carra_period_heights(self):
        self.carra_STN_fname="STN_CARRA_height_levels_20240324_20240410.nc"
        self.carra_STN_ds=xr.open_dataset(self.carra_path+self.carra_STN_fname,
                                          engine="netcdf4")
    def open_sea_ice_2024(self,calc_mean_of_campaign=True):
        self.carra_seaice_fname=\
            "CARRA_west_sea_ice_STN_campaign_March_April_2024.nc"
        self.seaice_ds=xr.open_dataset(self.carra_path+self.carra_seaice_fname)
        if calc_mean_of_campaign:
            STN_period_seaice=self.seaice_ds.sel(
                valid_time=slice("2024-03-24","2024-04-11"))
            self.mean_seaice=STN_period_seaice.mean(dim="valid_time")
            return self.mean_seaice
        else:
            return None
        
    def calculate_carra_theta_e(self,var_to_use="self"):
        """
        Updated CARRA dataset with new variables Theta and Theta_e.

        Raises
        ------
        Exception
            raises exception if variables needed to calculate Theta_e,
            such as T or RH, are missing.

        """
        if isinstance(var_to_use, str):
            if var_to_use=="self":
                carra_var=self.carra_STN_ds
        else:
            carra_var=var_to_use
        
        if not [i in list(carra_var.keys()) for i in ["r","t"]]:
            raise Exception("Some variables are not included in the dataset.",
                            "Re-download the correct Reanalysis data")
        print("Start calculating the different temperatures")
        T_profile   = carra_var["t"]
        RH_profile  = carra_var["r"]
        p_hPa       = carra_var["pres"]
        if float(p_hPa.max())>10000:
            p_hPa       /= 100
        
        p_hPa=p_hPa * units.hPa
        
        print("Calculate Potential Temperature")
        Theta_profile=mpcalc.potential_temperature(p_hPa, T_profile)
        carra_var["theta"]=xr.DataArray(np.array(Theta_profile),
            coords=carra_var["t"].coords,attrs={'units': "K",
                'long_name':"Potential Temperature"})
        
        print("Calculate Dewpoint")
        Td_profile=mpcalc.dewpoint_from_relative_humidity(T_profile,
                                                          RH_profile)
        
        print("Calculate Equivalent Potential Temperature -- takes a while")
        Te_profile=mpcalc.equivalent_potential_temperature(p_hPa,
                                                           T_profile,
                                                           Td_profile)
        carra_var["theta_e"]=xr.DataArray(
            np.array(Te_profile),coords=carra_var["t"].coords,
            attrs={'units': "K",
                  'long_name':"Equivalent Potential Temperature"})
        
        print("Theta_e calculated")
        return carra_var   
        
    def calculate_carra_q_from_rh(self,var_to_use="self"):
        """
        
        
        Returns
        -------
        None.
            
        """
        print("Calculate q from RH")
        #Using metpy functions to calculate specific humidity from RH
        if isinstance(var_to_use, str):
            if var_to_use=="self":
                carra_var=self.carra_STN_ds
        else:
            carra_var=var_to_use
        
        p_hPa       = carra_var["pres"]
        p_hPa=p_hPa * units.hPa
        
    
        rh=carra_var["r"].data/100
        temperature=carra_var["t"].data * units.K
        
        mixing_ratio=mpcalc.mixing_ratio_from_relative_humidity(
            p_hPa,temperature,rh)
        q=mpcalc.specific_humidity_from_mixing_ratio(mixing_ratio)
        specific_humidity=xr.DataArray(
            np.array(q),
                dims=["valid_time","heightAboveGround"])
    
        carra_var=carra_var.assign({"specific_humidity":specific_humidity})
        return carra_var
    
    def calculate_carra_wspeed_and_wdir(self,var_to_use="self"):
        """
        Calculates wspeed and wdir from u and v component.

        """
        #Wind speed
        if isinstance(var_to_use,str):
            if var_to_use=="self":
                carra_var=self.carra_ds
        else:
            carra_var=var_to_use
            carra_var["wspeed"]=np.sqrt(carra_var["u"]**2+\
                                       carra_var["v"]**2)
        print("Wspeed calculated")
        # Wind direction
        print("wind direction not yet calculated. Tbd.")
        return carra_var
