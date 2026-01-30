# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 12:40:08 2025

@author: u300737
"""

"""
This contains the classes of all relevant measurement platforms located at 
Station Noord, Greenland. In particular, the BELUGA platform is here defined as
class.
"""
import os 
import glob


from functools import reduce

import numpy as np
import pandas as pd
import xarray as xr

import math
from geopy.distance import geodesic


import metpy
import metpy.calc as mpcalc
from metpy.units import units

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates  as mdates
import seaborn as sns

import Performance
import logging_processing


 
class Measurement_Platforms_VRS():
    
    def __init__(self,main_path=os.getcwd()+"/../BELUGA_data/"):
        self.station_name       = "Villum_Research_Station"
        self.coordinates        = VRS_coords={"lat": 81.60250,"lon": -16.67,
                                              "height":31}
        self.main_path          = main_path
        
        print("Main path is:",main_path)
    def get_info(self):
        print("This is the class for ",self.station_name,
              "with coordinates",self.VRS_coords,
              ". Data can be found under:", self.main_path)

class BELUGA(Measurement_Platforms_VRS):
    def __init__(self, Measurement_Platforms_VRS):
        self.station_name       = Measurement_Platforms_VRS.station_name        
        self.coordinates        = Measurement_Platforms_VRS.coordinates
        self.main_path          = Measurement_Platforms_VRS.main_path
        self.flight_dates       = {
            "RF01":"20240324","RF02":"20240324","RF03":"20240325",
            "RF04":"20240325","RF05":"20240326","RF06":"20240326",
            "RF07":"20240327","RF08":"20240328","RF09":"20240328",
            "RF10":"20240329","RF11":"20240330","RF12":"20240401",
            "RF13":"20240401","RF14":"20240401","RF15":"20240402",
            "RF16":"20240402","RF17":"20240403","RF18":"20240404",
            "RF19":"20240404","RF20":"20240405","RF21":"20240406",
            "RF22":"20240407","RF23":"20240408","RF24":"20240408",
            "RF25":"20240409","RF26":"20240410","RF27":"20240411",
            "RF28":"20240412",
            #artificial merged dates
            "RF1312":"20240401","RF141312":"20240401"}
        self.temporary=True
    #-------------------------------------------------------------------------#
    # Barometric pressure calculation routines
    def get_measured_p_sfc(self,sensor="bme",data_df=pd.DataFrame()):
        if data_df.shape[0]==0:
            df=self.l1_df
        
        # find surface pressure by periods where mean gradients
        # are very low
        if sensor=="bme":
            sensor_pressure=sensor+"_pressure"
            hpa_factor=1
        elif sensor=="TMP_turb":
            sensor_pressure="p"
            hpa_factor=100
        else:    
            sensor_pressure="P_R"
            hpa_factor=100
        p=pd.DataFrame(data=df[sensor_pressure].values*hpa_factor,
                       columns=["pressure"],index=df.index)
        p["pressure"].loc[p["pressure"]==0]=np.nan
        p["int_idx"]=np.arange(p.shape[0])
        if sensor=="bme":
            t=df[sensor+"_temperature"]
        elif sensor=="TMP_turb":
            t=df["T"]
        
        else:
            t=df["T_R_fil"]
        # Minutely pressure gradients
        p_grads=p.diff().resample("30s").mean()
        low_p_grads=p_grads[abs(p_grads["pressure"]/hpa_factor)<0.05]
        # use this values for surface pressure
        intersect_index=low_p_grads.index.intersection(p.index)
        p_sfc_values_all=p.loc[intersect_index]
        #except:
        #    p_sfc_values_all=p.loc[low_p_grads.index[1:]]
            
        # 2nd condition
        quantile_value=0.95
        if self.flight=="RF14" or self.flight=="RF02" or self.flight=="RF09":
            quantile_value=0.9
        if self.flight=="RF07" and sensor=="bme":
            quantile_value=0.8
            
        p_95=p["pressure"].quantile(quantile_value) # changed from 0.9
        p_sfc_values=p_sfc_values_all.loc[p_sfc_values_all["pressure"]>p_95]
        t_sfc=t.loc[p_sfc_values.index]
        if t_sfc.isnull().all():
            p_q=p["pressure"].quantile(.85) # changed from 0.9
            p_q_sfc_values=p_sfc_values_all.loc[p_sfc_values_all["pressure"]>p_q]
            t_sfc=t.loc[p_q_sfc_values.index]
            t_sfc[t_sfc>-5]=np.nan
        return p,p_sfc_values,t_sfc
    
    def get_barometric_altitude(self,sensor="bme"):
        # Barometric altitude we need the measured surface pressure
        self.logger.info(f"Derive barometric altitude from {sensor}-sensor")
        p,p_sfc_meas,t_sfc=self.get_measured_p_sfc(
            sensor=sensor,data_df=pd.DataFrame())
        if not self.flight=="RF09":
            _,p_tendency=self.calc_sfc_pressure_tendency(p_sfc_meas)
            p_tendency_slope=float(p_tendency.slope)
            #if self.flight=="RF14" and sensor=="TMP_turb":
                #    p_tendency_slope=-1.8/36/50
                # apply regression coefficients to timesteps for continous time series
                # unused psfc
            p_sfc_ts=p_tendency_slope*p["int_idx"]+p_tendency.intercept
            self.logger.info(f"BELUGA-based surface pressure tendency of {p_tendency.slope*36}"+\
                         "hPa/h considered")
            p_sfc_ts_est=p_sfc_ts.copy()
                
        #p_sfc_ts=
        # Load ground-based BP BME-280 pressure
        if not self.flight=="RF27":
            gnd_bp_path=self.main_path+"/ground-based/"+self.flight+"/"
            bme_files=glob.glob(gnd_bp_path+"BME*")
            bme_columns=['yyyy', 'month', 'dd', 'hh',
            'mm', 'sssss', 'bme_temperature', 'bme_rh', 'bme_pressure']
        
            bme_gnd=pd.read_csv(bme_files[0], delim_whitespace=True,
                names=bme_columns, dtype=str, encoding='latin1', skiprows=1)
            bme_gnd['time'] = pd.to_datetime(
                bme_gnd[['yyyy', 'month', 'dd', 'hh', 'mm']].astype(str).agg(
                ' '.join, axis=1), format='%Y %m %d %H %M', errors='coerce') +\
                pd.to_timedelta(bme_gnd['sssss'].astype(float).fillna(0) , unit='s')
            
            bme_gnd[bme_columns] = bme_gnd[bme_columns].apply(
            pd.to_numeric, errors='coerce')
            
            bme_gnd.set_index('time', inplace=True)
        
            bme_gnd.drop(columns=['yyyy', 'month', 'dd', 'hh', 'mm', 'sssss'],
                         inplace=True)
            if not sensor=="TMP_turb":
                bme_gnd=bme_gnd.resample("1s").mean()
            else:
                bme_gnd=bme_gnd.loc[bme_gnd.index.dropna()]
                bme_gnd=bme_gnd.sort_index()
                bme_gnd=bme_gnd.resample("20ms").nearest()
                
            union_index=self.l1_df.index.union(bme_gnd.index)
            p_sfc_ts=pd.Series(data=np.nan,index=union_index)
            p_sfc_ts.loc[bme_gnd.index]=bme_gnd["bme_pressure"].values
            p_sfc_ts=p_sfc_ts.ffill().bfill()
            p_sfc_ts=p_sfc_ts.loc[self.l1_df.index]
            t_sfc   =pd.Series(data=np.nan, index=union_index)
            t_sfc.loc[bme_gnd.index]=bme_gnd["bme_temperature"].values
            t_sfc=t_sfc.ffill().bfill()
            t_sfc=t_sfc.loc[self.l1_df.index]
        if self.plot_processing:
            if sensor=="bme":
                temp_var="bme_temperature"
            elif sensor=="TMP_turb":
                temp_var="T"
            else:
                temp_var="T_R_fil"
            if not self.flight=="RF09":
                self.plot_p_sfc(p,p_sfc_meas, self.l1_df[temp_var], 
                            t_sfc,sensor=sensor)
        # Apply p_sfc series for barometric height calculation
        self.calc_barometric_altitude(p,sensor=sensor,T0=t_sfc.mean(),
                                      p0=p_sfc_ts)
        self.l1_df["z_b"][self.l1_df["z_b"]>1000]=np.nan
        if self.l1_df["z_b"].min()<0:
            offset=self.l1_df["z_b"].min()
            self.l1_df["z_b"]-=offset
        if self.plot_processing:
            self.plot_quicklook_zb(sensor=sensor)

                        
    def calc_barometric_altitude(self,p,
                sensor="bme",T0=288.15,p0=101325): 
        # find surface pressure by periods where mean gradients
        # are very low
        if T0<200:
            T0+=273.15
        R = 287.05 #fixed
        g = 9.81
        self.l1_df["z_b"]=(R * T0) / g * np.log(p0 / p["pressure"])
        if sensor=="bme":
            self.logger.info("Barometric height added to BP data frame")
        else:
            self.logger.info("Barometric height added to TMP data frame")
        if self.plot_processing:
            self.plot_quicklook_zb(sensor=sensor)
            return 

    def calc_sfc_pressure_tendency(self,p_sfc_values):
        from scipy.stats import linregress
        x = p_sfc_values["int_idx"]#.m .map(pd.Timestamp.toordinal)
        y = p_sfc_values["pressure"]
        
        p_tendency = linregress(x, y)
        #print(f'Slope: {p_tendency.slope}')
        #print(f'Intercept: {p_tendency.intercept}')
        # Fitted values
        p_fit = p_tendency.slope * x + p_tendency.intercept
        """
        weigthing approach
        # Convert datetime index to ordinal
        x = ts.index.map(pd.Timestamp.toordinal)
        y = ts.values

        # Create weights emphasizing more recent data
        # For example, linearly increasing weights towards the end
        weights = np.linspace(1, 2, len(ts))  # weights from 1 to 2

        # Perform weighted linear fit
        coefficients = np.polyfit(x, y, 1, w=weights)
        slope, intercept = coefficients

        # Calculate fitted values
        y_fit = np.polyval(coefficients, x)
        """
        return p_fit,p_tendency
    
    def add_VRS_georeferences(self,ds):
        
        ds["lat"]                           = self.coordinates["lat"]
        ds["lon"]                           = self.coordinates["lon"]
        ds["VRS_height"]                        = self.coordinates["height"]
        
        # add attributes based on CF-convention
        ds["lat"].attrs["standard_name"]    = "latitude"
        ds["lat"].attrs["long_name"]        = "latitude Villum Research Station"
        ds["lat"].attrs["units"]            = "degrees_north"
        
        ds["lon"].attrs["standard_name"]    = "longitude"
        ds["lon"].attrs["long_name"]        = "longitude Villum Research Station"
        ds["lon"].attrs["units"]            = "degrees_east"
        
        #ds["STN_height"].attrs["standard_name"] = "height"
        ds["VRS_height"].attrs["long_name"]     = "Station Nord height above mean sea level"
        ds["VRS_height"].attrs["units"]         = "m"
        #self.logger.info("Georeference data of Station Nord added to Dataset")
        return ds 
    
        
    def add_ABL_transition_type(self,ds):
        transition=self.flight_infos["type"].loc[ds.attrs["flight"]]
        ds.attrs["ABL_transition_type"]=transition
        return ds

    def add_rf_as_id_for_cf_role(self,ds):
        ds["flight_id"]=xr.DataArray(
            data=int(self.flight[2:4]))
        ds["flight_id"].attrs["long_name"]="research flight ID"
        ds["flight_id"].attrs["cf_role"]="timeseries_id"
        ds["flight_id"].attrs["units"]="1"
        return ds
    
    def add_L2_meta_data(self,ds, probe="BP"):
        
        import nc_attributes
        # Read dictionaries
        if probe=="BP":
            global_att_dict = nc_attributes.BP_global_attributes
            probe_att_dict     = nc_attributes.BP_var_attributes 
        elif probe=="TMP_met":
            global_att_dict = nc_attributes.TMP_met_global_attributes
            probe_att_dict     = nc_attributes.TMP_met_var_attributes
        elif probe=="TMP_turb":
            global_att_dict = nc_attributes.TMP_turb_global_attributes
            probe_att_dict  = nc_attributes.TMP_turb_var_attributes
        else:
            raise Exception("Caution! Wrong instrument probe name defined.")
        # Add Global Attributes from dictionary
        for att in [*global_att_dict.keys()]:
            ds.attrs[att]=global_att_dict[att]
        
        # Add variable attributes from dictionary
        #for coord in [*self.l2_ds.coords]:
        #    coord_attrs=probe_att_dict[coord]
        #    for att in [*coord_attrs.keys()]:
        #        self.l2_ds[coord].attrs[att]=coord_attrs[att]
        
        for var in [*ds.keys()]:
            var_attrs=probe_att_dict[var]
            for att in [*var_attrs.keys()]:
                ds[var].attrs[att]=var_attrs[att]
        ds=self.update_rf_specfic_attributes(ds)
        if self.flight=="RF08":
            ds["quality_flag"].attrs['flag_values'].append(5)
            ds["quality_flag"].attrs["flag_meaning"]+="bad_sonde_T_replaced_with_sonic_T "
        if probe=="TMP_met":
            if self.flight=="RF11" or self.flight=="RF25":
                ds["vv"].attrs["caution"]="for this RF, no OLA data was available. Wind extracted from radiosonde"
                ds["dd"].attrs["caution"]="for this RF, no OLA data was available. Wind extracted from radiosonde"
        ds.time.attrs['long_name']     = 'measurement time (UTC)'
        ds.time.attrs['standard_name'] = 'time'
        
        return ds
    
    def replace_inf(self,df):
        df.replace([np.inf,-np.inf],np.nan,inplace=True)
        return df
    
    def update_rf_specfic_attributes(self,ds):
        # Get current time 
        from datetime import datetime
        now = datetime.now()
        date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
        ds.attrs["flight_date"]         = self.flight_date
        ds.attrs["flight"]              = self.flight
        ds.attrs["processing_date"]     = date_time
        return ds
    
    def save_L2_data_as_nc_file(self,ds,version_number,fill_value=-9999.,
                                probe="BP"):
        if self.temporary:
            outp_path=self.raw_data_path+"/../temporary_data/" 
        
        else:
            outp_path=self.raw_data_path+"/../final_data/"
        
        os.makedirs(outp_path,exist_ok=True)
        l2_fname="BELUGA_VRS_L2_"+probe+"_"+self.flight+"_"+self.flight_date+\
            version_number+".nc"
        # Create fill values
        ds=ds.fillna(value=fill_value)
        nc_compression=dict(zlib=True,_FillValue=fill_value,
                            complevel=1,dtype=np.float64)
        #nc_encoding= {ds_var:nc_compression for ds_var in ds.variables}
        nc_encoding = {}
        for var_name, var in ds.variables.items():
            if np.issubdtype(var.dtype, np.number):
                nc_encoding[var_name] = nc_compression
            else:
                # Omit encoding for non-numeric variables 
                nc_encoding[var_name] = {}
        ds.time.encoding['units']      = "Seconds since 1970-01-01 00:00:00"
        ds.time.encoding["calendar"]   = "standard"
        ds.to_netcdf(outp_path+l2_fname,mode="w",engine="netcdf4",
                     encoding=nc_encoding)
        print("L2 data saved as:",outp_path+l2_fname)
        return outp_path+l2_fname
        
    def get_flight_information(self):
        self.flight_infos=pd.read_csv(
            self.main_path+"VRS_BELUGA_flight_infos.csv",
            index_col="Flight No")
        self.flight_infos["Start Time"]     = pd.to_datetime(
                self.flight_infos["Start Time"], dayfirst=True)
        self.flight_infos["End Time"]       = pd.to_datetime(
                self.flight_infos["End Time"], dayfirst=True)
        
        self.flight_infos["type"]           = "uncategorized"
        self.flight_infos["type"].loc[["RF07","RF08","RF10","RF15",
            "RF16","RF19","RF23","RF27"]]   = "clear-cloudy"
        self.flight_infos["type"].loc[["RF05","RF06"]]    = "day-night"
        self.flight_infos["type"].loc[["RF12","RF13","RF14","RF26"]]="low-level jet"

        #----------------------------------------------------------------------#
        self.LLJ_flights    = self.flight_infos.loc[\
                                self.flight_infos["type"]=="low-level jet"]
        self.cloud_flights  = self.flight_infos.loc[\
                                self.flight_infos["type"]=="clear-cloudy"]
        self.night_flights  = self.flight_infos.loc[\
                                self.flight_infos["type"]=="day-night"]
        #----------------------------------------------------------------------#
    
    def open_met_data(self,rf="RF12"):
        self.flight=rf
        print(self.main_path)
        files=glob.glob(self.main_path+self.flight+"_*")
        if not len(files)==1:
            raise FileExistsError("either no or several files label with",
                                  rf, "exist!")
        else:
            balloon_df=pd.read_csv(files[0],index_col=0)
            self.opened_file=files[0]
            balloon_df.index=pd.DatetimeIndex(balloon_df.index)
            self.met_df=balloon_df
            return balloon_df
   
    def save_met_data(self,rf="RF12", update_level="L1"):
        df=self.met_df
        file_path=self.opened_file
        file_name=file_path.split("\\")[-1]
        if update_level!="":
            file_name=update_level+"_"+file_name
        #Get the directory part of file_path
        dir_name = os.path.dirname(file_path)
        file_path=dir_name+"/"+file_name
        df.to_csv(file_path)
        print("Dataset saved as:",file_path)
    
    def open_all_met_data(self,level="L1"):
        """
    
        Returns
        -------
        cpgn_met_df : pd.DataFrame
            Merged df containing met data from balloon for all RFs.

        """
        
        self.flights=self.flight_infos.index
        data_list=[]
        if level=="L1":
            for rf in self.flights:
                try:
                    df=self.open_met_data(rf=rf)
                    print(rf," met data read")
                    data_list.append(df)
                except: 
                    print("no met data for",rf,
                      " available, continue to next flight")
            # merge data
            cpgn_met_df=pd.concat(data_list)
            
        elif level=="L2":
            self.raw_data_path      = self.main_path +\
           "/BELUGA_TMP_met_probe/raw_data/"
            l2_path=self.raw_data_path+"/../temporary_data/"
            version_number="v2.2.nc"
            l2_files=glob.glob(l2_path+"*"+version_number)
            l2_ds=xr.open_mfdataset(l2_files, combine="nested",concat_dim="time")
            l2_df=l2_ds[["T","rh","z_b","vv"]].to_dataframe()
            l2_df.rename(columns={"T":"TEMP","rh":"RH","z_b":"ALT",
                                  "vv":"hor_vv"},inplace=True)
            cpgn_met_df=l2_df.copy()
        else:
            raise Exception("Wrong format given. Only L1 and L2 possible")
        return cpgn_met_df
    
    def open_rad_data(self,level="l2"):
      if level=="l2": 
          fpath=self.main_path+\
             "/BELUGA_broadband_probe/temporary_data/"
          fname="BELUGA_VRS_L2_BP_"+self.flight+"_"+\
              self.flight_dates[self.flight]+".nc"
          file=fpath+fname   
          l2_ds=xr.open_dataset(file)
          return l2_ds
        
    def create_artificial_tnirs(self,sample_size):
        # Set the random seed for reproducibility
        mean_30m = -20
        min_value_30m = -100
        std_dev_30m = 20  
        # Adjusted to create a left-skewed distribution
        # Generate left-skewed data by reflecting a normal distribution
        data_30m = np.random.normal(loc=mean_30m,
            scale=std_dev_30m, size=sample_size)
        data_30m = np.clip(data_30m, min_value_30m, None) #Capping at -100 W/m²
        
        mean_hmax = -30
        min_value_hmax = -80
        std_dev_hmax = 20 
        # Adjusted to create a slightly right-skewed distribution
        
        # Generate slightly right-skewed data
        data_hmax = np.random.normal(loc=mean_hmax,
            scale=std_dev_hmax, size=int(sample_size*0.45))
        data_hmax = np.clip(data_hmax,
            min_value_hmax, None)  # Ensure minimum value is -80 W/m²
        tnir_30m  = pd.Series(data_30m)
        tnir_hmax = pd.Series(data_hmax)
        
        tnir_df=pd.DataFrame(data=np.nan, index=range(sample_size),
                             columns=["30m", "h$_{\mathrm{max}}$"])
        
        tnir_df["30m"]  = tnir_30m
        tnir_df["h$_{\mathrm{max}}$"] = tnir_hmax
        self.tnir_df    = tnir_df 
        return tnir_df
    
    def get_tnir_at_hmax_30m(self,artificial=True):
        """
        CAUTION, ath the moment this function creates artificial data.
        

        Returns
        -------
        None.

        """
        if artificial:    
            timesteps=3600*4*25 #assuming 25 RFs, recording 4 hours with 1Hz res 
            tnir_df=self.create_artifical_tnirs()
        else:
            tnir_df=self.tnir_df
        return tnir_df
    
    def open_turb_data(self):
        """
        This function is currently unused, 
        waiting for the processing of the turbulence data

        Returns
        -------
        None.

        """
        pass
    
    @staticmethod
    def calc_wspeed(df):
        df["hor_vv"]=np.sqrt(df["U_filtered"]**2+df["V_filtered"]**2)
        return df
    
    @staticmethod
    def calc_theta(df, temp_to_use="TEMP"):
        """
        This function calculates the potential temperature for the BELUGA
        temperature measurements

        Parameters
        ----------
        df : TYPE
            DESCRIPTION.
        temp_to_use : TYPE, optional
            DESCRIPTION. The default is "TEMP".

        Returns
        -------
        None.

        """
        
        df["Theta"]=np.nan
        # Prepare for metpy
        T_profile=df[temp_to_use]+273.15
        T_profile=T_profile * units.kelvin
        p_hPa=df["PRES"]
        p_hPa=p_hPa * units.hPa
        
        print("Calculate Potential Temperature")
        Theta_profile=mpcalc.potential_temperature(
            p_hPa.values * units.hPa,T_profile.values * units.kelvin)
        
        df["Theta"]=np.array(Theta_profile)
        return df
    
    @staticmethod
    def calc_theta_e(df,temp_to_use="TEMP"):
        df["Theta_e"]=np.nan
        # Prepare for metpy
        T_profile_temp=df[temp_to_use]+273.15
        T_profile=T_profile_temp.values * units.kelvin
        RH_profile  = df["RH"].values/100 #*units.percentage
        p_hPa=df["PRES"].values
        p_hPa=p_hPa * units.hPa
        print("Calculate Dewpoint")
        Td_profile=mpcalc.dewpoint_from_relative_humidity(T_profile,RH_profile)
        print("Calculate Equivalent Potential Temperature -- takes a while")
        Te_profile=mpcalc.equivalent_potential_temperature(p_hPa,
                                            T_profile,Td_profile)
        df["Theta_e"]=np.array(Te_profile)
        return df
    
    @staticmethod
    def calc_q_from_rh(df):
        """
        
        
        Returns
        -------
        None.
            
        """
        print("Calculate q from RH")
        #Using metpy functions to calculate specific humidity from RH
        df["q"]     = np.nan 
        p_hPa=np.array(df["PRES"])
        p_hPa=p_hPa * units.hPa
        
    
        RH_profile=np.array(df["RH"])/100
        T_profile=np.array(df["TEMP"])+273.15
        T_profile=T_profile * units.kelvin
        
        mr_profile=mpcalc.mixing_ratio_from_relative_humidity(
            p_hPa,T_profile,RH_profile)
        q_profile=mpcalc.specific_humidity_from_mixing_ratio(mr_profile)
        df["q"]=np.array(q_profile.to('g/kg'))
        return df
    #-------------------------------------------------------------------------#
    # Flight segmentation
    @staticmethod
    def find_profile_peaks(rf_df,min_height=200,alt_var="ALT",use_turb=False):
        from scipy.signal import find_peaks
        min_height=min_height
        if not use_turb:
            df=rf_df.copy()
        else:
            df=rf_df.resample("1s").mean()
        alt_max_idx,_ = find_peaks(df[alt_var],height=min_height,
            distance=500,width=2)
        alt_max_vals=df[alt_var].iloc[alt_max_idx]
        return alt_max_vals 
    
    def find_max_rf_height(rf_df,alt_var="ALT"):
        return rf_df[alt_var].max()
    
    def segment_flight_sections(self, rf_df,alt_var="ALT",
                                rate_threshold=0.7,use_turb=False):    
        def compute_height_change(group):
            # Calculate height change over each group
            if group['is_ascent_or_descent'].any():
                height_diff = group[alt_var].iloc[-1] - group[alt_var].iloc[0]
                segment_length = len(group)  
                return pd.Series({
                'height_change': height_diff,
                'segment_length': segment_length
                })
            else:
                return pd.Series({'height_change': 0, 
                                  'segment_length': len(group)})
        
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
            
        print("Perform flight segmentation:",self.flight)
        performance=Performance.performance()
        
        # Classes of segmentation:
            # ascent, descent, peak, near-ground
        segmentation_classes=["ascent","descent","peak","max"
                              "near_ground", "nearly_constant_altitude",
                              "uncategorised"]
        # find profile peaks as first separation
        profile_peaks=self.find_profile_peaks(rf_df,alt_var=alt_var,
                                              use_turb=use_turb)
        # get gradients of altitude
        if not use_turb:
            alt_grad=rf_df[alt_var].diff()
            alt_grad_30s=alt_grad.rolling("30s",center=True).mean()
        else:
            # resolution of turb is higher, regrid to 1 Hz
            alt_grad=rf_df[alt_var].resample("1s").mean().diff()
            alt_grad_30s=alt_grad.rolling("30s",center=True).mean()
            
        snd_derivative_alt=alt_grad_30s.diff().rolling(
            "30s",center=True).mean()
            
        rf_df["segments"]="none"
        alt_grad_threshold=rate_threshold
        # Assign indices
        asc_index=alt_grad_30s[alt_grad_30s>alt_grad_threshold].index
        dsc_index=alt_grad_30s[alt_grad_30s<-1*alt_grad_threshold].index
        con_alt_index=alt_grad_30s[alt_grad_30s.between(-1*alt_grad_threshold,
                                                    alt_grad_threshold)].index
        #0.005
        snd_dev_index=snd_derivative_alt[snd_derivative_alt.between(-.0075,.0075)].index
        con_index=con_alt_index.intersection(snd_dev_index)
        gnd_index=con_index.intersection(rf_df[rf_df[alt_var]<35].index)

        # Assign segment flags
        if not use_turb:    
            rf_df["segments"].loc[asc_index]="ascent"
            rf_df["segments"].loc[dsc_index]="descent"
        else:
            print("do the ascent flagging")
            for a,idx in enumerate(asc_index):
                performance.updt(len(asc_index),a)
                rf_df["segments"][rf_df.index.round("s")==idx]="ascent"
            print("do the descent flagging")
            for d,idx in enumerate(dsc_index):
                performance.updt(len(dsc_index),d)
                rf_df["segments"][rf_df.index.round("s")==idx]="descent"
        # Update flight segmentation with added profile information
        # Merge flag ascent/descent as boolean
        added_rf_df=pd.DataFrame(data=rf_df[[alt_var,"segments"]],
            columns=[alt_var,"segments"],index=rf_df.index)
        added_rf_df['is_ascent_or_descent'] = added_rf_df['segments'].isin(
                            ['ascent', 'descent'])
        added_rf_df["DATETIME"]=added_rf_df.index
        
        # Assign a group id to each connected ascent/descent segment
        added_rf_df['profile_group'] = \
            (~added_rf_df['is_ascent_or_descent']).cumsum()
        
        segment_stats = added_rf_df.groupby('profile_group').apply(
                        compute_height_change)

        # Step 4: Merge segment stats back into the original DataFrame
        added_rf_df = added_rf_df.merge(segment_stats, on='profile_group')

        # Step 5: Flag segments longer than 100 m in height change as 'profile'
        added_rf_df['segment_category'] = \
            added_rf_df.apply(assign_profile, axis=1)

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
        print("number the profiles")
        for i in range(len(splitted)):
            performance.updt(len(splitted),i)
            if len(splitted.iloc[i]) == 2:
                added_rf_df.at[i, "segment_category"] = \
                    splitted.iloc[i][0]+\
                    " " +str(int(added_rf_df["profile_number"].iloc[i])).zfill(2) +\
                        " "+splitted.iloc[i][1]
        added_rf_df=added_rf_df.dropna(subset="profile_number")
        added_rf_df.index=added_rf_df["DATETIME"]
        added_rf_df=added_rf_df.drop_duplicates(
            subset="DATETIME",keep="first")         
        rf_df["segments"].loc[added_rf_df.index]=added_rf_df["segment_category"]
        #---------------------------------------------------------#
        peak_period=15
        if use_turb:
            peak_period=750
        for peak in profile_peaks.index:
            
            if not use_turb:
                peak_idx=rf_df.index.get_loc(peak)
            else:
                peak_idx=rf_df.index.round("s").get_loc(peak).start
                
            rf_df["segments"].iloc[peak_idx-peak_period:\
                                   peak_idx+peak_period]="peak"
        
        if not use_turb:
            rf_df["segments"].loc[con_index]="nearly_constant_altitude"
            rf_df["segments"].loc[gnd_index]="near_ground"
        else:
            print("Do the constant index")
            for c,idx in enumerate(con_index):
                performance.updt(len(con_index),c)
                rf_df["segments"][rf_df.index.round("s")==idx]="nearly_constant_altitude"
            print("Do the gnd index")
            for g,idx in enumerate(gnd_index):
                performance.updt(len(gnd_index),g)
                rf_df["segments"][rf_df.index.round("s")==idx]="near_ground"
        rf_df["segments"].loc[profile_peaks.index]="max"
        # Sometimes peaks are interrupted by constant altitude segments due to
        # thresholds in vertical velocities. Change them to peaks. 
        for i in range(1, len(rf_df)):
            prev_segment = rf_df["segments"].iloc[i - 1]
            current_segment = rf_df["segments"].iloc[i]
            # If previous segment was 'peak' and current is 'nearly constant altitude', merge them
            if prev_segment=='peak' and current_segment=='nearly_constant_altitude':
                rf_df["segments"].iloc[i] = 'peak'
            # --> If previous was 'peak' and current is also 'peak',
            # it's a new 'peak' period, do nothing
            # If previous was 'peak' and current is something else, no merge
        return rf_df, profile_peaks
    
    def segments_cf_conform_IDs(self,rf_df):
        rf_df["flight_segments"]=0
        rf_df["flight_segments"].loc[rf_df["segments"]=="near ground"]=0
        rf_df["flight_segments"].loc[rf_df["segments"]=="const altitude"]=1
        rf_df["flight_segments"].loc[rf_df["segments"].str.contains("ascent")]=100
        rf_df["flight_segments"].loc[rf_df["segments"].str.contains("descent")]=-100
        rf_df["flight_segments"].loc[rf_df["segments"]=="peak"]=2
        rf_df["flight_segments"].loc[rf_df["segments"]=="max"]=3
        
        unique_ascents  = np.unique(rf_df["segments"].loc[rf_df["segments"].str.contains("ascent")])
        unique_descents = np.unique(rf_df["segments"].loc[rf_df["segments"].str.contains("descent")])  
        for asc in unique_ascents:
            if " " in asc:
                rf_df["flight_segments"].loc[rf_df["segments"]==asc]+=int(asc[8:10])
            else:
                continue
        for dsc in unique_descents:
            if " " in dsc:
                rf_df["flight_segments"].loc[rf_df["segments"]==dsc]-=int(dsc[8:10])
            else:
                continue
        del rf_df["segments"]
        return rf_df                                       
    def plot_p_sfc(self,p,p_sfc_values,t,t_sfc,sensor="bme"):
        mean_p_sfc=p_sfc_values.mean()
        mean_t_sfc=t_sfc.mean()
        #pressure quicklooks
        pres_fig=plt.figure(figsize=(18,16))
        ax1=pres_fig.add_subplot(311)
        ax1.plot(p.index,p["pressure"]/100,color="grey",lw=1,)
        ax1.scatter(p_sfc_values.index,
                    p_sfc_values["pressure"]/100,s=30,
        color="orange",label="$\overline{p_{\mathrm{sfc}}}$="+\
            str(round(mean_p_sfc["pressure"]/100,3))+" hPa")
        ax1.set_ylabel(sensor.upper()+" Pressure (hPa)")
        ax1.legend()
        ax1.invert_yaxis()
    
        ax2=pres_fig.add_subplot(312)
        ax2.scatter(p_sfc_values.index,p_sfc_values["pressure"]/100,s=100,
                color="orange",edgecolor="k",lw=1)
    
        p_fit,p_tendency=self.calc_sfc_pressure_tendency(p_sfc_values)
        p["p_sfc"]=p_tendency.slope*p["int_idx"]+p_tendency.intercept
        ax2.plot(p_sfc_values.index,p_fit/100,color="orange",ls="--",lw=2,
                 label="$dp/dt$"+\
                 str(round(p_tendency.slope*3600/100,1))+" hPa h"+"$^{-1}$")
        ax2.set_ylabel("Pressure")
        ax2.legend()
        ax3=pres_fig.add_subplot(313)
        ax3.plot(t,color="grey")
        ax3.scatter(t_sfc.index,t.loc[t_sfc.index],
            color="red",s=30,label="$\overline{T_{\mathrm{sfc}}}$="+\
            str(round(mean_t_sfc,1))+"$^{\circ}\mathrm{C}$")
        ax3.set_ylabel("Temperature (K)")
        ax3.legend()
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax3.set_xlabel(
            'UTC Time of '+\
                self.flight_date +\
                    " ("+self.flight+")")
        fig_name=sensor.upper()+"-based_p_sfc_tendency_"+self.flight+"_"+\
            self.flight_date+".png"
        fig_path=self.process_plot_path+"/"
        fname=fig_path+fig_name
        pres_fig.savefig(fname,dpi=300,bbox_inches="tight")
        self.logger.info("plotted p_sfc tendency method under: "+fname)
    
    def plot_quicklook_zb(self,sensor="bme"):
        #rough quicklook of Z_b
        matplotlib.rcParams.update({"font.size":18})
        z_b_fig=plt.figure(figsize=(12,9))
        ax1=z_b_fig.add_subplot(111)
        ax1.plot(self.l1_df["z_b"],color="orange",label="z_b")
        if sensor=="TMP_sonde":
            ax1.plot(self.l1_df["ALT_R"]-31,color="grey",lw=1,ls="--",
                     label="GPS")
        ax1.grid()
        ax1.set_ylabel("Barometric height in m")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.set_xlabel('UTC Time of '+self.flight_date +\
            " ("+self.flight+")")
        ax1.legend()
        fig_name=sensor.upper()+"-based_z_b_"+self.flight+"_"+\
            self.flight_date+".png"
        fig_path=self.process_plot_path+"/"
        fname=fig_path+fig_name
        z_b_fig.savefig(fname,dpi=300,bbox_inches="tight")
        self.logger.info("plotted z_b method under: "+fname)
    
#-----------------------------------------------------------------------------#
class Meteorological_Probe(BELUGA):
    # The code for processing the TMP_met measurement data is mainly created by 
    # Holger Siebert (L1 data) and Henning Dorff (L2 data) using knowledge and 
    # methods from e.g. Lonardi et al. (2022) & Pilz et al. (2023).
    #%% Basics
    def __init__(self, BELUGA,rf="RF12",
                 datatypes=["radiosonde","ultra-sonic anemometer"],
                 run_L1_processing=True,run_L2_processing=True,
                 plot_processing=True):
        
        
        self.BELUGA_cls         = BELUGA
        self.probe              = "TMP_met"
        self.station_name       = BELUGA.station_name        
        self.coordinates        = BELUGA.coordinates
        self.main_path          = BELUGA.main_path
        self.flight_dates       = BELUGA.flight_dates
        self.raw_data_path      = self.main_path +\
           "/BELUGA_TMP_met_probe/raw_data/"
        self.flight_infos       = BELUGA.flight_infos
        self.flight             = rf
        self.flight_date        = self.flight_dates[self.flight]
        self.datatypes          = datatypes
        self.plot_path          = self.main_path+"/plotting/"
            
        self.process_plot_path  = self.plot_path+"/test_processing/"+\
            self.flight+"/"
        os.makedirs(self.process_plot_path,exist_ok=True)
        # File names
        self.l1_fname="L1_TMP_met_"+self.flight+"_"+self.flight_date[0:4]+"-"+\
            self.flight_date[4:6]+"-"+self.flight_date[6:8]+"-1Hz.csv"
        
        # Switches
        self.run_L1_processing  = run_L1_processing
        self.run_L2_processing  = run_L2_processing
        self.plot_processing    = plot_processing
    
        # Configure logging once at the module level
        self.log_file_name='BELUGA_TMP_met_processing_'+self.flight+'.log'
        self.logger=logging_processing.setup_logging(self.log_file_name)
        # Constants
        self.C_P                = 1005.0  # J/(kg·K)
        self.DENSITY_AIR        = 1.2922  # kg/m³
        self.R_L                = 287
        self.R_v                = 461
        self.L                  = 2.5e6
        self.g                  = 9.81
        self.rho_s              = 1.2
        
    def read_raw_files(self):
        raw_fstart = self.raw_data_path+self.flight+"/"+\
            self.flight
        self.msr_file  = raw_fstart + '_hotwire_1s.txt'
        self.ola_file  = raw_fstart + '_ola_1s.txt'
        self.raso_file = raw_fstart + '_radiosonde.txt'

 
        # OpenLogArtemis (OLA) includes to different data sets    
        # "_O"  :: OLA data (Bosch sensor)
        # "_S" :: Trisonica data (included in OLA data set)

        # "_R" :: Radiosonde data incl XDATA (incl HMC 3D compass)
        # "_M" :: MSR-145 data logger (for hot-wire), here with 1 Hz  
 
        # defining col titles

        name_M  = ["Time_M"  ,"p_M","RH_M","T_RH_M","A2"] 
        col_M   = [0,1,3,4,6]
 
        col_O   = ["Time_O","S_O","T_SENS_O","T_AMB_SENS_O","P_O", "H_O","ALT_O","T_O",
           "S_S","D_S","U_S", "V_S","W_S","T_S","RH_S","P_S","PI_S","RO_S","MD_S"]
 
        name_R  = ["P_R","T_R","RH_R","ALT_R","HEAD_X","VELO_X","PI_X","RO_X", "Time_R"] 
        # VELO_X :: U_Sonic via XDATA "_X"
        col_R   = [3,4,5,7,11,12,13,14,20]

        # reading all 1 Hz data files 
        # OLA
        if os.path.exists(self.ola_file):
            self.ola_df = pd.read_csv(self.ola_file,
                sep=',',skiprows=0,header=0,names = col_O, 
                usecols = col_O, engine='python',encoding='latin1',
                index_col = "Time_O",parse_dates = True)
            self.ola_df = self.ola_df[~self.ola_df.index.duplicated(
                keep='first')] 
            # remove duplicated indexes
            # -->
            self.logger.info("OLA file read for "+self.flight)
        else:
            self.logger.info("No OLA file available for "+self.flight)
        #
        if os.path.exists(self.raso_file):
            self.raso_df = pd.read_csv(self.raso_file,
                sep=',',skiprows=0,header=0,names = name_R, 
                usecols = col_R, engine='python',encoding='latin1',
                index_col = "Time_R",parse_dates = True)
            self.raso_df = self.raso_df[~self.raso_df.index.duplicated(
                keep='first')] 
            # remove duplicated indexes
            # -->
            
            self.logger.info("Sonde file read for "+self.flight)
        else:
            self.logger.info("No sonde file available for "+self.flight)
        
        if os.path.exists(self.msr_file):
            self.msr_df = pd.read_csv(self.msr_file,  sep=',',skiprows=0,
                header=0,names = name_M, usecols = col_M, engine='python',
                encoding='latin1',index_col = "Time_M",
                parse_dates = True)
            self.msr_df = self.msr_df[~self.msr_df.index.duplicated(
                keep='first')] 
            # remove duplicated indexes
            # -->
            
            self.logger.info("MSR file read for "+self.flight)
            if self.msr_df.index.tz is not None:
                self.msr_df.index   = (self.msr_df.index   - \
                                       pd.Timedelta(seconds = 7200)).tz_localize(None)
                self.logger.info("MSR time zone adjusted")
            else:
                self.logger.info("MSR data is already in naive local time")   
        else:
            self.logger.info("No MSR file available for "+self.flight)   
        
    def read_l1_data(self):
        l1_fname=self.raw_data_path+"/../temporary_data/"+self.l1_fname
            #self.flight+"-"+self.flight_date[0:4]+"-"+\
            #self.flight_date[4:6]+"-"+self.flight_date[6:8]+"-1Hz.csv"
        self.l1_df=pd.read_csv(l1_fname)
        self.l1_df.rename(columns={"Date/Time in UTC":"time"},inplace=True)
        self.l1_df.index=pd.DatetimeIndex(self.l1_df["time"])
        del self.l1_df["time"]
        self.l1_df = self.l1_df.reset_index().drop_duplicates(
            subset='time', keep='first').set_index('time')
    
    def rename_columns(self):
        self.l1_df.rename(columns={"pp":"p","zb":"z_b","T_S":"sonic_T",
                                   "RH":"rh","UU":"vv","DD":"dd",
                                   "U_lon":"u_lon","U_lat":"u_lat",
                                   "PI":"pitch"},inplace=True)
    #%%Processing
    def run_processing(self,version_number="_v2.1",rf="RF12",
            datatypes=["radiosonde","ultrasonic"]):
        self.version_number=version_number
        self.logger.info(
            f"Processing flight {self.flight} on date {self.flight_date}\n")
        if self.run_L1_processing:
            self.read_raw_files()
            self.L0_to_L1_processing()
        if self.run_L2_processing:
            self.L1_to_L2_processing()
    
    #%% L0 to L1 processing
    def synchronise_times(self):    
        # estimate possible time shifts between different data sets and
        # correct for radiosonde (GPS based) is taken as reference
        # MSR RTC is slightly drifting and not GPS-based
        ##%% shift Sonic data inside the OLA data by dt_S seconds
        # shift Sonic data inside the OLA data by dt_S seconds
        self.logger.info("Synchronise L0 data files in time")
        if hasattr(self,"ola_df"):
            column_mask = self.ola_df.columns.str.contains("_O")
            ola         = self.ola_df.loc[:, column_mask]
            sonic       = self.ola_df.loc[:, ~column_mask]
            if hasattr(self,"raso_df"):
                if not self.flight=="RF12":
                    dt_S         = int(self.raso_df.P_R.idxmin().timestamp() -\
                                   self.ola_df.P_S.idxmin().timestamp() ) 
                else:
                    dt_S=0
                self.logger.info("Sonic and OLA data syncronized, shift by "+\
                                             str(dt_S)+" seconds")
                sonic.set_index((sonic.index + pd.Timedelta(seconds=dt_S)),
                            inplace=True)
            df_O = ola.merge(sonic, left_index=True, right_index=True,
                             how="outer").sort_index()
        if hasattr(self,"msr_df"):
            # # merge radiosonde and MSR data 
            if hasattr(self, "raso_df"):
                temporary_raso_df = self.raso_df.merge(
                self.msr_df, left_index=True, right_index=True, how='right') 
                if not self.flight=="RF12":
                    dt_M   = int(temporary_raso_df.p_M.idxmin().timestamp() - \
                         temporary_raso_df.P_R.idxmin().timestamp())  
                else:
                    dt_M=0
                self.msr_df.set_index((self.msr_df.index - \
                                       pd.Timedelta(seconds=dt_M)),inplace=True)
                self.logger.info("MSR is corrected for time shift to sonde by "+\
                             str(dt_M)+" seconds")

        if hasattr(self,"raso_df"):
            # define a period with all devices available
            s_R = self.raso_df.index.min()
            e_R = self.raso_df.index.max()

            if hasattr(self,"ola_df"):
                s_O = self.ola_df.index.min()
                e_O = self.ola_df.index.max()
    
                start = max(s_O.tz_localize(None),s_R)
                end   = min(e_O.tz_localize(None),e_R)
            else:
                start = s_R
                end   = e_R
        else:
            start = self.ola_df.index.min()
            end  =  self.ola_df.index.max()
            
        if hasattr(self,"msr_df"):
            s_M = self.msr_df.index.min()
            e_M = self.msr_df.index.max()

            start = max(start,s_M)
            end   = min(end,  e_M)
        
#        if hasattr(self,"ola_df"):
#            self.ola_df=self.ola_df.loc[start:end]
            
#        if hasattr(self,"raso_df"):
#            self.raso_df=self.raso_df.loc[start:end]
            
#        if hasattr(self,"msr_df"):
#            self.msr_df=self.msr_df.loc[start:end]
        
        self.logger.info("Period with all devices available: "+\
            str(start)+"  -- "+ str(end))
        
    def correct_p_offsets(self):
        #Pressure 
        #----------------------------------------------#
        # Correction using minimum value
        # with sonde data as reference"
        if hasattr(self,"msr_df"):
            # Correct for minimum value
            dp_M     = self.raso_df.P_R.min() - self.msr_df.p_M.min()
            self.msr_df.p_M +=  dp_M
            self.logger.info("p min MSR corrected: "+\
                str(self.msr_df.p_M.min())+" hPa")
        
        if hasattr(self,"ola_df"):        
            dp_O     = self.raso_df.P_R.min() - self.ola_df.P_O.min()
            self.ola_df.P_O += dp_O
            self.logger.info("p min OLA corrected: "+\
                             str(self.ola_df.P_O.min())+" hPa")
            #print("p max OLA:  ", df_O.P_O.max())
            #print("p max MSR:  ", df_M.p_M.max())
            #print("p max RASO: ", df_R.P_R.max())
    
    def remove_spikes(self):
    #----------------------------------------------#
        T_window  = 50
        rh_window = 5
        p_window  = 5
        # Temperature
        try:
            # create low-pass-filtered data to compare with original data set
            T_low        = self.raso_df.T_R.rolling(window = T_window, 
            center=True).median().fillna(0) 
            self.raso_df = self.raso_df.assign(
                T_R_fil = self.raso_df.T_R.where(
                    (self.raso_df.T_R - T_low).abs() < 10 , np.NaN).\
                    interpolate(method ='linear', limit_direction ='forward'))
            self.raso_df.T_R_fil = self.raso_df.T_R_fil[self.raso_df.T_R_fil<0] 
            self.logger.info("Sonde-based T despiked")
        except:
            self.logger.error("No spike removal for sonde-based T")
        
        # Relative humidity
        try:
            self.raso_df.RH_R = self.raso_df.RH_R.rolling(
                window = rh_window,center=True).median()
            self.logger.info("RH smoothed with rolling median filter using "+\
                             str(rh_window)+" elements.")
        except:
            self.logger.error("No smoothing for RH.")
        # Pressure
        try:
            self.raso_df.P_R=self.raso_df.P_R.rolling(
                window=p_window,center=True).median()
            self.logger.info("P smoothed with rolling median filter using " +\
                             str(p_window)+" elements.")
        except:
            self.logger.error("No smoothing for pressure.")

    def merge_all_TMP_met_datasets(self):
        
        if hasattr(self, "ola_df"): 
            self.ola_df = self.ola_df.tz_localize(None)
            self.l1_df   = self.raso_df.merge(self.ola_df, 
                left_index=True, right_index=True, how="outer").sort_index() 
            self.logger.info("OLA data with radiosonde data merged")
        else:
            self.logger.error("No OLA data available to merge,"+\
                              " l1_df includes only sonde data")
            self.l1_df = self.raso_df
        
        if hasattr(self,"msr_df"):    
            self.l1_df = self.l1_df.merge(self.msr_df,
                left_index=True, right_index=True, how="outer").sort_index() 
            self.logger.info("L1_df includes MSR data")
       
    def define_flight_start_end_and_cut(self):
        # Define Start and End of flight when ascent and descent => subrecord
        offset = 10 # Seconds

        id_a = (self.l1_df.where(
            self.l1_df.P_R.diff()< -0.1).first_valid_index()) -\
                pd.Timedelta(seconds=offset)
        id_e = (self.l1_df.where(
            self.l1_df.P_R.diff()  >  0.1).last_valid_index())  + \
                pd.Timedelta(seconds= offset)

        self.l1_df = self.l1_df[id_a:id_e]
        self.logger.info("Set flight data start to:"+ str(id_a))
        self.logger.info("Set flight data end to:"+ str(id_e))
        
    def correct_sonic_temperature(self):
        # estimate offset of sonic temperature and correct T_Sonic
        self.p_GND = self.l1_df.P_R.max()
        self.T_GND  =  self.l1_df.T_R_fil[(
            self.l1_df.P_R > (self.p_GND - 0.2))].mean()
        p_range = 2 # in hPa defines data 16 m above surface (g_GND - 2hPa)
        #           for removing surface effects
        if hasattr(self,"ola_df"):
            # get surface unaffected temperature median values
            Ts_m = self.l1_df.T_S[self.l1_df.P_O < 
                                  (self.p_GND - p_range)].median() 
            Tr_m = self.l1_df.T_R_fil[self.l1_df.P_R < 
                                  (self.p_GND - p_range)].median()
            dT= Tr_m-Ts_m
            # Offset correction
            self.l1_df.T_S += dT
            # Unaffected ground temperature
            # surface temperature defined between surface (p_max) and 8/5 m )
             
            self.logger.info("Sonic temperature offset of "+str(dT)+\
                             " K is corrected")
            
    def add_derived_parameters(self):
        def ES(T):    
            return 6.1 * np.exp(17.62 * T /(243.12 + T))  
        def e(ES,rH):
            return rH / 100 * (ES * 100)
        def absH(self, rH,ES,T):
            return (rH * ES) / self.R_v / (T+self.T_GND)
        def mixRatio(self,rH,ES,p):
            return rH/100 * ES / p * self.R_L/self.R_v  * 1e3
        def potT(T,p):
            return (T + self.T_GND) * (self.p_GND/p)**(self.R_L / self.C_P)   
        #def zb(self,p):
        #    return self.R_L * (self.T_GND + 273.15) / \
        #        self.g * np.log(self.p_GND / p)
        
        def rho(T,p):
            return(p*100/(self.R_L*(T+273.15)))
        
        self.l1_df  = self.l1_df.assign(
            potT = potT(self.l1_df.T_R_fil,self.l1_df.P_R))
        self.l1_df  = self.l1_df.assign(
            r_v  = mixRatio(self,self.l1_df.RH_R,
                            ES(self.l1_df.T_R_fil),self.l1_df.P_R))             
        
        # Get barometric pressure altitude
        self.BELUGA_cls.l1_df               = self.l1_df.copy()
        self.BELUGA_cls.logger              = self.logger
        self.BELUGA_cls.plot_processing     = self.plot_processing
        self.BELUGA_cls.process_plot_path   = self.process_plot_path
        self.BELUGA_cls.flight              = self.flight
        self.BELUGA_cls.flight_date         = self.flight_date
        
        self.BELUGA_cls.get_barometric_altitude(sensor="TMP_sonde")
        self.l1_df["zb"]=self.BELUGA_cls.l1_df["z_b"]
        
        if hasattr(self, "msr_df"):
            try:
                self.l1_df  = self.l1_df.assign(
                    rho  = rho(self.l1_df.T_RH_M,self.l1_df.P_R))
            except:
                self.l1_df  = self.l1_df.assign(
                    rho  = rho(self.l1_df.T_RH_M_x,self.l1_df.P_R))
    
    def include_radiosonde_based_windspeed(self):
        file_path=self.raw_data_path+"/../temporary_data/BELUGA_with_radiosonde_wind/"
        fname=glob.glob(file_path+self.flight+"*")
        wind_df=pd.read_csv(fname[0],index_col=0)
        wind_df.index=pd.DatetimeIndex(wind_df.index)
        wind_df = wind_df.loc[~wind_df.index.duplicated(keep='first'), :]
        intersect_index=self.l1_df.index.intersection(wind_df.index)
        self.l1_df["UU"].loc[intersect_index]=wind_df["UU"].loc[intersect_index]
        self.l1_df["DD"].loc[intersect_index]=wind_df["DD"].loc[intersect_index]
   
                    
    def average_cyclic_data(self):
        def cyclic_mean(phi):
            phi = phi / 180. * np.pi
            s = np.sin(phi).mean() # mean of unit vector coponents in x
            c = np.cos(phi).mean() # mean of unit vector coponents in y
            phi_mean = np.arctan(s/c) * 180 / np.pi #angle of mean unit vector 
                #(-90° to + 90° selecting valid quadrant according s & c
            if ((s>0) & (c>0)):
                return phi_mean 
            if ((s>0) & (c<0)):
                phi_mean = phi_mean + 180
                return phi_mean 
            if ((s<0) & (c<0)):
                phi_mean = phi_mean + 180
                return phi_mean
            if ((s<0) & (c>0)):
                phi_mean = phi_mean + 360
                return phi_mean
        
        def replace(value):
            if value < 180:
                return value + 360
            else:
                return value
        
        self.l1_df = self.l1_df.assign(
            HEAD_X_fil = self.l1_df.HEAD_X.rolling(
                5,center = True).apply(cyclic_mean))
        # 360° => 540° system for nicer plotting
        self.l1_df  = self.l1_df.assign(
             HEAD_X_fil_540 = self.l1_df.HEAD_X_fil.apply(replace))
        self.l1_df  = self.l1_df.assign(
             HEAD_X_540     = self.l1_df.HEAD_X.apply(replace))

    def thresholding_physical_ranges(self):
        # Physical thresholding
        # Temperature    
        try:
            bad_data_index=self.l1_df.loc[self.l1_df["T_R_fil"]>0].index
            self.l1_df["quality_flag"].loc[bad_data_index]=3
            self.l1_df[self.l1_df["T_R_fil"]>0] = np.nan
        except:
            self.l1_df[self.l1_df["T_R"]>0]     = np.nan
            
        # Relative humidity
        self.l1_df[self.l1_df["RH_R"]>120]        = 120
        # Wind speed
        try:
            self.l1_df[self.l1_df["S_S"]>20]         = np.nan
        except:
            self.logger.error("Windspeed could not be quality checked for this flight.")
        self.logger.info("checked data for physically plausible thresholds (T<0°C,RH<120,VV<20")
        
    def flag_for_quality(self):
        self.l1_df["quality_flag"]=0
        
        if "PI_S" in self.l1_df.columns:
            self.l1_df["quality_flag"].loc[abs(self.l1_df["PI_S"])>15]=1
            #"ok_wind_pitch_affected"
        self.l1_df["quality_flag"].loc[\
            self.l1_df["interpolation_flag"]==True]=2
            #"ok_interpolated"
        # 
        try:
            self.l1_df.loc[self.l1_df["T_R_fil"]==np.nan]=4
        except:            
            self.l1_df.loc[self.l1_df["T_R"]==np.nan]=4
        
        del self.l1_df["interpolation_flag"]   
        # call physical plausible thresholding
        self.thresholding_physical_ranges()
        self.logger.info("data flagged for quality")
    def interpolate_and_fillna(self):
        # Create a mask of original non-NaN points
        if "T_R_fil" in self.l1_df.columns:
            t_col="T_R_fil"
        else:
            t_col="T_R"
        original_mask = self.l1_df[t_col].notna()
        # interpolate for short periods frequently observed in the sonde data
        self.l1_df=self.l1_df.interpolate(method="linear",limit=10) 
        self.l1_df['interpolation_flag'] = \
        (~original_mask) & \
            self.l1_df[t_col].notna()
   
        #self.l1_df=self.l1_df.fillna(-999)
    
    def store_l1_data(self):
        original_cols=["P_R", "zb", "T_R_fil", "T_S", "RH_R", 
                       "S_S", "V_S", "U_S", "HEAD_X_fil", 
                       "PI_X","quality_flag"]
        output_cols = ["pp", "z_b", "T", "T_S", "RH",
                       "UU", "U_lon", "U_lat", "DD", "PI","quality_flag"]
        
        if hasattr(self,"ola_df"):
            self.logger.info("Save L1 data with wspeed and sonic T from OLA")
        else:
            self.l1_df = self.l1_df.assign(T_S  = np.nan)
            self.l1_df = self.l1_df.assign(S_S  = np.nan)
            self.l1_df = self.l1_df.assign(U_S  = np.nan)
            self.l1_df = self.l1_df.assign(V_S  = np.nan)
        
        self.l1_df_full=self.l1_df.copy()
        
        fpath=self.raw_data_path+"/../temporary_data/"
        # subsample relevant columns
        self.l1_df=self.l1_df_full[original_cols]
        # rename cols
        rename_cols=dict(zip(original_cols,output_cols))
        self.l1_df.rename(columns=rename_cols,inplace=True)
        # some research flights are filled with radiosonde-based data
        if self.flight in ["RF11","RF25"]:
            self.include_radiosonde_based_windspeed()
        
        # save df
        self.l1_df.to_csv(fpath+self.l1_fname, encoding='utf-8',
                          index=True, index_label = "Date/Time in UTC")
        store_string="TMP_met L1 data (1 Hz) is stored as "+\
                         fpath+self.l1_fname
        self.logger.info(store_string)
        print(store_string)        
        if self.plot_processing:
           self.plot_L1_processed_data()             
    def L0_to_L1_processing(self):
        self.logger.info(f"Starting L0 to L1 processing")
        self.synchronise_times()
        self.correct_p_offsets()
        self.remove_spikes()
        self.merge_all_TMP_met_datasets()
        self.define_flight_start_end_and_cut()
        self.correct_sonic_temperature()
        self.average_cyclic_data()
        self.interpolate_and_fillna()
            
        self.flag_for_quality()
        self.add_derived_parameters()
        self.store_l1_data()
    #%% L1 to L2 processing    
    def run_flight_segmentation(self,segm_alt_var="z_b"):
        # Use the automatised modules from the BELUGA cls, choosing the 
        # altitude variable from instrument probe that has to be considered
        # for flight segmentation
        self.seg_alt_var                 = segm_alt_var
        self.BELUGA_cls.flight           = self.flight
        # Adapt the height change rate for specific flights
        rate_threshold=.6
        if self.flight in ["RF03","RF04","RF06","RF08","RF10","RF11","RF13",
                           "RF14", "RF15","RF16","RF17","RF18",
                           "RF23","RF24","RF26"]:
            rate_threshold=.3
        if self.flight=="RF01":
            rate_threshold=.45
        segmented_df,self.profile_peaks  = \
            self.BELUGA_cls.segment_flight_sections(self.l2_df,
                rate_threshold=rate_threshold,
                alt_var=segm_alt_var)
        self.logger.info(
            "Data categorised in flight segments based on altitude")
        
        return segmented_df
 
    def L1_to_L2_processing(self):
        self.logger.info(f"Starting L1 to L2 processing")
        # Get L1 data
        try:
            self.read_l1_data()
            self.rename_columns()
        except:
            self.read_raw_files()
            self.L0_to_L1_processing()
            self.rename_columns()
        self.l1_df=self.BELUGA_cls.replace_inf(self.l1_df)
        if self.plot_processing:
            self.plot_thermodynamic_profiles()
            
        self.l2_df=self.l1_df.copy()
        # Flight segmentation
        self.l2_df=self.run_flight_segmentation(segm_alt_var="z_b")
        if self.plot_processing:
            self.plot_flight_segmentation(add_peak_infos=True)

        self.l2_df=self.segments_cf_conform_IDs(self.l2_df)
        
        if self.flight=="RF08":
            self.exchange_T_with_sonic_T()
                
            
        self.changeunits()
        self.df2netcdf()
        
        self.BELUGA_cls.flight_date=self.flight_date
        self.l2_ds=self.BELUGA_cls.add_L2_meta_data(self.l2_ds,probe=self.probe)
        self.logger.info("Meta data (global attributes and variable attributes)"+\
            "added to L2 dataset.")
            
        self.l2_ds=self.BELUGA_cls.add_VRS_georeferences(self.l2_ds)
        self.logger.info("Georeference data of Station Nord added to Dataset")
        self.l2_ds=self.BELUGA_cls.add_ABL_transition_type(self.l2_ds)
        self.logger.info("ABL transition type added to L2 dataset")
        # cf role
        self.l2_ds=self.BELUGA_cls.add_rf_as_id_for_cf_role(self.l2_ds)
        # save file
        self.BELUGA_cls.raw_data_path=self.raw_data_path
        self.l2_fname=self.BELUGA_cls.save_L2_data_as_nc_file(self.l2_ds,
            self.version_number,probe=self.probe)
        self.logger.info("temporary L2 data saved as: "+\
                         self.l2_fname)
   
    def df2netcdf(self):
        self.l2_ds=xr.Dataset.from_dataframe(self.l2_df)
        self.logger.info("Transformed dataframe to xr.Dataset")
    def exchange_T_with_sonic_T(self):
        high_values=self.l2_df[self.l1_df["T"]>-20].index
        self.l2_df["T"].loc[high_values]=self.l2_df["sonic_T"].loc[high_values].values
        self.l2_df["quality_flag"].loc[high_values]=5
        nan_index=self.l2_df["T"][self.l2_df["T"].isna()].index
        self.l2_df["T"].loc[nan_index]=self.l2_df["sonic_T"].loc[nan_index]
        self.l2_df["quality_flag"].loc[nan_index]=5
    def changeunits(self):
        self.l2_df["T"].loc[self.l2_df["T"]!=-999]              += 273.15
        self.l2_df["sonic_T"].loc[self.l2_df["sonic_T"]!=-999]  += 273.15
        self.logger.info(
            "Changed units for temperature to Kelvin")
    #%% Plot processing
    def plot_L1_processed_data(self):
        plt.close("all") 
        fig, (ax1,ax2,ax3,ax4,ax5 ) = \
            plt.subplots(figsize=(15,10),nrows = 5 )
        plt.rcParams['font.size'] = 14
        
        fig.suptitle(self.flight+" "+self.flight_date+\
                     "\n Overview of stored data", fontsize=20)
        
        ax1.set_ylim([0,1000])
        ax2.set_ylim([-35,0])
        ax3.set_ylim([0,360])
        ax4.set_ylim([0,15])
        ax5.set_ylim([-30,30])
        
        
        ax1.set_xlabel("Time")
        ax1.set_ylabel("$z_b$ / m")
        ax2.set_ylabel("$T$ / °C")
        ax3.set_ylabel("$D$ / °")
        ax4.set_ylabel("$U / \mathrm{m\,s}^{-1}$")
        ax5.set_ylabel("Pitch angle /   $^\circ$")
        ax3.set_yticks(np.arange(0, 360.1, 90))
        
        # Start plotting
        
        ax1.plot(self.l1_df.z_b)
        if hasattr(self, "ola_df"):
            ax2.plot(self.l1_df.T_S,label="$T_{\mathrm{sonic}}$")
        ax2.plot(self.l1_df["T"],label ="$T_{\mathrm{sonde}}$") 
        ax2.legend(loc='upper right') 
            
        # Angles
        ax3.scatter(self.l1_df_full.index,
                    self.l1_df_full.HEAD_X,
                    label=r"$D_{\mathrm{HMC}}$",s=1)
        ax3.scatter(self.l1_df.index,self.l1_df.DD,
                    label=r"$<D_{\mathrm{HMC}}>$",s=1)
        ax3.legend(loc='upper right') 
        
        # wind speeds
        if hasattr(self,"ola_df"):
            ax4.plot(abs(self.l1_df.UU),label    = r"$U_{\mathrm{S}}$")
            ax4.plot(abs(self.l1_df.U_lon),label = r"$U_{\mathrm{lon}}$")
            ax4.plot(abs(self.l1_df.U_lat),label = r"$U_{\mathrm{lat}}$")
        
            ax4.legend(loc='upper right')  
        
        # pitch angle
        ax5.plot(self.l1_df.PI,label = r"$\rm Pi $") 
        #cos = np.cos((df.PI_X / 180 * np.pi))
        # ax5.plot(cos,label = r"$\rm Pi $") 
        # fig.legend(loc='upper right')
        plot_path=self.process_plot_path
        fig_file="TMP_met_L1_variables_"+self.flight+"_"+\
            self.flight_date+".png"
        plt.savefig(plot_path+fig_file, dpi=300, bbox_inches="tight")
        plot_file_string="TMP_met L1 data plotted as: "+plot_path+fig_file
        self.logger.info(plot_file_string)
        print(plot_file_string)


    def plot_flight_segmentation(self,add_peak_infos=False, probe="TMP_met"):
        
        rf=self.flight
        
        df_initial = self.l2_df.copy()
        df         = df_initial.copy()
        only_asc_dsc_df=df_initial.loc["2024-03-24 13:00":"2024-03-24 13:50"]
        
        df["z_b"][df["z_b"]==np.Inf]=np.nan
        # Do not separate between profile numbers anymore
        df["segments"][df_initial["segments"].str.contains(
            'ascent', case=False, na=False)] = 'ascent'
        df["segments"][df_initial["segments"].str.contains(
            'descent', case=False, na=False)] = 'descent'
        #---------------------------------------------------------------------#
        matplotlib.rcParams.update({"font.size":28})
        #% BELUGA Plotting
        fig=plt.figure(figsize=(16,9))
        
        ax1=fig.add_subplot(111)
        ax1.plot(df["z_b"],lw=4,color="w")
        ax1.plot(df["z_b"],lw=2,color="grey")
        
        # Colour-code flight segmentation
        segmentation_classes=["ascent","descent","peak",
                              "near_ground", "nearly_constant_altitude"]
        segmentation_legend=["ascent","descent","peak",
                              "near ground", "nearly constant altitude"]
        segmentation_cls_colors=["blue","orange","red",
                                 "sienna","green"]
        
        for s, seg in enumerate(segmentation_classes):
            seg_df=df["z_b"][df["segments"]==seg]
            ax1.scatter(seg_df.index,seg_df,marker="s",s=40,
                        color=segmentation_cls_colors[s],
                        zorder=2,label=segmentation_legend[s])
            
        # Add extra marker for maximum heights    
        ax1.scatter(self.profile_peaks.index,
            self.profile_peaks, marker="o", s=100,
            color="red",edgecolor="k",zorder=3)
        # only ascent but not profile
        #if rf=="RF01":
            
        #    ax1.scatter(only_asc_dsc_df.loc[\
        #                    only_asc_dsc_df["segments"]=="ascent"].index,
        #                only_asc_dsc_df["z_b"].loc[\
        #                    only_asc_dsc_df["segments"]=="ascent"].values,
        #                s=50,color="lightblue",zorder=2)
            
        #    ax1.scatter(only_asc_dsc_df.loc[\
        #                    only_asc_dsc_df["segments"]=="descent"].index,
        #                only_asc_dsc_df["z_b"].loc[\
        #                    only_asc_dsc_df["segments"]=="descent"].values,
        #                s=50,color="yellow",zorder=2)
        max_max=df["z_b"].dropna().max()
        ax1.scatter(df["z_b"].idxmax(),max_max,
                    marker="o",s=200,color="darkred",
                    edgecolor="k",zorder=4)
        ax1.axhline(y=200, color="darkgrey", lw=3,ls="--")
        
        # Add profile peak infos
        if add_peak_infos:
            ax1.text(df["z_b"].idxmax(),max_max+20,
                 str(int(max_max))+" m")
            ax1.text(0.1,0.9,self.flight+": Number of profiles: "+\
                     str(2*len(self.profile_peaks)),
                 transform=ax1.transAxes)
        
        ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
        ax1.set_ylabel("Barometric height (m)")
        
        if max_max>600:
            ylim_max=900
        else:
            ylim_max=600
        
        ax1.set_ylim([0,ylim_max])
        ax1.legend(fontsize=24)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        
        sns.despine(offset=10)
        for axis in ['bottom','left']:
            ax1.spines[axis].set_linewidth(2)
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)

        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        plot_path=self.process_plot_path
        fig_name="BELUGA_"+probe+"_flight_segmentation_"+rf
        file_end=".png"
        fig_name+=file_end
        fname=plot_path+fig_name
        fig.savefig(fname, dpi=600, bbox_inches="tight")
        print("Figure saved as:", fname)
        self.logger.info("Figure saved as:"
                         +plot_path+fname)
        
        #self.get_barometric_altitude()
        #self.l2_df=self.l1_df.copy()
        #if self.plot_processing:
        #    self.plot_up_and_down_ward_terr_radiation()
        #    self.plot_net_terr_radiation()
        #self.l2_df=self.run_flight_segmentation(segm_alt_var="z_b")
        #if self.plot_processing:
        #    self.plot_flight_segmentation(add_peak_infos=True)
        #self.read_l1_data()
    def plot_thermodynamic_profiles(self):
        matplotlib.rcParams.update({"font.size":24})
        
        from matplotlib.lines import Line2D
        from matplotlib.gridspec import GridSpec

        df=self.l1_df.copy()
        df.name=self.flight
        moisture_label = "rh"
        moisture_unit  = " (%)"
        
        quick_scatter=plt.figure(figsize=(16,12))
        import beluga_plotting
        gs = GridSpec(2, 3, figure=quick_scatter,
                          height_ratios=[0.35,1],wspace=.2)    
        ax0=quick_scatter.add_subplot(gs[0,:])
        beluga_plotting.plot_RF_height(df,alt_var="z_b",ax_obj=ax0)      
        ax1=quick_scatter.add_subplot(gs[1,0])
        ax2=quick_scatter.add_subplot(gs[1,1])
        ax3=quick_scatter.add_subplot(gs[1,2])
        sub_fig_labels=["(b)","(c)","(d)"]
        # Create a proxy artist for the legend with a larger marker size
        legend_marker_size = 20  # Larger size for the legend
        beluga_lgd_handle = Line2D([0], [0], marker='o', color='salmon',
                markersize=legend_marker_size, label='BELUGA')
        # Scatter plots
        ax1.scatter(df["T"],df["z_b"],color="salmon",s=1) # BELUGA
        ax1.set_yticks([0,100,200,300,400,500,600,700])
        ax1.set_xlim([int(df["T"].dropna().min())-2,
                      int(df["T"].loc[lambda v: v<np.Inf].max())+2])
        ax1.set_ylabel("Barometric height (m)")
        ax1.set_xlabel("Temperature (K)")
        #ax1.legend(handles=[beluga_lgd_handle],
        #       loc="center left",fontsize=16)
        #-------------------------------------------------------------------------#
        # Wind speed
        ax2.scatter(df["vv"],df["z_b"],color="lightgreen",s=1)
        
        ax2.set_yticks([0,100,200,300,400,500,600,700])
        if not df["vv"].dropna().shape[0]==0:
            ax2.set_xlim([int(df["vv"].dropna().min())-2,
            int(df["vv"].loc[lambda v: v<np.Inf].max())+2])
            if df["vv"].max()<7.5:
                ax2.set_xlim([0,7.5])
            else:
                ax2.set_xlim([0,10])
        else:
            ax2.set_xlim([0,10])
        ax2.set_xlabel("Wind speed ($\mathrm{m\,s}^{-1}$)")
        #-------------------------------------------------------------------------#
        # Relative humidity
        ax3.scatter(df[moisture_label],df["z_b"],color="lightblue",s=1)
        ax3.set_yticks([0,100,200,300,400,500,600,700])
        
        ax3.set_xlim([0,100])
        ax3.set_xlabel(moisture_label.upper()+moisture_unit)
        
        ax1.set_ylim([0,700])
        ax2.set_ylim([0,700])
        ax3.set_ylim([0,700])
        
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
        
        ax1.text(0.03,0.95,sub_fig_labels[0],
                 transform=ax1.transAxes,fontsize=20)
        ax2.text(0.03,0.95,sub_fig_labels[1],
                 transform=ax2.transAxes,fontsize=20)
        ax3.text(0.03,0.95,sub_fig_labels[2],
                 transform=ax3.transAxes,fontsize=20)
        
        plt.subplots_adjust(wspace=0.5,hspace=0.3)
        sns.despine(offset=10)
        plot_path=self.process_plot_path
        os.makedirs(plot_path,exist_ok=True)
        major_name="Profiles_BELUGA_TMP_met_"
        fig_name=major_name+moisture_label+"_"+\
            self.flight+"_"+self.flight_date+".png"
        quick_scatter.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
        self.logger.info("Figure saved as:"+plot_path+fig_name)
        
#-----------------------------------------------------------------------------#
#%% Turbulence probe
class Turbulence_Probe(BELUGA):
    # The code for processing the TMP_met measurement data is mainly created by 
    # Holger Siebert (L1 data) and Henning Dorff (L2 data) using knowledge and 
    # methods from e.g. Lonardi et al. (2022) & Pilz et al. (2023).
    #%% Basics
    def __init__(self, BELUGA,rf="RF22",
                 datatypes=["hot wire anemometer","radiosonde"],
                 run_L1_processing=True,run_L2_processing=True,
                 plot_processing=True):
        
        
        self.BELUGA_cls         = BELUGA
        self.probe              = "TMP_turb"
        self.station_name       = BELUGA.station_name        
        self.coordinates        = BELUGA.coordinates
        self.main_path          = BELUGA.main_path
        self.flight_dates       = BELUGA.flight_dates
        self.raw_data_path      = self.main_path +\
           "/BELUGA_TMP_turb_probe/raw_data/"
        self.flight_infos       = BELUGA.flight_infos
        self.flight             = rf
        self.flight_date        = self.flight_dates[self.flight]
        self.datatypes          = datatypes
        self.plot_path          = self.main_path+"/plotting/"
            
        self.process_plot_path  = self.plot_path+"/test_processing/"+\
            self.flight+"/"
        os.makedirs(self.process_plot_path,exist_ok=True)
        # File names
        self.l1_fname="L1_TMP_turb_"+self.flight+"_"+self.flight_date[0:4]+"-"+\
            self.flight_date[4:6]+"-"+self.flight_date[6:8]+"-1Hz.csv"
        
        # Switches
        self.run_L1_processing  = run_L1_processing
        self.run_L2_processing  = run_L2_processing
        self.plot_processing    = plot_processing
    
        # Configure logging once at the module level
        self.log_file_name='BELUGA_TMP_turb_processing_'+self.flight+'.log'
        self.logger=logging_processing.setup_logging(self.log_file_name)
        # Constants
        self.C_P                = 1005.0  # J/(kg·K)
        self.DENSITY_AIR        = 1.2922  # kg/m³
        self.R_L                = 287
        self.R_v                = 461
        self.L                  = 2.5e6
        self.g                  = 9.81
        self.rho_s              = 1.2
        
    def read_raw_file(self):
        #raw_fstart = self.raw_data_path+self.flight+"/"+\
        #    self.flight
        #self.msr_file  = raw_fstart + '_hotwire_1s.txt'
        # --> to be changed
        pass
    
    def read_l1_data(self):
        #RF14-2024-04-01-50Hz
        self.l1_fname=self.flight+"-"+self.flight_date[0:4]+"-"+\
            self.flight_date[4:6]+"-"+self.flight_date[6:8]+"-50Hz.csv"
        self.l1_df=pd.read_csv(self.raw_data_path+"/../temporary_data/"+\
                               self.l1_fname,index_col=0)
        self.l1_df.index=pd.DatetimeIndex(self.l1_df.index)
        self.l1_df.index=self.l1_df.index.rename('time')
    # Processing
    def run_processing(self,rf,version_number):
        self.version_number=version_number
        self.logger.info(
            f"Processing flight {self.flight} on date {self.flight_date}\n")
        if self.run_L1_processing:
            self.read_raw_files()
            self.L0_to_L1_processing()
            #self.merge_L1_data()
        if self.run_L2_processing:
            self.L1_to_L2_processing()
    
    def L0_to_L1_processing(self):
        pass
    
    def update_z_b_for_pressure_tendency(self):
        self.BELUGA_cls.logger           = self.logger
        self.BELUGA_cls.l1_df            = self.l1_df.copy()
        self.BELUGA_cls.flight_date      = self.flight_date
        self.BELUGA_cls.flight           = self.flight
        self.BELUGA_cls.plot_processing  = self.plot_processing
        self.BELUGA_cls.process_plot_path= self.process_plot_path
        self.BELUGA_cls.get_barometric_altitude(sensor="TMP_turb")
        self.l1_df["z_b"]=self.BELUGA_cls.l1_df["z_b"]
        del self.l1_df["zb"]
    def flag_for_quality(self):
        self.l1_df["quality_flag"]=0
        self.l1_df["quality_flag"].loc[self.l1_df["U"]==np.nan]=3
    
    #%% L1 to L2 processing
    def L1_to_L2_processing(self):
        self.logger.info(f"Starting L1 to L2 processing")
        # Get L1 data
        try:
            self.read_l1_data()
            #self.rename_columns()
        except:
            self.read_raw_files()
            self.L0_to_L1_processing()
            self.rename_columns()
        self.l1_df=self.BELUGA_cls.replace_inf(self.l1_df)
        self.update_z_b_for_pressure_tendency()
        if self.plot_processing:
           pass
            #self.plot_highres_windspeeds()
            # self.plot_thermodynamic_profiles()
        self.l2_df=self.l1_df.copy()
        # Flight segmentation
        self.l2_df=self.run_flight_segmentation(segm_alt_var="z_b")
        if self.plot_processing:
            self.plot_flight_segmentation(add_peak_infos=True)
        self.l2_df=self.segments_cf_conform_IDs(self.l2_df)
        
        self.changeunits()
        self.df2netcdf()
        
        self.BELUGA_cls.flight_date=self.flight_date
        self.BELUGA_cls.flight=self.flight
        self.l2_ds=self.BELUGA_cls.add_L2_meta_data(self.l2_ds,probe=self.probe)
        self.logger.info("Meta data (global attributes and variable attributes)"+\
            "added to L2 dataset.")
            
        self.l2_ds=self.BELUGA_cls.add_VRS_georeferences(self.l2_ds)
        self.logger.info("Georeference data of Station Nord added to Dataset")
        self.l2_ds=self.BELUGA_cls.add_ABL_transition_type(self.l2_ds)
        self.logger.info("ABL transition type added to L2 dataset")
        # cf role
        self.l2_ds=self.BELUGA_cls.add_rf_as_id_for_cf_role(self.l2_ds)
        
        # save file
        self.BELUGA_cls.raw_data_path=self.raw_data_path
        self.l2_fname=self.BELUGA_cls.save_L2_data_as_nc_file(self.l2_ds,
            self.version_number,probe=self.probe)
        self.logger.info("temporary L2 data saved as: "+\
                         self.l2_fname)
    def df2netcdf(self):
        self.l2_ds=xr.Dataset.from_dataframe(self.l2_df)
        self.logger.info("Transformed dataframe to xr.Dataset")
    
    def changeunits(self):
        if self.l2_df["T"].mean()<200:
            self.l2_df["T"].loc[self.l2_df["T"]!=-999] += 273.15
            self.logger.info(
                "Changed units for temperature to Kelvin")

    def run_flight_segmentation(self,segm_alt_var="z_b"):
        # Use the automatised modules from the BELUGA cls, choosing the 
        # altitude variable from instrument probe that has to be considered
        # for flight segmentation
        self.seg_alt_var                 = segm_alt_var
        self.BELUGA_cls.flight           = self.flight
        # Adapt the height change rate for specific flights
        rate_threshold=.6
        if self.flight in ["RF03","RF04","RF06","RF08","RF10","RF13",
                           "RF14", "RF15","RF16","RF17","RF18",
                           "RF23","RF24","RF26"]:
            rate_threshold=.3
        segmented_df,self.profile_peaks  = \
            self.BELUGA_cls.segment_flight_sections(self.l2_df,
                rate_threshold=rate_threshold,
                alt_var=segm_alt_var,use_turb=True)
        self.logger.info(
            "Data categorised in flight segments based on altitude")
        
        return segmented_df

    #%% Plot processing
    def plot_L2_processed_data(self):
        pass

    def plot_flight_segmentation(self,add_peak_infos=False, probe="TMP_turb"):
        
        rf=self.flight
        df_initial = self.l2_df.copy()
        df         = df_initial.copy()
        df["z_b"][df["z_b"]==np.Inf]=np.nan
        # Do not separate between profile numbers anymore
        df["segments"][df_initial["segments"].str.contains(
            'ascent', case=False, na=False)] = 'ascent'
        df["segments"][df_initial["segments"].str.contains(
            'descent', case=False, na=False)] = 'descent'
        #---------------------------------------------------------------------#
        matplotlib.rcParams.update({"font.size":28})
        #% BELUGA Plotting
        fig=plt.figure(figsize=(16,9))
        
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
        ax1.scatter(self.profile_peaks.index,
            self.profile_peaks, marker="o", s=100,
            color="red",edgecolor="k",zorder=3)
        
        max_max=df["z_b"].dropna().max()
        ax1.scatter(df["z_b"].idxmax(),max_max,
                    marker="o",s=200,color="darkred",
                    edgecolor="k",zorder=4)
        ax1.axhline(y=200, color="darkgrey", lw=3,ls="--")
        
        # Add profile peak infos
        if add_peak_infos:
            ax1.text(df["z_b"].idxmax(),max_max+20,
                 str(int(max_max))+" m")
            ax1.text(0.1,0.9,self.flight+": Number of profiles: "+\
                     str(2*len(self.profile_peaks)),
                 transform=ax1.transAxes)
        
        ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
        ax1.set_ylabel("Barometric height (m)")
        
        if max_max>600:
            ylim_max=900
        else:
            ylim_max=600
        
        ax1.set_ylim([0,ylim_max])
        ax1.legend(fontsize=24)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        
        sns.despine(offset=10)
        for axis in ['bottom','left']:
            ax1.spines[axis].set_linewidth(2)
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)

        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        plot_path=self.process_plot_path
        fig_name="BELUGA_"+probe+"_flight_segmentation_"+rf
        file_end=".png"
        fig_name+=file_end
        fname=plot_path+fig_name
        fig.savefig(fname, bbox_inches="tight")
        print("Figure saved as:", fname)
        self.logger.info("Figure saved as:"
                         +plot_path+fname)


#%% Broadband probe
class Broadband_Probe(BELUGA):
    # The code for processing the Broadband Probe is mainly created by 
    # Henning Dorff, Fan Wu and Joshua Müller using knowledge and methods from
    # e.g. Lonardi et al. (2022) & Pilz et al. (2023).
    
    #%% Basics
    def __init__(self, BELUGA,rf="RF12",datatypes=["ADC","ADS","BME","DOF"],
                 run_L1_processing=True,run_L2_processing=True,
                 plot_processing=True):
        
        
        self.BELUGA_cls         = BELUGA
        self.probe              = "BP"
        self.station_name       = BELUGA.station_name        
        self.coordinates        = BELUGA.coordinates
        self.main_path          = BELUGA.main_path
        self.flight_dates       = BELUGA.flight_dates
        self.raw_data_path      = self.main_path +\
           "/BELUGA_broadband_probe/raw_data/"
        self.flight_infos       = BELUGA.flight_infos
        self.flight             = rf
        self.flight_date        = self.flight_dates[self.flight]
        self.datatypes          = datatypes
        self.plot_path          = self.main_path+"/plotting/"
        self.process_plot_path  = self.plot_path+"/test_processing/"+\
            self.flight+"/"
        os.makedirs(self.process_plot_path,exist_ok=True)
        
        # Switches
        self.run_L1_processing  = run_L1_processing
        self.run_L2_processing  = run_L2_processing
        self.plot_processing    = plot_processing
        
        # Configure logging once at the module level
        self.log_file_name='BELUGA_BP_processing_'+self.flight+'.log'
        self.logger=logging_processing.setup_logging(self.log_file_name)
        # Create a new file handler in write mode
        # Constants
        self.C_P                = 1005.0  # J/(kg·K)
        self.DENSITY_AIR        = 1.2922  # kg/m³
        self.SIGMA              = 5.67e-8  # W/m²K⁴
        self.R0                 = 40000  # Ohm
        self.GAIN               = 1
        # Calibration coefficients / Sensitivites for up-downward instruments
        # + 2 values for solar irradiance (redundant here)
        self.sensBPL            = np.array([111,111,10.17,10.57])
        self.sens               = self.sensBPL*10**(-6) # in V / (W/m²)
        self.alpha              = 1.0295e-3
        self.beta               = 2.391e-4
        self.gamma              = 1.568e-7
    
       
    
    def get_all_files_dict(self): # create dictionary of filepaths
        import re    
        file_dates_dict = {} # create empty dictionary
        for fname in os.listdir(self.raw_data_path): # loc relevant filenames
            if 'logfile' in fname:
                continue
            if 'GPS' in fname:
                continue
            if 'am2315' in fname:
                continue
            # search digits in filenames to locate days
            match1 = re.search(r'(\d{8})', fname) 
            match2 = re.search(r'(\d{6})', fname)
            if match1:
                day_key = match1.group(1) # define day_key within filename
                # join datapath with filenames
                file_dates_dict.setdefault(day_key, []).append(
                    os.path.join(self.raw_data_path, fname))
                
            elif match2:
                day_key = '20' + match2.group(1)
                file_dates_dict.setdefault(day_key, []).append(
                    os.path.join(self.raw_data_path, fname))
            self.file_dates_dict=file_dates_dict

    def check_for_relevant_start_time(self,files,ftype="ADC"):
        # get start time of research flight
        start_time=self.BELUGA_cls.flight_infos.loc[self.flight]["Start Time"]
        file_times=[]
        if ftype=="ADC":
            for file in files:
                split_file=file.split("\\")
                file_times.append(
                pd.to_datetime(
                    "20"+split_file[-1][0:2]+"-"+split_file[-1][2:4]+"-"+\
                        split_file[-1][4:6]+" "+\
                        split_file[-1][7:9]+":"+split_file[-1][9:11]))
        elif ftype=="ADS":
            for file in files:
                split_file=file.split("\\")
                file_times.append(
                    pd.to_datetime(
                    split_file[-1][10:14]+"-"+split_file[-1][14:16]+"-"+\
                        split_file[-1][16:18]+" "+\
                        split_file[-1][19:21]+":"+split_file[-1][21:23]))
        elif ftype=="BME":
            #BME_3_20240325_102043
            for file in files:
                split_file=file.split("\\")
                file_times.append(
                pd.to_datetime(
                    split_file[-1][6:10]+"-"+split_file[-1][10:12]+"-"+\
                    split_file[-1][12:14]+" "+split_file[-1][15:17]+\
                        ":"+split_file[-1][17:19]))
        elif ftype=="DOF":
            for file in files:
                split_file=file.split("\\")
                file_times.append(
                pd.to_datetime(split_file[-1][8:12]+"-"+split_file[-1][12:14]+\
                "-"+split_file[-1][14:16]+" "+split_file[-1][17:19]+":"+\
                    split_file[-1][19:21]))
                    
        timedelta_series=[abs(pd.Timedelta(file_time-start_time)) \
                          for file_time in file_times]
        min_t=np.argmin(timedelta_series)
        file=files[min_t]
        return file
    
    def write_out_temporary_files(self,fileformat=".csv"):
        """
        This function writes out the data that is constrained on
        the selection of variables. The output file format is either csv or nc.

        Returns
        -------
        None.

        """
        output_path=self.raw_data_path+"/../temporary_data/"
        os.makedirs(output_path,exist_ok=True)
        for ftype in self.datatypes:
            data_cls_variable = f"{ftype}_df"
            data_df = getattr(self, data_cls_variable)
            #data_df.index=pd.DatetimeIndex(data_df["DATETIME"]) maybe unneeded
            output_file_name="BP_"+ftype+"_L0_var_selection_"+self.flight+"_"+\
                self.flight_date
            
            if fileformat=="nc":
                # convert to xr, save as nc
                xr.Dataset.from_dataframe(data_df).to_netcdf(
                    path=output_path+output_file_name,
                    mode="w",engine="NETCDF4") 
            else:
                data_df.to_csv(path_or_buf=output_path+output_file_name)
    
    def resample_and_interpolate_data(self):
         for ftype in self.datatypes:
             data_cls_variable = f"{ftype.lower()}_df"
                 
             
             data_df = getattr(self, data_cls_variable)
             
             #----------------------------------------------------------------#
             # Old version to be checked
             ## Resample to 1-second intervals with mean values, but only if 
             #  more than three values are there
             def custom_mean(x):
                 # Only compute mean if more than 3 non-NaN values are present
                 if x.dropna().size > 3:
                     return x.mean()
                 else:
                     return pd.NA
             for col in data_df.columns:
                 pd.to_numeric(data_df[col], errors='coerce')

             df_resampled=data_df.resample("1s").agg(custom_mean)
             
             # Create a mask of original non-NaN points
             original_mask = df_resampled.notna()
             # Step 3: Interpolate with a limit of 5 seconds
             
             try:
                 df_resampled = df_resampled.apply(
                     lambda col: pd.to_numeric(col, errors='coerce') \
                         if col.dtypes == 'object' else col)

                 df_resampled = df_resampled.interpolate(
                     method='time', limit=5)
                 self.logger.info(
                     f"{ftype} data was gap-filled by interpolation.")
             except:
                 self.logger.error(f"{ftype} data could not be gap-filled.")
             #----------------------------------------------------------------#
             # Step 4: Create 'interpolated' flag
             # Flag points that are now non-NaN but were originally NaN
             try:
                 df_resampled['interpolation_flag'] = \
                 (~original_mask[original_mask.columns[-1]]) & \
                     df_resampled[df_resampled.columns[-1]].notna()
                 self.logger.info(f"1Hz {ftype} data with interpolation flag.")
             except:
                 self.logger.error(f"no interpolation flag could be provided.") 
             #df_resampled=df_resampled.fillna(-999)
             #self.logger.info("Filled nans with -999 as missing value")
             setattr(self, data_cls_variable, df_resampled)
    
    #def merge_and_interpolate_netcdfs(day_key, file_paths): 
        # create nc files merged by day
    #    datasets = [xr.open_dataset(fp).sortby('time') for fp in file_paths] 
    # open ftype_day files and sort by time
    #    resampled_datasets = [ds.resample(time='1s').interpolate('linear') \
    #    for ds in datasets] 
    # resamples to common time grid and fills in missing values
    
    #    merged_ds = xr.merge(resampled_datasets) # merge files by day
    #    merged_ds.to_netcdf(f'D:/Documents/arcticamplification/beluga/all_broadband/all/merged_{day_key}.nc')
    #    # no resampling - only T shows up
    
    def get_files(self,ftype="ADC"):
        path=self.raw_data_path
        if      ftype =="ADC":
            fnames=self.flight_date[2:]+"*_Data_ADC.txt"
        elif ftype=="ADS":
            fnames="ads1115_3_"+self.flight_date+"*.txt"
        elif    ftype=="BME":    
            fnames="BME_3_"+self.flight_date+"*.txt"
        elif ftype=="DOF":
            fnames="dof10_3_"+self.flight_date+"*.txt"
        else:
            raise Exception("Wrong file type (ftype) chosen.")
        fpath=path+fnames
        files=glob.glob(fpath)
        self.raw_file=self.check_for_relevant_start_time(files,ftype=ftype)
        
    #%% --- Raw Data (Level 0) read Functions ---
    def read_adc(self): # radiation in Vs to correct for instrument T
        ftype="ADC"    
        columns = ['Not_timest', 'hh', 'mm', 'sssss', 'adc_Vch1V',
                   'adc_vswd', 'adc_h3V', 'adc_vswu', 'VarName9',
                   'adc_vlwd', 'VarName11', 'adc_vlwu', 'VarName13']
        self.get_files(ftype=ftype)
        try:
            self.adc_df = pd.read_csv(self.raw_file,
                delim_whitespace=True, names=columns,
                dtype=str, encoding='latin1', skiprows=1)
        except FileNotFoundError as e:
            self.logger.error(f"File not found for {self.flight_date}: {e}")
            FileNotFoundError("No ADC files for given day found:"+\
                              self.flight_date)
        
        self.adc_df[columns] = self.adc_df[columns].apply(
                                    pd.to_numeric, errors='coerce')
        self.adc_df['time'] = pd.to_datetime(self.adc_df['Not_timest'],
                                             unit='s')
        self.adc_df.set_index('time', inplace=True)
        self.adc_df.drop(columns=['Not_timest', 'hh', 'mm', 'sssss', 
                    'VarName9', 'VarName11', 'VarName13'], inplace=True)
        self.adc_df.name = 'adc_'+self.flight
        self.logger.info("ADC data was read")
        
    def read_ads(self): # upward and downward radiation in W m-2
        ftype="ADS"    
        self.get_files(ftype=ftype)
        columns = ['yyyy', 'month', 'dd', 'hh', 'mm', 'sssss', 'ads_VarName7',
                   'ads_digital_u_u', 'ads_digital_u_d', 'ads_digital_u_ref',
                   'ads_conversion_factor'] 
        #, 'a', 'b', 'c', 'd', 'e', 'ads_upward', 'ads_downward']
        #, (58, 68), (68, 78), (78, 88), (88, 97),
        #(97, 106), (106, 114), (114, 122)]
        
        colspecs = [(0, 4), (4, 7), (7, 10), (10, 13), (13, 16), (16, 23),
                    (23, 29), (29, 35), (35, 41), (41, 47), (47, 58)] 
        try:
            self.ads_df= pd.read_fwf(self.raw_file,
                colspecs=colspecs, names=columns, dtype=str,
                encoding='latin1', skiprows=1)
        except FileNotFoundError as e:
            self.logger.error(f"File not found for {self.flight_date}: {e}")
            FileNotFoundError("No ADS files for given day found:"+\
                              self.flight_date)
        if not self.flight=="RF02":    
            self.ads_df[columns] = self.ads_df[columns].apply(
                            pd.to_numeric, errors='coerce')
        date_df_str=self.ads_df[['yyyy', 'month', 'dd', 'hh', 'mm']].astype(str)
        agg_date   = date_df_str.agg(" ".join,axis=1)
        self.ads_df['time'] = pd.to_datetime(
            agg_date,format='%Y %m %d %H %M',errors='coerce') +\
            pd.to_timedelta(self.ads_df['sssss'].fillna(0).astype(float), unit='s')
        self.ads_df.set_index('time', inplace=True)
        self.ads_df.drop(columns=['yyyy', 'month', 'dd', 'hh', 'mm', 'sssss'],
                         inplace=True)
        if self.flight=="RF02":
            columns_new = ['ads_VarName7','ads_digital_u_u', 'ads_digital_u_d',
                           'ads_digital_u_ref','ads_conversion_factor'] 
            
            self.ads_df[columns_new] = self.ads_df[columns_new].apply(
                            pd.to_numeric, errors='coerce')    
        self.ads_df.name = 'ads'+self.flight
        self.logger.info("ADS data was read")
        
    
    def read_bme(self): # thermodynamic parameters
        ftype="BME"
        self.get_files(ftype)
        
        columns = ['yyyy', 'month', 'dd', 'hh', 'mm', 'sssss', 
                   'bme_temperature', 'bme_rh', 'bme_pressure']
        try:
            self.bme_df = pd.read_csv(self.raw_file, delim_whitespace=True,
                names=columns, dtype=str, encoding='latin1', skiprows=1)
        
        except FileNotFoundError as e:
            self.logger.error(f"File not found for {self.flight_date}: {e}")
            FileNotFoundError("No BME files for given day found:"+\
                                  self.flight_date)
            
        self.bme_df[columns] = self.bme_df[columns].apply(
            pd.to_numeric, errors='coerce')
        self.bme_df['time'] = pd.to_datetime(
            self.bme_df[['yyyy', 'month', 'dd', 'hh', 'mm']].astype(str).agg(
                ' '.join, axis=1), format='%Y %m %d %H %M', errors='coerce') +\
            pd.to_timedelta(self.bme_df['sssss'].fillna(0) , unit='s')
        self.bme_df.set_index('time', inplace=True)
        
        self.bme_df.drop(columns=['yyyy', 'month', 'dd', 'hh', 'mm', 'sssss'],
                         inplace=True)
        self.logger.info(f"BME data was read")
        
        
    def read_dof(self): # altitude sensor for correcting tilt and for pressure
        ftype="DOF"    
        columns = ['yyyy', 'month', 'dd', 'hh', 'mm','sssss',
                   'dof_roll', 'dof_pitch', 'dof_yaw', 'dof_pressure',
                   'dof_temperature', 'dof_height']
        self.get_files(ftype)
        try:
            self.dof_df = pd.read_csv(self.raw_file, delim_whitespace=True,
                names=columns, dtype=str, encoding='latin1', skiprows=1)
        except FileNotFoundError as e:
            self.logger.error(f"File not found for {self.flight_date}: {e}")
            FileNotFoundError("No DOF files for given day found:"+\
                              self.flight_date)
        
        if not self.flight=="RF20":    
            self.dof_df[columns] = self.dof_df[columns].apply(
            pd.to_numeric, errors='coerce')
        date_df_str=self.dof_df[['yyyy', 'month', 'dd', 'hh', 'mm']].astype(str)
        agg_date   = date_df_str.agg(" ".join,axis=1)
        
        self.dof_df['time'] = pd.to_datetime(
            agg_date, format='%Y %m %d %H %M', errors='coerce') +\
            pd.to_timedelta(self.dof_df['sssss'].fillna(0).astype(float), unit='s')
        self.dof_df.set_index('time', inplace=True)
        self.dof_df.drop(columns=['yyyy', 'month', 'dd', 'hh', 'mm', 'sssss'],
                         inplace=True)
        if self.flight=="RF20":
            columns_new = ['dof_roll', 'dof_pitch', 'dof_yaw',
                           'dof_pressure','dof_temperature', 'dof_height']
            self.dof_df[columns_new] = self.dof_df[columns_new].apply(
                            pd.to_numeric, errors='coerce')
        self.logger.info(f"DOF data was read")
        
    
    def read_raw_files(self):
        self.logger.info(f"Read raw files (Level 0) for {self.flight_date}")
        read_funcs={"ADC":self.read_adc,
                    "ADS":self.read_ads,
                    "BME":self.read_bme,
                    "DOF":self.read_dof}
        [read_funcs[data_type]() for data_type in self.datatypes]
        
    #%% Run processing
    def run_processing(self,version_number="_v2.1",rf="RF12",
            datatypes=["ADC","ADS","BME","DOF"]):
        #with open(self.log_file_name,"w") as log_file:
        #    log_file.write(
        self.version_number=version_number
        self.logger.info(
            f"Processing flight {self.flight} on date {self.flight_date}\n")
        if self.run_L1_processing:
            self.read_raw_files()
            self.L0_to_L1_processing(rf)
            self.merge_L1_data()
            self.flag_for_quality()
        if self.run_L2_processing:
            self.L1_to_L2_processing()
            
    #%%L0 to L1 processing    
    def L0_to_L1_processing(self,rf):
        """
        This routine takes the raw files and processes them to the relevant
        irradiances by clipping irrelevant variables.
        

        Returns
        -------
        None.

        """
        self.is_temperature_corrected = False
        self.is_inertia_corrected     = False
        # Logging-------------------------------------------------------------#
        #self.logger.info(f"Starting L0 to L1 processing")
        #---------------------------------------------------------------------#
        self.flight=rf
        self.flight_date=self.flight_dates[self.flight]
        self.logger.info(f"Starting L0 to L1 processing")
        self.resample_and_interpolate_data()
        try:
            self.correct_irradiances_for_temp_and_compare()
            self.is_temperature_corrected=True
        except:
            self.logger.error("Error in temperature correction")
        try:
            self.correct_irradiances_for_inertia_and_compare()
            self.is_inertia_corrected=True
        except:
            self.logger.error("No inertia correction applied")
            self.is_inertia_corrected=False
        self.logger.info("Completed L0 to L1 processing successfully.")
        
    
    def u_to_resistence(self,u_th): # Vs to resistance
        resistance= self.R0*u_th/(self.u_ref - u_th)
        return resistance
    def resistence_to_temp(self,resist_to_use): # resistance to instrument T
        radiometer_temp=1/(self.alpha+self.beta*np.log(resist_to_use)+\
                                self.gamma * np.log(resist_to_use) ** 3)
        return radiometer_temp
    
    def u_to_irr(self,u_to_use,sensitivity_to_use,
                 spectral="terrestial"): # for u_to_terr_irr
        irradiance=u_to_use/(sensitivity_to_use * self.GAIN)
        return irradiance
    
    def u_to_terr_irr(self,u_to_use,sensitivity_to_use,
                      temp=""):
        if isinstance(temp, pd.Series):
            if temp.shape[0]==0:
                temp=self.radiometer_temp
        elif temp=="":
             temp=self.radiometer_temp
        else:
            pass
        # corrects irradiance for instrument T, for tilt not needed
        terr_irradiance=self.u_to_irr(u_to_use,sensitivity_to_use)+\
            self.SIGMA*(temp**4)
        return terr_irradiance
    
    def calc_pure_irradiances(self):
        # get pure uncorrected irradiances
        self.L0_lw_d=self.u_to_terr_irr(self.adc_df["adc_vlwd"],
                                        self.sens[2],0)
        self.L0_lw_u=self.u_to_terr_irr(self.adc_df["adc_vlwu"],
                                        self.sens[3],0)
        
        self.logger.info(
            f"Calculated raw terrestrial up- and downward irradiances")
        
    def correct_irradiances_for_temp(self):
        self.logger.info("Convert ADS digital counts to voltages")
        
        u_upw = self.ads_df['ads_digital_u_u'] *\
            self.ads_df['ads_conversion_factor'] # upward Vs
        u_dwn = self.ads_df['ads_digital_u_d'] * \
            self.ads_df['ads_conversion_factor'] # downward Vs
        
        self.u_ref = self.ads_df['ads_digital_u_ref'] *\
            self.ads_df['ads_conversion_factor'] # reference Vs
            
        self.logger.info("Calculate radiometer body temperature from ADS")
        resistence_upw=self.u_to_resistence(u_upw)
        resistence_dwn=self.u_to_resistence(u_dwn)
        temp_lw_u = self.resistence_to_temp(resistence_upw) # upward T
        temp_lw_d = self.resistence_to_temp(resistence_dwn) # downward T
        
        #Terrestrial irradiance from ADC voltages with temp correction (body T)
        self.lw_d_tcorr = self.u_to_terr_irr(self.adc_df['adc_vlwd'], 
            self.sens[2], temp_lw_d) # Dwnward terrestrial irr
        self.lw_u_tcorr = self.u_to_terr_irr(self.adc_df['adc_vlwu'], 
            self.sens[3], temp_lw_u) # Upward terrestrial irr
        
        self.logger.info(
            "Corrected terrestrial irradiances for radiometer temperature")

    def correct_irradiances_for_temp_and_compare(self):
        self.logger.info("Correct radiometer values for body temperature")
        # Terrestrial irradiance from ADS file without correction
        if not hasattr(self,"L0_lw_d"):
            self.calc_pure_irradiances()
        
        self.correct_irradiances_for_temp()
        if self.plot_processing:
            self.plot_temperature_correction_comparison()
            
    def correct_irradiances_for_interia(self,add_simple_deconvolution=True):
        from scipy.fft import fft, ifft
        # Given parameters
        self.pyrgeo_tau = 6 #3.6# findout the right tau for ballon
        self.fcut       = .2 # cut-off frequency
        self.window     = 3 #1 # smoothing window
        self.xdt        = 1
        
        self.rolling_window="7s"

        # check inertia times
        # function for inertia correction (see Ehrlich and Wendisch, 2015)                
        def decon_rt_fast(inp_arr,tau,fcut,window, xdt):
            
            uneven=0
            #time_arr=inp_arr.index
            inp_arr.index=pd.DatetimeIndex(inp_arr.index)
            time_arr=inp_arr.index.hour*3600 +\
                inp_arr.index.minute*60 + inp_arr.index.second + \
                    inp_arr.index.microsecond/1e6
            if (len(time_arr) % 2 == 1):
                uneven=1
        
                last_inp = inp_arr[-1]
                last_time = inp_arr.index[-1]
        
                inp_arr = inp_arr[0:-1]
                time_arr = inp_arr.index
        
            xnt = len(time_arr)
        
            # Convolution Function
            xfcon_all = 1./tau*np.exp(-(np.arange(xnt)*xdt)/tau)
            xt99 = tau*np.log(1./0.00001)
            ximax = int(xt99/xdt)
            if (ximax % 2 == 1):
                ximax=ximax+1
                
            inp = np.empty(len(inp_arr))*0.
        
            for i,val in enumerate(time_arr):
                #print(i)
                stopID=0
        
                ia = 0+i*ximax*6
                ie = ximax*12 + (i+1)*ximax*6
                if (ie >= xnt):
                    ie=xnt
                    stopID=1
        
                p_time = time_arr[ia:ie]
                p_inp = inp_arr[ia:ie]
        
                pnt = len(p_time)
        
                if (pnt % 2 == 1):
                    p_inp = p_inp[1:]
                    p_time = p_time[1:]
                    pnt=pnt-1
        
                pfcon_all = 1./tau*np.exp(-(np.arange(pnt)*xdt)/tau)
                pt99 = tau*np.log(1./0.00001)
                pimax = int(pt99/xdt)
        
                if (pimax % 2 == 1):
                    pimax=pimax+1
                    
                zero = np.zeros(pimax)           #np.empty(pimax)*0.
                fcon = np.hstack((pfcon_all,zero))
        
                zero = np.zeros(pimax)+ p_inp[0] #np.empty(pimax)*0.+ p_inp[0]
                datacon = np.hstack((p_inp,zero))
        
                ftdatacon = fft(datacon)
                ftfcon = fft(fcon)
                ftfcon = ftfcon/ftfcon[0]
        
                pntzero = pnt+pimax #pnt #pnt+pimax 
                freq = np.arange(pntzero/2.)/xdt/pntzero
        
                ftrm = window/window*np.sin(np.pi*window*freq)/\
                    (np.pi*window*freq)
                ftrm[0]=1
                ftrm = np.hstack((ftrm,np.flip(ftrm)))
                
                ftdata0 = ftdatacon/ftfcon
                ftdata1 = ftdata0
        
                ifcut = np.argwhere(freq > fcut)
                ifcut = ifcut[0][0]
                ftdata1[ifcut:-ifcut] = complex(0,0)
                ftdata = ftdata1*ftrm
        
                datazero = np.real(ifft(ftdata))
                pdata = datazero[0:pnt]
        
        
        
                if (i != 0):
                    ibb = 3*ximax
                if (i == 0):
                    ibb = 0
                icc = 9*ximax
                ib = ia+ibb
                ic = ia+icc
        
                if (stopID != 1):
                    inp[ib:ic] = pdata[ibb:icc]
                if (stopID == 1):
                    try:
                        inp[ib:] = pdata[ibb:]
                    except:
                        max_idx=len(inp)+ib
                        inp[ib:] = pdata[ibb:max_idx]
                    break
        
        
            if (uneven == 1):
                inp_arr = np.hstack((inp_arr,last_inp))
                inp = np.hstack((inp,last_inp))
                time_arr = np.hstack((time_arr,last_time))
        
            return inp


        # Fill nans for FFT
        self.lw_d_tcorr = self.lw_d_tcorr.interpolate().ffill().bfill()
        self.lw_u_tcorr = self.lw_u_tcorr.interpolate().ffill().bfill()
        
        try:
            self.lw_u_L1 = decon_rt_fast(self.lw_u_tcorr,self.pyrgeo_tau,
                                  self.fcut,self.window,self.xdt)
            self.logger.info("Deconvoluted upward irradiance ")
        except:
            self.logger.error(
                "Error in deconvolution of upward irradiance (not performed)")
        try:
            self.lw_d_L1 = decon_rt_fast(self.lw_d_tcorr,self.pyrgeo_tau,
                                         self.fcut,self.window,self.xdt)
            self.logger.info("Deconvoluted downward irradiance")
        except:
            self.logger.error("Error in deconvolution of downward irradiance (not performed)")
        
        # Convert L1 Irradiances to pd.Series
        try:
            self.lw_u_L1=pd.Series(data=self.lw_u_L1,
                                   index=self.lw_u_tcorr.index)
            self.lw_d_L1=pd.Series(data=self.lw_d_L1,
                                   index=self.lw_d_tcorr.index)
        except:
            self.lw_u_L1=pd.Series(data=self.lw_u_tcorr,
                                   index=self.lw_u_tcorr.index)
            self.lw_d_L1=pd.Series(data=self.lw_d_tcorr,
                                   index=self.lw_d_tcorr.index)
        
        # Add 3s running mean after deconvolution to correct for
        # Gibbs Theorem as suggested by Pilz et al. (2023)
        
        self.lw_u_L1 = self.lw_u_L1.rolling(self.rolling_window).mean()
        self.lw_d_L1 = self.lw_d_L1.rolling(self.rolling_window).mean()
        
        
        if add_simple_deconvolution:
            self.simple_deconvolution()
            if self.plot_processing:
                if self.flight=="RF12":
                    self.quicklook_inertia_correction_profile()
                else: 
                    pass
    def simple_deconvolution(self,add_running_mean=False):
        self.lw_u_dt=self.lw_u_tcorr.diff().dropna()
        self.lw_d_dt=self.lw_d_tcorr.diff().dropna()
        self.lw_u_dt.index=self.lw_u_tcorr.index[0:-1]
        self.lw_d_dt.index=self.lw_d_tcorr.index[0:-1]
        self.lw_u_simple=self.lw_u_tcorr[0:-1]+self.pyrgeo_tau*self.lw_u_dt
        self.lw_d_simple=self.lw_d_tcorr[0:-1]+self.pyrgeo_tau*self.lw_d_dt
        if add_running_mean:
            self.lw_u_simple=self.lw_u_simple.rolling("3s").mean()
            self.lw_d_simple=self.lw_d_simple.rolling("3s").mean()
        
    def merge_L1_data(self):
            
        #Prepare irradiance dataset
        if self.is_inertia_corrected:
            f_up_values   = np.array(self.lw_u_L1)
            f_down_values = np.array(self.lw_d_L1)
                
        else:
            f_up_values   = np.array(self.lw_u_tcorr)
            f_down_values = np.array(self.lw_d_tcorr)
            
        self.irr_df=pd.DataFrame(data=np.nan,columns=["F_up","F_down","F_net"],
                                 index=self.lw_u_tcorr.index)
        self.irr_df["interpolation_flag"]=self.adc_df["interpolation_flag"]
        
        self.irr_df["F_up"]   = f_up_values
        self.irr_df["F_down"] = f_down_values
        self.irr_df["F_net"]  = self.irr_df["F_down"]-self.irr_df["F_up"]
        
        # data dfs
        data_dfs=[self.irr_df,self.bme_df,self.dof_df]
        # Get íntersection index
        # Initialize with the index of the first DataFrame
        intersection_idx = data_dfs[0].index
        # Loop through remaining DataFrames and find their interseciton indices
        for df in data_dfs[1:]:
            intersection_idx = intersection_idx.intersection(df.index)
        print(self.irr_df.shape)
        print(intersection_idx)
        # Allocate extra columns
        self.irr_df["bme_temperature"] = np.nan
        self.irr_df["bme_pressure"]    = np.nan
        self.irr_df["bme_rh"]          = np.nan
        
        self.irr_df["roll"]            = np.nan
        self.irr_df["pitch"]           = np.nan
        self.irr_df["yaw"]             = np.nan
        self.irr_df["dof_pressure"]    = np.nan
        
        self.irr_df["bme_temperature"].loc[intersection_idx]=\
            self.bme_df["bme_temperature"].loc[intersection_idx].values
        self.irr_df["bme_pressure"].loc[intersection_idx]=\
            self.bme_df["bme_pressure"].loc[intersection_idx].values
        self.irr_df["bme_rh"].loc[intersection_idx]=\
            self.bme_df["bme_rh"].loc[intersection_idx].values
        
        self.irr_df["roll"].loc[intersection_idx]=\
            self.dof_df["dof_roll"].loc[intersection_idx].values
        self.irr_df["pitch"].loc[intersection_idx]=\
            self.dof_df["dof_pitch"].loc[intersection_idx].values
        self.irr_df["yaw"].loc[intersection_idx]=\
            self.dof_df["dof_yaw"].loc[intersection_idx].values
        self.irr_df["dof_pressure"].loc[intersection_idx]=\
            self.dof_df["dof_pressure"].loc[intersection_idx].values
            
        irr_csvname=self.raw_data_path+"/../temporary_data/"+\
            "L1_BP_irradiances_"+self.flight+"_"+self.flight_date+".csv"
        self.irr_df.to_csv(path_or_buf=irr_csvname,index=True)
        print("irradiances saved as:",irr_csvname)
        self.logger.info("L1 irradiances temporarily saved as: "+\
                         irr_csvname)
    
    def correct_irradiances_for_inertia_and_compare(self):
        self.logger.info(
            "Correct for inertia, based on Ehrlich & Wendisch (2025)")
        self.correct_irradiances_for_interia()
        if self.plot_processing:
            self.plot_inertia_correction_comparison()
            self.plot_differences_inertia_correction()
    def flag_for_quality(self, ok_pitch=2.5,bad_pitch=5):
        self.irr_df["quality_flag"]=0
        self.irr_df["quality_flag"].loc[\
            self.irr_df["interpolation_flag"]==True]=1
        self.irr_df["quality_flag"].loc[abs(self.irr_df["pitch"])>ok_pitch]=2
        self.irr_df["quality_flag"].loc[abs(self.irr_df["pitch"])>bad_pitch]=3
        self.irr_df["quality_flag"].loc[\
            self.irr_df["interpolation_flag"]==True]=1
        self.irr_df["quality_flag"].loc[self.irr_df["F_net"]==np.nan]=4
        del self.irr_df["interpolation_flag"]    
    
    def read_l1_data(self):
        l1_fname=self.raw_data_path+"/../temporary_data/"+\
            "L1_BP_irradiances_"+self.flight+"_"+self.flight_date+".csv"
        self.l1_df=pd.read_csv(l1_fname,index_col=0)
        self.l1_df.index=pd.DatetimeIndex(self.l1_df.index)
        
    #%% L1 to L2 processing
    def L1_to_L2_processing(self):
        self.logger.info(f"Starting L1 to L2 processing")
        # Get L1 data
        if hasattr(self, "irr_df"):
            self.l1_df=self.irr_df.copy()
            del self.irr_df
        else:
            self.read_l1_data()
            
        self.l1_df=self.BELUGA_cls.replace_inf(self.l1_df)
        # prepare for z_b calculation
        # Get barometric pressure altitude
        self.BELUGA_cls.l1_df               = self.l1_df.copy()
        self.BELUGA_cls.logger              = self.logger
        self.BELUGA_cls.plot_processing     = self.plot_processing
        self.BELUGA_cls.process_plot_path   = self.process_plot_path
        self.BELUGA_cls.flight              = self.flight
        self.BELUGA_cls.flight_date         = self.flight_date
        
        self.BELUGA_cls.get_barometric_altitude()
        self.l1_df["z_b"]=self.BELUGA_cls.l1_df["z_b"]
        
        self.l2_df=self.l1_df.copy()
        if self.plot_processing:
            self.plot_up_and_down_ward_terr_radiation()
            self.plot_net_terr_radiation()
        self.l2_df=self.run_flight_segmentation(segm_alt_var="z_b")
        if self.plot_processing:
            self.plot_flight_segmentation(add_peak_infos=True)
        # CF conform flight segmentation
        self.l2_df=self.segments_cf_conform_IDs(self.l2_df)
        
        # 
        self.df2netcdf()
        self.changeunits()
        self.l2_ds=self.add_L2_meta_data(self.l2_ds,probe=self.probe)
        self.logger.info("Meta data (global attributes and variable attributes)"+\
            "added to L2 dataset.")
        if self.flight=="RF09":
            self.l2_ds.attrs["warning"]="Significant loss of data from this flight."+\
                "Sonde stopped working during first profile. "+\
                    "Micro-SD Card was broken because of unkown reason. Repaired in the evening."

        self.l2_ds=self.BELUGA_cls.add_VRS_georeferences(self.l2_ds)
        self.logger.info("Georeference data of Station Nord added to Dataset")
        self.l2_ds=self.BELUGA_cls.add_ABL_transition_type(self.l2_ds)
        self.logger.info("ABL transition type added to L2 dataset")
        # cf role
        self.l2_ds=self.BELUGA_cls.add_rf_as_id_for_cf_role(self.l2_ds)
        # save file
        self.BELUGA_cls.raw_data_path = self.raw_data_path
        self.BELUGA_cls.flight_date   = self.flight_date
        self.l2_fname=self.BELUGA_cls.save_L2_data_as_nc_file(self.l2_ds,
            self.version_number,probe=self.probe)
        self.logger.info("temporary L2 data saved as: "+\
                         self.l2_fname)
        
    def get_measured_p_sfc(self,sensor="bme",data_df=pd.DataFrame()):
        if data_df.shape[0]==0:
            irr_df=self.l1_df
        
        # find surface pressure by periods where mean gradients
        # are very low
        p=pd.DataFrame(data=irr_df[sensor+"_pressure"].values,
                       columns=["pressure"],index=irr_df.index)
        p["int_idx"]=np.arange(p.shape[0])
        t=irr_df[sensor+"_temperature"]
        # Minutely pressure gradients
        p_grads=p.diff().resample("1min").mean()
        low_p_grads=p_grads[abs(p_grads)<0.5]
        # use this values for surface pressure
        p_sfc_values=p.loc[low_p_grads.index[1:]]
        # 2nd condition
        quantile_value=0.95
        if self.flight=="RF14":
            quantile_value=0.9
        p_90=p["pressure"].quantile() # changed from 0.9
        p_sfc_values=p_sfc_values.loc[p_sfc_values["pressure"]>p_90]
        t_sfc=t.loc[p_sfc_values.index]
        return p,p_sfc_values,t_sfc
    
    def get_barometric_altitude(self,sensor="bme"):
        # Barometric altitude we need the measured surface pressure
        self.logger.info(f"Derive barometric altitude from {sensor}-sensor")
        
        p,p_sfc_meas,t_sfc=self.get_measured_p_sfc(
            sensor=sensor,data_df=pd.DataFrame())
        _,p_tendency=self.calc_sfc_pressure_tendency(p_sfc_meas)
        # apply regression coefficients to timesteps for continous time series
        p_sfc_ts=p_tendency.slope*p["int_idx"]+p_tendency.intercept
        self.logger.info(f"Surface pressure tendency of {p_tendency.slope*36}"+\
                         "hPa/h considered")
        if self.plot_processing:
            self.plot_p_sfc(p,p_sfc_meas, self.l1_df["bme_temperature"], 
                            t_sfc,sensor=sensor)
        # Apply p_sfc series for barometric height calculation
        self.calc_barometric_altitude(p,sensor=sensor,T0=t_sfc.mean(),
                                      p0=p_sfc_ts)
                    
    def calc_barometric_altitude(self,p,
        sensor="bme",T0=288.15,p0=101325): 
        # find surface pressure by periods where mean gradients
        # are very low
        if T0<200:
            T0+=273.15
        R = 287.05 #fixed
        g = 9.81
        self.l1_df["z_b"]=(R * T0) / g * np.log(p0 / p["pressure"])
        self.logger.info("Barometric height added to BP data frame")

    def calc_sfc_pressure_tendency(self,p_sfc_values):
        from scipy.stats import linregress
        x = p_sfc_values["int_idx"]#.m .map(pd.Timestamp.toordinal)
        y = p_sfc_values["pressure"]
        
        p_tendency = linregress(x, y)
        #print(f'Slope: {p_tendency.slope}')
        #print(f'Intercept: {p_tendency.intercept}')
        # Fitted values
        p_fit = p_tendency.slope * x + p_tendency.intercept
        """
        weigthing approach
        # Convert datetime index to ordinal
        x = ts.index.map(pd.Timestamp.toordinal)
        y = ts.values

        # Create weights emphasizing more recent data
        # For example, linearly increasing weights towards the end
        weights = np.linspace(1, 2, len(ts))  # weights from 1 to 2

        # Perform weighted linear fit
        coefficients = np.polyfit(x, y, 1, w=weights)
        slope, intercept = coefficients

        # Calculate fitted values
        y_fit = np.polyval(coefficients, x)
        """
        return p_fit,p_tendency
    
    def df2netcdf(self):
        self.l2_df.rename(columns={"F_up": "F_upw_terr", "F_down": "F_dnw_terr"})
        self.l2_ds=xr.Dataset.from_dataframe(self.l2_df)
        #self.l2_ds.assign_coords({"time":})
        self.logger.info("Transformed dataframe to xr.Dataset")
    
    def changeunits(self):
        self.l2_ds["bme_temperature"] += 273.15
        self.l2_ds["bme_pressure"]    /= 100
        self.l2_ds["dof_pressure"]    /= 100
        self.logger.info(
            "Changed units for temperature (to K) and pressure (to hPa)")        
    
    def run_flight_segmentation(self,segm_alt_var="z_b"):
        # Use the automatised modules from the BELUGA cls, choosing the 
        # altitude variable from instrument probe that has to be considered
        # for flight segmentation
        self.seg_alt_var                 = segm_alt_var
        self.BELUGA_cls.flight           = self.flight
        rate_threshold=.6
        if self.flight in ["RF03","RF04","RF06","RF08","RF10","RF13",
                           "RF14", "RF15","RF16","RF17","RF18",
                           "RF23","RF24","RF26"]:
            rate_threshold=.3
        elif self.flight=="RF11":
            rate_threshold=.1
        segmented_df,self.profile_peaks  = \
            self.BELUGA_cls.segment_flight_sections(self.l2_df,
                rate_threshold=rate_threshold, alt_var=segm_alt_var)
        self.logger.info(
            "Data categorised in flight segments based on altitude")
        
        return segmented_df
    
    
    
    
    
    # def add_meta_data(self,ds):
    #     import nc_attributes
    #     # Read dictionaries
    #     global_att_dict = nc_attributes.global_attributes
    #     bp_att_dict     = nc_attributes.BP_var_attributes 
    #     # Add Global Attributes from dictionary
    #     for att in [*global_att_dict.keys()]:
    #         self.l2_ds.attrs[att]=global_att_dict[att]
        
    #     # Add variable attributes from dictionary
    #     #for coord in [*self.l2_ds.coords]:
    #     #    coord_attrs=bp_att_dict[coord]
    #     #    for att in [*coord_attrs.keys()]:
    #     #        self.l2_ds[coord].attrs[att]=coord_attrs[att]
        
    #     for var in [*self.l2_ds.keys()]:
    #         var_attrs=bp_att_dict[var]
    #         for att in [*var_attrs.keys()]:
    #             self.l2_ds[var].attrs[att]=var_attrs[att]
    #     self.update_rf_specfic_attributes()
    #     self.logger.info(
    #         "Meta data (global attributes and variable attributes) \
    #             added to L2 dataset.")
                
    
    
    # --- Save Results ---
    #result_ds = xr.Dataset({
    #    'datetime': ds['time'],
    #    'lw_down': lw_d,
    #    'lw_up': lw_u,
    #    'lw_net': ( lw_net),
    #    'bme_height' : (bme_height),
    #    'dof_height' : (dof_height),
    #    'temperature' : (temperature),
    #    'rh' : (rh),
    #    'pressure' : (bme_pressure)
    #}, coords={'time': ds['time']})

    #%% Plot processing
    def plot_temperature_correction_comparison(self):
        temp_corr_fig=plt.figure(figsize=(10, 6))
        ax1=temp_corr_fig.add_subplot(111)
        matplotlib.rcParams.update({"font.size":18})
        fig_name=self.flight+"_"+"BP_L1_processing_Temp_correction.png"
        # L0 data (uncorrected for radiometer temperature)
        ax1.plot(self.L0_lw_u, c='lightblue',lw=2, label='F_up, no T corr')
        ax1.plot(self.L0_lw_d, c='orange', lw=2, label='F_down, no T corr')
        # L1 data (correct for temperature)
        ax1.plot(self.lw_u_tcorr, c='darkblue', label='F_up, T corr')
        ax1.plot(self.lw_d_tcorr, c='darkred', label='F_down, T corr')
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        #ax1.text(0.1,0.9,self.flight,transform=ax1.transAxes)
        ax1.set_xlabel('UTC Time of '+ self.flight_date +" ("+self.flight+")")
        ax1.set_ylabel('Irradiance ($\mathrm{W\,m}^{-2}$)')
        ax1.grid(True)
        ax1.legend(fontsize=14)
        sns.despine(offset=10)
        temp_corr_fig.savefig(self.process_plot_path+fig_name,
                    dpi=300,bbox_inches="tight")
        self.logger.info("Temperature correction comparison plotted under: "+\
                         self.process_plot_path+fig_name)
    
    def plot_inertia_correction_comparison(self,add_simple_deconvolution=True):
        iner_corr_fig=plt.figure(figsize=(12,9))
        matplotlib.rcParams.update({"font.size":20})
        fig_name=self.flight+"_"+"BP_L1_processing_Inertia_correction.png"
        ax1=iner_corr_fig.add_subplot(211)
        ax2=iner_corr_fig.add_subplot(212)
        # Irradiance F_upward 
        ax1.plot(self.lw_u_L1,c="darkblue",lw=3,label="inertia corrected")
        ax1.plot(self.lw_u_tcorr,c="lightblue",lw=2,label="uncorrected",
                 zorder=3)
        
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.set_xlabel("")
        ax1.set_ylabel("$F_{\mathrm{up}} (\mathrm{W\,m}^{-2})$")
        ax1.grid(True)
        # Irradiance F_downward
        ax2.plot(self.lw_d_L1,c="darkred",lw=3,label="inertia corrected")
        ax2.plot(self.lw_d_tcorr,c="salmon",lw=1,label="uncorrected",zorder=3)
        if add_simple_deconvolution:
            ax1.plot(self.lw_u_simple,c="dimgrey",lw=1,ls="--",
                     label="simple inertia corr",zorder=5)
            ax2.plot(self.lw_d_simple,c="dimgrey",lw=1,ls="--",
                     label="simple inertia corr",zorder=5)
        
        ax2.set_ylabel("$F_{\mathrm{down}} (\mathrm{W\,m}^{-2})$")
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.set_xlabel('UTC Time of '+ self.flight_date +" ("+self.flight+")")
        ax2.grid(True)
        
        ax1.legend()
        ax2.legend()
        sns.despine(offset=10)
        iner_corr_fig.savefig(self.process_plot_path+fig_name, 
                              dpi=300, bbox_inches="tight")
        self.logger.info("Inertia correction comparison plotted under: "+\
                         self.process_plot_path+fig_name)
    def plot_differences_inertia_correction(self):
        iner_diff_fig=plt.figure(figsize=(12,9))
        matplotlib.rcParams.update({"font.size":20})
        fig_name=self.flight+"_"+\
            "BP_L1_processing_Inertia_correction_difference.png"
        ax1=iner_diff_fig.add_subplot(211)
        ax2=iner_diff_fig.add_subplot(212)
        # Irradiance F_upward 
        ax1.plot(self.lw_u_L1-self.lw_u_tcorr,c="darkblue",lw=3,label="inertia corrected-raw")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.set_xlabel("")
        ax1.set_ylabel("$F_{\mathrm{up}} (\mathrm{W\,m}^{-2})$")
        ax1.grid(True)
        # Irradiance F_downward
        ax2.plot(self.lw_d_L1-self.lw_d_tcorr,c="darkred",
                 lw=3,label="inertia corrected")
        
        ax2.set_ylabel("$F_{\mathrm{down}} (\mathrm{W\,m}^{-2})$")
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax2.set_xlabel('UTC Time of '+ self.flight_date +" ("+self.flight+")")
        ax2.grid(True)
        
        ax1.legend()
        ax2.legend()
        sns.despine(offset=10)
        iner_diff_fig.savefig(self.process_plot_path+fig_name, 
                              dpi=300, bbox_inches="tight")
        self.logger.info("Inertia correction difference plotted under: "+\
                         self.process_plot_path+fig_name)
    def quicklook_inertia_correction_profile(self):
        if self.flight=="RF12":
            time_period=["2024-04-01 09:45",
                         "2024-04-01 10:04"]
        merge_raw_index=self.lw_u_tcorr.index.intersection(
            self.lw_d_tcorr.index)
        merge_raw_index=merge_raw_index.intersection(self.bme_df.index)
        
        merge_simple_index=self.lw_u_simple.index.intersection(
            self.lw_d_simple.index)
        merge_simple_index=merge_simple_index.intersection(
            self.bme_df.index)
        
        merge_corr_index=self.lw_u_L1.index.intersection(
            self.lw_d_L1.index)
        merge_corr_index=merge_corr_index.intersection(
            self.bme_df.index)
        
        raw_df    = pd.DataFrame(data=np.nan,index=merge_raw_index,
                                 columns=["F_up","F_down","pressure"])
        simple_df = pd.DataFrame(data=np.nan,index=merge_simple_index,
                                 columns=["F_up","F_down","pressure"])
        corr_df   = pd.DataFrame(data=np.nan,index=merge_corr_index,
                                 columns=["F_up","F_down","pressure"]) 
        
        raw_df["F_up"]            = self.lw_u_tcorr.loc[merge_raw_index].values
        raw_df["F_down"]          = self.lw_d_tcorr.loc[merge_raw_index].values
        raw_df["pressure"]    = self.bme_df["bme_pressure"].loc[merge_raw_index].values
        
        simple_df["F_up"]         = self.lw_u_simple.loc[merge_simple_index].values
        simple_df["F_down"]       = self.lw_d_simple.loc[merge_simple_index].values
        simple_df["pressure"] = self.bme_df["bme_pressure"].loc[merge_simple_index].values
        
        corr_df["F_up"]           = self.lw_u_L1.loc[merge_corr_index].values
        corr_df["F_down"]         = self.lw_d_L1.loc[merge_corr_index].values
        corr_df["pressure"]   = self.bme_df["bme_pressure"].loc[merge_corr_index].values
        
        raw_cut_df    = raw_df.loc[time_period[0]:time_period[1]]
        simple_cut_df = simple_df.loc[time_period[0]:time_period[1]]
        corr_cut_df   = corr_df.loc[time_period[0]:time_period[1]]
        
        profile_fig=plt.figure(figsize=(16,16))
        
        ax1=profile_fig.add_subplot(121)
        ax1.plot(simple_cut_df["F_up"].values,simple_cut_df["pressure"].values/100,
                 color="grey",ls="--", lw=2, label="simple corr")
        ax1.plot(raw_cut_df["F_up"].values,raw_cut_df["pressure"].values/100,
                 color="lightblue",lw=2,label="uncorrected",zorder=2)
        ax1.plot(corr_cut_df["F_up"].values,corr_cut_df["pressure"].values/100,
                 color="darkblue",lw=3,label="inertia corrected")
        ax1.set_ylabel("Pressure (hPa)")
        ax1.invert_yaxis()
        
        ax2=profile_fig.add_subplot(122)
        ax2.plot(simple_cut_df["F_down"],simple_cut_df["pressure"]/100,
                 color="grey",ls="--", lw=2, label="simple corr")
        ax2.plot(raw_cut_df["F_down"],raw_cut_df["pressure"]/100,
                 color="salmon",lw=2,label="uncorrected",zorder=2)
        ax2.plot(corr_cut_df["F_down"],corr_cut_df["pressure"]/100,
                 color="darkred",lw=3,label="inertia corrected")
        ax2.set_ylabel("Pressure (hPa)")
        ax2.invert_yaxis()
        ax1.legend()
        ax2.legend()
        ax1.set_xlabel("$F_{\mathrm{up}} (\mathrm{W\,m}^{-2})$")
        ax2.set_xlabel("$F_{\mathrm{down}} (\mathrm{W\,m}^{-2})$")
        
        fig_name=self.flight+"_"+"BP_quicklook_profile_Inertia_correction.png"
        
        profile_fig.savefig(self.process_plot_path+fig_name, 
                              dpi=300, bbox_inches="tight")
        self.logger.info("Quicklook profile inertia correction comparison under: "+\
                         self.process_plot_path+fig_name)
    
    
    
    def plot_quicklook_zb(self,sensor="bme"):
        #rough quicklook of Z_b
        matplotlib.rcParams.update({"font.size":18})
        z_b_fig=plt.figure(figsize=(12,9))
        ax1=z_b_fig.add_subplot(111)
        ax1.plot(self.l1_df["z_b"],color="orange",label="z_b")
        if sensor=="TMP_sonde":
            ax1.plot(self.l1_df["ALT_R"],color="grey",lw=1,ls="--",
                     label="GPS")
        ax1.grid()
        ax1.set_ylabel("Barometric height in m")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.set_xlabel('UTC Time of '+self.flight_date +\
            " ("+self.flight+")")
        ax1.legend()
        fig_name=sensor.upper()+"-based_z_b_"+self.flight+"_"+\
            self.flight_date+".png"
        fig_path=self.process_plot_path+"/"
        fname=fig_path+fig_name
        z_b_fig.savefig(fname,dpi=300,bbox_inches="tight")
        self.logger.info("plotted z_b method under: "+fname)
    def plot_up_and_down_ward_terr_radiation(self,):
        import beluga_plotting
        import matplotlib.ticker as mticker
        import matplotlib.colors as mcolors
        from matplotlib.gridspec import GridSpec

        matplotlib.rcParams.update({"font.size":24})
        df=self.l2_df.copy()
        df["ALT"]=df["z_b"]
        df.name=self.flight
        max_height=df["z_b"].max()    
        quick_scatter=plt.figure(figsize=(16,16))
        import beluga_plotting
        gs = GridSpec(2, 2, figure=quick_scatter,
            height_ratios=[0.35,1],wspace=.2)    
        ax0=quick_scatter.add_subplot(gs[0,:])
        beluga_plotting.plot_RF_height(df,ax_obj=ax0)      
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
                       norm=norm_main,s=3,label=self.flight) # BELUGA
        ax1.set_ylabel("Barometric height AGL (m)")
        ax2.scatter(df["F_down"],df["z_b"],c=df.index,
                    cmap=sns_cmap,norm=norm_main,s=3,label=self.flight)
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
        ax1.set_xlabel('Upward Thermal-\nInfrared Irradiance ($\mathrm{W\,m}^{-2})$')
        ax2.set_xlabel('Downward Thermal-\nInfrared Irradiance ($\mathrm{W\,m}^{-2})$')
        f_up_p5,f_up_p95=df["F_up"].quantile([0.05,0.95])
        f_down_p5,f_down_p95=df["F_down"].quantile([0.05,0.95])
        
        ax1.set_xlim([f_up_p5//10*10-10,f_up_p95//10*10+10])
        ax2.set_xlim([f_down_p5//10*10-10,f_down_p95//10*10+10])
        
        plt.subplots_adjust(wspace=0.5,hspace=0.3)
        sns.despine(offset=10)
        plot_path=self.process_plot_path+"/"
        fig_name="Vertical_profiles_F_up_F_down_"+self.flight+"_"+\
            self.flight_date+".png"
        fname=plot_path+fig_name
        quick_scatter.savefig(fname,dpi=300,bbox_inches="tight")
        self.logger.info("Figure saved as:"+fname)
        
    def plot_net_terr_radiation(self,):
        import beluga_plotting
        import matplotlib.ticker as mticker
        import matplotlib.colors as mcolors
        from matplotlib.gridspec import GridSpec

        matplotlib.rcParams.update({"font.size":24})
        df=self.l2_df.copy()
        df["ALT"]=df["z_b"]
        df.name=self.flight
        max_height=df["z_b"].max()    
        quick_scatter=plt.figure(figsize=(16,16))
        import beluga_plotting
        gs = GridSpec(2, 1, figure=quick_scatter,
            height_ratios=[0.35,1],wspace=.2)    
        ax0=quick_scatter.add_subplot(gs[0])
        beluga_plotting.plot_RF_height(df,ax_obj=ax0)      
        ax1=quick_scatter.add_subplot(gs[1])
        sub_fig_labels=["(b)"]
        # Create a proxy artist for the legend with a larger marker size
        # Scatter plots
        sns_cmap=sns.color_palette("icefire", as_cmap=True)
        ## Normalize time index for color mapping
        time=df.index.values
        norm_main = mcolors.Normalize(vmin=time.min(), vmax=time.max())
        # Plot 3: Balloon    
        ax1.scatter(df["F_net"],df["z_b"],c=df.index,cmap=sns_cmap,
            norm=norm_main,s=3,label=self.flight) 
        ax1.set_ylabel("Barometric height AGL (m)")
        #-------------------------------------------------------------------------#
        if max_height<500:
           ax1.set_ylim([0,500])
           ax1.set_yticks([0,250,500])   
        elif 500<max_height<700:
            ax1.set_ylim([0,700])
            ax1.set_yticks([0,100,200,300,400,500,600,700])
        elif 700<max_height<800:
            ax1.set_ylim([0,800])
            ax1.set_yticks([0,200,400,600,800])
        else:
            ax1.set_ylim([0,1000])
            ax1.set_yticks([0,200,400,600,800,1000])
        
        for axis in ['bottom','left']:
            ax1.spines[axis].set_linewidth(2)
            ax1.yaxis.set_tick_params(width=2,length=6)
            ax1.xaxis.set_tick_params(width=2,length=6)
            
        ax1.text(0.03,0.95,sub_fig_labels[0],transform=ax1.transAxes,
                 fontsize=20)
        ax1.set_xlabel('Net Terrestrial \nRadiation ($\mathrm{W\,m}^{-2})$')
        f_net_p5,f_net_p95=df["F_net"].quantile([0.05,0.95])
        ax1.set_xlim([f_net_p5//10*10,f_net_p95//10*10+10])
        plt.subplots_adjust(wspace=0.5,hspace=0.3)
        sns.despine(offset=10)
        plot_path=self.process_plot_path+"/"
        fig_name="Vertical_profiles_F_net_"+self.flight+"_"+\
            self.flight_date+".png"
        fname=plot_path+fig_name
        quick_scatter.savefig(fname,dpi=300,bbox_inches="tight")
        self.logger.info("Figure saved as:"+fname)
        
        
    def plot_flight_segmentation(self,add_peak_infos=False):
        
        rf=self.flight
        
        df_initial = self.l2_df.copy()
        df         = df_initial.copy()
        
        # Do not separate between profile numbers anymore
        df["segments"][df_initial["segments"].str.contains(
            'ascent', case=False, na=False)] = 'ascent'
        df["segments"][df_initial["segments"].str.contains(
            'descent', case=False, na=False)] = 'descent'
        #---------------------------------------------------------------------#
        #% BELUGA Plotting
        fig=plt.figure(figsize=(16,12))
        
        ax1=fig.add_subplot(111)
        ax1.plot(df["z_b"],lw=4,color="w")
        ax1.plot(df["z_b"],lw=2,color="grey")
        
        # Colour-code flight segmentation
        segmentation_classes=["ascent","descent","peak",
                              "near_ground", "nearly_constant_altitude"]
        
        segmentation_cls_colors=["blue","orange","red",
                                 "sienna","green"]
        
        for seg,clr in zip(segmentation_classes,segmentation_cls_colors):
            seg_df=df["z_b"][df["segments"]==seg]
            ax1.scatter(seg_df.index,seg_df,marker="s",s=40,
                        color=clr,zorder=2,label=seg)
        
        # Add extra marker for maximum heights    
        ax1.scatter(self.profile_peaks.index,
            self.profile_peaks, marker="o", s=100,
            color="red",edgecolor="k",zorder=3)
        max_max=df["z_b"].dropna().max()
        ax1.scatter(df["z_b"].idxmax(),max_max,
                    marker="o",s=200,color="darkred",
                    edgecolor="k",zorder=4)
        ax1.axhline(y=200, color="darkgrey", lw=3,ls="--")
        
        # Add profile peak infos
        if add_peak_infos:
            ax1.text(df["z_b"].idxmax(),max_max+20,
                 str(int(max_max))+" m")
            ax1.text(0.1,0.9,"Number of profiles: "+\
                     str(2*len(self.profile_peaks)),
                 transform=ax1.transAxes)
        
        ax1.set_xlabel(str(df.index.date[0])+" Time (UTC)")
        ax1.set_ylabel("Barometric height (m)")
        ax1.set_title(rf)
        
        if max_max>600:
            ylim_max=900
        else:
            ylim_max=600
        
        ax1.set_ylim([0,ylim_max])
        ax1.legend()
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        
        sns.despine(offset=10)
        for axis in ['bottom','left']:
            ax1.spines[axis].set_linewidth(2)
        ax1.yaxis.set_tick_params(width=2,length=6)
        ax1.xaxis.set_tick_params(width=2,length=6)

        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        plot_path=self.process_plot_path
        fig_name="BELUGA_BP_flight_segmentation_"+rf
        file_end=".png"
        fig_name+=file_end
        fname=plot_path+fig_name
        fig.savefig(fname,dpi=300, bbox_inches="tight")
        print("Figure saved as:", fname)
        self.logger.info("Figure saved as:"
                         +plot_path+fname)
#%% Radiosondes   
class Radiosondes(Measurement_Platforms_VRS):
    def __init__(self, Measurement_Platforms_VRS,rf,
            run_L1_processing=True,run_L2_processing=True,
            plot_processing=True):
        self.station_name       = Measurement_Platforms_VRS.station_name        
        self.coordinates        = Measurement_Platforms_VRS.coordinates
        self.main_path          = Measurement_Platforms_VRS.main_path
        self.flight_dates       = {
            "RFS1":"20240321","RFS2":"20240322","RFS3":"20240323",
            "RFS4":"20240331","RF01":"20240324","RF02":"20240324","RF03":"20240325",
            "RF04":"20240325","RF05":"20240326","RF06":"20240326",
            "RF07":"20240327","RF08":"20240328","RF09":"20240328",
            "RF10":"20240329","RF11":"20240330","RF12":"20240401",
            "RF13":"20240401","RF14":"20240401","RF15":"20240402",
            "RF16":"20240402","RF17":"20240403","RF18":"20240404",
            "RF19":"20240404","RF20":"20240405","RF21":"20240406",
            "RF22":"20240407","RF23":"20240408","RF24":"20240408",
            "RF25":"20240409","RF26":"20240410","RF27":"20240411",
            "RF28":"20240412",
            #artificial merged dates
            "RF1312":"20240401","RF141312":"20240401"}
        
        self.launches       = {
            "20240321":1,"20240322":2,"20240323":3,
            "20240324":4,"20240325":5,"20240326":6,
            "20240327":7,"20240328":8,"20240329":9,
            "20240330":10,"20240331":11,"20240401":12,
            "20240402":13,"20240403":14,"20240404":15,
            "20240405":16,"20240406":17,"20240407":18,
            "20240408":19,"20240409":20,"20240410":21,
            "20240411":22}
        
        self.flight             = rf
        self.flight_date        = self.flight_dates[self.flight]
        
        self.raw_data_path      = self.main_path +\
           "/radiosondes/"
        self.plot_path          = self.main_path+"/plotting/"
        self.process_plot_path  = self.plot_path+"/test_processing/"+\
            self.flight+"/"
        os.makedirs(self.process_plot_path,exist_ok=True)
        # Switches
        self.run_L1_processing  = run_L1_processing
        self.run_L2_processing  = run_L2_processing
        self.plot_processing    = plot_processing
    
        # Configure logging once at the module level
        self.log_file_name='Radiosonde_processing_'+self.flight_date+'.log'
        self.logger=logging_processing.setup_logging(self.log_file_name)
    
    def get_balloon_flight_information(self):
        self.flight_infos=pd.read_csv(self.main_path+"flight_schedule.csv",
                                      index_col="Flight No")
        self.flight_infos["Start Time"]     = pd.to_datetime(
                self.flight_infos["Start Time"], dayfirst=True)
        self.flight_infos["End Time"]       = pd.to_datetime(
                self.flight_infos["End Time"], dayfirst=True)
        
        self.flight_infos["type"]           = "uncategorized"
        self.flight_infos["type"].loc[["RF07","RF08","RF10","RF15",
            "RF16","RF19","RF23","RF27"]]   = "cloudfree-cloudy"
        self.flight_infos["type"].loc[["RF05","RF06"]]    = "day-to-night"
        self.flight_infos["type"].loc[["RF12","RF13","RF14","RF26"]]="LLJ"

        #----------------------------------------------------------------------#
        self.LLJ_flights    = self.flight_infos.loc[\
                                self.flight_infos["type"]=="LLJ"]
        self.cloud_flights  = self.flight_infos.loc[\
                                self.flight_infos["type"]=="clear to cloudy"]
        self.night_flights  = self.flight_infos.loc[\
                                self.flight_infos["type"]=="Polar night to day"]    
    def add_VRS_georeferences(self,ds):
        
        ds["VRS_lat"]                           = self.coordinates["lat"]
        ds["VRS_lon"]                           = self.coordinates["lon"]
        ds["VRS_height"]                        = self.coordinates["height"]
        #Attributes:
        # add attributes based on CF-convention
        ds["VRS_lat"].attrs["standard_name"]    = "latitude"
        ds["VRS_lat"].attrs["long_name"]        = "latitude Villum Research Station"
        ds["VRS_lat"].attrs["units"]            = "degrees_north"
        
        ds["VRS_lon"].attrs["standard_name"]    = "longitude"
        ds["VRS_lon"].attrs["long_name"]        = "longitude Villum Research Station"
        ds["VRS_lon"].attrs["units"]            = "degrees_east"
        
        #ds["STN_height"].attrs["standard_name"] = "height"
        ds["VRS_height"].attrs["long_name"]     = "Villum Research Station height above mean sea level"
        ds["VRS_height"].attrs["units"]         = "m"
        
        # add attributes based on CF-convention
        ds.attrs["Location"] = "Villum Research Station: "+str(self.coordinates["lat"])+\
            "°N, "+str(self.coordinates["lon"])+"°E."+\
            " Height above ground level: "+str(self.coordinates["height"])+" m"
        return ds 
    
    def add_ABL_transition_type(self,ds):
        try:
            self.get_balloon_flight_information()
            transition=self.flight_infos["type"].loc[self.flight]
        except:
            transition="uncategorised"
            #self.flight_infos["type"].loc[self.flight]
            
        ds.attrs["ABL_transition_type"]=transition
        return ds
    
    def add_L2_meta_data(self):
        
        import nc_attributes
        global_att_dict = nc_attributes.Radiosonde_global_attributes
        sonde_att_dict  = nc_attributes.Radiosonde_var_attributes 
        
        # Add Global Attributes from dictionary
        for att in [*global_att_dict.keys()]:
            self.l2_ds.attrs[att]=global_att_dict[att]
        
        for var in [*self.l2_ds.keys()]:
            #if var=="time_since_launch":
        #
            var_attrs=sonde_att_dict[var]
            for att in [*var_attrs.keys()]:
                self.l2_ds[var].attrs[att]=var_attrs[att]
        # Set attributes for the time and time_since_launch coordinates
        self.l2_ds.time.attrs['long_name']     = 'measurement time (UTC)'
        self.l2_ds.time.attrs['standard_name'] = 'time'
        self.l2_ds.time.attrs['units']         = "seconds since 1970-01-01 00:00:00 UTC"
        self.l2_ds.time.attrs["calendar"]      = "standard",

        if "time_since_launch" in [*self.l2_ds.variables]:
            self.l2_ds.time_since_launch.attrs['long_name'] = 'Time since launch'
            self.l2_ds.time_since_launch.attrs['standard_name'] = 'time' 
        # Not a perfect fit, but best available
        self.l2_ds=self.update_rf_specfic_attributes(self.l2_ds)
        
    
    def replace_inf(self,df):
        df.replace([np.inf,-np.inf],np.nan,inplace=True)
        return df
    
    def update_rf_specfic_attributes(self,ds):
        # Get current time 
        from datetime import datetime
        now = datetime.now()
        date_time = now.strftime("%Y-%m-%d, %H:%M:%S")
        ds.attrs["flight_date"]         = self.flight_date
        ds.attrs["processing_date"]     = date_time
        return ds
    
    def save_L2_data_as_nc_file(self,ds,version_number,fill_value=-9999.,
                                temporary=True):
        if temporary:
            outp_path=self.raw_data_path+"/temporary_data/" 
        else:
            outp_path=self.raw_data_path+"/final_data/"
        
        os.makedirs(outp_path,exist_ok=True)
        l2_fname="Radiosonde_VRS_L2_"+self.flight_date+".nc"
        # Create fill values
        ds=ds.fillna(value=fill_value)
        nc_compression=dict(zlib=True,_FillValue=fill_value,
                            complevel=1,dtype=np.float64)
        
        nc_encoding = {}
        for var_name, var in ds.variables.items():
            if np.issubdtype(var.dtype, np.number):
                nc_encoding[var_name] = nc_compression
            else:
                # Omit encoding for non-numeric variables 
                nc_encoding[var_name] = {}
        
        ds.to_netcdf(outp_path+l2_fname,mode="w",engine="netcdf4",
                     encoding=nc_encoding)
        self.logger.info("Radiosonde L2-data saved as: "+outp_path+l2_fname)
        print("radiosonde time",ds.time[0],ds.time[-1])
        return outp_path+l2_fname
    
    def run_processing(self, version_number,
                       run_L1_processing=True,
                       run_L2_processing=False):
        self.version_number=version_number
        self.logger.info(
            f"Processing radiosonde da for launch date {self.flight_date}\n")
        if self.run_L1_processing:
            self.read_raw_profile()
            self.read_GPS_file()
            try:
                self.read_profile_file()
                self.profile_file_exists=True
            except:
                self.profile_file_exists=False
                self.logger.info("No internally-processed profile file available.")
            if self.raw_file_exists:
                self.combine_raw_and_gps_file()
            self.L0_to_L1_processing()
            
            if self.plot_processing:
                if self.profile_file_exists and self.raw_file_exists:
                    self.plot_compare_raw_and_profile_datasets(
                        with_wind_components=True)
            #self.merge_L1_data()
        if self.run_L2_processing:
            self.L1_to_L2_processing()
    
    def get_trajectory_id(self):
        self.l2_ds["launch"]=xr.DataArray(
            data=self.launches[self.flight_date],
            dims      = ["time"],
            coords    = {"time":self.l2_ds.time.values})
        self.l2_ds["launch"].attrs["long_name"]="radiosonde launch ID"
        self.l2_ds["launch"].attrs["cf_role"]="trajectory_id"
        self.l2_ds["launch"].attrs["units"]="1"
    def read_raw_profile(self):
        raw_cols = ['time_since_launch', 'time', 'T', 'RH', 'Tu_Cap']
        date_fpath=self.main_path+"/radiosondes/"+\
            self.flight_date+"/"
        if os.path.exists(date_fpath):
            raw_files=glob.glob(date_fpath+self.flight_date+"*-raw*")
            if len(raw_files)==0:
                self.logger.info("no raw profile data available for flight date ",
                                 self.flight_date)
                self.raw_file_exists=False
            else:
                df_raw = pd.read_csv(raw_files[0],
                    encoding='latin1', names=raw_cols, 
                    skiprows=1, sep='\t', on_bad_lines='skip',
                    na_values=['--', '---','--    ',
                               '--          '])
                datetime_str_raw = self.flight_date + 'T' + df_raw['time'].str.strip()
                dt_pattern = r'(\d{4})(\d{2})(\d{2})T(\d{2}):(\d{2}):(\d{2})\.(\d{1})'
                td_pattern = r'(\d{2}):(\d{2})\.(\d{1})'
    
                raw_valid_times_mask = pd.Series(np.ones(
                    len(datetime_str_raw),dtype=bool))
                
                raw_seconds_since_launch = pd.Series(np.zeros(
                    len(datetime_str_raw), dtype=float))
                
                for i, time_str in enumerate(datetime_str_raw):
                    timedelta_str = df_raw['time_since_launch'].iloc[i]
                    try:
                        timedelta_seconds_raw = sum(x * float(t) for x, 
                                                    t in zip([60, 1],
                                                             timedelta_str.split(':')))
                    except ValueError:
                        print(f"Invalid timedelta format: {timedelta_str}")
                        timedelta_seconds_raw = np.nan
                    raw_seconds_since_launch[i] = timedelta_seconds_raw
                # ocassionally strange values exist
                df_raw=df_raw[df_raw["time"]!='                   ']
                datetime_str_raw=datetime_str_raw.loc[df_raw.index]
                raw_valid_times_mask=raw_valid_times_mask.loc[df_raw.index]
                raw_seconds_since_launch=raw_seconds_since_launch.loc[df_raw.index]
                df_raw['time']              = pd.to_datetime(
                    datetime_str_raw[raw_valid_times_mask],
                    format='%Y%m%dT%H:%M:%S.%f')
                #df_raw["time"]              = pd.to_datetime(
                #    datetime_str_raw,format=='%Y%m%dT%H:%M:%S.%f')
                df_raw['time_since_launch'] = pd.to_timedelta(
                    raw_seconds_since_launch, unit='s')
                self.df_raw=df_raw
                self.df_raw["T"] = self.df_raw["T"].apply(
                    pd.to_numeric,errors='coerce')
                #
                #self.df_raw["T"]=self.df_raw["T"].astype(float)
                self.df_raw=self.df_raw.loc[self.df_raw["T"]<0]
                self.raw_file_exists=True
        else:
            self.logger.info("No radiosonde launch for flight date ",
                             self.flight_date)
            self.raw_file_exists=False
            
    def read_GPS_file(self):
        GPS_cols = ['time_since_launch', 'time', 'distance',
                    'vertical_speed', 'horizontal_speed', 'longitude',
                    'latitude', 'altitude']
        
        date_fpath=self.main_path+"/radiosondes/"+\
            self.flight_date+"/"
        
        if os.path.exists(date_fpath):
            gps_files=glob.glob(date_fpath+self.flight_date+"*-GPS*")
            if len(gps_files)==0:
                self.logger.info(
                    "no GPS profile data available for flight date ",
                    self.flight_date)
            else:
                df_gps = pd.read_csv(gps_files[0], encoding='latin1',
                            names=GPS_cols, skiprows=1,sep='\t', 
                            on_bad_lines='skip', na_values=['--', '---'])
    
                datetime_str_gps = self.flight_date+'T'+\
                    df_gps['time'].str.strip()
    
                ### iterate thorugh the datetime_str array to
                # ensure that we only take values, where the time is correctly recorded, 
                ### i.e. it matches the format
    
                dt_pattern = r'(\d{4})(\d{2})(\d{2})T(\d{2}):(\d{2}):(\d{2})\.(\d{1})'
                td_pattern = r'(\d{2}):(\d{2})\.(\d{1})'
    
                gps_valid_times_mask = np.ones(len(datetime_str_gps), dtype=bool)
    
                gps_seconds_since_launch = np.zeros(len(datetime_str_gps),
                                                    dtype=float)
    
                for i, time_str in enumerate(datetime_str_gps):
                    timedelta_str = df_gps['time_since_launch'].iloc[i]
                    try:
                        timedelta_seconds_gps = sum(x * float(t) for x, t \
                                        in zip([60, 1], timedelta_str.split(':')))
                    except ValueError:
                        print(f"Invalid timedelta format: {timedelta_str}")
                        timedelta_seconds_gps = np.nan
                        gps_valid_times_mask[i] = False
                    gps_seconds_since_launch[i] = timedelta_seconds_gps
                # assign new time values to df       
                df_gps['time'] = pd.to_datetime(datetime_str_gps[\
                        gps_valid_times_mask],format='%Y%m%dT%H:%M:%S.%f')
                
                df_gps['time_since_launch'] = pd.to_timedelta(
                    gps_seconds_since_launch, unit='s')
                
                self.df_gps=df_gps
        else:
            self.logger.info("No radiosonde launch for flight date ",
                     self.flight_date)
        
    def read_profile_file(self,remove_descent_from_data=True):
        """
        Caution! Files might be affected by erroneous data transfer. 
        
        Columns: 'Zeit [sec]', 'P [hPa]', 'T [°C]', 'Hu [%]', 
        'Ws [m/s]', 'Wd [°]', 'Ws U [m/s]', 'Ws V [m/s]',
        'Geopot [m]', 'Dew [°C]']
        """
        date_fpath=self.main_path+"/radiosondes/"+\
            self.flight_date+"/"
        
        if os.path.exists(date_fpath):
            profile_files=glob.glob(date_fpath+self.flight_date+"*-Profil*")
            if len(profile_files)==0:
                self.logger.info(
                    "no GPS profile data available for flight date ",
                    self.flight_date)
            else:
                df_profile = pd.read_csv(profile_files[0], encoding='latin1',
                    sep='\t', on_bad_lines='skip',index_col=1, 
                    na_values=['--', '---'])
                
                self.profile_df=df_profile
                
                self.profile_df.index=self.flight_date[0:4]+"-"+\
                    self.flight_date[4:6]+"-"+\
                        self.flight_date[6:8]+" "+self.profile_df.index
                self.profile_df.index=pd.DatetimeIndex(self.profile_df.index)
           
            if remove_descent_from_data:
                ascent_end_time=self.profile_df["Geopot [m]"].idxmax()
                start_time=self.profile_df.index[0]
                self.profile_df=self.profile_df.loc[start_time:ascent_end_time]
            
            # remove outliers
            self.profile_df[self.profile_df["Geopot [m]"] < 5] = np.nan 
                    
    def combine_raw_and_gps_file(self):
        start_time  = max(self.df_raw['time'].min(), 
                     self.df_gps['time'].min())
        end_time    = min(self.df_raw['time'].max(),
                   self.df_gps['time'].max())
        # Cut to periods 
        self.df_raw = self.df_raw[(self.df_raw['time'] >= start_time) &\
                              (self.df_raw['time'] <= end_time)]
        self.df_gps = self.df_gps[(self.df_gps['time'] >= start_time) &\
                              (self.df_gps['time'] <= end_time)]
        
        # Reset index based on synchronised periods
        self.df_raw.set_index('time', inplace=True)
        self.df_gps.set_index('time', inplace=True)
    
        ### convert all variables except time and time_since_launch to float
        for df in [self.df_raw, self.df_gps]:
            for col in df.columns:
                if col not in ['time', 'time_since_launch']:
                    df[col] = pd.to_numeric(df[col], 
                        errors='coerce', downcast='float')
                    #if not re.match(dt_pattern, time_str):
                    #    print(f"Time format mismatch: {time_str}")
                    #    raw_valid_times_mask[i] = False
    
        ds_raw = self.df_raw.to_xarray()
        ds_gps = self.df_gps.to_xarray()
    
        self.combined_ds = xr.merge([ds_raw, ds_gps], compat="override")
        
    def check_for_outliers(self):
        l1_ds = self.combined_ds.copy()
        l1_ds = l1_ds.where(l1_ds['altitude'] >= 5)#
        if "longitude" in [*l1_ds.variables]:
            l1_ds = l1_ds.where(l1_ds['longitude'] != 0)
        if "vertical_speed" in [*l1_ds.variables]:
            l1_ds = l1_ds.where(l1_ds['vertical_speed'] >= 0)
        self.logger.info("Outliers removed from data.")
        self.l1_ds=l1_ds
    def remove_descent(self):
        # find the max value of the profile altitude 
        # and remove all data afterwards
        ascent_end_time=self.l1_ds["altitude"].idxmax()
        start_time=self.l1_ds["time"][0]
        self.l1_ds=self.l1_ds.sel(time=slice(start_time,ascent_end_time))
        self.logger.info("Descent data is erased.")
    def resample_ds_resolution(self,resolution="2s"):
        if self.raw_file_exists:
            l1_df=self.l1_ds.to_dataframe()
            time_since_launch    = l1_df["time_since_launch"]
            if self.profile_file_exists:
                l1_df=l1_df.loc[self.profile_df.index[0]:\
                                self.profile_df.index[-1]]
                
        else:
            l1_df=self.profile_df    
        
        l1_res_df=l1_df.resample(resolution).mean()
        #Interpolate gaps
        l1_res_df=l1_res_df.interpolate(method="time",limit=5)
        if self.raw_file_exists:
            time_since_launch_res = time_since_launch.reindex(
                l1_res_df.index, method='nearest')
            l1_res_df["time_since_launch"]=time_since_launch_res.values
        
        ds=l1_res_df.to_xarray()
        
        return ds
    
    #def calculate_wind_components_from_wspeed_wdir(self):
    #    from metpy.calc import wind_components
    #    ds=self.l1_ds.copy()
    #    u,v,
    #    u_da = xr.DataArray(u_vals, coords=[ds_proc.time], dims=["time"])
    #    v_da = xr.DataArray(v_vals, coords=[ds_proc.time], dims=["time"])
    def calculate_wind_components_from_displacement(self):
        """
        Calculates u and v wind components from time-series position data.

        This function takes an xarray.Dataset containing time, latitude, and 
        longitude data, computes the eastward (u) and northward (v) wind 
        components based on the balloon's drift, and returns the dataset 
        with 'u' and 'v' variables populated.

        """
    
        def _calculate_bearing(lat1, lon1, lat2, lon2):
            """Calculates the initial bearing in radians between two points."""
            lat1_rad, lon1_rad = math.radians(lat1), math.radians(lon1)
            lat2_rad, lon2_rad = math.radians(lat2), math.radians(lon2)
            
            dLon = lon2_rad - lon1_rad
            y = math.sin(dLon) * math.cos(lat2_rad)
            x = math.cos(lat1_rad) * math.sin(lat2_rad) - \
                math.sin(lat1_rad) * math.cos(lat2_rad) * math.cos(dLon)
        
            return math.atan2(y, x)
        
        ds=self.l1_ds.copy()
        #Pre-process the data to handle NaNs for a continuous track
        # I think it's ok to interpolate the lat/lon
        ds['latitude'] = ds['latitude'].interpolate_na(dim='time')
        ds['longitude'] = ds['longitude'].interpolate_na(dim='time')
    
        ds_proc = ds.ffill(dim='time').dropna(
            subset=['latitude', 'longitude'], dim='time', how='any'
        )
    
        # 2. Extract cleaned data as numpy arrays for efficient iteration
        lat = ds_proc['latitude'].values
        lon = ds_proc['longitude'].values
        times = ds_proc['time'].values
    
        # 3. Initialize arrays to store the calculated values
        u_vals = np.full(lat.shape, np.nan)
        v_vals = np.full(lat.shape, np.nan)
    
        # 4. Iterate through data points to calculate velocity between each point
        for i in range(1, len(lat)):
            pos1 = (lat[i-1], lon[i-1])
            pos2 = (lat[i], lon[i])
    
            # Time difference in seconds
            dt = (times[i] - times[i-1]) / np.timedelta64(1, 's')
            if dt == 0:
                continue
    
            # Calculate horizontal distance (m) and speed (m/s)
            distance_m = geodesic(pos1, pos2).meters
            h_speed = distance_m / dt
            
            # Calculate the bearing (direction of travel) in radians
            bearing_rad = _calculate_bearing(lat[i-1], lon[i-1], lat[i], lon[i])
            
            # Convert speed and bearing to u and v components
            u_vals[i] = h_speed * math.sin(bearing_rad) # Eastward component
            v_vals[i] = h_speed * math.cos(bearing_rad) # Northward component
    
        # 5. Create new DataArrays and align them with the original dataset's time index
        u_da = xr.DataArray(u_vals, coords=[ds_proc.time], dims=["time"])
        v_da = xr.DataArray(v_vals, coords=[ds_proc.time], dims=["time"])
        
        ds['u'] = u_da.reindex_like(ds)
        ds['v'] = v_da.reindex_like(ds)
        
        self.l1_ds=ds
        
    def calc_MAE_to_profile_file(self):
        
        raw_1min=self.l1_ds[["T","RH","u","v","altitude"]].to_dataframe().\
                          resample("1min").mean()
        
        profile_1min=self.profile_df.resample("1min").mean()
        #, , 'Ws [m/s]', 'Wd [°]',
        #       'Ws U [m/s]', 'Ws V [m/s]', 'Geopot [m]', 'Dew [°C]']
        profile_1min=profile_1min.rename(columns={'T [°C]':"T",'Hu [%]':"RH",
                        'Ws U [m/s]':"u",'Ws V [m/s]':"v",
                        'Geopot [m]':"altitude"})
        
        mae_1min=abs(raw_1min-profile_1min)
        mae_all=mae_1min.mean()
        
        
        
        
        # check for maximum errors
        self.profile_ok = {"T":False,"RH":False,
                           "u":False,"v":False}
        if mae_all["T"]<1:
            self.profile_ok["T"]  = True
        if mae_all["RH"]<2:
            self.profile_ok["RH"] = True
        if mae_all["u"]<1:
            self.profile_ok["u"]  = True
        if mae_all["v"]<1:
            self.profile_ok["v"]  = True 
    def rename_profile_vars(self,profile_df):
        # Sometimes p, T have different namings. thus doubling rename 
        
        profile_df=profile_df.rename(columns={'P [hPa]':"p",
                                              'P [hPa] ':"p",
                                              'T [°C]':"T",
                                              'T [°C]  ':"T",
                                              'Hu [%]':"RH",
                    'Ws U [m/s]':"u",'Ws V [m/s]':"v",
                    'Ws [m/s]':"wspeed", 'Wd [°]':"wdir",
                    'Geopot [m]':"altitude"})
        # Reset index based on synchronised periods
        profile_df.index.name="time"
        return profile_df
    
    def add_pressure_to_l1_data(self):
        l1_2s_ds=self.l1_2s_ds.copy()
        l1_2s_df=l1_2s_ds.to_dataframe()
        l1_2s_df["time"]=l1_2s_df.index
        
        start_time  = max(l1_2s_df['time'].min(), 
                     self.profile_df.index[0])
        end_time    = min(l1_2s_df['time'].max(),
                   self.profile_df.index[-1])
        
        # Cut to periods 
        l1_2s_df = l1_2s_df.loc[start_time:end_time]
        
        profile_df = self.profile_df.loc[start_time:end_time]
        profile_df = self.rename_profile_vars(profile_df)
    
        vars_to_add=[*self.l1_2s_ds.variables]
        
        l2_ds=self.l1_2s_ds.copy()
        
        if "p" in profile_df.columns:
            l2_ds["p"]=xr.DataArray(
            data      = profile_df["p"].resample("2s").mean().values,
            dims      = ["time"],
            coords    = {"time":l1_2s_df.index})
            self.logger.info("pressure added from profile file")
        self.l2_ds=l2_ds
        
    def exchange_l1_data_with_profile_data(self):
        l1_2s_ds=self.l1_2s_ds.copy()
        l1_2s_df=l1_2s_ds.to_dataframe()
        l1_2s_df["time"]=l1_2s_df.index
        
        start_time  = max(l1_2s_df['time'].min(), 
                     self.profile_df.index[0])
        end_time    = min(l1_2s_df['time'].max(),
                   self.profile_df.index[-1])
        # Cut to periods 
        l1_2s_df = l1_2s_df.loc[start_time:end_time]#[(l1_2s_df['time'] >= start_time) &\
                              #(l1_2s_df['time'] <= end_time)]
        
        profile_df = self.profile_df.loc[start_time:end_time]
        profile_df = self.rename_profile_vars(profile_df)
        l1_2s_df.set_index('time', inplace=True)
        
        vars_to_add=[*self.l1_2s_ds.variables]
        
        times=l1_2s_df.index
        
        l2_ds=xr.Dataset({
            "time": xr.DataArray(
                data      = times,
                dims      = ['time'],
                coords    = {"time":times})
                }
            )
        
        for var in vars_to_add:
            if var == "time":
                continue
            var_da=xr.DataArray(
                data      = l1_2s_df[var].values.astype(float),
                dims      = ["time"],
                coords    = {"time":times})
            # overwrite if var is in profile_ok flag and true
            if var in [*self.profile_ok.keys()]:
                if self.profile_ok[var]:
                    var_da=xr.DataArray(
                        data      = profile_df[var].values,
                        dims      = ["time"],
                        coords    = {"time":times})
                    self.logger.info(
                        var+" exchanged with appropriate data of profile file") 
            l2_ds[var]=var_da
        
        l2_ds["p"]=xr.DataArray(
            data      = profile_df["p"].values,
            dims      = ["time"],
            coords    = {"time":times})
        self.logger.info("pressure added from profile file")
        self.l2_ds=l2_ds
        

    def L0_to_L1_processing(self):
        # mask some of the very bad data...

        if not hasattr(self,"combined_ds"):
            if self.raw_file_exists:
                self.read_GPS_file()
                self.read_raw_profile()
                self.combine_raw_and_gps_file() # creates the self.combined_ds
            else:
                self.profile_df=self.rename_profile_vars(self.profile_df)
                self.combined_ds=self.profile_df.to_xarray()
            
        self.check_for_outliers()
        self.remove_descent()
        if self.raw_file_exists:
            self.calculate_wind_components_from_displacement()
        
        
    def L1_to_L2_processing(self):
        self.l1_2s_ds=self.resample_ds_resolution()
        if self.profile_file_exists and self.raw_file_exists:
            self.calc_MAE_to_profile_file()
            self.add_pressure_to_l1_data()
            #self.exchange_l1_data_with_profile_data()
        else:
            self.l2_ds=self.l1_2s_ds.copy()
        
        if "time_since_launch" in [*self.l2_ds.variables]:
            self.l2_ds["time_since_launch"]=\
                    self.l2_ds["time_since_launch"].astype(float)
            self.logger.info("No raw vars exchanged. Treat data with caution")
        if not self.raw_file_exists:
            del self.l2_ds["Zeit [sec]"]
            del self.l2_ds["Dew [°C]"]              
        
        self.add_L2_meta_data()
        self.l2_ds=self.add_VRS_georeferences(self.l2_ds)
        self.l2_ds=self.add_ABL_transition_type(self.l2_ds)
        self.get_trajectory_id()
        if self.plot_processing:
            self.plot_L2_radiosonde()
        self.save_L2_data_as_nc_file(self.l2_ds, self.version_number,
                                     temporary=self.temporary)
    def plot_L2_radiosonde(self):
        """
        this routine plots a quicklook showing the vertical profiles of the 
        radiosonde measurements for single launches
        """
        ds=self.l2_ds
        matplotlib.rcParams.update({"font.size":12})
        fig, ax = plt.subplots(1, 8, figsize=(24, 9), sharey=True)
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
        fig_name="Quicklook_Radiosonde_L2_data_"+\
            self.flight_date+".png"
        fig_path=self.process_plot_path+"/"
        fname=fig_path+fig_name
        plt.savefig(fname,dpi=300,bbox_inches="tight")

        self.logger.info("Radiosonde L2file wind comparison" +\
                         " (Profile and Raw files)"+\
                         "plotted under: "+fname)
        
        
    def plot_compare_raw_and_profile_datasets(self,with_wind_components=True):
        matplotlib.rcParams.update({"font.size":20})
        profile_fig=plt.figure(figsize=(18,16))
        ax1=profile_fig.add_subplot(131)
        # Pressure
        try:
            ax1.plot(self.profile_df["P [hPa]"],
                 self.profile_df["Geopot [m]"]/1000,
                 color="k",lw=2,label="profile")
        except:
            ax1.plot(self.profile_df["P [hPa] "],
                 self.profile_df["Geopot [m]"]/1000,
                 color="k",lw=2,label="profile")
            
        ax1.set_ylabel("Altitude (km)")
        ax1.set_xlabel("Pressure (hPa)")
        ax1.set_ylim([0,15])
        
        ax2=profile_fig.add_subplot(132)
        # Temperature
        try:
            
            ax2.plot(self.profile_df["T [°C]"]+273.15, 
                 self.profile_df["Geopot [m]"]/1000,
                 color="red",lw=2,label="profile file")
        except:
            ax2.plot(self.profile_df["T [°C]  "]+273.15, 
                 self.profile_df["Geopot [m]"]/1000,
                 color="red",lw=2,label="profile file")
        #  raw data
        t_raw=self.l1_ds[["T","altitude"]].to_dataframe()
        t_raw.dropna(inplace=True)
        ax2.plot(t_raw["T"]+273.15,
                 t_raw["altitude"]/1000,
                 color="grey",lw=1,ls="--",label="raw file",zorder=2)
        ax2.set_xlabel("Temperature (K)")
        ax2.set_ylim([0,15])
        ax2.legend()
        
        ax3=profile_fig.add_subplot(133)
        # RH
        rh_raw=self.l1_ds[["RH","altitude"]].to_dataframe()
        rh_raw.dropna(inplace=True)
                          
        ax3.plot(self.profile_df["Hu [%]"], 
                 self.profile_df["Geopot [m]"]/1000,
                 color="blue",lw=2,label="profile file")
        
        ax3.plot(rh_raw["RH"],rh_raw["altitude"]/1000,
                 color="grey",lw=1,ls="--",label="raw file",zorder=2)
        
        ax3.set_xlabel("Relative humidity")
        ax3.set_ylim([0,15])
        ax3.legend()
        fig_name="Radiosonde_file_comparison_Profile_and_Raw_"+\
            self.flight_date+".png"
        fig_path=self.process_plot_path+"/"
        fname=fig_path+fig_name
        profile_fig.savefig(fname,dpi=300,bbox_inches="tight")
        self.logger.info("Radiosonde file comparison (Profile and Raw files)"+\
                         "plotted under: "+fname)
        
        # Add additional wind plot if desired by boolean
        if with_wind_components:
            matplotlib.rcParams.update({"font.size":20})
            wind_fig=plt.figure(figsize=(16,16))
            ax1=wind_fig.add_subplot(121)
            
            u_raw=self.l1_ds[["u","altitude"]].to_dataframe()
            u_raw.dropna(inplace=True)
            
            # u-wind
            u_raw_2s=u_raw.resample("2s").mean()
            
            ax1.plot(u_raw["u"],u_raw["altitude"]/1000,
                color="grey",lw=1,ls="--",label="raw file",zorder=0)
            ax1.plot(u_raw_2s["u"],u_raw_2s["altitude"]/1000,
                color="sandybrown",lw=1,ls="-",label="raw file (2s) mean",
                zorder=1)
            
            ax1.plot(self.profile_df['Ws U [m/s]'],
                     self.profile_df["Geopot [m]"]/1000,
                     color="brown",lw=2,label="profile file",zorder=2)
           

            ax1.set_ylabel("Altitude (km)")
            ax1.set_xlabel("u-wind (m/s)")
            ax1.set_ylim([0,15])
            ax1.legend()
            #-----------------------------------------------------------------#
            # v-wind
            ax2=wind_fig.add_subplot(122)
            #  raw data
            v_raw=self.l1_ds[["v","altitude"]].to_dataframe()
            v_raw.dropna(inplace=True)
            v_raw_2s=v_raw.resample("2s").mean()
            
            ax2.plot(v_raw["v"],v_raw["altitude"]/1000,
                     color="grey",lw=1,ls="--",
                     label="raw file",zorder=0)
            
            ax2.plot(v_raw_2s["v"],v_raw_2s["altitude"]/1000,
                     color="salmon",lw=1,ls="-",
                     label="raw file (2s) mean",zorder=1)
            
            ax2.plot(self.profile_df['Ws V [m/s]'], 
                     self.profile_df["Geopot [m]"]/1000,
                     color="red",lw=2,label="profile file",zorder=1)
            ax2.set_xlabel("v wind (m/s)")
            ax2.set_ylim([0,15])
            ax2.legend()
            fig_name="Radiosonde_file_wind_comparison_Profile_and_Raw_"+\
                self.flight_date+".png"
            fig_path=self.process_plot_path+"/"
            fname=fig_path+fig_name
            wind_fig.savefig(fname,dpi=300,bbox_inches="tight")
            self.logger.info("Radiosonde file wind comparison" +\
                             " (Profile and Raw files)"+\
                             "plotted under: "+fname)
            
            
#def remove_outliers(df,variable):
#    # height threshold =50
#    height threshold = 50
    
#    # define quantiles for height above surface
#    quantiles=df[df["ALT"]>height_threshold]].quantile([.1,.9])
#    sigma=df[variable].between(quantiles[0],quantiles[1]).std()
#    lower_threshold=quantiles[0]-sigma
#    upper_threshold=quantiles[1]+sigma
#    df[df[variable]<lower_threshold]=np.nan
#    df[df[variable]>upper_threshold]=np.nan
#    return df

#-----------------------------------------------------------------------------#
#def compute_height(pressure):
#return 44330.8 * (1 - (pressure / max(pressure)) ** 0.190263)

#def potT(T,p):
#    return (T + T0) * (p0/p)**(RL / cp)

#def ES(T):    
#    return 6.1 * np.exp(17.62 * T /(243.12 + T))  

#def e(ES,rH):
#    return rH / 100 * (ES * 100)

#def absH(rH,ES,T):
#    return (rH * ES) / Rv / (T+T0)

#def mixRatio(rH,ES,p):
#    return rH/100 * ES / p * RL/Rv * 1e-3
#-----------------------------------------------------------------------------#