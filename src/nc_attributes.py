# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 13:49:45 2025

@author: u300737
"""
#%% Global attributes
TMP_met_global_attributes={
    "title": "Tethered balloon-borne meteorological profile measurement data"+\
       " characterising the Arctic atmopsheric boundary layer at Station Nord", 
    
    "authors":"Dorff, Henning; Siebert, Holger; Schaefer, Michael;" +\
        " Mueller, Joshua; Navale, Komal; Wu, Fan; Ehrlich, Andre;" +\
            " Wendisch, Manfred",
    
    "institutes":" (1) Leipzig Institute for Meteorology," +\
        " Leipzig University; "+"(2) Leibniz Institute for" +\
            " Tropospheric Research e.V., Leipzig, Germany",
    
    "contact":"henning.dorff@uni-leipzig.de,"+\
        " a.ehrlich@uni-leipzig.de, siebert@tropos.de",
    
    "abstract": "The tethered balloon system BELUGA (Balloon-bornE moduLar"+\
        " Utility for profilinG the lower Atmosphere) was operated during a"+\
        " collaborative measurement campaign at Station Nord (Greenland) in"+\
        " spring 2024 (2024-03-19 to 2024-04-18). With the BELUGA system, frequent"+\
        " profiling of the Arctic atmospheric boundary layer was performed to examine"+\
        " transitions between typical states of the Arctic atmospheric" +\
        " boundary layer. The balloon payload included a combined turbulence"+\
        " meteorological probe (TMP), containing a radiosonde-based extended" +\
        " meteorological package and a hot-wire package for" +\
        " turbulence, as well as a broadband radiation package (BP)." +\
        " The present dataset contains the processed Level-2" +\
        " meteorological data measured by the TMP"+\
        " specified for each of the 28 flights from 24 March to 12 April 2024."+\
        " This TMP data is published together with the BP measurement data.",
        
    "campaign":"Collaborative measurement campaign at Station Nord (Greenland) in"+\
        " spring 2024 (2024-03-19 to 2024-04-18)",
    
    "project": "Third Phase of ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes,"+\
        " and Feedback Mechanisms (AC)³"+\
        " (URI: https://gepris.dfg.de/gepris/projekt/268020496)."+\
        " Sub project A02: Balloon-borne observations and dedicated simulations"+\
        " of the transitions between typical states of the Arctic atmospheric"+\
        " boundary layer.", 
    
    "funding": "German Research Foundation (DFG) (URI: http://www.dfg.de/en/),"+\
        " grant/award no. 268020496: TRR 172: ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes, and Feedback"+\
        " Mechanisms (URI: https://gepris.dfg.de/gepris/projekt/268020496)",
    
    "source":"balloon-borne measurement data (BELUGA)" ,
    
    "flight_date":"change for your flight date",
    
    "comments": "Spike-removed meteorological data, with sonic temperature bias-corrected with"+\
        " sonde-based temperature, wind speed values corrected for pitch angle."+\
        " Flight segmentation categorising periods of specific flight"+\
        " maneuvers (e.g. profile ascents and descents) added."+\
        " The flight altitude is referenced to barometric altitude.",
    
    "processing_date" :"change to actual date of processing.",
    
    "processing_level": "PANGAEA data processing level 2",
    
    "further_details": "This dataset complenents previous Arctic ABL BELUGA"+\
        " observations described in Lonardi et al., 2022;"+\
        " https://doi.org/10.1525/elementa.2021.000120 within the scope of the"+\
        " MOSAIC campaign (Shupe et al., 2022;"+\
        " https://doi.org/10.1525/elementa.2021.00060)",
    
    "Conventions": "CF-1.12",
    
    "licence": "Creative Commons Attribution NonCommercialShareAlike 4.0"+\
                "International (CC BY-NC-SA 4.0)",
    
    "featureType":"TimeSeries",
    
    "keywords":"Atmospheric boundary layer (ABL); Arctic; BELUGA;"+\
        " Thermodynamics; Lapse-rate; Station Nord; Tethered balloon;"+\
        " vertical profile;" " in situ"}
    
TMP_turb_global_attributes={
    "title": "Tethered balloon-borne turbulence profile measurement data"+\
       " characterising the Arctic atmospheric boundary layer at Station Nord", 
    
    "authors":"Dorff, Henning; Schaefer, Michael; Siebert, Holger;" +\
        " Mueller, Joshua; Navale, Komal; Ehrlich, Andre;" +\
            " Wendisch, Manfred",
    
    "institutes":"(1) Leipzig Institute for Meteorology, Leipzig, Germany" +\
        " Leipzig University; "+"(2) Leibniz Institute for" +\
            " Tropospheric Research e.V., Leipzig, Germany",
    
    "contact":"henning.dorff@uni-leipzig.de,"+\
        " a.ehrlich@uni-leipzig.de, siebert@tropos.de",
    
    "abstract": "The tethered balloon system BELUGA (Balloon-bornE moduLar"+\
        " Utility for profilinG the lower Atmosphere) was operated during a"+\
        " collaborative measurement campaign at Station Nord (Greenland) in"+\
        " spring 2024 (2024-03-19 to 2024-04-18). With the BELUGA system, frequent"+\
        " profiling of the Arctic atmospheric boundary layer was performed to examine"+\
        " transitions between typical states of the Arctic atmospheric" +\
        " boundary layer. The balloon payload included a combined turbulence"+\
        " meteorological probe (TMP), containing a radiosonde-based extended" +\
        " meteorological package and a hot-wire anemometer package for" +\
        " turbulence, as well as a broadband radiation package (BP)." +\
        " The present dataset contains the processed" +\
        " Level-2 turbulence data measured by the TMP"+\
        " specified for each of the available flights from 24 March to 12 April 2024."+\
        " This TMP data is published together with the BP measurement data.",
        
    "campaign":"Collaborative measurement campaign at Station Nord (Greenland) in"+\
        " spring 2024 (2024-03-19 to 2024-04-18)",
    
    "project": "Third Phase of ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes,"+\
        " and Feedback Mechanisms (AC)³"+\
        " (URI: https://gepris.dfg.de/gepris/projekt/268020496)."+\
        " Sub project A02: Balloon-borne observations and dedicated simulations"+\
        " of the transitions between typical states of the Arctic atmospheric"+\
        " boundary layer.", 
    
    "funding": "German Research Foundation (DFG) (URI: http://www.dfg.de/en/),"+\
        " grant/award no. 268020496: TRR 172: ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes, and Feedback"+\
        " Mechanisms (URI: https://gepris.dfg.de/gepris/projekt/268020496)",
    
    "source":"balloon-borne measurement data (BELUGA)" ,
    
    "flight_date":"change for your flight date",
    
    "comments": "Hot-wire post-calibration.."+\
        " Flight segmentation categorising periods of specific flight"+\
        " maneuvers (e.g. profile ascents and descents) added."+\
        " The flight altitude is referenced to barometric altitude.",
    
    "processing_date" :"change to actual date of processing.",
    
    "processing_level": "PANGAEA data processing level 2",
    
    "further_details": "This dataset complenents previous BELUGA turbulence observations"+\
        " of the Arctic ABL described in Akansu et al., 2023;"+\
        " https://doi.org/10.1038/s41597-023-02582-5) within the scope of the"+\
        " MOSAIC campaign (Shupe et al., 2022; https://doi.org/10.1525/elementa.2021.00060)",
    
    "Conventions": "CF-1.12",
    
    "licence": "Creative Commons Attribution NonCommercialShareAlike 4.0"+\
                " International (CC BY-NC-SA 4.0)",
    
    "featureType":"TimeSeries",
    
    "keywords":"Atmospheric boundary layer (ABL); Arctic; BELUGA;"+\
        " Thermodynamics; Turbulence; Wind speed; Station Nord; Tethered balloon;"+\
        " vertical profiles; in situ; (AC)³"}
    
BP_global_attributes={
    "title": "Tethered balloon-borne broadband thermal-infrared"+\
        " irradiance profile measurement data"+\
       " characterising the Arctic atmospheric boundary layer at Station Nord", 
    
    "authors":"Dorff, Henning; Schaefer, Michael; Siebert, Holger;"+\
        " Mueller, Joshua; Navale, Komal; Wu, Fan; Ehrlich, Andre;"+\
            " Wendisch, Manfred",
    
    "institutes":" (1) Leipzig Institute for Meteorology, Leipzig, Germany" +\
        " Leipzig University; "+"(2) Leibniz Institute for" +\
            " Tropospheric Research e.V., Leipzig, Germany",
    
    "contact":"henning.dorff@uni-leipzig.de,"+\
        " a.ehrlich@uni-leipzig.de, siebert@tropos.de",
    
    "abstract": "The tethered balloon system BELUGA (Balloon-bornE moduLar"+\
        " Utility for profilinG the lower Atmosphere) was operated during a"+\
        " collaborative measurement campaign at Station Nord (Greenland) in"+\
        " spring 2024 (2024-03-18 to 2024-04-18). BELUGA performed frequent"+\
        " profiling of the Arctic atmospheric boundary layer to examine"+\
        " transitions between typical states of the Arctic atmospheric"+\
        " boundary layer. The balloon payload included a combined turbulence"+\
        " meteorological probe (TMP), containing a radiosonde-based " +\
        " meteorological package and a hot-wire anemometer package for" +\
        " turbulence, as well as a broadband radiation probe (BP)." +\
        " The present dataset contains the processed Level-2 up- and downward"+\
        " thermal-infrared irradiances measured by the BP specified for"+\
        " each of the 28 flights from 24 March to 12 April 2024. The BP data is"+\
        " published with the TMP measurement data.",
        
    "campaign":"Measurement campaign at Station Nord (Greenland) in"+\
        " spring 2024 (2024-03-19 to 2024-04-18)",
    
    "project": "Third Phase of ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes,"+\
        " and Feedback Mechanisms (AC)3" +\
        " (URI: https://gepris.dfg.de/gepris/projekt/268020496)."+\
        " Sub project A02: Balloon-borne observations and dedicated simulations"+\
        " of the transitions between typical states of the Arctic atmospheric"+\
        " boundary layer.", 
    
    "funding": "German Research Foundation (DFG) (URI: http://www.dfg.de/en/),"+\
        " grant/award no. 268020496: TRR 172: ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes, and Feedback"+\
        " Mechanisms (URI: https://gepris.dfg.de/gepris/projekt/268020496)",
    
    "source" : "balloon-borne measurement data (BELUGA)" ,
    
    "flight_date" : "change for your flight date",
    
    
    "comments": "Processing steps: Calibrated radiation data corrected for"+\
        " pyrgeometer"+" sensor temperature. Sensor inertia correction"+\
        " following (Ehrlich and"+" Wendisch, 2015; 10.5194/amt-8-3671-2015)."+\
        " Flight segmentation categorising periods of specific flight"+\
        " maneuvers (e.g. profile ascents and descents),"+\
        " with flight altitude referenced to barometric altitude.",
    
    "processing_date" :"change to actual date of processing.",
    
    "processing_level": "PANGAEA data processing level 2",
    "further_details": "This dataset complenents previous Arctic ABL BELUGA"+\
        " observations described in Pilz et al., 2023;"+\
        " https://doi.org/10.1038/s41597-023-02423-5 within the scope of the"+\
        " MOSAIC campaign (Shupe et al., 2022;"+\
        " https://doi.org/10.1525/elementa.2021.00060)",
    
    
    "Conventions" : "CF-1.12",
    
    "licence" : "Creative Commons Attribution NonCommercialShareAlike 4.0"+\
                "International (CC BY-NC-SA 4.0)",
    
    "featureType":"TimeSeries",
    
    "keywords":"Atmospheric boundary layer (ABL); Arctic; BELUGA;"+\
        " Broadband radiation; Station Nord; Radiation fluxes; Tethered balloon;"+\
        " vertical profile",
    ##Domain attribute(s):O2A Registry URI: \
    ##https://hdl.handle.net/10013/sensor.936bd700-3456-4d3b-a9fa-674ff87e7735
    }
    
Radiosonde_global_attributes={
    "title": "Daily radiosonde profile measurement data complementing "+\
    " balloon-borne measurements at Station Nord for characterising"+\
    " the Arctic atmospheric boundary layer", 
    
    "authors":"Dorff, Henning; Mueller, Joshua; Schaefer, Michael;"+\
        " Siebert, Holger; Ehrlich, Andre; Wendisch, Manfred",
    
    "institutes":" (1) Leipzig Institute for Meteorology, Leipzig, Germany" +\
        " Leipzig University; (2) Leibniz Institute for"+\
            " Tropospheric Research e.V., Leipzig, Germany",
    
    "contact":"henning.dorff@uni-leipzig.de,"+\
        " joshua.mueller@uni-leipzig.de, siebert@tropos.de",
    
    "abstract": "During the collaborative measurement campaign conducted at"+\
         " Station Nord (Greenland) in spring 2024 (2024-03-18 to 2024-04-18),"+\
         " tethered balloon-borne measurement flights were supported by daily"+\
         " radiosonde launches. The Level-2 (L2) radiosonde profiles provide"+\
         " supplementary information regarding the meteorological conditions"+\
         " during the measurement campaign. In particular, they sampled the"+\
         " vertical thermodynamic properties throughout the entire troposphere"+\
         " above the heights reachable for the balloon, whose purpose was"+\
         " profiling of the Arctic atmospheric boundary layer to examine"+\
         " transitions between its typical states .",
         
        
    "campaign":"Measurement campaign at Station Nord (Greenland) in "+\
        "spring 2024 (2024-03-19 to 2024-04-18)",
    
    "project": "Third Phase of ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes,"+\
        " and Feedback Mechanisms (AC)3" +\
        " (URI: https://gepris.dfg.de/gepris/projekt/268020496)."+\
        " Sub project A02: Balloon-borne observations and dedicated simulations"+\
        " of the transitions between typical states of the Arctic atmospheric"+\
        " boundary layer.", 
    
    "funding": "German Research Foundation (DFG) (URI: http://www.dfg.de/en/),"+\
        " grant/award no. 268020496: TRR 172: ArctiC Amplification:"+\
        " Climate Relevant Atmospheric and SurfaCe Processes, and Feedback"+\
        " Mechanisms (URI: https://gepris.dfg.de/gepris/projekt/268020496)",
    
    "source": "Radiosonde measurement data, sonde type: DFM-17 (GRAW)",
    
    "flight_date": "change for your flight date",
    
    
    "comments": "Processing steps: Spikes removed (cleaned of outliers)."+\
        " Wind speed values were corrected for radiosonde trajectory.",
    
    "processing_date":"change to actual date of processing.",
    
    "processing_level": "PANGAEA data processing level 2",
    
    "Conventions" : "CF-1.12",
    
    "licence" : "Creative Commons Attribution NonCommercialShareAlike 4.0"+\
                "International (CC BY-NC-SA 4.0)",
    
    "featureType":"trajectory",
    
    "keywords":"Arctic troposphere; Arctic; BELUGA; Radiosonde;"+\
        "Thermodynamics; Station Nord; Vertical Profile; In situ; Troposphere"}
    
#%% Turbulence_meteorological_probe_met variable attributes
TMP_met_var_attributes={
    #'p', 'z_b', 'T', 'sonic_T', 'rh', 'vv', 'dd'
    # Date/Time in UTC: based on seconds of day from radiosonde GPS
    "time":{
        "standard_name" :"time",
        "units":"seconds since 1970-01-01 00:00:00 UTC",
        "calendar":"standard",
        "long_name":"measurement time",
        "source":"GPS sensor"},
    #T: static temperature in °C  from the radiosonde (smoothed with a rolling median with 5 elements
    "T":{
        "standard_name":"air_temperature",
        "long_name":"air_temperature_radiosonde",
        "units":"K",
        "source":"Balloon-borne radiosonde"},
    # Caution for this
    "sonic_T":{
        "standard_name":"air_temperature",
        "long_name":"sensor_temperature_sonic_anemometer",
        "units":"K",
        "source":"Ultra-sonic anemometer",
        "comments": "corrected for a mean offset compared to sonde-based T"
        },
    #pp: static pressure in hPa from radiosonde
    "p":{
        "standard_name":"air_pressure",
        "long_name":"measured static air pressure",
        "units":"hPa",
        "source":"Balloon-borne radiosonde"},
	#RH: relative humidity in % from the radiosonde
    "rh":{
        "standard_name":"relative_humidity",
        "long_name":"measured relative humidity",
        "units":"%",
        "source":"Balloon-borne radiosonde"},
    
    # Zb: barometric height in m above ground level:
    # RL * (T_GND + 273.15) / g * np.log(p_GND / p)
    # Surface temperature T_GND is estimated as a mean value 
    # between ground and p_GND -0.2hPa (about 1.6 m)
    "z_b":{
        "standard_name" : "barometric_altitude",
        "long_name"     : "Pressure measurement based barometric altitude",
        "units"         : "m",
        "source"        : "Based on radiosonde pressure",
        "comments"  : "z_b=RL * (T_GND + 273.15) / g * np.log(p_GND / p)"
        },
    "vv":{
        "standard_name" : "wind_speed",
        "long_name"     : "horizontal wind speed",
        "units"         : "m s-1",
        "source"        : "ultrasonic anemometer",
        "comments"      : "vv=sqrt(u^2+v^2)"},
    
    "dd":{
        "standard_name" : "wind_from_direction",
        "long_name"     : "horizontal wind direction",
        "units"         : "degree",
        "source"        : "Honeywell 3 axis digital compass HMC 5883L",
        "comments"      : "Actually just orientation of the probe: "+\
            "When perfectly aligned into the wind, these values should be "+\
            "equal to wind direction. However, due to high latitude with "+\
                "predominant high magnetic declination and inclination, "+\
                    "this should be not over interpreted!"},
    
    "u_lat":{
        "standard_name":"y_wind",
        "long_name":"latitudinal sensor displacement",
        "units":"m s-1",},
        #"comments":"Should I add something?"},
    
    "u_lon":{
        "standard_name":"x_wind",
        "long_name":"longitudinal sensor displacement",
        "units":"m s-1"},
    
    "pitch":{
        "standard_name":"platform_pitch_angle",
        "long_name":"pitch angle of instrument probe",
        "units":"degrees"},
    
    "flight_segments":{
        "long_name"     : "flight segmentation categories",
        "units": "1",
        "comments":"defined categories:"+\
        " 0 = near ground, 2 = constant altitude,"+\
        " 100+ x = ascent + profile_no,"+\
        " -100- x = descent - profile_no,"+\
        " 3 = peak, 4 = max"
        },
        
    "quality_flag":{
        'long_name': 'Quality flag',
        'standard_name': 'status_flag',  # CF standard name
        'units': '1',  # No units for flags
        'flag_values': [0,1,2,3,4],  # Required for CF
        'flag_meaning': 'good_data '  # Use spaces to separate meanings
                        'ok_wind_pitch_affected '
                        'ok_interpolated '
                        'bad_data '
                        'missing_data '}
            }
#%% Turbulence_meteorological_probe_turb variable attributes
TMP_turb_var_attributes={
    #'p', 'z_b', 'T', 'sonic_T', 'rh', 'vv', 'dd'
    # Date/Time in UTC: based on seconds of day from radiosonde GPS
    "time":{
        "standard_name" :"time",
        "units":"seconds since 1970-01-01 00:00:00 UTC",
        "calendar":"standard",
        "long_name":"measurement time",
        "source":"GPS sensor"},
    #pp: static pressure in hPa from radiosonde
    "p":{
        "standard_name":"air_pressure",
        "long_name":"measured static air pressure",
        "units":"hPa",
        "source":"Balloon-borne MSR"},
	
    #T: static temperature in °C from the radiosonde
    "T":{
        "standard_name":"air_temperature",
        "long_name":"air_temperature",
        "units":"K",
        "source":"Balloon-borne MSR"},
    "E_U":{
        "long_name":"raw sensor voltage",
        "units": "V",
        "source":"Balloon-borne hot wire anemometer"},
    # Zb: barometric height in m above ground level:
    # RL * (T_GND + 273.15) / g * np.log(p_GND / p)
    # Surface temperature T_GND is estimated as a mean value 
    # between ground and p_GND -0.2hPa (about 1.6 m)
    "z_b":{
        "standard_name" : "barometric_altitude",
        "long_name"     : "Pressure measurement based barometric altitude",
        "units"         : "m",
        "source"        : "Based on MSR pressure",
        "comments"      : "z_b=RL * (T_GND + 273.15) / g * np.log(p_GND / p)"
        },
    "rho":{
        "standard_name" : "air_density",
        "long_name"     : "air density",
        "units"         : "kg m-3",
        "source"        : "MSR",
        },
    "U":{
        "standard_name" :"wind_speed",
        "long_name"     :"wind speed (50 Hz)",
        "units"         :"m s-1"},
    
    "flight_segments":{
        "long_name"     : "flight segmentation categories",
        "units": "1",
        "comments":"defined categories are: profile_ascent_xy, "+\
            "profile_descent_xy, peak, constant altitude, near ground"},
    "quality_flag":{
        'long_name': 'Quality flag',
        'standard_name': 'status_flag',  # CF standard name
        'units': '1',  # No units for flags
        'flag_values': [0,1,2,3],  # Required for CF
        'flag_meaning': 'good_data '  # Use spaces to separate meanings
                        'ok_interpolated_data '
                        'bad_data '
                        'missing_data '}
        }
    #"quality_flag":{
    #    "long_name":"data quality flag",
    #    "units": "1",
    #    "comments":"major quality classes are good, ok and bad, 
        #with specifications given e.g. ok_interpolated.",}
    

    #%% Broadband probe variable attributes
BP_var_attributes={
    "time":{
        "standard_name" :"time",
        "units":"seconds since 1970-01-01 00:00:00 UTC",
        "calendar":"standard",
        "long_name":"measurement time"},
    "F_up":{
        "standard_name":"upwelling_longwave_flux_in_air",
        "long_name":"upward terrestrial (long-wave) radiation (irradiance)",
        "units":"W m-2",
        "accuracy": "5 W m-2",
        "source":"Kipp & Zonen CGR4 pyrgeometer"},
    
    "F_down":{
        "standard_name":"downwelling_longwave_flux_in_air",
        "long_name":"downward terrestrial (long-wave) radiation (irradiance)",
        "units":"W m-2",
        "source":"Kipp & Zonen CGR4 pyrgeometer",
        "accuracy": "5 W m-2"},
    
    "F_net":{
        "standard_name":"net_downward_longwave_flux_in_air",
        "long_name":"net terrestrial (long-wave) radiation. ",
        "units":"W m-2",
        "comments":"F_net=F_down-F_up, positive if there is downward"+\
            " longwave flux as F_down and F_up values are positively defined.",
        "accuracy": " 7 W m-2"},
	
    "bme_temperature":{
        "standard_name":"air_temperature",
        "long_name":"measured air temperature",
        "units":"K",
        "source":"BME sensor"},
        
    "bme_pressure":{
        "standard_name":"air_pressure",
        "long_name":"measured static air pressure",
        "units":"hPa",
        "source":"BME sensor"},
	
    "bme_rh":{
        "standard_name":"relative_humidity",
        "long_name":"measured relative humidity",
        "units":"%",
        "source":"BME sensor"},
    
    "roll":{
        "standard_name":"platform_roll",
        "long_name":"roll angle sensor platform",
        "units":"degree",
        "source":"DOF altitude sensor"},
	
    "pitch":{
        "standard_name":"platform_pitch",
        "long_name":"pitch angle sensor platform",
        "units": "degree",
        "source":"DOF altitude sensor"},
    
    "yaw":{
        "standard_name":"platform_yaw",
        "long_name":"yaw angle sensor platform",
        "units":"degree",
        "source":"DOF altitude sensor"},
    
    "dof_pressure":{
        "standard_name" : "air_pressure",
        "long_name"     : "measured static air pressure",
        "source"        : "DOF altitude sensor",
        "units"         : "hPa"},
	
    "z_b":{
        "standard_name" : "barometric_altitude",
        "long_name"     : "Pressure measurement based barometric altitude",
        "units"         : "m",
        "instrument"    : "BME sensor"
        },
    
    "flight_segments":{
        "long_name"     : "flight segmentation categories",
        "units": "1",
        "comments":"defined categories:"+\
        " 0 = near ground, 1 = constant altitude,"+\
        " 100+ x = ascent + profile_no,"+\
        " -100- x = descent - profile_no,"+\
        " 3 = peak, 4 = max"},
        #"comments":"defined categories are: profile_ascent_xy, "+\
        #    "profile_descent_xy peak, constant altitude, near ground"},
    
    "quality_flag":{
        'long_name': 'Quality flag',
        'standard_name': 'status_flag',  # CF standard name
        'units': '1',  # No units for flags
        'flag_values': [0,1,2,3,4],  # Required for CF
        'flag_meaning': 'good_data '  # Use spaces to separate meanings
                        'ok_interpolated_data '
                        'ok_data_absolute_pitch_high_but_below_20_deg '
                        'bad_data_absolute_pitch_above_20_deg '
                        'missing_data '}
        }
#%% Radiosonde variable attributes
Radiosonde_var_attributes={
        "time":{
            "standard_name" :"time",
            "units":"seconds since 1970-01-01 00:00:00 UTC",
            "calendar":"standard",
            "long_name":"measurement time"},
                
        'T': {"standard_name":"air_temperature",
              "long_name":'Air Temperature',
              "units": 'K',
              },
        "time_since_launch":{
            "long_name":"time duration since radiosonde launch",
            "units":"seconds"},
        
        'RH': {'standard_name':'relative_humidity',
               'long_name':'Relative Humidity',
               'units':'%'},
        
        "p": {"standard_name":"air_pressure",
              "long_name":"Static Air Pressure",
              "units":"hPa"},
        
        'Tu_Cap': {'long_name':'Capacitor Temperature Sensor',
                   'units':'K'},
        
        'altitude': {"standard_name":'altitude',
                     "long_name":'Geometric altitude above mean sea level',
                     "units":'m'},
        
        'latitude':{'standard_name':'latitude',
                    'long_name': 'Latitude',
                    'units':'degrees_north'},
        
        'longitude':{'standard_name': 'longitude',
                     'long_name': 'Longitude',
                     'units':'degrees_east'},
        
        'horizontal_speed':{'standard_name':"wind_speed",
                'long_name':'Horizontal Speed of the Radiosonde',
                'units':'m s-1'},
        
        "wspeed":{"standard_name":"wind_speed",
                  "long_name":"Horizontal wind speed",
                  "units":'m s-1'},
        
        "wdir":{"standard_name":"wind_from_direction",
                "long_name":"Wind direction",
                "units": "degrees"},
        
        'vertical_speed':{
            'long_name':'Vertical Speed of the Radiosonde',
            'units':'m s-1'},
        
        'distance':{'long_name':'Slant range distance from launch point',
                    'units':'m'},
        
        'u':{'standard_name':"eastward_wind",
             'long_name': "zonal wind component",
             'units':"m s-1"},
        
        'v':{'standard_name':"northward_wind",
             'long_name':'meridional wind component',
             'units':"m s-1"}
        }