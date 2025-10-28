# -*- coding: utf-8 -*-
"""
This is the running script to perform the Level-0 to Level-2 
processing of the BELUGA and radiosonde measurement data gathered
during the comprehensive measurement campaign at Station Nord in 
spring 2024.

Four instrument specfic L2- datasets can be obtained for the 
BELUGA-borne Turbulent Meteorological Probe (with two components for 
meteorological and turbulent data), for the BELUGA-borne Broadband 
Radiation package (BP), and for the tropospheric profiles recorded by 
the daily radiosondes.
The research flights and instrument components to be processed are 
definded based on switches in the processing config file.
    
Created on Tue Jun 24 18:14:49 2025

@author: Henning Dorff

"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

main_path=os.getcwd()
sys.path.insert(1,main_path+"/src/")
sys.path.insert(2,main_path+"/plotting/")

import Measurement_Platforms as Platforms

# Open and allocate classes
Platforms_STN_cls = Platforms.Measurement_Platforms_STN(
                        main_path=main_path+"/BELUGA_data/")
BELUGA_cls        = Platforms.BELUGA(Platforms_STN_cls)
BELUGA_cls.get_flight_information()

# get the switches and processing definitions from 
# the processing_config file
import STN_processing_config
config_dict=STN_processing_config.config_dict
if config_dict["final_processing"]:
    BELUGA_cls.temporary=False

# Run processing
if config_dict["process_BP"]:
    for rf in config_dict["BP_flights_to_process"]:
        print("BP Process ",rf)
        # Allocate class and run processing
        BELUGA_BP         = Platforms.Broadband_Probe(BELUGA_cls,rf=rf,
            run_L1_processing = config_dict["run_l1_processing"],
            run_L2_processing = config_dict["run_l2_processing"],
            plot_processing   = config_dict["plot_processing"])
        BELUGA_BP.run_processing(rf=rf,
            version_number=config_dict["version_number"])
        
if config_dict["process_TMP_met"]:
    for rf in config_dict["TMP_met_flights_to_process"]:
        print("TMP_met processing of:", rf)
        # Allocate class and run processing
        BELUGA_TMP_met    = Platforms.Meteorological_Probe(BELUGA_cls,rf=rf,
            run_L1_processing = config_dict["run_l1_processing"],
            run_L2_processing = config_dict["run_l2_processing"],
            plot_processing   = config_dict["plot_processing"])

        BELUGA_TMP_met.run_processing(rf=rf,
            version_number = config_dict["version_number"])

if config_dict["process_TMP_turb"]:
    for rf in config_dict["TMP_turb_flights_to_process"]:
        print("TMP_turb processing of:",rf)
        BELUGA_TMP_turb = Platforms.Turbulence_Probe(
            BELUGA_cls,rf=rf,
            run_L1_processing = False,
            run_L2_processing = config_dict["run_l2_processing"],
            plot_processing   = config_dict["plot_processing"])
        
        BELUGA_TMP_turb.run_processing(rf=rf,
            version_number=config_dict["version_number"])

if config_dict["process_radiosonde"]:
    for rf in config_dict["radiosonde_flights_to_process"]:
        print(rf)
        Radiosonde_STN    = Platforms.Radiosondes(
            Platforms_STN_cls,rf=rf,
            run_L1_processing = config_dict["run_l1_processing"],
            run_L2_processing = config_dict["run_l2_processing"],
            plot_processing   = config_dict["plot_processing"])
        Radiosonde_STN.temporary=False
        Radiosonde_STN.run_processing(
            config_dict["version_number"],
            run_L1_processing=config_dict["run_l1_processing"],
            run_L2_processing=config_dict["run_l2_processing"])
        