# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 18:14:49 2025

@author: u300737
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

main_path=os.getcwd()#+"/BELUGA_data/"
sys.path.insert(1,main_path+"/src/analysing/")
sys.path.insert(2,main_path+"/plotting/")

import Measurement_Platforms as Platforms

# Open and allocate classes
Platforms_STN_cls = Platforms.Measurement_Platforms_STN(
                        main_path=main_path+"/BELUGA_data/")
BELUGA_cls        = Platforms.BELUGA(Platforms_STN_cls)
BELUGA_cls.get_flight_information()

# Run processing
BP_flights_to_process=[#"RF01","RF02","RF03","RF04",
                       #"RF05","RF06",
                       #"RF07","RF08",
                       #"RF09",
                       #"RF10","RF11","RF12","RF13","RF14",
                       #"RF15","RF16","RF17","RF18",
                       "RF19",
                       "RF20","RF21","RF22","RF23","RF24",
                       "RF25","RF26","RF27",
                       "RF28",
                    ]
                    
TMP_met_flights_to_process=["RF01","RF02","RF03","RF04",
                            "RF05","RF06",
                            "RF07","RF08",
                            "RF10","RF11","RF12","RF13",
                            "RF14","RF15","RF16","RF17",
                            "RF18","RF19","RF20","RF21",
                            "RF22","RF23","RF24","RF25",
                            "RF26","RF27",
                            # not able to process
                            #"RF28",
                            #"RF09",
                            ]

TMP_turb_flights_to_process=["RF13","RF14","RF15",
                             "RF16","RF17","RF18",
                             "RF20","RF21","RF22",
                             "RF23","RF24","RF26",
                             "RF27"
                             ]

radiosonde_flights_to_process=[
    "RFS1","RFS2","RFS3",
    "RFS4","RF01","RF02","RF03",
    "RF05","RF06", 
    "RF07",
    "RF08", "RF09", "RF10","RF11","RF12",
    "RF13","RF14","RF15","RF16","RF17","RF18",
    "RF19","RF20","RF21","RF22","RF23","RF24",
    "RF25","RF26","RF27"
    ]


process_BP          = True
process_TMP_met     = False
process_TMP_turb    = False
process_radiosonde  = False

run_l1_processing   = True
run_l2_processing   = True
plot_processing     = True

version_number="_v2.5"
final_processing    = True
if final_processing:
    BELUGA_cls.temporary=False
if process_BP:
    for rf in BP_flights_to_process:
        print("BP Process ",rf)
        # Allocate class and run processing
        BELUGA_BP         = Platforms.Broadband_Probe(BELUGA_cls,rf=rf,
                                run_L1_processing = run_l1_processing,
                                run_L2_processing = run_l2_processing,
                                plot_processing   = plot_processing)
        
        BELUGA_BP.run_processing(rf=rf,version_number=version_number)
        
if process_TMP_met:
    for rf in TMP_met_flights_to_process:
        print("TMP_met processing of:", rf)
        # Allocate class and run processing
        BELUGA_TMP_met    = Platforms.Meteorological_Probe(BELUGA_cls,rf=rf,
                                run_L1_processing = run_l1_processing,
                                run_L2_processing = run_l2_processing,
                                plot_processing   = plot_processing)

        BELUGA_TMP_met.run_processing(rf=rf,version_number=version_number)

if process_TMP_turb:
    for rf in TMP_turb_flights_to_process:
        print("TMP_turb processing of:",rf)
        BELUGA_TMP_turb = Platforms.Turbulence_Probe(BELUGA_cls,rf=rf,
                            run_L1_processing = False,
                            run_L2_processing = run_l2_processing,
                            plot_processing   = plot_processing)
        BELUGA_TMP_turb.run_processing(rf=rf,version_number=version_number)

if process_radiosonde:
    for rf in radiosonde_flights_to_process:
        print(rf)
        Radiosonde_STN    = Platforms.Radiosondes(Platforms_STN_cls,rf=rf,
                            run_L1_processing = run_l1_processing,
                            run_L2_processing = run_l2_processing,
                            plot_processing   = plot_processing)
        Radiosonde_STN.temporary=False
        Radiosonde_STN.run_processing(version_number,
                           run_L1_processing=True,
                           run_L2_processing=True)
        