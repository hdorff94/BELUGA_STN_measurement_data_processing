# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 17:04:50 2025

@author: u300737
"""
config_dict={
    "version_number":"_v2.5",
    "BP_flights_to_process":["RF01","RF02",
        #"RF03","RF04","RF05","RF06",
        #"RF07","RF08","RF10","RF11",
        #"RF12","RF13","RF14","RF15",
        #"RF16","RF17","RF18","RF19",
        #"RF20","RF21","RF22","RF23",
        #"RF24",
        #"RF25","RF26","RF27",#"RF28",
        #"RF09" no available data
        ],
    "TMP_met_flights_to_process":["RF01","RF02",
        #"RF03",#"RF04","RF05","RF06","RF07","RF08",
        #"RF10","RF11","RF12","RF13",
        #"RF14","RF15","RF16","RF17",
        #"RF18","RF19","RF20","RF21",
        #"RF22","RF23","RF24","RF25",
        #"RF26","RF27"
        #"RF09" no available data 
        ],
    "TMP_turb_flights_to_process":["RF13",
        #"RF14","RF15","RF16","RF17","RF18",
        #"RF20","RF21","RF22","RF23","RF24",
        #"RF26","RF27",
        ],
    "radiosonde_flights_to_process":[
        "RFS1",
        #"RFS2","RFS3","RFS4",
        "RF01","RF02","RF03",
        "RF05","RF06", 
        "RF07",
        "RF08", "RF09", "RF10","RF11","RF12",
        "RF13","RF14","RF15","RF16","RF17","RF18",
        "RF19","RF20","RF21","RF22","RF23","RF24",
        "RF25","RF26","RF27"
        ],
    
    # Instrument packages to be processed
    "process_BP"        : False,
    "process_TMP_met"   : False,
    "process_TMP_turb"  : True,
    "process_radiosonde": False,
    # processing level (L1 in csv, L2 in nc)
    "run_l1_processing" : True,
    "run_l2_processing" : True,
    # store data in directory for publication
    "final_processing"  : True,
    # Quicklooks
    "plot_processing"   : True,
    }



