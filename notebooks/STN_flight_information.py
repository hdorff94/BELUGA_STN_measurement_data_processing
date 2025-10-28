import os
import pandas as pd

def get_flight_information(path=os.getcwd()):
    """
    Returns the flight infos for all BELUGA STN flights as dataframe.
    Information is equivalent to Table 1 of ESSD manuscript.
    """
    STN_flight_table_path=path+"/../BELUGA_data/"
    flight_infos=pd.read_csv(STN_flight_table_path+"STN_BELUGA_flight_infos.csv",index_col="Flight No")
    flight_infos["Start Time"]     = pd.to_datetime(flight_infos["Start Time"], dayfirst=True)
    flight_infos["End Time"]       = pd.to_datetime(flight_infos["End Time"], dayfirst=True)
    flight_infos["type"]           = "uncategorized"
    flight_infos["type"].loc[["RF07","RF08","RF10","RF15",
            "RF16","RF19","RF23","RF27"]]   = "clear to cloudy"
    flight_infos["type"].loc[["RF05","RF06"]]    = "Polar night to day"
    flight_infos["type"].loc[["RF12","RF13","RF14","RF26"]]="LLJ"
    #----------------------------------------------------------------------#
    # flight
    # self.LLJ_flights    = self.flight_infos.loc[\
    #                            self.flight_infos["type"]=="LLJ"]
    # self.cloud_flights  = self.flight_infos.loc[\
    #                            self.flight_infos["type"]=="clear to cloudy"]
    # self.night_flights  = self.flight_infos.loc[\
    #                            self.flight_infos["type"]=="Polar night to day"]
    return flight_infos    