# BELUGA_STN_measurement_data_processing
This repository contains all routines for the post-processing of the Level-0 (L0) to Level-2 (L2) data measured by the Balloon-bornE moduLar Utility for profilinG the lower Atmosphere (BELUGA) during spring 2024 at a dedicated measurement campaign at Station Nord (Greenland).
In specific, this is a python-based object-oriented processing of the four instrument packages, which are the two components of the Turbulent meteorological probe (TMP_met & TMP_turb), the broadband radiation package, and supplementary daily radiosonde profiles.
The emerging L2-data is published and accessible via the data collection: .

For each instrument component, inhereting subclasses consist of all relevant routines. Level-0 (L0) to Level-1 (L1) processing comprises a sub-selection of the relevant variables and
first quality checks/post-calibrations. L1-data is provided as csv. The L2-data is given as netCDF files, with additional processing with a quality flag,
 the barometric altitude as common altitude reference, a flight segmentation, and additional variables and attributes improving the emerging scientific analysis and data handling.
During processing, logger files are created to deliver information of all conducted processing steps.

### Meteorological data of the TMP_met

### Thermal-infrared irradiances based on the broadband radiation probe (BP)

## Executing the processing
The processing can be executed by the script
```python run_processing.py ```
The instruments and research flights to be chosen for post-processing have to be set manually in the config-file "STN_processing_config"
and is then read by the processing routine

This version of the processing code is a prelimary version.