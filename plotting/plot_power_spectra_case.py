# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 21:37:50 2025

@author: u300737
"""
import os
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib import ticker
def PSD(x,fs):
    xF = np.fft.fft(x)       # estimates the fft of array a
    N = len(xF)
    time_step = 1/fs
    f = np.fft.fftfreq(N, time_step)
    f = f[1:int(N/2)]
    df = f[2] - f[1]
    PSD     = 2 * (abs(xF[1:int(N/2)])/ N)**2 /df
    return(f,PSD)
    

def PSD_Smooth(f1,t_bins,fs,f,PSD):    
    start_log = np.log10(f1)
    stop_log  = np.log10(fs/2)
    X     = f
    Y     = PSD
     
    bins = np.logspace(start_log,stop_log, 
                       num = t_bins, endpoint = True)   
    # define bin boundaries
    idx  = np.digitize(X,bins)                                                  
    #    sort X-values into bins
    Y  = [np.average(Y[idx==k]) for k in range(t_bins)]                   
    # average Y-values for each bin
    bins = bins[1:t_bins]                                                   
    # remove 1st point - it's a NaN
    Y        = Y[1:t_bins]                                                      
    # remove 1st point - it's a NaN  
    return(bins,Y)


flight="RF22"
date="20240407"
fname="BELUGA_STN_L2_TMP_turb_"+flight+"_"+date+"_v2.5.nc"
data_path=os.getcwd()+"/../BELUGA_data/BELUGA_TMP_turb_probe/"+\
    "final_data//"
file=data_path+fname
# load level 2 turbulence data from RF
turb_ds=xr.open_dataset(file)
#%% Plotting power spectrum
MSR=turb_ds["U"].to_series()
MSR_s = MSR["2024-04-07 15:02:01.000" :"2024-04-07 15:22:01.000"] # for RF22 only - select 20 minutes at almost constant height of 10 m above surface

fs = 50             # sampling frequency in Hz
start_freq = 0.001  # start frequency for smoothing
n_bins = 10         # bins per decade

MSR_s  = MSR_s.dropna()#
ps     = PSD(MSR_s,fs)  # calculatzing the power spectral density (normalized such as the integral will provide the variance)
sps    = PSD_Smooth(start_freq,n_bins,fs,ps[0],ps[1]) # smooth spectral density over logarithmic equidistant bins with n_bins per decade
 
fig, ax = plt.subplots(figsize=(12,8),nrows = 1)
plt.rcParams['font.size'] = 22

title_str =  (flight + "  " +  date + " -- " + "15:02 - 15:22 UTC in 10 m AGL")

ax.loglog(ps [0] ,ps[1],label=" raw spectrum",lw=2,color="mediumaquamarine")
ax.loglog(sps[0],sps[1] ,'-',lw=7,color="white")
ax.loglog(sps[0],sps[1] ,'o--',label=" smoothed spectrum",lw=4,color="darkgreen",zorder=2)
ax.loglog(sps[0],sps[0]**(-5/3.)/1e3,label=" -5/3 - inertial subrange model",lw=3,color="saddlebrown")
ax.set_ylabel(r'Power spectral density ($\rm m^2 s^{-2} Hz^{-1}$)')
ax.set_xlabel(r'Frequency  [$\rm   Hz$]')
#ax.set_title("Power spectrum of calibrated hot-wire data " + "\n" + title_str )
# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Make axes bolder
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.xaxis.set_tick_params(width=2)
ax.yaxis.set_tick_params(width=2)

# Enable minor ticks on log scale
ax.minorticks_on()
# Optional: set minor tick formatter for y-axis if needed
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs='auto', numticks=10))
ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))
ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs='auto', numticks=10))

# Set custom length and linewidth for major and minor ticks
ax.tick_params(which='major', length=10, width=2)
ax.tick_params(which='minor', length=5, width=2)

plt.legend()
plt.show()
plot_path=os.getcwd()+"/../plots/"
fig_name="fig07_power_spectra_TMP_turb_"+flight+".pdf"
fig.savefig(plot_path+fig_name,dpi=300,bbox_inches="tight")
print("Figure saved as:",plot_path+fig_name)