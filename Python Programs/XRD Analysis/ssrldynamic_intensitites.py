#%% The program imports a stacked .csv of dynamic scattering data and outputs 
# a list of integrated intensities

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from xrdfunctions import *

# %% Create dictionary for storing analyzed results - to be run at the start of a new session
samplelist = []
results = {"Sample":"Results"}

#%% Import Data
samplename = 'B1_lighton'
perov_import = csv_to_np('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/' + samplename + '.csv')
light_on = 1 # if light on set to 1 else set to 0
interval = 2 # time interval between frames in minutes

# Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

# Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.92 #lower limit in q
q_2 = 2.18 # upper limit in q 
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

# Initialize parameters of interest
files = num_files(perov)
time = np.arange(0, files*interval, interval) + interval/2 #time interval of the experiment

#  Do Curve Fitting of Initial Frame (assumes if light is on, new composition)
    if light_on:
    image = 0 #index of file you want to look at
    p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
    upper_limit = [1, 3000, q_2, 5]
    lower_limit = [0, 0, q_1, 0]
    popt,pcov = curve_fit(pvoigt, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
    plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
    plt.plot(q_sub,pvoigt(q_sub, *popt),'c--', label='Model') #plot best fit
    init_q = popt[2] #this is the inital peak position that will be used as a reference
    FWHM = popt[3] 
    

# %% Dynamics of Process
#  Plot Integrated Intensity over time
totalintensity = np.zeros(files)
for frame in range(files):
    totalintensity[frame] = sum(perov_fit[:,frame])
    
plt.plot(time, totalintensity/totalintensity[0], 'r--')
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Integrated Intensity [Normalized]',size=14)#Define y-axis label

# Plot Binned Intenisty
iod_bin = np.zeros(files)
orig_bin = np.zeros(files)
brom_bin = np.zeros(files)

iod_lim = find_nearest(q_sub, init_q-FWHM)
brom_lim = find_nearest(q_sub, init_q+FWHM)

for frame in range(files):
    iod_bin[frame] = sum(perov_fit[0:iod_lim,frame])
    orig_bin[frame] = sum(perov_fit[iod_lim:brom_lim,frame])
    brom_bin[frame] = sum(perov_fit[brom_lim:,frame])

plt.plot(time, iod_bin, 'r--')
plt.plot(time, orig_bin, 'k--')
plt.plot(time, brom_bin, 'b--')
plt.plot(time, brom_bin+iod_bin+orig_bin, 'k-')

plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Binned Intensity [a.u.]',size=14)#Define y-axis label

# %% Save Variables
samplelist.append(samplename) #creates a list of the samples you've looked at this session
datatostore = [time, totalintensity, iod_bin, orig_bin, brom_bin] #creates a list of relevant data
results[samplename] = datatostore #adds an item to a dictionary where the sample name is linked to these results
print(samplelist)