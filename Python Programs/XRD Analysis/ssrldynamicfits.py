#this is a base file for doing static fits (i.e. perovskite structure at the end of an experiment)
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd

#%% Import Data
from xrdfunctions import *
perov_import = csv_to_np('/Users/rbelisle/Desktop/SSRL_data/B1_xraydeg.csv')

#%% Covert to Q
q = two_to_q(perov_import[:,0],0.9763)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 


#%% #Trim and Remove Background
miller = [1, 0, 0] #peak you are trying to capture
q_1 = 0.84
q_2 = .96
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)

#%% #Loop through Data
files = num_files(perov)
time = np.zeros(files)
lattice = np.zeros(files)
intensity = np.zeros(files)
files = num_files(perov)
p0 = [10,0.9, 0.01] #best guess for number of peaks an initial values
for frame in range(files): 
    perov_to_fit = back_subtract(q_sub,perov_sub[:,frame],10)
    popt, pcov = curve_fit(gaussian, q_sub, perov_to_fit, p0)
    intensity[frame] = popt[0]
    lattice[frame] = q_to_a(popt[1],miller)
    time[frame] = frame*20+20
    print('Intensity:', popt[0])
    print('Lattice Spacing:', q_to_a(popt[1],miller))
    p0=popt

#%% Calculate and share lattice spacging
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(time,lattice,'ro--') #plot subfield of data
plt.xlabel('Time [s]',size=14) #Define x-axis label
plt.ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.axis([0, 600, 6.14, 6.145])
plt.savefig('MAPbI2Br_Xray_Lattice.png')

# %%