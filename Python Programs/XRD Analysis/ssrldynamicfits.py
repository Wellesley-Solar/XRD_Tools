#this is a base file for doing dynamics c fits (i.e. perovskite structure at the end of an experiment)
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from xrdfunctions import *

#%% Import Data
perov_import = csv_to_np('/Users/rbelisle/Desktop/5050onoff/lighton.csv')

#%% Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

#%% #Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.94
q_2 = 2.18
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

#%% Do Curve Fitting For One Peak
image = 0 #index of file you want to look at
p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [1, 3000, q_2, 5]
lower_limit = [0, 0, q_1, 0]
popt,pcov = curve_fit(pvoigt, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,pvoigt(q_sub, *popt),'c--', label='Model') #plot best fit
init_q = popt[2]

#%% Do Curve Fitting for Three Peaks TEST CELL (run this before running full loop to make sure nothing weird is happening)
image = 1 #index of file you want to look at
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
popt,pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
# Plot everything
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,three_gaussians(q_sub, *popt),'c--', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'b-.', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],popt[4],popt[5]), 'b--', label='Peak 2')
plt.plot(q_sub, gaussian(q_sub, popt[6],popt[7],popt[8]), 'b:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
print('Intensity:', popt[0], popt[3], popt[6])
print('Lattice Spacing:', q_to_a(popt[1],miller),q_to_a(popt[4],miller), q_to_a(popt[7],miller))


#%% #Loop through light on Data
#min_lattice_change = 1 
#offset = 2*math.pi/(min_lattice_change) #sets a maximum change in q in 2 minute interval
files = num_files(perov)-1 #determine the number of images I have
time = np.zeros(files)
lattice1 = np.zeros(files)
lattice2 = np.zeros(files)
lattice3 = np.zeros(files)
phase1 = np.zeros(files) #generate empty arrays for the parameters I'll be capturing
phase2 = np.zeros(files)
phase3 = np.zeros(files)
chem1 = np.zeros(files)
chem2 = np.zeros(files)
chem3 = np.zeros(files)
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
for frame in range(files): 
    popt, pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,frame+1], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
    intensity[frame] = popt[0]
    lattice1[frame] = q_to_a(popt[1],miller)
    lattice2[frame] = q_to_a(popt[4],miller)
    lattice3[frame] = q_to_a(popt[7],miller)
    phase1[frame] = sum(gaussian(q_sub,*popt[0:3]))
    phase2[frame] = sum(gaussian(q_sub,*popt[3:6]))
    phase3[frame] = sum(gaussian(q_sub,*popt[6:]))
    chem1[frame] = q_to_chem(popt[1],miller)
    chem2[frame] = q_to_chem(popt[4],miller)
    chem3[frame] = q_to_chem(popt[7],miller)
    time[frame] = frame*20+20
    print('Intensity:', popt[0])
    print('Lattice Spacing:', q_to_a(popt[1],miller))
    p0=popt
    upper_limit = [30, popt[1], .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
    lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, popt[7], 0]

#%% Loop through Dark
#p0 = popt 
#bounds assume opposite direction change (lower q peak gets larger, higher peak gets smaller)
#%% Plot lattice spacings over time
plt.figure(figsize=(5,4)) #make plot larger
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.plot(time,lattice1,'r.')
plt.plot(time,lattice2, 'k.')
plt.plot(time, lattice3, 'b.')

# %%cell to look at intensity changes over time
# have area of phases already
#plt.plot(time, (phase1+phase2+phase3)/(phase1[0]+phase2[0]+phase3[0]))
plt.plot(time, -1*(chem1*phase1+chem2*phase2+chem3*phase3)/(phase1+phase2+phase3))



# %%
