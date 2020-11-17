#this is a base file for doing dynamics c fits (i.e. perovskite structure at the end of an experiment)
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd

#%% Import Data
from xrdfunctions import *
perov_import = csv_to_np('/Users/rbelisle/Desktop/Xraydeg/C1_MAPbBr50_Xraydeg/c1_deg.csv')

#%% Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

#%% #Trim and Remove Background
miller = [1, 0, 0] #peak you are trying to capture
q_1 = 0.8
q_2 = 1.15
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


#%% #Loop through Data
files = num_files(perov) #determine the number of images I have
time = np.zeros(files)
lattice1 = np.zeros(files)
lattice2 = np.zeros(files)
lattice3 = np.zeros(files)
intensity = np.zeros(files) #generate empty arrays for the parameters I'll be capturing
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
for frame in range(files): 
    popt, pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,frame+1], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
    intensity[frame] = popt[0]
    lattice1[frame] = q_to_a(popt[1],miller)
    lattice2[frame] = q_to_a(popt[4],miller)
    lattice3[frame] = q_to_a(popt[7],miller)
    time[frame] = frame*20+20
    print('Intensity:', popt[0])
    print('Lattice Spacing:', q_to_a(popt[1],miller))
    p0=popt
    upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
    lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]

#%% Plot lattice spacings over time
plt.figure(figsize=(5,4)) #make plot larger
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.plot(time,lattice1,'r.')
plt.plot(time,lattice2, 'k.')
plt.plot(time, lattice3, 'b.')
# %% Cell for XRay degredation analysis should make this a function 
files = num_files(perov)
lead = time = np.zeros(files)
for frame in range(files): 
    lead[frame] = perov_fit[80,frame]/perov_fit[35,frame]
    time[frame] = frame*20
plt.figure(figsize=(5,4)) #make plot larger
plt.xlabel('X-ray Exposure [s]',size=14) #Define x-axis label
plt.ylabel('Perovskite to Lead Iodide Ratio',size=14) #Define x-axis label
plt.plot(time, lead,'ko')
plt.ylim(10, 15)

# %% Test Cell: Not for permanent 
def perovtopbi2(q, intensity):
    #array is a 1D vector of  q values
    #intensity is a 1D vector of perovskite intensities
    leadiodide_q = 0.9 #rough q of lead iodide peak 
    pad = 10 #number of points around leadiodide_q to look for true max
    peak = find_nearest(q,leadiodide_q)
    leadiodide_inensity = max(intensity[peak-pad:peak+pad])
    ratio = max(intensity)/leadiodide_inensity
    return ratio



# %%
perovtopbi2(q_sub,perov_fit[:,0])
# %%
files = num_files(perov)
time = np.zeros(files)
lead = np.zeros(files)
# time = np.zeros(files)
print(lead)
for frame in range(files): 
    time[frame] = frame*20
    lead[frame] = perovtopbi2(q_sub, perov_fit[:,frame])

print(lead)
print(time)

# %%
