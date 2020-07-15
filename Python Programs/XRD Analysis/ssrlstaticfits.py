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
perov_import = csv_to_np('/Users/rbelisle/Desktop/SSRL_data/startend_full.csv')
#%% Covert to Q
q = two_to_q(perov_import[:,0],0.9763)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 
#%% #Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.92
q_2 = 2.16
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_substract(q_sub,perov_sub[:,file],10) #remove background from that file
#%% Do Curve Fitting For One Peak
image = 0 #index of file you want to look at
p0 = [100, 2.04, .001] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [500, 2.2, 5]
lower_limit = [0, 1.8, 0]
popt,pcov = curve_fit(gaussian, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,gaussian(q_sub, *popt),'c--', label='Model') #plot best fit
#%% Do Curve Fitting for Three Peaks
image = 0 #index of file you want to look at
init_q = 2.08389035 #q value from initial structure
p0 = [5, 2.0, .01, 10, init_q, .001, 10, 2.10, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [20, init_q, .3, 20, init_q+.0001, .3, 20, 2.15, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, 2.00, 0, 0, init_q-0.0001, 0, 0, init_q, 0]
popt,pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
# %% Plot everything
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
