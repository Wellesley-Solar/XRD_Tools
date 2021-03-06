#this is a base file for doing static fits (i.e. perovskite structure at the end of an experiment)
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
import math as math

#%% Import Data
from xrdfunctions import csv_to_np, two_to_q, trim_data, num_files, back_subtract, pvoigt, gaussian, three_gaussians, q_to_a
#perov_import = csv_to_np('/Volumes/GoogleDrive/Shared drives/Wellesley Solar/Current Projects/ Hoke Effect/Inital_all_compositions/initial_full_2.csv')
perov_import = csv_to_np('/Users/rbelisle/Desktop/SSRL_Data_113/contacts_surface_wedge_750bins.csv')
#%% Covert to Q
wavelength = 0.982225 # 0.982381 #value from qsas with optimized energy
q = two_to_q(perov_import[:,0],wavelength)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 
#%% #Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.94
q_2 = 2.20
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file
#%% Do Curve Fitting For One Peak
sample = 'Bromine67'
name = 'x = 0.67'
image = 0 #index of file you want to look at
p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [1, 3000, q_2, 5]
lower_limit = [0, 0, q_1, 0]
popt,pcov = curve_fit(pvoigt, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,pvoigt(q_sub, *popt),'c--', label='Model') #plot best fit
init_q = popt[2]
#%% Do Curve Fitting for Three Peaks
image = 5 #index of file you want to look at
#init_q = 2.08389035 #q value from initial structure
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
popt,pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
# Plot everything
fig = plt.figure(figsize=(8, 5)) 
plt.plot(q_sub,perov_fit[:,image],'m-', label=name, linewidth = 2.0) #plot subfield of data
plt.plot(q_sub,three_gaussians(q_sub, *popt),'k-', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'k-.', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],popt[4],popt[5]), 'k--', label='Peak 2')
plt.plot(q_sub, gaussian(q_sub, popt[6],popt[7],popt[8]), 'k:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
print('Intensity:', popt[0], popt[3], popt[6])
print('Lattice Spacing:', q_to_a(popt[1],miller),q_to_a(popt[4],miller), q_to_a(popt[7],miller))


fig2 = plt.figure(figsize=(8, 5))
plt.plot(q_sub,perov_fit[:,image]-three_gaussians(q_sub, *popt),'m.', label=name+'Residual') #plot subfield of data
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Residuals [a.u.]',size=14)#Define y-axis label
# %%
fig.savefig('/Users/rbelisle/Desktop/SSRL_Data_113/'+sample,bbox_inches='tight', dpi=900, transparent=True)
fig2.savefig('/Users/rbelisle/Desktop/SSRL_Data_113/'+sample+'Residual',bbox_inches='tight', dpi=900, transparent=True)

# %% Temporary stacked plot
fig = plt.figure(figsize=(6, 7)) 
image = 4
plt.plot(q_sub,perov_fit[:,image]/max(perov_fit[:,image]), 'k--')
image = 5
plt.plot(q_sub,perov_fit[:,image]/max(perov_fit[:,image]), 'k-', label = 'Control')


offset = 1.05
image = 2
plt.plot(q_sub,perov_fit[:,image]/max(perov_fit[:,image])+offset, 'b--')
image = 3
plt.plot(q_sub,perov_fit[:,image]/max(perov_fit[:,image])+offset, 'b-', label = 'PTAA')

offset = 2.1
image = 0
plt.plot(q_sub,perov_fit[:,image]/max(perov_fit[:,image])+offset, 'r--')
image = 1
plt.plot(q_sub,perov_fit[:,image]/max(perov_fit[:,image])+offset, 'r-', label = 'C60')


plt.legend(loc="upper left")#Put legend in upper left hand corner
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Normalized Intensity [a.u.]',size=14)#Define y-axis label

# %%
fig.savefig('/Users/rbelisle/Desktop/SSRL_Data_113/'+'surfacecontacts.png',bbox_inches='tight', dpi=900, transparent=True)

# %%
