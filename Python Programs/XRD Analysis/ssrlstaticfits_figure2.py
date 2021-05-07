#using SSRL static fits to make figure for direct observations of the halide segregation
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
import math as math
from xrdfunctions import csv_to_np, two_to_q, trim_data, num_files, back_subtract, pvoigt, gaussian, three_gaussians, q_to_a

#%% Import Data
perov_import = csv_to_np('/Users/rbelisle/Desktop/SSRL_Data_113/all_wedge_onoff.csv')

#   Covert to Q
wavelength = 0.982381 #value from qsas with optimized energy
q = two_to_q(perov_import[:,0],wavelength)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

#   Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.95
q_2 = 2.18
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

#%% DETERMENING PHASES
sample = 'Bromine67'
name = 'x = 0.67'
image_start = 9 #index of file before illumination
image_sun = 0 #index of file at one sun


#   Do Curve Fitting For One Peak
p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [1, 3000, q_2, 5]
lower_limit = [0, 0, q_1, 0]
popt,pcov = curve_fit(pvoigt, q_sub, perov_fit[:,image_start], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
plt.plot(q_sub,perov_fit[:,image_start],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,pvoigt(q_sub, *popt),'c--', label='Model') #plot best fit
init_q = popt[2]
init = popt
fwhm.append(init[-1])

#%%
#   Do Curve Fitting for Light
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
popt,pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,image_sun], p0, bounds=(lower_limit, upper_limit), maxfev=8000)

#   Plot everything
fig = plt.figure(figsize=(8, 5)) 
plt.plot(q_sub,perov_fit[:,image_sun],'b-', label=name, linewidth = 2.0) #plot subfield of data
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
plt.plot(q_sub,perov_fit[:,image]-three_gaussians(q_sub, *popt),'b.', label=name+'Residual') #plot subfield of data
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Residuals [a.u.]',size=14)#Define y-axis label
# %% Save Figures
#fig.savefig('/Users/rbelisle/Desktop/SSRL_Data_113/'+sample,bbox_inches='tight', dpi=900, transparent=True)
#fig2.savefig('/Users/rbelisle/Desktop/SSRL_Data_113/'+sample+'Residual',bbox_inches='tight', dpi=900, transparent=True)

# %% INTENSITY CALCULATIONS
#   from Laura to do phase calculations should integrate peak area
#   will use fits to do calculations
#   still having issues with chemistry calculations
initial_crystal = sum(pvoigt(q_sub,*init))
print("Starting Composition:", q_to_chem(init_q,miller), initial_crystal)

final_crystal = sum(three_gaussians(q_sub,*popt))
print("Percentage of Original Crystallinity:", final_crystal/initial_crystal)

iodide_crystal = sum(gaussian(q_sub,*popt[0:3]))
original_crystal = sum(gaussian(q_sub,*popt[3:6]))
bromide_crystal = sum(gaussian(q_sub,*popt[6:]))

print("Iodine Rich:", q_to_chem(popt[1],miller), iodide_crystal/final_crystal)
print("Original Phase:", q_to_chem(popt[4],miller), original_crystal/final_crystal)
print("Bromide Rich:", q_to_chem(popt[7],miller), bromide_crystal/final_crystal)
print("Final Chemistry:", 1/(iodide_crystal+original_crystal+bromide_crystal)*(q_to_chem(popt[1],miller)*iodide_crystal+q_to_chem(popt[4],miller)*original_crystal+q_to_chem(popt[7],miller)*bromide_crystal))

# %% DO CHEMISTRY ANALYSIS
perov_init = np.column_stack((perov_fit[:,1],perov_fit[:,3],perov_fit[:,6],perov_fit[:,10]))

#   for mixed compositions
files = num_files(perov_init)
lattice = np.zeros(files)
p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [1, 3000, q_2, 5]
lower_limit = [0, 0, q_1, 0]
for frame in range(files):
    popt,pcov = curve_fit(pvoigt, q_sub, perov_init[:,frame], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
    lattice[frame] = q_to_a(popt[2], miller)

#   for bromine composition
popt,pcov = curve_fit(pvoigt, qchem, Br, p0, bounds=(lower_limit, upper_limit), maxfev=6000)
lattice_br = q_to_a(popt[2], miller)

#   for iodine composition
p0 = [300, 2, .01,300, 2, .01]
upper_limit = [3000, q_2, 5, 3000, q_2, 5]
lower_limit = [0, q_1, 0, 0, q_1, 0]
popt,pcov = curve_fit(two_gaussians, qchem, I, p0, bounds = (lower_limit, upper_limit), maxfev=6000)
lattice_I = .5*q_to_a(popt[4], [0,0,4])/2 + .5*q_to_a(popt[1], [2,2,0])/np.sqrt(2)

chem = np.array([0, .33, .5, .67, .75, 1]) #percent bromine for high confidence fits
fits = np.concatenate(([lattice_I], lattice, [lattice_br])) 
def linear(x,m,b):
    return m*x + b
popt_chem, pcov = curve_fit(linear, chem, fits)

plt.figure(num = 1, figsize=(8,6))
fig1,ax1 = plt.subplots()
ax1.set_xlabel('Bromine Fraction [x]',size=14) #Define x-axis label
ax1.set_ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.plot(chem,fits, 'ko')
plt.plot(chem,linear(chem,*popt_chem),'k--')
#plt.savefig('Initial Lattice.png')

def q_to_chem(center,plane):
    #center is the center of an xrd peak
    #plane is a list of the formal [h,k,l]
    #takes peak position in q and converts it to bromine percentage
    #using linear fit for bromine fraction on lattice spacing
    slope = popt_chem[0]
    intercept = popt_chem[1]
    br_frac = 1/slope*(q_to_a(center,plane)-intercept)
    return br_frac

# %% Make compond plot of initial, hour of light, and two hours of dark
fig3, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, figsize=(5,10))
ax3.plot(q_sub,norm(perov_fit[:,1]), color = 'b')
ax3.plot(q_sub,norm(perov_fit[:,3]), color = 'r')
ax3.plot(q_sub,norm(perov_fit[:,6]), color = 'g')
ax3.plot(q_sub,norm(perov_fit[:,10]), color = 'm')
ax3.plot(qchem,norm(Br), color = 'y')
ax3.plot(qchem,norm(I), color = 'k')
ax3.set_xlabel('Q [$\AA^{-1}$]',size=12)
ax3.set_ylabel('Norm. Intensity [a.u.]', size=12)

ax2.plot(q_sub,norm(perov_fit[:,2]), color = 'b')
ax2.plot(q_sub,norm(perov_fit[:,4]), color = 'r')
ax2.plot(q_sub,norm(perov_fit[:,7]), color = 'g')
ax2.plot(q_sub,norm(perov_fit[:,11]), color = 'm')
ax2.set_ylabel('Norm. Intensity [a.u.]', size=12)

ax1.plot(q_sub,norm(perov_fit[:,0]), color = 'b')
ax1.plot(q_sub,norm(perov_fit[:,5]), color = 'r')
ax1.plot(q_sub,norm(perov_fit[:,8]), color = 'g')
ax1.plot(q_sub,norm(perov_fit[:,9]), color = 'm')
ax1.set_ylabel('Norm. Intensity [a.u.]', size=12)

# %% Save compound figure
fig3.savefig('/Users/rbelisle/Desktop/SSRL_Data_113/'+'Figure2',bbox_inches='tight', dpi=900, transparent=True)

    

