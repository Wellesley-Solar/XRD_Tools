#this program takes in XRD data across our compositional space and predicts
#the peak position of a new composition
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from xrdfunctions import csv_to_np, two_to_q, trim_data, num_files, back_subtract, pvoigt, gaussian, three_gaussians, q_to_a
wavelength = 0.982381 #value from qsas with optimized energy

#Import Data
perov_import = csv_to_np('/Volumes/GoogleDrive/Shared drives/Wellesley Solar/Current Projects/ Hoke Effect/Inital_all_compositions/initial_full.csv')
# Covert to Q
q = two_to_q(perov_import[:,0],wavelength)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.95
q_2 = 2.18
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)

#%% #Loop through Data
files = num_files(perov)
lattice = np.zeros(files)
perov_fit = perov_sub
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

p0 = [0.5, 50, 2.02, 0.01] #best guess for number of peaks an initial values
upper_limit = [1, 100, q_2, 5]
lower_limit = [0, 0, q_1, 0]
for frame in range(files): 
    popt, pcov = curve_fit(pvoigt, q_sub, perov_fit[:,frame], p0, bounds = [lower_limit, upper_limit], maxfev = 4000)
    lattice[frame] = q_to_a(popt[2],miller)
    plt.plot(q_sub, perov_fit[:,frame])
    plt.plot(q_sub, pvoigt(q_sub, *popt))
    error = np.sqrt(np.diag(pcov))
    print(error)
    print('Intensity:', popt[1], error[1])
    print('Lattice Spacing:', q_to_a(popt[2],miller), error[2])
# %% #Assign chemistry and fit
chem = np.array([0, .33, .5, .67, .75, 1]) #percent bromine for high confidence fits
fits = lattice #lattice spacings for cubic perovskites
fits[0] = 6.272 #psuedocubic fit from additional analysis 
def linear(x,m,b):
    return m*x + b
popt, pcov = curve_fit(linear, chem, fits)

plt.figure(num = 1, figsize=(8,6))
fig1,ax1 = plt.subplots()
ax1.set_xlabel('Bromine Fraction [x]',size=14) #Define x-axis label
ax1.set_ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.plot(chem,fits, 'ko')
plt.plot(chem,linear(chem,*popt),'k--')
plt.savefig('Initial Lattice.png')
# %% predict q for specific bromine compositions
def predict_lattice(m,b,br):
    return m*br+b

lattice_20 = predict_lattice(popt[0],popt[1],0.2)
print(popt)

def predict_q(br,index):
    lattice = predict_lattice(popt[0],popt[1],br)
    return 2*np.pi*np.sqrt(index[0]**2+index[1]**2+index[2]**2)/lattice

def predict_chem(q,index):
    lattice = q_to_a(q,index)
    return (lattice - popt[1])/popt[0]

q_20 = predict_q(0.2,[1,0,0])
br_20 = predict_q(q_20, [1,0,0])
print(q_20, br_20)

# %%
