#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
import pandas as pd
import math
from scipy.optimize import curve_fit #for some reason I've found I have to call this out explicity to get curve_fit to run.

#%%
from xrdfunctions import * #imports all functions in file named xrdfunctions
#xrdfunctions can be found on our github

#%%
#import csv file
csv_to_np(put_filename_here)
perov=np.array(data)

#%%
#Plotting initial frame of data
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0],perov[:,1],'r-', label='$MAPbIBr_2$ initial') #plot two-theta versus XRD intensity
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)
plt.title('initial')
plt.legend(loc="upper right")
plt.show()

#%%
#Define wavelength and convert two-theta to inverse angstroms
#Note that our initial 2-theta values are in degrees, not radians
q = two_to_q(perov[:,0],0.9763)
perov_intensity=perov[:,1] 
plt.figure(figsize=(8,6)) 
plt.plot(Q,perov_intensity, marker='.',color='r')
plt.title('Initial:')
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.show()

#%%
#choose a peak and find its limits 
trim_data(q, perov_intensity, limit1, limit2)
q_sub = x[set1:set2]
perov_sub = data[set1:set2,:]
plt.plot(q_sub,perov_sub[:,-1])
plt.show()

#%%
#size
size = perov_sub.shape
print (size)
q_bins = size[0]
num_frames = size[1]

#%%
#remove the background
no_back=back_subtract(x, data, length)
int_correct=np.array([no_back[1]])
    #x is a 1D array of two theta or q values 
    #data is an array of x-ray intensities 
    #length is the number of values on the edges of the data you want to use to create a linear background (limit1 to limit2)

print(no_back)
plt.plot(no_back)
plt.show()

#%%
#single gaussian fit
p0 = [a, b, c] #p0 = [height, center, width] guesses
intensity_1 = np.zeros((num_frames)) #create correct size arrays for running in the loop
lattice_1= np.zeros((num_frames)) 
for j in range(num_frames):
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0) #fit with gaussian function after giving guesses
    intensity_1[j] = popt[0] #peak height
    plane = (math.sqrt(h**2+k**2+l**2))
    lattice_1[j] = (2*math.pi/popt[1]) #center/lattice spacing #h, k, and l are miller indices
    p0 = popt 
    
#%%
#if two gaussian fit is needed
popt, pcov = curve_fit(two_gaussians, np.array(q_sub), int_correct, p0=[a1, b1, c1, a2, b2, c2]) #p0 is guesses
p0 = popt 
pars_1 = popt[0:3]
pars_2 = popt[3:6]
gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
plt.figure(figsize=(12,10)) 
plt.plot(np.array(q_sub), int_correct, color = 'red', label='perov')
plt.plot(np.array(q_sub), two_gaussians(np.array(q_sub),*popt), color = 'orange')
print(popt [0:3])
print(popt [3:6])
plt.plot(np.array(q_sub), gauss_peak_1,color='blue')
plt.plot(np.array(q_sub),gauss_peak_2,color='purple') 
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.legend(loc="upper right")
print(max(int_correct))

#%%
#if three gaussian fit is needed
popt, pcov = curve_fit(three_gaussians, np.array(q_sub), int_correct, p0=[a1, b1, c1, a2, b2, c2, a3, b3, c3]) #p0 is guesses
p0 = popt 
pars_1 = popt[0:3]
pars_2 = popt[3:6]
pars_3 = popt[6:9]
gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
gauss_peak_3 = gaussian(np.array(q_sub), *pars_3)
plt.figure(figsize=(12,10)) 
plt.plot(np.array(q_sub), int_correct, color = 'red', label='perov')
plt.plot(np.array(q_sub), three_gaussians(np.array(q_sub),*popt), color = 'orange')
print(popt [0:3])
print(popt [3:6])
print(popt [6:9])
plt.plot(np.array(q_sub), gauss_peak_1,color='blue')
plt.plot(np.array(q_sub),gauss_peak_2,color='purple') 
plt.plot(np.array(q_sub),gauss_peak_3,color='green') 
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.legend(loc="upper right")
print(max(int_correct))

#%%
#convert from frames to time
time=frames_to_time(num_frames,speed,start_lag)
#find speed and start_lag in tif text file


#%%
#plot lspace
plt.figure(figsize=(12,10))
plt.plot(time, lattice_1, marker='bo') #change from lattice_1 if multi gauss
plt.xlabel('time(s)')
plt.ylabel('Lattice Spacing(angstrom)')
plt.title('Lattice Spacing vs. Time')
plt.xlim()
plt.ylim()
plt.legend(frameon=True, fancybox=True,framealpha=1, shadow=False, borderpad=1, 
           title="figure_info", loc='lower right', fontsize='14')
plt.show()


#%%
#plot intensity
plt.figure(figsize=(12,10))
plt.plot(time, intensity_1, marker='mo') #change from intensity_1 if multi gauss
plt.xlabel('time(s)')
plt.ylabel('Relative Intensity')
plt.title('Intensity vs. Time')
plt.xlim()
plt.ylim()
plt.legend(frameon=True, fancybox=True,framealpha=1, shadow=False, borderpad=1, 
           title="figure_info", loc='lower right', fontsize='14')
plt.show()

#%%
#chemistry guess 
index = [h,k,l] #miller indices --> fill in
halide_frac=q_to_chem(popt[1],index)
plt.figure(figsize=(12,10))
plt.plot(time, halide_frac, marker='ko') 
plt.xlabel('time(s)')
plt.ylabel('Halide Fraction')
plt.title('Halide Fraction vs. Time')
plt.xlim()
plt.ylim()
plt.legend(frameon=True, fancybox=True,framealpha=1, shadow=False, borderpad=1, 
           title="figure_info", loc='lower right', fontsize='14')
plt.show()
