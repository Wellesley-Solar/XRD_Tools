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
#Define wavelength and convert two-theta to Angstroms
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
#remove background
no_back=back_subtract(x, data, length)
    #x is a 1D array of two theta or q values (ie q_sub in old code)
    #data is an array of x-ray intensities (ie int_correct in old code)
    #length is the number of values on the edges of the data you want to use to create a linear background (limit1 to limit2)
print(no_back)
plt.plot(no_back)
plt.show()

#%%
p0 = [200, q, .01] #p0 = [height, center, width] guesses
intensity_1 = np.zeros((num_frames)) #create correct size arrays for running in the loop
lattice_1= np.zeros((num_frames)) 
for j in range(num_frames):
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0) #fit with gaussian function after giving guesses
    intensity_1[j] = popt[0] #peak height
    lattice_1[j] = 2*math.pi/popt[1] #center/lattice spacing #need to fit for correct index
    p0 = popt 

#%%
#convert from frames to time
time=frames_to_time(num_frames,speed,start_lag)
#find speed and start_lag in tif text file


#%%
plt.figure(figsize=(12,10))
plt.plot(time, lattice_1, marker='bo') 
plt.xlabel('time(s)')
plt.ylabel('Lattice Spacing(angstrom)')
plt.title('Lattice Spacing vs. Time')
plt.xlim()
plt.ylim()
plt.legend(frameon=True, fancybox=True,framealpha=1, shadow=False, borderpad=1, 
           title="figure_info", loc='lower right', fontsize='14')
plt.show()


#%%
plt.figure(figsize=(12,10))
plt.plot(time, intensity_1, marker='mo') 
plt.xlabel('time(s)')
plt.ylabel('Relative Intensity')
plt.title('Intensity vs. Time')
plt.xlim()
plt.ylim()
plt.legend(frameon=True, fancybox=True,framealpha=1, shadow=False, borderpad=1, 
           title="figure_info", loc='lower right', fontsize='14')
plt.show()

    
