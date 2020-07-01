#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
import pandas as pd
import math



#Read csv file:
def readcsv(filename):
    data = pd.read_csv(filename)
    return(np.array(data))


#Limiting to a Peak
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin() #finds the index of the value closest to the target
    return array[idx]

#Single gaussian fit to plot
def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 

#Turns two theta in Q
def two_to_q(two_theta, wave):
    #two_theta is a 1D array of two_theta angles
    #wave is the X-ray energy in angstroms
    rad_theta = two_theta/2*np.pi/180
    q = 4*np.pi*np.sin(rad_theta)/wave
    return q

#Plotting initial frame of data
perov = readcsv('put_filename_here')
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0],perov[:,1],'r-', label='$MAPbIBr_2$ initial') #plot two-theta versus XRD intensity
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)


plt.title('initial')
plt.legend(loc="upper right")

#Define wavelength and convert two-theta to Angstroms
#Note that our initial 2-theta values are in degrees, not radians
#Two theta to q

q = two_to_q(perov[:,0],0.9763)
Q= q.tolist()
y=perov[:,1]
plt.figure(figsize=(8,6)) 
plt.plot(Q,y, marker='.',color='r')
plt.title('Initial:')
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.show()

#choose a peak and find its limits 
q_1 = .98
q_2 = 1.15
limit1 = q.index(find_nearest(q, q_1))
limit2 = q.index(find_nearest(q, q_2))
q_sub = q[limit1:limit2]
perov_sub = perov[limit1:limit2,1:-1]# range will depend on the number of frames you have 
plt.plot(q_sub,perov_sub[:,-1])

#remove background
size = perov_sub.shape
print (size)
q_bins = size[0]
#print (q_bins)
num_frames = size[1]
slope = np.zeros((num_frames,1))
intercept = np.zeros((num_frames,1))
back = np.zeros((q_bins,num_frames))
int_correct = np.zeros((q_bins,num_frames))

#accept a data file and range and returns average values at start and end of range
for i in range(num_frames): 
    slope[i] = ((np.mean(perov_sub[-10:-1,i])-np.mean(perov_sub[0:10,i]))/(np.mean(q[limit2-10:limit2])-np.mean(q[limit1:limit1+10])))
    intercept[i]=perov_sub[0,i]-slope[i]*q_sub[0]
    back[:,i] = [slope[i]*element+intercept[i] for element in q_sub]
    int_correct[:,i] = perov_sub[:,i]-np.array(back[:,i])
plt.plot(np.array(q_sub),int_correct)

p0 = [height, center, width] 
intensity_1 = np.zeros((num_frames)) #create correct size arrays for running in the loop
lattice_1= np.zeros((num_frames))
for j in range(num_frames):
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0)
    intensity_1[j] = popt[0]
    lattice_1[j] = 2*math.pi/popt[1] #need to fit for correct index
    p0 = popt #once you have the initial fit use that as your next guess, we expect things to be continuous so this helps 

#plot
time = np.zeros(num_frames) #create empty array for time of correct length
for x in range(num_frames):
    time[x] = x*20 + 20 #xray dose is 20s per frame

plt.plot(time, lattice_1)
plt.xlabel('time(s)')
plt.ylabel('Lattice Spacing')


# In[ ]:




