#%%
#open file
#determine compositions/sample
#convert to Q
#crop data
#remove background
#fit guassians
#determine if should do 1,2,3 gaussians
#save lattice/intensity/strain and time/composition
#ideally have a function that does the fitting so we can call that function independent if doing
#in situ or start end

import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy
from scipy.optimize import curve_fit
from math import pi,sin

#%%
#Defining all functions that will be used in program
def readcsv(filename):
    data = pd.read_csv(filename)
    return(np.array(data))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 

def two_gaussians(x, h1, c1, w1, h2, c2, w2):
        return (gaussian(x, h1, c1, w1) +
            gaussian(x, h2, c2, w2))

def three_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3):
        return (gaussian(x, h1, c1, w1) +
            gaussian(x, h2, c2, w2)+gaussian(x, h3, c3, w3))
#%%
perov = readcsv('/Users/rbelisle/Desktop/startendxrd/initial.csv') #openfile of interest
plt.plot(perov[:,0],perov[:,1]) #plot inititial data
plt.title('Initial:')
plt.xlabel('2-theta')
plt.ylabel('Intensity')
plt.show()

# %%
# Convert two theta to Q and set data limits
#should turn this into a function that accepts an array and limits
q =[4*math.pi*math.sin(math.pi/180*row/2)/0.9763 for row in perov[:,0]]
#plt.plot(q,y, marker='.',color='blue')
plt.title('Initial:')
plt.xlabel('Q')
plt.ylabel('Intensity')
q_1 = 1.94
q_2 = 2.18
limit1 = q.index(find_nearest(q, q_1))
limit2 = q.index(find_nearest(q, q_2))
q_sub = q[limit1:limit2]

#q_sub = q[limit1:limit2]
#int_sub = y[limit1:limit2]
perov_sub = perov[limit1:limit2,1:]
plt.plot(q_sub,perov_sub[:,-1])

#%%
#remove background
size = perov_sub.shape
q_bins = size[0]
num_frames = size[1]
slope = np.zeros((num_frames,1))
b = np.zeros((num_frames,1))
back = np.zeros((q_bins,num_frames))
int_correct = np.zeros((q_bins,num_frames))

#accept a data file and range and returns average values at start and end of range
for j in range(num_frames): 
    slope[j] = ((np.mean(perov_sub[-10:-1,j])-np.mean(perov_sub[0:10,j]))/(np.mean(q[limit2-10:limit2])-np.mean(q[limit1:limit1+10])))
    b[j]=perov_sub[0,j]-slope[j]*q_sub[0]
    back[:,j] = [slope[j]*element+b[j] for element in q_sub]
    int_correct[:,j] = perov_sub[:,j]-np.array(back[:,j])

plt.plot(np.array(q_sub),int_correct)

# %%
#iterative guassian fitting
p0 = [50, 2.04, 0.01]#, 20, 2.04, 0.01] 
intensity_1 = np.zeros((num_frames)) #create correct size arrays for running in the loop
intensity_2= np.zeros((num_frames))
lattice_1= np.zeros((num_frames))
lattice_2 = np.zeros((num_frames))
print(lattice_1)
for j in range(num_frames):
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0)
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    #intensity_2[j] = popt[3]
    #lattice_2[j] = 4*math.pi/popt[4]
    p0 = popt #once you have the initial fit use that as your next guess, we expect things to be continuous so this helps with that
    print(lattice_1)
#%%
intensity_1 = np.zeros((num_frames))
intensity_2= np.zeros((num_frames))
intensity_3= np.zeros((num_frames))
lattice_1= np.zeros((num_frames))
lattice_2 = np.zeros((num_frames))
lattice_3 = np.zeros((num_frames))
# %%
#guassian fits with plotting
#pick fitting and initial guess
# it would be good to iteratively do these guesses, or use info from qlimit
#real issue in fitting - maybe this is four peaks? Perhaps I can look at inflection point instead
i = 3
j = 16
p0 =  [20, 2.08, 0.02, 10,2.08, 0.03,5,2.01, 0.05 ]
plt.plot(np.array(q_sub), int_correct[:,j],color='black')
if i == 2:
    popt, pcov = curve_fit(two_gaussians, np.array(q_sub), int_correct[:,j], p0[0:6])
    #plt.plot(np.array(q_sub), two_gaussians(np.array(q_sub),*popt))
    pars_1 = popt[0:3]
    pars_2 = popt[3:6]
    gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
    gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
    plt.plot(np.array(q_sub), gauss_peak_1, color='red')
    plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
    plt.plot(np.array(q_sub), gauss_peak_1+gauss_peak_2,color='green')
    print('lattice spacing:', [4*math.pi/popt[1], 4*math.pi/popt[4]])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    intensity_2[j] = popt[3]
    lattice_2[j] = 4*math.pi/popt[4]

if i == 3:
    popt, pcov = curve_fit(three_gaussians, np.array(q_sub), int_correct[:,j], p0)
    #plt.plot(np.array(q_sub), three_gaussians(np.array(q_sub),*popt))
    pars_1 = popt[0:3]
    pars_2 = popt[3:6]
    pars_3 = popt[6:9]
    gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
    gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
    gauss_peak_3 = gaussian(np.array(q_sub), *pars_3)
    plt.plot(np.array(q_sub), gauss_peak_1, color='magenta')
    plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
    plt.plot(np.array(q_sub),gauss_peak_3,color='green')
    plt.plot(np.array(q_sub),gauss_peak_3+gauss_peak_1+gauss_peak_2,color='red',linestyle='dashed')
    print('lattice spacing:', [4*math.pi/popt[1], 4*math.pi/popt[4],4*math.pi/popt[7]])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    intensity_2[j] = popt[3]
    lattice_2[j] = 4*math.pi/popt[4]
    intensity_2[j] = popt[6]
    lattice_3[j] = 4*math.pi/popt[7]

else: 
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0[0:3])
    plt.plot(np.array(q_sub), gaussian(np.array(q_sub),*popt))
    print('lattice spacing:', [4*math.pi/popt[1]])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 

print(popt)
#plt.savefig('time_42fit.png')
#popt_uncertainties = np.sqrt(np.diag(pcov))
#uncertainty = sum(popt_uncertainties)
#print('uncertainties:', popt_uncertainties)
#print('uncertainty:', uncertainty)

plt.xlabel('Q')
plt.ylabel('Intensity')

# %%
time = np.zeros(num_frames) #create empty array for time of correct length
for x in range(num_frames):
    time[x] = x*20 + 20 #xray dose is 20s per frame


# %%
#this all to better visualize what the peak fitting is doing to the data shape
gauss_peak_1 = gaussian(np.array(q_sub),11.3, 2.105, .0199)
gauss_peak_2 = gaussian(np.array(q_sub), 10.3, 2.1322, .00914)
gauss_peak_3 = gaussian(np.array(q_sub), 16.25, 2.051, .0232)
gauss_peak_4 = gaussian(np.array(q_sub), 1.6, 2.01,.005)
plt.plot(np.array(q_sub), int_correct[:,j],color='black')
plt.plot(np.array(q_sub), gauss_peak_1, color='magenta')
plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
plt.plot(np.array(q_sub),gauss_peak_3,color='green')
plt.plot(np.array(q_sub),gauss_peak_4,color='cyan')
plt.plot(np.array(q_sub),gauss_peak_4+gauss_peak_3+gauss_peak_1+gauss_peak_2,color='red',linestyle='dashed')

# %%
x = 0
for j in range(6): 
    plt.plot(np.array(q_sub),int_correct[:,x])
    x = x+5
plt.savefig('67Br33I.png')

# %%
plt.imshow(perov_sub)
fig, ax = plt.subplots(figsize=(6,6))
ax.imshow(perov_sub, interpolation='none', extent=[0,60,q_sub[0],q_sub[-1]])
ax.set_aspect(2)


# %%
