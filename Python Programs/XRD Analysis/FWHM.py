#%% 
import numpy as np
import csv
import os
import decimal
import matplotlib.pyplot as plt
from matplotlib import pyplot as mp
import pandas as pd
import math
import scipy
from scipy.optimize import curve_fit
from math import pi,sin
from heapq import nsmallest


def fwhm(y_values_q, x_values):
    y_values, q_l, q_r,q,inten = [], [], [],[],[]

    # To make 'y_values_temp', a numpy array, into a python list
    for x in range(0,len(y_values_q)):
        y_values.append(y_values_q[x])
    peak_height = max(y_values)
    print(peak_height)
    half_peak_height = max(y_values)/2
    print(half_peak_height)
    # Splitting the y_values data into before and after x_value at peak height
    y_l_q = y_values[0:y_values.index(peak_height)]
    #print(y_l_temp)
    y_r_q = y_values[y_values.index(peak_height):len(y_values)]
    #print(y_r_temp)
    # Finds 1st closest value to half_peak_height in y_l and y_r
    y_l = nsmallest(1, y_l_q, key=lambda x: abs(x-half_peak_height))
    print (y_l)
    y_r = nsmallest(1, y_r_q, key=lambda x: abs(x-half_peak_height))
    print (y_r)
    inten.append(y_l)
    inten.append(y_r)
    #print(inten)
    # Gets x_value pairs for y_l and y_r
    q_l.append(x_values[y_values.index(y_l)])#[0]
    q_r.append(x_values[y_values.index(y_r)])#[0]+ len(y_l) -1
    fwhm_n = ((abs(q_l[0] - q_r[0])))
    q.append(q_l)
    q.append(q_r)
    #print (temp_l[0])
    print(q_l)
    print(q_r)
    #print ( fwhm_n)
    plt.plot(q,inten)
    return(float(fwhm_n))
#%%
#Defining all functions that will be used in program
def readcsv(filename):
    data = pd.read_csv(filename)
    return(np.array(data))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_index(array, value):
    #array is a 1D vector of two_theta or q values
    #value is the specific angle or q for which you want the index
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
 
def back_subtract(x, data, length):
    #x is a 1D array of two theta or q values
    #data is an array of x-ray intensities
    #length is the number of values on the edges of the data you want to use to create a linear background 
    x_linear = np.hstack((x[0:length], x[-length:-1])) #I'm taking the starting and ending values
    data_linear = np.hstack((data[0:length], data[-length:-1])) #We'll use these to fit a straight line
    slope, intercept = np.polyfit(x_linear, data_linear, 1) #Do linear fit
    back = slope*x+intercept 
    data_correct=(data-back)
    return data_correct

           
#%%
perov = readcsv('/Users/rbelisle/Desktop/SSRL_Data_113/Xraydeg/B1_MAPbI2Br_xraydeg/b1_deg.csv') #openfile of interest
plt.plot(perov[:,0],perov[:,1]) #plot inititial data
plt.title('Initial:')
plt.xlabel('2-theta')
plt.ylabel('Intensity')
plt.show()

# %%
# Convert two theta to Q and set data limits
#should turn this into a function that accepts an array and limits
wave = 0.982381 #   wavelength from GSAS-II
q =[4*math.pi*math.sin(math.pi/180*row/2)/wave for row in perov[:,0]]
#plt.plot(q,y, marker='.',color='blue')
plt.title('Initial:')
plt.xlabel('Q')
plt.ylabel('Intensity')
q_1 = 1.94
q_2 = 2.18
limit1 = q.index(find_nearest(q, q_1))
limit2 = q.index(find_nearest(q, q_2))
q_sub = q[limit1:limit2]


#%%
#remove background
size = perov_sub.shape
q_bins = size[0]
num_frames = size[1]
perov_fit = perov_sub
for file in range(num_frames): 
    perov_fit[:,file] = back_subtract(np.array(q_sub),perov_sub[:,file],10) #remove background from that file
plt.plot(q_sub,perov_sub[:,0]) #    plot to ensure it worked

# %%

y_values = perov_sub[:,0]
x_values = np.array(q_sub)
#for x in range(0,len(y_values_temp)):
#        y_values.append(y_values_temp[x]) I don't think you need to do this, since you already have an arrat of y-values
peak_height = max(y_values)
print(peak_height)

    
half_peak_height = max(y_values)/2
print(half_peak_height)
# Splitting the y_values data into before and after x_value at peak height
y_l_temp = y_values[0:find_index(y_values, peak_height)]
print(y_l_temp)
y_r_temp = y_values[find_index(y_values, peak_height):] #adjusted so it works with nparray, see new find index function
print(y_r_temp)


# Finds 1st closest value to half_peak_height in y_l and y_r
y_l = find_nearest(y_l_temp,half_peak_height) #can use find nearest again here
print (y_l)

y_r = find_nearest(y_r_temp,half_peak_height) #can use find nearest again here
print (y_r)

#   I think find index will work here
temp_l = find_index(y_l_temp,y_l)
temp_r = find_index(y_r_temp, y_r)
fwhm = abs(q_sub[temp_r]-q_sub[temp_l])
print(temp_l)
print(temp_r)
print('FWHM', fwhm)
# %% Perhaps a cleaner way to do this is pull the FWHM from the pseudo voight fit
def normal_gaussian(x, a, b, c): 
    #nomralized gaussian curve for XRD analysis:
    #x is a 1D array of two theta or q values
    #a is the instensity 
    #b is the peak position and 
    #c is the variance (FWHM = sqrt(2ln2)*c)
    return a/(c*np.sqrt(2*math.pi))*np.exp(-(x - b)**2/(2*c**2))

def lorentz(x, a, b, c):
    #generic lorentzian curve, for xrd analysis
    #x is a 1D array of two theta or q values
    #a is the max intensity of the peak and representative of crystalling
    #b is the peak position and 
    # c is the FWHM
    return a/np.pi*((c/2)/((x-b)**2+(c/2)**2))

def pvoigt(x, e, a, b, c):
    #pseudovoigt curve common in xrd analysis
    #linear combination of lorentzian and gaussian curves
    #e is the fraction that is lorentzian
    # a is the intensity
    # b is the centroid
    # c is the FWHM
    c_g = c/(2*np.sqrt(2*np.log(2)))
    return e*lorentz(x, a, b, c) + (1-e)*normal_gaussian(x,a,b,c_g)

p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [1, 3000, q_2, 5]
lower_limit = [0, 0, q_1, 0]
popt,pcov = curve_fit(pvoigt, q_sub, perov_fit[:,0], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
plt.plot(q_sub,perov_fit[:,0],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,pvoigt(q_sub, *popt),'c--', label='Model') #plot best fit

# %%
