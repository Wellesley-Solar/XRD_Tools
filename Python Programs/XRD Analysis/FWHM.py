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


def fwhm(y_values_temp, x_values):
    y_values, temp_l, temp_r = [], [], []

    # To make 'y_values_temp', a numpy array, into a python list
    for x in range(0,len(y_values_temp)):
        y_values.append(y_values_temp[x])
    peak_height = max(y_values)
    print(peak_height)
    half_peak_height = max(y_values)/2
    print(half_peak_height)
    # Splitting the y_values data into before and after x_value at peak height
    y_l_temp = y_values[0:y_values.index(peak_height)]
    print(y_l_temp)
    y_r_temp = y_values[y_values.index(peak_height):len(y_values)]
    print(y_r_temp)
    # Finds 1st closest value to half_peak_height in y_l and y_r
    y_l = nsmallest(1, y_l_temp, key=lambda x: abs(x-half_peak_height))
    print (y_l)
    y_r = nsmallest(1, y_r_temp, key=lambda x: abs(x-half_peak_height))
    print (y_r)
    # Gets x_value pairs for y_l and y_r
    temp_l.append(x_values[y_l.index(y_l[0])])
    temp_r.append(x_values[y_r.index(y_r[0]) + len(y_l) -1])
    fwhm_n = ((temp_l[0] - temp_r[0]))
    print(temp_l)
    print(temp_r)
    print ( fwhm_n)
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
 


           
#%%
perov = readcsv('B1_MaPbI2Br_Light60min.csv') #openfile of interest
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

