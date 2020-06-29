#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
import pandas as pd

# %%
#Defining all functions that will be used in program
def csv_to_np(filename):
    data = pd.read_csv(filename)
    return(np.array(data))

def two_to_q(two_theta, wave):
    rad_theta = two_theta/2*np.pi/180
    q = 4*np.pi*np.sin(rad_theta)/wave
    return q

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def trim_data(x, data, limit1, limit2):
    set1 = find_nearest(x,limit1)
    set2 = find_nearest(x,limit2)
    return x[set1:set2], data[set1:set2,:]

def back_substract(x, data, length):
    x_linear = np.hstack(x[0:length], x[-length:-1]) #I'm taking the starting and ending values
    data_linear = np.hstack(data[0:length], data[-length:-1]) #We'll use these to fit a straight line
    slope, intercept = np.polyfit(x_linear, data_linear, 1) #Do linear fit
    back = slope*x+intercept 
    return data-back

def gaussian(x, a, b, c): 
    return a*np.exp(-(x - b)**2/(2*c**2))

def two_gaussians(x, a1, b1, c1, a2, b2, c2):
        return (gaussian(x, a1, b1, c1) +
            gaussian(x, a2, b2, c2))

def multi_gaussian(x, trips):
    # x is 1D array of 2-theta or q values for our fitting
    # trips is an array of fits i.e. [[200, 1, .01], [150, 1.05. .02]]
    peaks = [gaussian(x, fit[0], fit[1], fit[2]) for fit in trips]
    return np.sum(peaks, axis=0)


