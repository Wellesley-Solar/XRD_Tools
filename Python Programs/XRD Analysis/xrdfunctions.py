#%%
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
import pandas as pd

# %%
#Defining all functions that will be used in program
def csv_to_np(filename):
    #filemane is a csv of xrd data with the first column being two-theta 
    # and subsequent columns being xray intensities
    data = pd.read_csv(filename)
    return(np.array(data))

def two_to_q(two_theta, wave):
    #two_theta is a 1D array of two_theta angles
    #wave is the X-ray energy in angstroms
    rad_theta = two_theta/2*np.pi/180
    q = 4*np.pi*np.sin(rad_theta)/wave
    return q

def find_nearest(array, value):
    #array is a 1D vector of two_theta or q values
    #value is the specific angle or q for which you want the index
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def trim_data(x, data, limit1, limit2):
    #x is a 1D array of two theta or q values
    #data is an array of x-ray intensities
    #limit1 and limit2 are what you'd like to trime your data to 
    set1 = find_nearest(x,limit1)
    set2 = find_nearest(x,limit2)
    return x[set1:set2], data[set1:set2,:]

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

def gaussian(x, a, b, c): 
    #generic gaussian curve, for XRD analysis:
    #x is a 1D array of two theta or q values
    #a is the max intensity of the peak and representative of crystalling
    #b is the peak position and 
    # c is the FWHM
    return a*np.exp(-(x - b)**2/(2*c**2))

def two_gaussians(x, a1, b1, c1, a2, b2, c2):
        return (gaussian(x, a1, b1, c1) +
            gaussian(x, a2, b2, c2))

def multi_gaussian(x, guesses):
    #NOTE This function does not work with curve fitting yet TBD
    # x is 1D array of 2-theta or q values for our fitting
    # trips is an array of fits i.e. [[200, 1, .01], [150, 1.05. .02]]
    peaks = [gaussian(x, fit[0], fit[1], fit[2]) for fit in guesses]
    return np.sum(peaks, axis=0)

#%% Define three peak fitting
def three_gaussians(x, a1, b1, c1, a2, b2, c2, a3, b3, c3):
    return (gaussian(x, a1, b1, c1) +
            gaussian(x, a2, b2, c2)+ #this would be your initial peak center in Q
            gaussian(x, a3, b3, c3))

def lorentz(x, a, b, c):
    #generic lorentzian curve, for xrd analysis
    #x is a 1D array of two theta or q values
    #a is the max intensity of the peak and representative of crystalling
    #b is the peak position and 
    # c is the FWHM
    return a/np.pi*(c/((x-b)**2+c**2))

def pvoigt(x, e, a, b, c):
    #pseudovoigt curve common in xrd analysis
    #linear combination of lorentzian and gaussian curves
    #e is the fraction that is lorentzian
    return e*lorentz(x, a, b, c) + (1-e)*gaussian(x,a,b,c)

def mult_pvoigt(x, e, a, b, c, a2, b2, c2, a3, b3, c3):
    return pvoigt(x, e, a, b, c) + pvoigt(x, e, a2, b2, c2) +  pvoigt(x, e, a3, b3, c3) 

def q_to_a(center,plane):
    #center is the center of an xrd peak
    #plane is a list of the formal [h,k,l]
    a = 2*math.pi*math.sqrt(plane[0]**2+plane[1]**2+plane[2]**2)/center
    return a

def num_files(data):
    #data is the file you are currently analyzing
    #returns columns of data (i.e. frames for XRD data from SSRL)
    size = data.shape
    return size[1]

def q_to_chem(center,plane):
    #center is the center of an xrd peak
    #plane is a list of the formal [h,k,l]
    #takes peak position in q and converts it to bromine percentage
    #using linear fit for bromine fraction on lattice spacing
    slope = -0.34683009554825994
    intercept = 6.255901161926431
    br_frac = -1/slope*(q_to_a(center,plane)-intercept)
    return br_frac

def frames_to_time(x,speed,start_lag):
    #x=num_frames
    #speed=shutter speed in s
    #start_lag=t at x=0
    seconds=np.array([(x*speed)+ start_lag])
    return seconds
