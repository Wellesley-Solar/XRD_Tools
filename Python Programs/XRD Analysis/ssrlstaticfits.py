#this is a base file for doing static fits (i.e. perovskite structure at the end of an experiment)
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd

#%% Import Data
from xrdfunctions import *
perov_import = csv_to_np('/Users/rbelisle/Desktop/startendxrd/xrdstartend2.csv')

#%% Covert to Q
q = two_to_q(perov_import[:,0],0.9763)
perov = perov_import[:,1:-1] #seperate XRD intensities from rest of data 


#%% #Trim and Remove Background
q_1 = 1.98
q_2 = 2.2
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)

file_to_fit = 2 #choose the index of the file you want ot look at
perov_to_fit = back_substract(q_sub,perov_sub[:,file_to_fit],10)

#%% #Do Curve Fitting
p0 = [20,2.05, 0.01, 20, 2.05,.01] #best guess for number of peaks an initial values
popt, pcov = curve_fit(two_gaussians, q_sub, perov_to_fit, p0)
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q_sub,perov_to_fit,'r-', label='60 Min Sun') #plot subfield of data
plt.plot(q_sub,two_gaussians(q_sub, *popt),'b--', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'g--', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],popt[4],popt[5]), 'c--', label='Peak 2')
plt.xlabel('Q [$\AA^{-1}$]',size=12) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=12)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner

