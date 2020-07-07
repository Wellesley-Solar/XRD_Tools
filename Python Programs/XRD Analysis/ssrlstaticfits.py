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
perov_import = csv_to_np('/Users/rbelisle/Desktop/integrationtest.csv')
#%% Covert to Q
q = two_to_q(perov_import[:,0],0.9763)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 


#%% #Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.95
q_2 = 2.18
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)

file_to_fit = -1 #choose the index of the file you want ot look at
data_to_fit = back_substract(q_sub,perov_sub[:,file_to_fit],10)

#%% Define three peak fitting
def three_gaussians(x, a1, b1, c1, a2, c2, a3, b3, c3):
    return (gaussian(x, a1, b1, c1) +
            gaussian(x, a2, 2.0456, c2)+ #this would be your initial peak center in Q
            gaussian(x, a3, b3, c3))

#%% Do Curve Fitting For three peaks
init_q = 2.0456 #q value from initial structure
p0 = [5, 2.02, .001, 5, .01, 6, 2.07, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
popt,pcov = curve_fit(three_gaussians, q_sub, data_to_fit, p0, maxfev=5000)

# %% Plot everything
plt.plot(q_sub,data_to_fit,'b-', label='Data') #plot subfield of data
plt.plot(q_sub,three_gaussians(q_sub, *popt),'c--', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'b-.', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],init_q,popt[4]), 'b--', label='Peak 2')
plt.plot(q_sub, gaussian(q_sub, popt[5],popt[6],popt[7]), 'b:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
plt.savefig('MAPbI2Br_60minDark_three.png')
print('Intensity:', popt[0], popt[3], popt[5])
print('Lattice Spacing:', q_to_a(popt[1],miller),q_to_a(popt[6],miller))

# %%
