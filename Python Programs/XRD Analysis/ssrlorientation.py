#This program analyzes data taken in a grazing incicdence geometry at 
#beam line 11-3 at SSRl to plot intensity of the [100] peak as a fucntion of chi
#code developed by Becky Belisle 12/2020
#%% Import the required packages
import numpy as np
import pylab 
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
import scipy.optimize 
from scipy.optimize import curve_fit
import pyFAI
import pygix
import fabio
import pandas as pd
from xrdfunctions import gaussian, lorentz, pvoigt

#%% Import the image you would like to analyze
directory = '/Users/rbelisle/Desktop/OrientationAnalysis/'
file_interest = 'C1_MAPbBr50_L60Off_30s_01230231_0060.tif'
dataFile = directory+file_interest
data = fabio.open(dataFile).data # Use fabio to open file as a np.array

# plot the imported file (this will not be calibrated)
fig1, ax1 = plt.subplots()
ax1.set_xlabel("x-pixel (#)")
ax1.set_ylabel("y-pixel (#)") 
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.imshow(data, vmin = 3, vmax = 100, origin='lower') #use imshow to rasterize np.array

# %% Load calibration file and set sample geomtetry
pg = pygix.Transform() # set up the transformation file
pg.load('/Users/rbelisle/Desktop/lab6_calib_315_0122.poni') # load the callibration file
pg.sample_orientation = 3    # set sample orientation: 1 is horizontal, 2 is vertical, 3 is horizontal rotated by 180 deg, 4 is vertical rotated by 180 deg (for more details see: https://github.com/tgdane/pygix/blob/master/pygix/transform.py)
pg.incident_angle = 3     # set indicent x-ray angle in deg 
pg.tilt_angle = 0             # tilt angle of sample in deg (misalignment in "chi")
pg # optionally print geometry

#%% Run calibration accounting for solid angle correction
ii_2d, qxy_2d, qz_2d = pg.transform_reciprocal(data, correctSolidAngle=True, method='bbox', npt=3000)

#plot the corrected file in the correct units
fig2, ax2 = plt.subplots()
ax2.set_xlabel("q$_\mathregular{xy}$ ($\mathring{A}^\mathregular{-1}$)")
ax2.set_ylabel("q$_\mathregular{z}$ ($\mathring{A}^\mathregular{-1}$)") 
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.imshow(ii_2d, extent=(np.min(qxy_2d)/10,np.max(qxy_2d)/10,np.min(qz_2d)/10,np.max(qz_2d)/10), vmin = 3, vmax = 100, origin='lower') 
plt.ylim(0,3.5)
plt.xlim(-2,2)

#%% Save corrected figure
temp = file_interest.split(sep='.')
filename = temp[0]+'_corrected.png'
fig2.savefig(directory+filename, bbox_inches='tight', dpi=900, transparent=True)
#%% Do 1D plot to identify peak position and width
Itest, q = pg.profile_sector(data, correctSolidAngle=True, radial_range = (9.5, 11), method='bbox', npt=3000) #look for peak between 1 and 1.1 Ang^-1
p0 = [.1, 140, 10.4, .2] #initial guess at fit paramenters
popt,pcov = curve_fit(pvoigt, q, Itest-Itest[0], p0, maxfev=6000) #fit peak using psuedovoigt
plt.plot(q,pvoigt(q, *popt),'c--', label='Model') #plot best fit
plt.plot(q,Itest-Itest[0])
center = popt[2]
width = popt[3]
#%% Using parameters from 1D plot make plot of intensity versus scattering angle
back, chi = I, chi = pg.profile_chi(data, correctSolidAngle=True, radial_pos = center+1, radial_width = 2*width, chi_range = (0,90), method='bbox', npt=3000) #looks at range from 0 to 90 degrees
I, chi = pg.profile_chi(data, correctSolidAngle=True, radial_pos = center, radial_width = 2*width, chi_range = (0,90), method='bbox', npt=3000)
Iback = I-back #correct for background intensity
Icorrect = Iback*np.sin(chi*math.pi/180) #correct for amount of scattered material
fig3, ax3 = plt.subplots()
ax3.set_xlabel('$\chi$ ($^\circ$)')
ax3.set_ylabel('Intensity (a.u.)')
plt.plot(chi,Icorrect)

#%% Save corrected figure
temp = file_interest.split(sep='.')
filename = temp[0]+'_orientation.png'
fig3.savefig(directory+filename, bbox_inches='tight', dpi=900, transparent=True)
