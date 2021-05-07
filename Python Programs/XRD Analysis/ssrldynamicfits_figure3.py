#this is a base file for doing dynamics c fits (i.e. perovskite structure at the end of an experiment)
#assumes three peaks
#if residual too large does not plot
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from xrdfunctions import *

# %% Create dictionary for storing analyzed results - to be run at the start of a new session
samplelist = []
results = {"Sample":"Results"}

#%% Import Data
samplename = 'D2_lighton'
perov_import = csv_to_np('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/' + samplename + '.csv')
light_on = 1 # if light on set to 1 else set to 0
interval = 2 # time interval between frames in minutes

# Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 


# Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.92
q_2 = 2.18
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

# Initialize parameters of interest
files = num_files(perov)
time = np.arange(0, files*interval, interval) + interval/2 #time interval of the experiment
lattice1 = np.zeros(files) #lattice spacing for iodine rich phase
lattice2 = np.zeros(files) #lattice spacing for original composition
lattice3 = np.zeros(files) #lattice spacing for bromide composition
phase1 = np.zeros(files) #phase gives area under the curve
phase2 = np.zeros(files)
phase3 = np.zeros(files)

# %%  Do Curve Fitting For One Peak
image = 0 #index of file you want to look at
p0 = [0.2, 300, 2, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [1, 3000, q_2, 5]
lower_limit = [0, 0, q_1, 0]
popt,pcov = curve_fit(pvoigt, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=6000)
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,pvoigt(q_sub, *popt),'c--', label='Model') #plot best fit
init_q = popt[2] #this is the inital peak position that will be used as a reference
FWHM = popt[3] 
#   save initial data
lattice1[image] = q_to_a(popt[2],miller) #lattice spacing for iodine rich phase
lattice2[image] = q_to_a(popt[2],miller) #lattice spacing for original composition
lattice3[image] = q_to_a(popt[2],miller) #lattice spacing for bromide composition
phase1[image] =  0 #phase gives area under the curve
phase2[image] = sum(pvoigt(q_sub,*popt))
phase3[image] = 0 

# %% Dynamics of Process
#  Plot Max Intensity over time
maxintensity = np.zeros(files)
for frame in range(files):
    maxintensity[frame] = max(perov_fit[:,frame])
    
plt.plot(time, maxintensity/maxintensity[0], 'k--')
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Max Intensity [Normalized]',size=14)#Define y-axis label

#  Plot Integrated Intensity over time
totalintensity = np.zeros(files)
for frame in range(files):
    totalintensity[frame] = sum(perov_fit[:,frame])
    
plt.plot(time, totalintensity/totalintensity[0], 'r--')
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Integrated Intensity [Normalized]',size=14)#Define y-axis label

# Plot Binned Intenisty
iod_bin = np.zeros(files)
orig_bin = np.zeros(files)
brom_bin = np.zeros(files)

iod_lim = find_nearest(q_sub, init_q-FWHM)
brom_lim = find_nearest(q_sub, init_q+FWHM)

for frame in range(files):
    iod_bin[frame] = sum(perov_fit[0:iod_lim,frame])
    orig_bin[frame] = sum(perov_fit[iod_lim:brom_lim,frame])
    brom_bin[frame] = sum(perov_fit[brom_lim:,frame])

plt.plot(time, iod_bin, 'r--')
plt.plot(time, orig_bin, 'k--')
plt.plot(time, brom_bin, 'b--')
plt.plot(time, brom_bin+iod_bin+orig_bin, 'k-')

plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Binned Intensity [a.u.]',size=14)#Define y-axis label

#  Plot Binned Q
iod_q = np.zeros(files)
orig_q = np.zeros(files)
brom_q = np.zeros(files)
av_q = np.zeros(files)

iod_lim = find_nearest(q_sub, init_q-FWHM)
brom_lim = find_nearest(q_sub, init_q+FWHM)

for frame in range(files):
    iod_q[frame] = sum(perov_fit[0:iod_lim,frame]*q_sub[0:iod_lim])/sum(perov_fit[0:iod_lim,frame])
    orig_q[frame] = sum(perov_fit[iod_lim:brom_lim,frame]*q_sub[iod_lim:brom_lim])/sum(perov_fit[iod_lim:brom_lim,frame])
    brom_q[frame] = sum(perov_fit[brom_lim:,frame]*q_sub[brom_lim:])/sum(perov_fit[brom_lim:,frame])
    av_q[frame] = sum(perov_fit[:,frame]*q_sub)/sum(perov_fit[:,frame])



plt.plot(time, iod_q, 'r--')
plt.plot(time, orig_q, 'k--')
plt.plot(time, brom_q, 'b--')
plt.plot(time, av_q, 'k-')

plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Binned Average Q [$\AA^{-1}$]',size=14)#Define y-axis label

# %% Save Variables
samplelist.append(samplename) 
datatostore = [time, maxintensity, iod_bin, orig_bin, brom_bin, iod_q, orig_q, brom_q, av_q]
results[samplename] = datatostore
print(samplelist)

#%% PEAK FITTING OPTION 1: Start from intiial frame and work forward in time
# Do Curve Fitting for Three Peaks 

image = 1 #index of file you want to look at
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
popt,pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=8000)

#   Plot everything
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,three_gaussians(q_sub, *popt),'c--', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'b-.', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],popt[4],popt[5]), 'b--', label='Peak 2')
plt.plot(q_sub, gaussian(q_sub, popt[6],popt[7],popt[8]), 'b:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
print('Intensity:', popt[0], popt[3], popt[6])
print('Lattice Spacing:', q_to_a(popt[1],miller),q_to_a(popt[4],miller), q_to_a(popt[7],miller))

#   Loop through Data from initial time until final frame
p0 = [5, init_q-.01, .01, 5, init_q, .001, 10, init_q+.01, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]

if light_on:    
    for frame in range(files-1): 
        popt, pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,frame+1], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
        residual = sum(np.abs(perov_fit[:,frame+1]-three_gaussians(q_sub,*popt)))
        
        if residual <= 30:
            lattice1[frame+1] = q_to_a(popt[1],miller)
            lattice2[frame+1] = q_to_a(popt[4],miller)
            lattice3[frame+1] = q_to_a(popt[7],miller)
            phase1[frame+1] = sum(gaussian(q_sub,*popt[0:3]))
            phase2[frame+1] = sum(gaussian(q_sub,*popt[3:6]))
            phase3[frame+1] = sum(gaussian(q_sub,*popt[6:]))
        
        else:
            print(frame, "bad fit", residual)
        
        p0=popt
        upper_limit = [30, popt[1], .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
        lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, popt[7], 0]

else: 
    for frame in range(files): 
        popt, pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,frame+1], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
        lattice1[frame] = q_to_a(popt[1],miller)
        lattice2[frame] = q_to_a(popt[4],miller)
        lattice3[frame] = q_to_a(popt[7],miller)
        phase1[frame] = sum(gaussian(q_sub,*popt[0:3]))
        phase2[frame] = sum(gaussian(q_sub,*popt[3:6]))
        phase3[frame] = sum(gaussian(q_sub,*popt[6:]))
        chem1[frame] = q_to_chem(popt[1],miller)
        chem2[frame] = q_to_chem(popt[4],miller)
        chem3[frame] = q_to_chem(popt[7],miller)
        time[frame] = frame*20+20
        print('Intensity:', popt[0])
        print('Lattice Spacing:', q_to_a(popt[1],miller))
        p0=popt
        upper_limit = [30, init_q+.001, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
        lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q-0.001, 0]

#   Plot lattice spacings over time
plt.figure(figsize=(5,4)) #make plot larger
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.plot(time,lattice1,'r.')
plt.plot(time,lattice2, 'k.')
plt.plot(time, lattice3, 'b.')
plt.ylim((q_to_a(q_2,miller),q_to_a(q_1,miller)))

#%% OPTION 2: Start from final and work backward in time

#   Do Curve Fitting for final three peaks TEST CELL (run this before running full loop to make sure nothing weird is happening)
image = -1 #index of file you want to look at
p0 = [5, init_q-.05, .01, 5, init_q, .001, 10, init_q+.05, .01] #best guess for initial values in format [a1, b1, c1, a2, c2, a3, b3, c3]
upper_limit = [30, init_q, .3, 30, init_q+.001, .3, 30, q_2, .1] #upper limits fot the range where one peak is iodine rich and the other is bromine
lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]
popt,pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,image], p0, bounds=(lower_limit, upper_limit), maxfev=8000)

#   save final data
lattice1[image] = q_to_a(popt[1],miller) #lattice spacing for iodine rich phase
lattice2[image] = q_to_a(popt[4],miller) #lattice spacing for original composition
lattice3[image] = q_to_a(popt[7],miller) #lattice spacing for bromide composition
phase1[image] =  sum(gaussian(q_sub,*popt[0:3]))
phase2[image] = sum(gaussian(q_sub,*popt[3:6]))
phase3[image] = sum(gaussian(q_sub,*popt[6:]))

# Plot everything
plt.plot(q_sub,perov_fit[:,image],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,three_gaussians(q_sub, *popt),'c--', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'b-.', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],popt[4],popt[5]), 'b--', label='Peak 2')
plt.plot(q_sub, gaussian(q_sub, popt[6],popt[7],popt[8]), 'b:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
print('Intensity:', popt[0], popt[3], popt[6])
print('Lattice Spacing:', q_to_a(popt[1],miller),q_to_a(popt[4],miller), q_to_a(popt[7],miller))
print('Residual:', sum(np.abs(perov_fit[:,image]-three_gaussians(q_sub,*popt))) )
p0 = popt
#    Do reversed in time fits of data if the light it son
if light_on:   
    for frame in reversed(range(files)): 
        popt, pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,frame], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
        residual = sum(np.abs(perov_fit[:,frame]-three_gaussians(q_sub,*popt)))
        
        if residual <= 30:
            lattice1[frame] = q_to_a(popt[1],miller)
            lattice2[frame] = q_to_a(popt[4],miller)
            lattice3[frame] = q_to_a(popt[7],miller)
            phase1[frame] = sum(gaussian(q_sub,*popt[0:3]))
            phase2[frame] = sum(gaussian(q_sub,*popt[3:6]))
            phase3[frame] = sum(gaussian(q_sub,*popt[6:]))
        
        else:
            print(frame, "bad fit", residual)
        
        p0=popt
        upper_limit = [popt[0], init_q+.001, .3, 30, init_q+.001, .3, popt[6], popt[7], .3] #upper limits fot the range where one peak is iodine rich and the other is bromine
        lower_limit = [0, popt[1], 0, 0, init_q-0.001, 0, 0, init_q-.001, 0]

else: 
    fits = []
    for frame in range(files): 
        popt, pcov = curve_fit(three_gaussians, q_sub, perov_fit[:,frame+1], p0, bounds=(lower_limit, upper_limit), maxfev=8000)
        lattice1[frame] = q_to_a(popt[1],miller)
        lattice2[frame] = q_to_a(popt[4],miller)
        lattice3[frame] = q_to_a(popt[7],miller)
        phase1[frame] = sum(gaussian(q_sub,*popt[0:3]))
        phase2[frame] = sum(gaussian(q_sub,*popt[3:6]))
        phase3[frame] = sum(gaussian(q_sub,*popt[6:]))
        chem1[frame] = q_to_chem(popt[1],miller)
        chem2[frame] = q_to_chem(popt[4],miller)
        chem3[frame] = q_to_chem(popt[7],miller)
        time[frame] = frame*2
        fits.append(popt)
        print('Intensity:', popt[0])
        print('Lattice Spacing:', q_to_a(popt[1],miller))
        p0=popt
        upper_limit = [popt[0], init_q, .3, 30, init_q+.001, .3, popt[6], popt[7], .3] #upper limits fot the range where one peak is iodine rich and the other is bromine
        lower_limit = [0, q_1, 0, 0, init_q-0.001, 0, 0, init_q, 0]

#   Plot lattice spacings over time
plt.figure(figsize=(5,4)) #make plot larger
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Lattice Spacing [$\AA$]',size=14)#Define y-axis label
plt.plot(time,lattice1,'ro')
plt.plot(time,lattice2, 'ko')
plt.plot(time, lattice3, 'bo')
plt.ylim((5.8,6.3))

#%% plot arbitraty frame
choice = 0
# Plot everything
popt = fits[choice]
plt.plot(q_sub,perov_fit[:,choice],'b-', label='Data') #plot subfield of data
plt.plot(q_sub,three_gaussians(q_sub, *popt),'c--', label='Model') #plot best fit
plt.plot(q_sub, gaussian(q_sub, popt[0],popt[1],popt[2]), 'b-.', label='Peak 1')
plt.plot(q_sub, gaussian(q_sub, popt[3],popt[4],popt[5]), 'b--', label='Peak 2')
plt.plot(q_sub, gaussian(q_sub, popt[6],popt[7],popt[8]), 'b:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner

#%% plot series of frames
figtime = plt.figure(figsize=(6, 4))
times = np.array([0, 10, 30, 60, 90, 118]) #list of times you would like to look at
for i in range(len(times)):
    plt.plot(q_sub,perov_fit[:,int(times[i]/interval)])


plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner



# %% PLOTTING OF DYNAMICS  
# Generic Comparative Plot 
index = 3 # chose which of the following [time, maxintensity, iod_bin, orig_bin, brom_bin, iod_q, orig_q, brom_q, av_q]
normalized = 1 #if want data normalized chose 1, else chose 0

if normalized:
    plt.plot(results[samplelist[0]][0],results[samplelist[0]][index]/max(results[samplelist[0]][index]), 'b.--', label = 'x=0.33')
    plt.plot(results[samplelist[1]][0],results[samplelist[1]][index]/max(results[samplelist[1]][index]), 'r.--',  label = 'x=0.50')
    plt.plot(results[samplelist[2]][0],results[samplelist[2]][index]/max(results[samplelist[2]][index]), 'g.--',  label = 'x=0.67')
    plt.plot(results[samplelist[3]][0],results[samplelist[3]][index]/max(results[samplelist[3]][index]), 'm.--',  label = 'x=0.75')
else: 
    plt.plot(results[samplelist[0]][0],results[samplelist[0]][index], 'b.--', label = 'x=0.33')
    plt.plot(results[samplelist[1]][0],results[samplelist[1]][index], 'r.--',  label = 'x=0.50')
    plt.plot(results[samplelist[2]][0],results[samplelist[2]][index], 'g.--',  label = 'x=0.67')
    plt.plot(results[samplelist[3]][0],results[samplelist[3]][index], 'm.--',  label = 'x=0.75')

plt.xlabel('Time [min]',size=14) #Define x-axis label
#plt.ylabel('Max Intensity [Normalized]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner

# %% Generic Comparative Plot Ligh On and Off
index = 4 # chose which of the following [time, maxintensity, iod_bin, orig_bin, brom_bin, iod_q, orig_q, brom_q, av_q]
normalized = 1 #if want data normalized chose 1, else chose 0
dark_shift = 1
dark_start_time = 60

if normalized:
    chem = 0
    plt.plot(results[samplelist[chem]][0],results[samplelist[chem]][index]/results[samplelist[chem]][index][0], 'b.-', label = 'x=0.33')
    plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,results[samplelist[chem+dark_shift]][index]/results[samplelist[chem]][index][0], 'b.-')
    
    chem = 2
    plt.plot(results[samplelist[chem]][0],results[samplelist[chem]][index]/results[samplelist[chem]][index][0], 'r.-', label = 'x=0.5')
    plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,results[samplelist[chem+dark_shift]][index]/results[samplelist[chem]][index][0], 'r.-')
    
    chem = 4
    plt.plot(results[samplelist[chem]][0],results[samplelist[chem]][index]/results[samplelist[chem]][index][0], 'g.-', label = 'x=0.67')
    plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,results[samplelist[chem+dark_shift]][index]/results[samplelist[chem]][index][0], 'g.-')
    
    chem = 6
    plt.plot(results[samplelist[chem]][0],results[samplelist[chem]][index]/results[samplelist[chem]][index][0], 'm.-', label = 'x=0.75')
    plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,results[samplelist[chem+dark_shift]][index]/results[samplelist[chem]][index][0], 'm.-')

else: 
    plt.plot(results[samplelist[0]][0],results[samplelist[0]][index], 'b.--', label = 'x=0.33')
    plt.plot(results[samplelist[1]][0],results[samplelist[1]][index], 'r.--',  label = 'x=0.50')
    plt.plot(results[samplelist[2]][0],results[samplelist[2]][index], 'g.--',  label = 'x=0.67')
    plt.plot(results[samplelist[3]][0],results[samplelist[3]][index], 'm.--',  label = 'x=0.75')

plt.xlabel('Time [min]',size=14) #Define x-axis label
#plt.ylabel('Max Intensity [Normalized]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Original Peak [Normalized]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
plt.xlim(0, 180)
#plt.ylim(0, 1.1)
#fig.savefig('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/'+'originalpeak.png', bbox_inches='tight', dpi=900, transparent=True)

# %% Total Intensity plot 
def totalcounts(chem):
    return results[samplelist[chem]][2]+results[samplelist[chem]][3]+results[samplelist[chem]][4]
dark_shift = 1
dark_start_time = 60
fig = plt.figure(figsize=(8, 5))
plt.tick_params(right=True, labelright=True)


chem = 0
plt.plot(results[samplelist[chem]][0],totalcounts(chem)/totalcounts(chem)[0], 'b.-', label = 'x=0.33')
plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,totalcounts(chem+dark_shift)/totalcounts(chem)[0], 'b.-')
chem = 2
plt.plot(results[samplelist[chem]][0],totalcounts(chem)/totalcounts(chem)[0], 'r.-', label = 'x=0.5')
plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,totalcounts(chem+dark_shift)/totalcounts(chem)[0], 'r.-')
chem = 4
plt.plot(results[samplelist[chem]][0],totalcounts(chem)/totalcounts(chem)[0], 'g.-', label = 'x=0.67')
plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,totalcounts(chem+dark_shift)/totalcounts(chem)[0], 'g.-')
chem = 6
plt.plot(results[samplelist[chem]][0],totalcounts(chem)/totalcounts(chem)[0], 'm.-', label = 'x=0.75')
plt.plot(results[samplelist[chem+dark_shift]][0]+dark_start_time,totalcounts(chem+dark_shift)/totalcounts(chem)[0], 'm.-')


plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Total Counts [Normalized]',size=14)#Define y-axis label
plt.legend(loc="lower right")#Put legend in upper left hand corner
plt.xlim(0, 180)
plt.ylim(0, 1.1)
#fig.savefig('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/'+'totalintensity.png', bbox_inches='tight', dpi=900, transparent=True)
# %% Original Peak Intensity - perhaps the best for looking at the decay of the original peak
#   Though likely complicated by how closely the peaks overlap?
plt.plot(results[samplelist[0]][0],results[samplelist[0]][3]/max(results[samplelist[0]][3]), 'b.--', label = 'x=0.33')
plt.plot(results[samplelist[1]][0],results[samplelist[1]][3]/max(results[samplelist[1]][3]), 'r.--',  label = 'x=0.50')
plt.plot(results[samplelist[2]][0],results[samplelist[2]][3]/max(results[samplelist[2]][3]), 'g.--',  label = 'x=0.67')
plt.plot(results[samplelist[3]][0],results[samplelist[3]][3]/max(results[samplelist[3]][3]), 'm.--',  label = 'x=0.75')
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Original Peak Intensity [Normalized]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner


# %%    Bromine Peak Intensity
plt.plot(results[samplelist[0]][0],results[samplelist[0]][4]/max(results[samplelist[0]][4]), 'b.--', label = 'x=0.33')
plt.plot(results[samplelist[1]][0],results[samplelist[1]][4]/max(results[samplelist[1]][4]), 'r.--',  label = 'x=0.50')
plt.plot(results[samplelist[2]][0],results[samplelist[2]][4]/max(results[samplelist[2]][4]), 'g.--',  label = 'x=0.67')
plt.plot(results[samplelist[3]][0],results[samplelist[3]][4]/max(results[samplelist[3]][4]), 'm.--',  label = 'x=0.75')
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Bromine Peak Intensity [Normalized]',size=14)#Define y-axis label
plt.legend(loc="lower right")#Put legend in upper left hand corner

# %%
# %% Iodide Peak Intensity - perhaps the best for looking at the decay of the original peak
#   Though likely complicated by how closely the peaks overlap?
plt.plot(results[samplelist[0]][0],results[samplelist[0]][3]/max(results[samplelist[0]][3]), 'b.--', label = 'x=0.33')
plt.plot(results[samplelist[1]][0],results[samplelist[1]][3]/max(results[samplelist[1]][3]), 'r.--',  label = 'x=0.50')
plt.plot(results[samplelist[2]][0],results[samplelist[2]][3]/max(results[samplelist[2]][3]), 'g.--',  label = 'x=0.67')
plt.plot(results[samplelist[3]][0],results[samplelist[3]][3]/max(results[samplelist[3]][3]), 'm.--',  label = 'x=0.75')
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Original Peak Intensity [Normalized]',size=14)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner
# %% Dynamics of single composition
chemistry = 6
dark_shift = 1
fig2 = plt.figure(figsize=(6, 4))
#   Light on 
original_intensity = (results[samplelist[chemistry]][4][0]+results[samplelist[chemistry]][3][0]+results[samplelist[chemistry]][2][0])
plt.plot(results[samplelist[chemistry]][0],(results[samplelist[chemistry]][3]+results[samplelist[chemistry]][2]+results[samplelist[chemistry]][4])/(original_intensity), 'k-', label = 'Total Intensity')
plt.plot(results[samplelist[chemistry]][0],(results[samplelist[chemistry]][3])/original_intensity, 'k.--', label = 'Original Peak')
plt.plot(results[samplelist[chemistry]][0],(results[samplelist[chemistry]][2])/original_intensity, 'r.--', label = 'Iodide Peak')
plt.plot(results[samplelist[chemistry]][0],(results[samplelist[chemistry]][4])/original_intensity, 'b.--', label = 'Bromide Peak')

#   Light off
plt.plot(results[samplelist[chemistry+dark_shift]][0]+60,(results[samplelist[chemistry+dark_shift]][3]+results[samplelist[chemistry+dark_shift]][2]+results[samplelist[chemistry+dark_shift]][4])/(original_intensity), 'k-', label = 'Total Intensity')
plt.plot(results[samplelist[chemistry+dark_shift]][0]+60,(results[samplelist[chemistry+dark_shift]][3])/original_intensity, 'k.--', label = 'Original Peak')
plt.plot(results[samplelist[chemistry+dark_shift]][0]+60,(results[samplelist[chemistry+dark_shift]][2])/original_intensity, 'r.--', label = 'Iodide Peak')
plt.plot(results[samplelist[chemistry+dark_shift]][0]+60,(results[samplelist[chemistry+dark_shift]][4])/original_intensity, 'b.--', label = 'Bromide Peak')


plt.title(samplelist[chemistry], size=16)
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Peak Intensities [Normalized]',size=14)#Define y-axis label
fig2.savefig('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/'+samplelist[chemistry]+'.png', bbox_inches='tight', dpi=900, transparent=True)
# %% fit decay and rise
bin_scale = 1
chemistry = 0
dark_shift = 1
long_time = np.arange(0, 300, 1)
def exp_decay(x, a, b, c):
    return a*np.exp(-x/b)+c

def exp_rise(x, a, b, c):
    return a*(1-np.exp(-x/b))+c

#decay
p0 = [1.5, 20, .1]
popt, pcov = curve_fit(exp_decay, results[samplelist[chemistry]][0], results[samplelist[chemistry]][3]/max(results[samplelist[chemistry]][3]), p0=p0, maxfev=8000)

fig1 = plt.figure(figsize=(8, 5))
plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][3]/max(results[samplelist[chemistry]][3]), 'k.-', label = 'Max Intensity')
plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*popt))
decay = popt

#rise
p0 = [.7, 50, .3]
popt, pcov = curve_fit(exp_rise, results[samplelist[chemistry+dark_shift]][0], results[samplelist[chemistry+dark_shift]][3]*bin_scale/max(results[samplelist[chemistry]][3]), p0 = p0, maxfev=8000)
fig2 = plt.figure(figsize=(8, 5))
plt.plot(results[samplelist[chemistry+dark_shift]][0],results[samplelist[chemistry+dark_shift]][3]*bin_scale/max(results[samplelist[chemistry]][3]), 'k.-', label = 'Max Intensity')
plt.plot(results[samplelist[chemistry+dark_shift]][0], exp_rise(results[samplelist[chemistry+dark_shift]][0],*popt))
rise = popt
#combine plot
fig3 = plt.figure(figsize=(6, 4))
plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][3]/max(results[samplelist[chemistry]][3]), 'k.-', label = 'Max Intensity')
plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*decay))
plt.plot(results[samplelist[chemistry+dark_shift]][0]+60,results[samplelist[chemistry+dark_shift]][3]*bin_scale/max(results[samplelist[chemistry]][3]), 'k.-', label = 'Max Intensity')
plt.plot(long_time+60, exp_rise(long_time,*rise))
# %%
# %% fit decay and rise for all

long_time = np.arange(0, 180, 1)
rise_times = []
rise_error = []
decay_times = []
decay_error = []
chemistry = 0
i = 0

while i < 8:

    if (i%2)==0:
        #decay
        p0 = [1, 15, .1]
        popt, pcov = curve_fit(exp_decay, results[samplelist[chemistry]][0], results[samplelist[chemistry]][3]/max(results[samplelist[chemistry]][3]), p0=p0, maxfev=8000)

        fig1 = plt.figure(figsize=(8, 5))
        plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][3]/max(results[samplelist[chemistry]][3]), 'k.-', label = 'Max Intensity')
        plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*popt))
        decay = popt
        decay_times.append(decay)
        decay_error.append(np.sqrt(np.diag(pcov)))

    else:
        #rise
        p0 = [.7, 100, .3]
        popt, pcov = curve_fit(exp_rise, results[samplelist[chemistry]][0], results[samplelist[chemistry]][3]*bin_scale/max(results[samplelist[chemistry-1]][3]), p0 = p0, maxfev=8000)
        fig2 = plt.figure(figsize=(8, 5))
        plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][3]*bin_scale/max(results[samplelist[chemistry-1]][3]), 'k.-', label = 'Max Intensity')
        plt.plot(results[samplelist[chemistry]][0], exp_rise(results[samplelist[chemistry]][0],*popt))
        rise = popt
        rise_times.append(rise)
        rise_error.append(np.sqrt(np.diag(pcov)))

    print(samplelist[chemistry])
    print(decay_error)
    print(rise_error)
    chemistry += 1
    i += 1
    

# %% Figure with all fits or original peak intensity
fig4 = plt.figure(figsize=(7, 5))
plt.xlabel('Time [min]',size=14) #Define x-axis label
plt.ylabel('Original Peak Intensity [Normalized]',size=14)#Define y-axis label

plot_light = 1#if 1 plot light on
plot_dark = 1 #if 1 plot dark on
dark_shift = 1
index = 3

if plot_light:
    chemistry = 0
    plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][index]/max(results[samplelist[chemistry]][index]), 'b.', label = 'x=0.33')
    plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*decay_times[0]), 'b--')

    chemistry = 2
    plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][index]/max(results[samplelist[chemistry]][index]), 'r.',  label = 'x=0.50')
    plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*decay_times[1]), 'r--')

    chemistry = 4
    plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][index]/max(results[samplelist[chemistry]][index]), 'g.',  label = 'x=0.67')
    plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*decay_times[2]), 'g--')

    chemistry = 6      
    plt.plot(results[samplelist[chemistry]][0],results[samplelist[chemistry]][index]/max(results[samplelist[chemistry]][index]), 'm.',  label = 'x=0.75')
    plt.plot(results[samplelist[chemistry]][0], exp_decay(results[samplelist[chemistry]][0],*decay_times[3]), 'm--')
    
    shift = 60 

if plot_dark:
    chemistry = 0
    plt.plot(results[samplelist[chemistry+dark_shift]][0]+shift,results[samplelist[chemistry+dark_shift]][index]/max(results[samplelist[chemistry]][index]), 'b.')
    plt.plot(long_time+shift, exp_rise(long_time,*rise_times[0]),'b--')

    chemistry = 2
    plt.plot(results[samplelist[chemistry+dark_shift]][0]+shift,results[samplelist[chemistry+dark_shift]][index]*bin_scale/max(results[samplelist[chemistry]][index]), 'r.')
    plt.plot(long_time+shift, exp_rise(long_time,*rise_times[1]),'r--')

    chemistry = 4
    plt.plot(results[samplelist[chemistry+dark_shift]][0]+shift,results[samplelist[chemistry+dark_shift]][index]/max(results[samplelist[chemistry]][index]), 'g.')
    plt.plot(long_time+shift, exp_rise(long_time,*rise_times[2]),'g--')

    chemistry = 6      
    plt.plot(results[samplelist[chemistry+dark_shift]][0]+shift,results[samplelist[chemistry+dark_shift]][index]/max(results[samplelist[chemistry]][index]), 'm.')
    plt.plot(long_time+shift, exp_rise(long_time,*rise_times[3]),'m--')
plt.ylim(0, 1.1)
plt.xlim(0, 240)
fig4.savefig('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/'+'originalpeal_fits.png', bbox_inches='tight', dpi=900, transparent=True)

# %%
fig5 = plt.figure(figsize=(7, 5))
plt.xlabel('Bromine Francion',size=14) #Define x-axis label
plt.ylabel('Rate Constant [min]',size=14)#Define y-axis label
chem = [.33, .5, .6, .75]

rise_time = [rise_times[0][1], rise_times[1][1], rise_times[2][1], rise_times[3][1]]
decay_time = [decay_times[0][1], decay_times[1][1], decay_times[2][1], decay_times[3][1]]
rise_errors = [rise_error[0][1], rise_error[1][1], rise_error[2][1], rise_error[3][1]]
decay_errors = [decay_error[0][1], decay_error[1][1], decay_error[2][1], decay_error[3][1]]

#plt.plot(chem,rise_time, 'bo--')
#plt.plot(chem, decay_time, 'ro--')
plt.errorbar(chem, decay_time, decay_errors, fmt = 'bo--', label = 'Light')
plt.errorbar(chem, rise_time, rise_errors, fmt = 'ro--', label = 'Dark')
plt.legend(loc="upper right")#Put legend in upper left hand corner
fig5.savefig('/Users/rbelisle/Desktop/dynamics_on_off_750bins_wedge/'+'rate.png', bbox_inches='tight', dpi=900, transparent=True)


# %%
