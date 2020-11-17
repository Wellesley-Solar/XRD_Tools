#%%This is a program for calculating the lead iodide to perovskite (100) peak ratio over time
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
from xrdfunctions import *
#%% Import Data
perov_import = csv_to_np('/Users/rbelisle/Desktop/Xraydeg/CSV/e1_deg.csv')

#%% Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

#%% #Trim to isolate first peaks and Remove Background
miller = [1, 0, 0] #peak you are trying to capture
q_1 = 0.8
q_2 = 1.15
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

# %%
files = num_files(perov)
time = np.zeros(files)
lead = np.zeros(files)
for frame in range(files): 
    time[frame] = frame*20
    lead[frame] = perovtopbi2(q_sub, perov_fit[:,frame])

# %% Plot result over time
plt.figure(num = 1, figsize=(5,4))
fig1,ax1 = plt.subplots()
ax1.set_xlabel('Time [s]',size=14) #Define x-axis label
ax1.set_ylabel('$Perovskite$ (100)/$PbI_2$ (001)',size=14)#Define y-axis label
plt.plot(time,lead/(lead[0]), 'ko')

