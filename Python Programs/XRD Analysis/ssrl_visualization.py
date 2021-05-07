#this is a file for making plots of XRD over time
#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colorbar
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from xrdfunctions import *

#%% Import Data
perov_import = csv_to_np('/Users/rbelisle/Desktop/5050onoff/lighton.csv')

#%% Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 

#%% #Trim and Remove Background
miller = [2, 0, 0] #peak you are trying to capture
q_1 = 1.96
q_2 = 2.15
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

# %% Option 1 
fig1, ax1 = plt.subplots()
ax1.set_ylabel("q ($\mathring{A}^\mathregular{-1}$)")
ax1.set_xlabel("time [s]") 
#ax1.imshow(perov_fit, vmin = 5, vmax = 50, origin='lower') 

ax1.imshow(perov_fit, cmap = 'jet', extent=(0,60,min(q_sub),max(q_sub)), aspect = 'auto', vmin = 3, vmax = 40, origin='lower') 
# %% Option 2
perov_2 = np.transpose(perov_fit)

fig1, ax1 = plt.subplots()
ax1.set_xlabel("q ($\mathring{A}^\mathregular{-1}$)")
ax1.set_ylabel("time [s]") 
ax1.imshow(perov_2, cmap = 'jet', extent=(min(q_sub),max(q_sub),0,60), aspect = 'auto', vmin = 1, vmax = 30, origin='lower') 
# %%
perov_import = csv_to_np('/Users/rbelisle/Desktop/5050onoff/darkon.csv')
#   Covert to Q
q = two_to_q(perov_import[:,0],0.982381)
perov = perov_import[:,1:] #seperate XRD intensities from rest of data 
#   Trim and Remove Background
q_sub, perov_sub = trim_data(q,perov,q_1,q_2)
perov_fit2 = perov_sub
files = num_files(perov_sub)
for file in range(files): 
    perov_fit2[:,file] = back_subtract(q_sub,perov_sub[:,file],10) #remove background from that file

# %%
test = np.concatenate((perov_fit,perov_fit2),axis=1) 

# %% Option 2
perov_plot = np.transpose(test)

fig1, ax1 = plt.subplots()
ax1.set_xlabel("Q [$\mathring{A}^\mathregular{-1}$]")
ax1.set_ylabel("Time [s]") 
ax1.imshow(perov_plot, cmap = 'jet', extent=(min(q_sub),max(q_sub),0,240), aspect = 'auto', vmin = 1, vmax = 40, origin='lower') 
# %%
# %% Option 2
perov_2 = np.transpose(perov_fit2)

fig1, ax1 = plt.subplots()
ax1.set_xlabel("q ($\mathring{A}^\mathregular{-1}$)")
ax1.set_ylabel("time [s]") 
ax1.imshow(perov_2, cmap = 'jet', extent=(min(q_sub),max(q_sub),60,240), aspect = 'auto', vmin = 1, vmax = 30, origin='lower') 
# %%
fig3, ax3 = plt.subplots(nrows=2, sharex=True, figsize=(5, 5))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

ax3[1].imshow(np.transpose(perov_fit), cmap = 'jet', extent=(min(q_sub),max(q_sub),0,60), aspect = 'auto', vmin = 1, vmax = 40, origin='lower') 

ax3[0].imshow(np.transpose(perov_fit2), cmap = 'jet', extent=(min(q_sub),max(q_sub),60,240), aspect = 'auto', vmin = 1, vmax = 40, origin='lower') 
ax3[1].set_xlabel("Q ($\mathring{A}^\mathregular{-1}$)")
plt.show()
# need to make it so that in dark is scaled 3X fatter
# %%
color = 'jet'
fig = plt.figure(figsize=(6, 8)) 
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
ax0 = plt.subplot(gs[0])
ax0.imshow(np.transpose(perov_fit2), cmap = color, extent=(min(q_sub),max(q_sub),60,240), aspect = 'auto', vmin = 1, vmax = 35, origin='lower') 
ax0.set_ylabel("Light Off [min]")
ax0.tick_params(labelbottom=False)   
ax1 = plt.subplot(gs[1])
ax1.imshow(np.transpose(perov_fit), cmap = color, extent=(min(q_sub),max(q_sub),0,60), aspect = 'auto', vmin = 1, vmax = 35, origin='lower') 
ax1.set_xlabel("Q [$\mathring{A}^\mathregular{-1}$]")
ax1.set_ylabel("Light On [min]")
ax1.xaxis.set_ticks_position('both')
# %%
fig.savefig('/Users/rbelisle/Desktop/5050onoff/overtime.png',bbox_inches='tight', dpi=900, transparent=True)
# %%
# %%
color = 'jet'
fig = plt.figure(figsize=(8, 5)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3]) 
ax0 = plt.subplot(gs[1])
ax0.imshow(perov_fit2, cmap = color, extent=(60,240, min(q_sub),max(q_sub)), aspect = 'auto', vmin = 1, vmax = 35, origin='lower') 
ax0.set_xlabel("Light Off [min]")
ax0.tick_params(labelleft=False)   
ax1 = plt.subplot(gs[0])
ax1.imshow(perov_fit, cmap = color, extent=(0,60, min(q_sub),max(q_sub)), aspect = 'auto', vmin = 1, vmax = 35, origin='lower') 
ax1.set_ylabel("Q [$\mathring{A}^\mathregular{-1}$]")
ax1.set_xlabel("Light On [min]")
ax1.yaxis.set_ticks_position('both')
# %%
fig.savefig('/Users/rbelisle/Desktop/5050onoff/overtime2.png',bbox_inches='tight', dpi=900, transparent=True)

# %%
