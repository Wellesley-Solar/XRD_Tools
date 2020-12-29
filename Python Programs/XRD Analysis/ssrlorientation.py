#%% numpy and plotting
import numpy as np
import pylab 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 

# pyFAI
import pyFAI
# pygix
import pygix
import fabio
# pandas
import pandas as pd
#%% Import default file 
dataFile = '/Users/rbelisle/Desktop/5050onoff/D7a_MAPbIBr2_PTAA_L60On_30s_01231947_0001.tif'
data = fabio.open(dataFile).data
fig1, ax1 = plt.subplots()
ax1.set_xlabel("x-pixel (#)")
ax1.set_ylabel("y-pixel (#)") 
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.imshow(data, vmin = 3, vmax = 100, origin='lower')
# %%
