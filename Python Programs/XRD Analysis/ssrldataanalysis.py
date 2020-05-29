#%%
#open file
#determine compositions/sample
#convert to Q
#crop data
#remove background
#fit guassians
#determine if should do 1,2,3 gaussians
#save lattice/intensity/strain and time/composition
#ideally have a function that does the fitting so we can call that function independent if doing
#in situ or start end

import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import pandas as pd
import math
import scipy
# import scipy.signal
#%%
#Open the csv file and plot the first and last frames of the data. Index -1 plots last frame in a data set'''
# this could be improved such that whole data file is imported insteaof specific rows
with open('/Users/rbelisle/Desktop/startendxrd/xrdstartend2.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    x=[]
    y=[]
    z=[]
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            x.append(float(row[0]))
            y.append(float(row[10]))
            z.append(float(row[-1]))
            line_count += 1
plt.plot(x,y, marker='o', color='blue')
plt.plot(x,z, marker='o', color='red')
plt.title('Initial:')
plt.xlabel('2-theta')
plt.ylabel('Intensity')
plt.show()

# %%
# Convert two theta to Q and set data limits
import math
print(math.pi)
q =[4*math.pi*math.sin(math.pi/180*row/2)/0.9763 for row in x]
plt.plot(q,y, marker='.',color='blue')
plt.title('Initial:')
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.show()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

q_1 = 1.9
q_2 = 2.2
limit1 = q.index(find_nearest(q, q_1))
limit2 = q.index(find_nearest(q, q_2))

q_sub = q[limit1:limit2]
int_sub = y[limit1:limit2]
plt.plot(q_sub,int_sub)
#%%
#remove background
slope = (np.mean(y[limit2-10:limit2])-np.mean(y[limit1:limit1+10]))/(np.mean(q[limit2-10:limit2])-np.mean(q[limit1:limit1+10]))

b = y[limit1]-slope*q[limit1]
back = [slope*element+b for element in q]

int_correct = np.array(int_sub)-np.array(back[limit1:limit2])

plt.plot(np.array(q_sub),int_correct)

# %%
#guassian fits
from scipy.optimize import curve_fit

def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 

def two_gaussians(x, h1, c1, w1, h2, c2, w2):
    return (gaussian(x, h1, c1, w1) +
        gaussian(x, h2, c2, w2))

#p0 is the orignial guessess for the files, if we guess that both peaks are in the center still seems to fit. 
#Can likely use this to get guess values - fit original peak for chemistry adnd then use that as guess

popt, pcov = curve_fit(two_gaussians, np.array(q_sub), int_correct, p0=[50, 2.07, 0.03, 70,2.10, 0.01])

plt.plot(np.array(q_sub), int_correct)
plt.plot(np.array(q_sub), two_gaussians(np.array(q_sub),*popt))

pars_1 = popt[0:3]
pars_2 = popt[3:6]
gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)

popt_uncertainties = np.sqrt(np.diag(pcov))
uncertainty = sum(popt_uncertainties)
print('lattice spacing:', [4*math.pi/popt[1], 4*math.pi/popt[4]])
print('uncertainties:', popt_uncertainties)
print('uncertainty:', uncertainty)

plt.plot(np.array(q_sub), gauss_peak_1)
plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
plt.plot(np.array(q_sub),gauss_peak_1+gauss_peak_2,color='green')
plt.xlabel('Q')
plt.ylabel('Intensity')

# %%
#guassian fit (single)
from scipy.optimize import curve_fit

def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 


#p0 is the orignial guessess for the files, if we guess that both peaks are in the center still seems to fit. 
#Can likely use this to get guess values - fit original peak for chemistry adnd then use that as guess

popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct, p0=[80, 2.11, 0.01])

plt.plot(np.array(q_sub), int_correct)
plt.plot(np.array(q_sub), gaussian(np.array(q_sub),*popt))

pars_1 = popt[0:3]
gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)

popt_uncertainties = np.sqrt(np.diag(pcov))
uncertainty = sum(popt_uncertainties)

print('lattice spacing:', 4*math.pi/popt[1])
print('uncertainties:', popt_uncertainties)
print('uncertainty:', uncertainty)

plt.plot(np.array(q_sub), gauss_peak_1)

plt.xlabel('Q')
plt.ylabel('Intensity')

# %%
#three guassian fits
from scipy.optimize import curve_fit

def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 

def three_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3):
    return (gaussian(x, h1, c1, w1) +
        gaussian(x, h2, c2, w2)+gaussian(x, h3, c3, w3))

#p0 is the orignial guessess for the files, if we guess that both peaks are in the center still seems to fit. 
#Can likely use this to get guess values - fit original peak for chemistry adnd then use that as guess

popt, pcov = curve_fit(three_gaussians, np.array(q_sub), int_correct, p0=[5, 2.03, 0.05, 10,2.07, 0.05, 35, 2.14, 0.01])

plt.plot(np.array(q_sub), int_correct)
#plt.plot(np.array(q_sub), three_gaussians(np.array(q_sub),*popt))

pars_1 = popt[0:3]
pars_2 = popt[3:6]
pars_3 = popt[6:9]
gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
gauss_peak_3 = gaussian(np.array(q_sub), *pars_3)

popt_uncertainties = np.sqrt(np.diag(pcov))
uncertainty = sum(popt_uncertainties)
print('lattice spacing:', [4*math.pi/popt[1], 4*math.pi/popt[4],4*math.pi/popt[7]])
print('uncertainties:', popt_uncertainties)
print('uncertainty:', uncertainty)

plt.plot(np.array(q_sub), gauss_peak_1, color='yellow')
plt.plot(np.array(q_sub),gauss_peak_2,color='blue')
plt.plot(np.array(q_sub),gauss_peak_3,color='red')
#plt.plot(np.array(q_sub),gauss_peak_1+gauss_peak_2,color='green')
plt.xlabel('Q [1/Å]')
plt.ylabel('Intensity  [a.u.]')
plt.savefig('50% Bromine after Light.png')



# %%
