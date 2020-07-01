#!/usr/bin/env python
# coding: utf-8

# In[36]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.optimize import curve_fit
import pandas as pd


# In[37]:


def readcsv(filename):
    data = pd.read_csv(filename)
    return(np.array(data))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 

def two_gaussians(x, h1, c1, w1, h2, c2, w2):
        return (gaussian(x, h1, c1, w1) +
            gaussian(x, h2, c2, w2))

def three_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3):
        return (gaussian(x, h1, c1, w1) +
            gaussian(x, h2, c2, w2)+gaussian(x, h3, c3, w3))

def two_to_q(two_theta, wave):
    #two_theta is a 1D array of two_theta angles
    #wave is the X-ray energy in angstroms
    rad_theta = two_theta/2*np.pi/180
    q = 4*np.pi*np.sin(rad_theta)/wave
    return q
def trim_data(x, data, limit1, limit2):
    #x is a 1D array of two theta or q values
    #data is an array of x-ray intensities
    #limit1 and limit2 are what you'd like to trime your data to 
    set1 = find_nearest(x,limit1)
    set2 = find_nearest(x,limit2)
    return x[set1:set2], data[set1:set2,:]


# In[310]:


perov = readcsv('E1MAPbBr75L60On.csv') #openfile of interest
plt.plot(perov[:,0],perov[:,1]) #plot inititial data
plt.title('Initial:')
plt.xlabel('2-theta')
plt.ylabel('Intensity')
plt.show()
size = perov.shape
print (size)


# In[311]:


q = two_to_q(perov[:,0],0.9763)
Q= q.tolist()
y=perov[:,1]
plt.figure(figsize=(8,6)) 
plt.plot(Q,y, marker='.',color='r')
plt.title('Initial:')
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.show()


# In[312]:


#choose a peak and find its limits 
#plt.plot(Q,y, marker='.',color='r')
q_1 = 1
q_2 = 1.15
limit1 = Q.index(find_nearest(q, q_1))
limit2 = Q.index(find_nearest(q, q_2))
q_sub = q[limit1:limit2]
perov_sub = perov[limit1:limit2,1:-1]# range will depend on the number of frames you have 
plt.plot(q_sub,perov_sub[:,-1])


# In[313]:


#%%
#remove background
size = perov_sub.shape
print (size)
q_bins = size[0]
num_frames = size[1]
slope = np.zeros((num_frames,1))
b = np.zeros((num_frames,1))
back = np.zeros((q_bins,num_frames))
int_correct = np.zeros((q_bins,num_frames))


# In[314]:


#accept a data file and range and returns average values at start and end of range
for j in range(num_frames): 
    slope[j] = ((np.mean(perov_sub[-10:-1,j])-np.mean(perov_sub[0:10,j]))/(np.mean(q[limit2-10:limit2])-np.mean(q[limit1:limit1+10])))
    b[j]=perov_sub[0,j]-slope[j]*q_sub[0]
    back[:,j] = [slope[j]*element+b[j] for element in q_sub]
    int_correct[:,j] = perov_sub[:,j]-np.array(back[:,j])

plt.plot(np.array(q_sub),int_correct)


# In[315]:


import math
p0 = [400, 1.06, 0.02]   #third parameter length of half peak width/standard deviation
intensity_1 = np.zeros((num_frames)) #create correct size arrays for running in the loop
intensity_2= np.zeros((num_frames))
lattice_1= np.zeros((num_frames))
lattice_2 = np.zeros((num_frames))
for j in range(num_frames):
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0)
    intensity_1[j] = popt[0]
    lattice_1[j] = 2*math.pi/popt[1] 
    #intensity_2[j] = popt[3]
    #lattice_2[j] = 4*math.pi/popt[4]
    p0 = popt #once you have the initial fit use that as your next guess, we expect things to be continuous so this helps with that
    #print (j)
#print(lattice_1[1:-1])
print (lattice_1)


# In[316]:


#%%
intensity_1 = np.zeros((num_frames))
intensity_2= np.zeros((num_frames))
intensity_3= np.zeros((num_frames))
lattice_1= np.zeros((num_frames))
lattice_2 = np.zeros((num_frames))
lattice_3 = np.zeros((num_frames))
# %%


# In[317]:


#guassian fits with plotting
#pick fitting and initial guess
# it would be good to iteratively do these guesses, or use info from qlimit
#real issue in fitting - maybe this is four peaks? Perhaps I can look at inflection point instead
i = 2
j = 0
p0 =  [400, 1.06, 0.02, 200,1.061, 0.03]
plt.plot(np.array(q_sub), int_correct[:,j],color='black')
if i == 2:
    popt, pcov = curve_fit(two_gaussians, np.array(q_sub), int_correct[:,j], p0[0:6])
    #plt.plot(np.array(q_sub), two_gaussians(np.array(q_sub),*popt))
    pars_1 = popt[0:3]
    pars_2 = popt[3:6]
    gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
    gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
    plt.plot(np.array(q_sub), gauss_peak_1, color='red')
    plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
    plt.plot(np.array(q_sub), gauss_peak_1+gauss_peak_2,color='green')
    print('lattice spacing:', [4*math.pi/popt[1]/2, 4*math.pi/popt[4]/2])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    intensity_2[j] = popt[3]
    lattice_2[j] = 4*math.pi/popt[4]

if i == 3:
    popt, pcov = curve_fit(three_gaussians, np.array(q_sub), int_correct[:,j], p0)
    #plt.plot(np.array(q_sub), three_gaussians(np.array(q_sub),*popt))
    pars_1 = popt[0:3]
    pars_2 = popt[3:6]
    pars_3 = popt[6:9]
    gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
    gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
    gauss_peak_3 = gaussian(np.array(q_sub), *pars_3)
    plt.plot(np.array(q_sub), gauss_peak_1, color='magenta')
    plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
    plt.plot(np.array(q_sub),gauss_peak_3,color='green')
    plt.plot(np.array(q_sub),gauss_peak_3+gauss_peak_1+gauss_peak_2,color='red',linestyle='dashed')
    print('lattice spacing:', [4*math.pi/popt[1]/2, 4*math.pi/popt[4]/2,4*math.pi/popt[7]/2])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    intensity_2[j] = popt[3]
    lattice_2[j] = 4*math.pi/popt[4]
    intensity_2[j] = popt[6]
    lattice_3[j] = 4*math.pi/popt[7]

else: 
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0[0:3])
    plt.plot(np.array(q_sub), gaussian(np.array(q_sub),*popt))
    print('lattice spacing:', [4*math.pi/popt[1]])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1]

print(popt)
#plt.savefig('time_42fit.png')
#popt_uncertainties = np.sqrt(np.diag(pcov))
#uncertainty = sum(popt_uncertainties)
#print('uncertainties:', popt_uncertainties)
#print('uncertainty:', uncertainty)

plt.xlabel('Q')
plt.ylabel('Intensity')




# In[318]:


plt.plot(np.array(q_sub), int_correct[:,j],label = "initial",color='red')
plt.plot(np.array(q_sub), gauss_peak_1+gauss_peak_2,color='green')
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.title("MaPbBr75")
plt.legend() 


# In[319]:


print('Intensity:', popt[0])

#Caculate and pring d-spacing
d = 2*np.pi/popt[1] #Applying d = 2*pi/Q
print('d-Spacing: ', d) 

#Print lattice constant
miller = [1, 0, 0] #need to guess miller indices of peak
a = d/np.sqrt(miller[0]**2+miller[1]**2+miller[2]**2) #calculate a using a = d/sqrt(h^2+k^2+l^2) for a cubic lattice
print('Lattice Spacing:', a)


# In[320]:


i = 3
j = -1
p0 =  [122, 1.06, 0.02, 61,1.059, 0.03,30.5,1.058,0.05]
plt.plot(np.array(q_sub), int_correct[:,j],color='black')


# 
# 
# 

# In[336]:


## guassian fits with plotting
#pick fitting and initial guess
# it would be good to iteratively do these guesses, or use info from qlimit
#real issue in fitting - maybe this is four peaks? Perhaps I can look at inflection point instead
i = 3
j = 28
p0 = [120, 1.06, 0.02,62,1.06, 0.03, 30.5, 1.06, 0.05]
plt.plot(np.array(q_sub), int_correct[:,j],color='black')
if i == 2:
    popt, pcov = curve_fit(two_gaussians, np.array(q_sub), int_correct[:,j], p0[0:6])
    #plt.plot(np.array(q_sub), two_gaussians(np.array(q_sub),*popt))
    pars_1 = popt[0:3]
    pars_2 = popt[3:6]
    gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
    gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
    plt.plot(np.array(q_sub), gauss_peak_1, color='red')
    plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
    plt.plot(np.array(q_sub), gauss_peak_1+gauss_peak_2,color='green')
    print('lattice spacing:', [4*math.pi/popt[1]/2, 4*math.pi/popt[4]/2])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    intensity_2[j] = popt[3]
    lattice_2[j] = 4*math.pi/popt[4]

if i == 3:
    popt, pcov = curve_fit(three_gaussians, np.array(q_sub), int_correct[:,j], p0)
    #plt.plot(np.array(q_sub), three_gaussians(np.array(q_sub),*popt))
    pars_1 = popt[0:3]
    pars_2 = popt[3:6]
    pars_3 = popt[6:9]
    gauss_peak_1 = gaussian(np.array(q_sub), *pars_1)
    gauss_peak_2 = gaussian(np.array(q_sub), *pars_2)
    gauss_peak_3 = gaussian(np.array(q_sub), *pars_3)
    plt.plot(np.array(q_sub), gauss_peak_1, color='magenta')
    plt.plot(np.array(q_sub),gauss_peak_2,color='blue') 
    plt.plot(np.array(q_sub),gauss_peak_3,color='green')
    plt.plot(np.array(q_sub),gauss_peak_3+gauss_peak_1+gauss_peak_2,color='red',linestyle='dashed')
    print('lattice spacing:', [4*math.pi/popt[1]/2, 4*math.pi/popt[4]/2,4*math.pi/popt[7]/2])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1] 
    intensity_2[j] = popt[3]
    lattice_2[j] = 4*math.pi/popt[4]
    intensity_2[j] = popt[6]
    lattice_3[j] = 4*math.pi/popt[7]

else: 
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0[0:3])
    plt.plot(np.array(q_sub), gaussian(np.array(q_sub),*popt))
    print('lattice spacing:', [4*math.pi/popt[1]/2])
    intensity_1[j] = popt[0]
    lattice_1[j] = 4*math.pi/popt[1]

print(popt)
#plt.savefig('time_42fit.png')
#popt_uncertainties = np.sqrt(np.diag(pcov))
#uncertainty = sum(popt_uncertainties)
#print('uncertainties:', popt_uncertainties)
#print('uncertainty:', uncertainty)

plt.xlabel('Q')
plt.ylabel('Intensity')


# In[337]:


plt.plot(np.array(q_sub), int_correct[:,j],label = "final",color='black')
#plt.plot(np.array(q_sub), gauss_peak_1+gauss_peak_2,color='green')
plt.plot(np.array(q_sub),gauss_peak_3+gauss_peak_1+gauss_peak_2,color='red',linestyle='dashed')
plt.xlabel('Q')
plt.ylabel('Intensity')
plt.title("MaPbBr75")
plt.legend() 


# In[338]:


print('Intensity:', popt[0])

#Caculate and pring d-spacing
d = 2*np.pi/popt[1] #Applying d = 2*pi/Q
print('d-Spacing: ', d) 

#Print lattice constant
miller = [1, 0, 0] #need to guess miller indices of peak
a = d/np.sqrt(miller[0]**2+miller[1]**2+miller[2]**2) #calculate a using a = d/sqrt(h^2+k^2+l^2) for a cubic lattice
print('Lattice Spacing:', a)


# In[ ]:




