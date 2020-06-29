#!/usr/bin/env python
# coding: utf-8

# In[259]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
import pandas as pd


# In[260]:


def csv_to_np(filename): #changing excel to numpy
    data = pd.read_csv(filename)
    return(np.array(data))
perov = csv_to_np(r'C:\Users\tiann\OneDrive\Desktop\Summer Research 2020\D1_MaPbBlBr2_XrayDegGSASExcel.csv')


# In[324]:


plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0],perov[:,1],'r-', label='$MAPbIBr_2$ initial') #plot two-theta versus XRD intesntiy
plt.plot(perov[:,0], perov[:, -1], 'g-', label = '$MAPbIBr_2$ final')
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)
plt.legend(loc="upper right")
plt.title('first and last frames xraydeg')


# In[325]:


plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0],perov[:,1],'r-', label='$MAPbIBr_2$ initial') #plot two-theta versus XRD intesntiy
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)
plt.title('initial')
plt.legend(loc="upper right")


# ### Seeing if I can remove the background of this graph

# In[326]:


plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0], perov[:, 30], 'g-', label = '$MAPbIBr_2$ final')
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)
plt.title('Final')
plt.legend(loc="upper right")


# In[327]:


#Note that our initial 2-theta values are in degrees, not radians
#We also need to define the X-ray wavelength we're using
wave = 0.9763 #wavelength used at SSRL in Angstroms
theta = perov[:,0]/2*np.pi/180 #convert from 2-theta to theta in radians
q = 4*np.pi*np.sin(theta)/wave #convert to q


# In[328]:


#We can now check that our conversion went as expected by ploting our data versus q
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q,perov[:,1],'r-', label='$MAPbIBr_2 Initial') #plot Q versus XRD intesntiy
plt.xlabel('Q [$\AA^{-1}$]',size=12)
plt.ylabel('Intensity [a.u.]',size=12)
plt.legend(loc="upper right")
plt.title('Initial')


# In[329]:


def find_nearest(array, target):
    array = np.asarray(array) # In case the input is not an array type
    idx = (np.abs(array - target)).argmin() # Finds the index of the value closest to the target
    return idx
#We use our new function to find the index of real value closest to our desired limit
q = 4*np.pi*np.sin(theta)/wave
q_1 = 0.98 #lower limit for Q we'll consider
q_2 = 1.3 # upper limit - ideally there is only one peak between these values

limit1 = find_nearest(q, q_1) #First our lower limit
limit2 = find_nearest(q, q_2) #And of our higher limit
print('limit1=', limit1)
print('limit2=', limit2)


# In[330]:


# Having extablished our limits we can now trim the data

q_sub = q[limit1:limit2] # We'll reduce the domain of Q
perov_sub = perov[limit1:limit2, 1:-1] # And correspondingly shrink the size of our diffraction data set. 
#YOU NEED THE () ABOVE FOR (1,-1) OR ELSE THE CODE WILL NOT FUNCTION

# You'll also notice I'm dropping the first column of data in the perov array, since that's now accounted for in q_sub
#Now let's go ahead and plot again to see if our trimming worked
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q_sub,perov_sub[:,1],'g-', label='$MAPbIBr_2$') #plot subfield of data
plt.xlabel('Q [$\AA^{-1}$]',size=12) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=12)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper right hand corner
plt.title('Initial')

q_linear = np.hstack((q_sub[0:10], q_sub[-11:-1])) #I'm taking the starting and ending values
perov_linear = np.hstack((perov_sub[0:10,0], perov_sub[-11:-1,0])) #We'll use these to fit a straight line
slope, intercept = np.polyfit(q_linear, perov_linear, 1) #Do linear fit
back = slope*q_sub+intercept #Create background array of the form Background = Ax+B
plt.plot(q_sub, back)


# 

# In[331]:


#Let's begin by getting our data ready to analyze
perov_fit = perov_sub[:,0]-back #We'll begin by subtracting the background we calculated for this piece of data

#Now let's define a function we'll want to fit to - this is analagous to the "straight-line-model" from tutorial 03
#We'll call our function gaussian and it will calculate the expression described above
def gaussian(x, a, b, c): 
    return a*np.exp(-(x - b)**2/(2*c**2))

#We'll also give an initial guess for our fits based off of a visual interpretaion of our data
p0 = [300, 1.05, 0.1]

#Use scipy.optimize.curve_fit to fit our desired data
popt, pcov = scipy.optimize.curve_fit(gaussian, q_sub, perov_fit, p0)


# In[332]:


#To confirm our fits it's always nice to plot our model versus our data.
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q_sub,perov_fit,'g-', label='$MAPbIBr_2$ Initial') #plot subfield of data
plt.plot(q_sub,gaussian(q_sub, *popt),'b--', label='Model') #plot best fit
plt.xlabel('Q [$\AA^{-1}$]',size=12) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=12)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner

#NOTE: subtracting background moves plot to zero


# In[333]:


#With confidence in our fit we can now go ahead and print/make note of/tabulate our parameters of interest. 
#Print peak intensity
print('Final Intensity:', popt[0])

#Caculate and pring d-spacing
d = 2*np.pi/popt[1] #Applying d = 2*pi/Q
print('d-Spacing: ', d) 

#Print lattice constant
miller = [1, 0, 0] #need to guess miller indices of peak
a = d/np.sqrt(miller[0]**2+miller[1]**2+miller[2]**2) #calculate a using a = d/sqrt(h^2+k^2+l^2) for a cubic lattice
print('Lattice Spacing:', a)


# ### With Github Code

# In[334]:


def csv_to_np(filename): #changing excel to numpy
    data = pd.read_csv(filename)
    return(np.array(data))
perov = csv_to_np(r'C:\Users\tiann\OneDrive\Desktop\Summer Research 2020\D1_MaPbBlBr2_XrayDegGSASExcel.csv')


# In[335]:


plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0],perov[:,1],'r-', label='$MAPbIBr_2$ initial') #plot two-theta versus XRD intesntiy
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)
plt.title('initial')
plt.legend(loc="upper right")


# In[336]:


#Note that our initial 2-theta values are in degrees, not radians
#We also need to define the X-ray wavelength we're using
wave = 0.9763 #wavelength used at SSRL in Angstroms
theta = perov[:,0]/2*np.pi/180 #convert from 2-theta to theta in radians
q = 4*np.pi*np.sin(theta)/wave #convert to q


# In[337]:


#We can now check that our conversion went as expected by ploting our data versus q
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q,perov[:,1],'r-', label='$MAPbIBr_2 Initial') #plot Q versus XRD intesntiy
plt.xlabel('Q [$\AA^{-1}$]',size=12)
plt.ylabel('Intensity [a.u.]',size=12)
plt.legend(loc="upper right")
plt.title('Initial')


# In[338]:


def find_nearest(array, target):
    array = np.asarray(array) # In case the input is not an array type
    idx = (np.abs(array - target)).argmin() # Finds the index of the value closest to the target
    return idx
#We use our new function to find the index of real value closest to our desired limit
q = 4*np.pi*np.sin(theta)/wave
q_1 = 0.98 #lower limit for Q we'll consider
q_2 = 1.4 # upper limit - ideally there is only one peak between these values

limit1 = find_nearest(q, q_1) #First our lower limit
limit2 = find_nearest(q, q_2) #And of our higher limit
print('limit1=', limit1)
print('limit2=', limit2)


# In[382]:


q_sub = q[limit1:limit2] # We'll reduce the domain of Q
perov_sub = perov[limit1:limit2, 1:-1] 
plt.plot(q_sub, perov_sub[:,-1])


# In[403]:


size = perov_sub.shape
q_bins = size[0]
num_frames = size[1]
print('number of frames=', num_frames)
slope = np.zeros((num_frames,1))
b = np.zeros((num_frames,1))
back = np.zeros((q_bins,num_frames))
int_correct = np.zeros((q_bins,num_frames))
print('size=',size)
print('q bins=', q_bins)


# In[384]:


#accept a data file and range and returns average values at start and end of range
for j in range(num_frames): 
    slope[j] = ((np.mean(perov_sub[-10:-1,j])-np.mean(perov_sub[0:10,j]))/(np.mean(q[limit2-10:limit2])-np.mean(q[limit1:limit1+10])))
    b[j]=perov_sub[0,j]-slope[j]*q_sub[0]
    back[:,j] = [slope[j]*element+b[j] for element in q_sub]
    int_correct[:,j] = perov_sub[:,j]-np.array(back[:,j]) #y-axis
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(np.array(q_sub),int_correct)


# In[385]:


def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2)) 
#trying to find intensity


# In[386]:


import scipy
from scipy.optimize import curve_fit
import math
from math import pi,sin
p0 = [300, 1.05, .01] 
intensity_1 = np.zeros((num_frames)) #create correct size arrays for running in the loop
lattice_1= np.zeros((num_frames))
#print(lattice_1)
for j in range(num_frames):
    popt, pcov = curve_fit(gaussian, np.array(q_sub), int_correct[:,j], p0)
    intensity_1[j] = popt[0]
    lattice_1[j] = 2*math.pi/popt[1] #need to fit for correct index
    p0 = popt #once you have the initial fit use that as your next guess, we expect things to be continuous so this helps with that
print(lattice_1)
print('time=', time)
print('intensity=', intensity_1)
plt.plot(time, intensity_1)
plt.xlabel('time(s)')
plt.ylabel('Intensity [a.u]')
plt.title('Effect of Xray Exporsure on Intensity')


# In[381]:


time = np.zeros(num_frames) #create empty array for time of correct length
for x in range(num_frames):
    time[x] = x*20 + 20 #xray dose is 20s per frame
print(time)
print(num_frames)


# In[401]:


print(lattice_1)
print(time)

plt.plot(time, lattice_1, 'go')
plt.xlabel('time(s)')
plt.ylabel('Lattice Spacing')
plt.title('Effect of Xray Exposure on Lattice Spacing')


# In[397]:


print(lattice_1)
print(time)
plt.figure(figsize=(20,12))
plt.ylim(0,10)
plt.plot(time, lattice_1, 'go')
plt.xlabel('time(s)', fontsize=18)
plt.ylabel('Lattice Spacing [$\AA$]', fontsize = 18)
plt.title('Effect of Xray Exposure on Lattice Spacing', fontsize=18)


# In[370]:


# Exercise 4
import math

planes = [[1,0,0],[1,1,0], [1,1,1], [2,0,0], [3,0,0]]
#two_d = np.array([planes])
#print(two_d)
#print(two_d[:,1]**2)
for plane in planes: #hones in on one miller index
    a = 6.03 #FEEDBACK Better would be to link this directly to previous fit for maximum accuracy
    sumplane = 0
    for n in plane: #hones in on each value within the miller index #FEEDBACK NESTED FOR LOOPS! NICE!!
        millersqr = n**2
        sumplane += millersqr
    d = a/((sumplane)**.5)
    print('d=', d)
    print('miller index sum squared=', sumplane) #FEEDBACK Might want to print specific miller index for readability? personally I find it easier to interpret [110] than 2.
    lambbda = .9763
    Q = (2*np.pi)/d
    print('Q=', Q)
    theta1 = (math.asin((Q*lambbda)/(4*np.pi)))
    print('theta1', theta1)
    answer = theta1*2
    print('radian=', answer)
    print('angle=', answer*180/np.pi)
    print()

#FEEDBACK This is really good! Slight improvement would be to generate an array of q values to that they're easier to store and reference later for plotting/comparisson/etc. How could you build that array as you step through your outer loop? How could you then plot those values on top of the actual XRD data/compare them tp the XRD data. 


# In[371]:


plt.figure(figsize=(8,6)) #make plot larger
plt.plot(perov[:,0], perov[:, 30], 'g-', label = '$MAPbIBr_2$ final')
plt.xlabel('2-theta [$^o$]',size=12) 
plt.ylabel('Intensity [a.u.]',size=12)
plt.title('Final')
plt.legend(loc="upper right")


# In[185]:


#Note that our initial 2-theta values are in degrees, not radians
#We also need to define the X-ray wavelength we're using
wave = 0.9763 #wavelength used at SSRL in Angstroms
theta = perov[:,0]/2*np.pi/180 #convert from 2-theta to theta in radians
q = 4*np.pi*np.sin(theta)/wave #convert to q


# In[186]:


#We can now check that our conversion went as expected by ploting our data versus q
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q,perov[:,30],'g-', label='$MAPbIBr_2 Final') #plot Q versus XRD intesntiy
plt.xlabel('Q [$\AA^{-1}$]',size=12)
plt.ylabel('Intensity [a.u.]',size=12)
plt.legend(loc="upper right")
plt.title('Final')


# In[31]:


def find_nearest(array, target):
    array = np.asarray(array) # In case the input is not an array type
    idx = (np.abs(array - target)).argmin() # Finds the index of the value closest to the target
    return idx
#We use our new function to find the index of real value closest to our desired limit
q = 4*np.pi*np.sin(theta)/wave
q_1 = 0.98 #lower limit for Q we'll consider
q_2 = 1.2 # upper limit - ideally there is only one peak between these values

limit1 = find_nearest(q, q_1) #First our lower limit
limit2 = find_nearest(q, q_2) #And of our higher limit
print('limit1=', limit1)
print('limit2=', limit2)


# In[135]:


# Having extablished our limits we can now trim the data

q_sub = q[limit1:limit2] # We'll reduce the domain of Q
perov_sub = perov[limit1:limit2, (1,-1)] # And correspondingly shrink the size of our diffraction data set. 
#YOU NEED THE () ABOVE FOR (1,-1) OR ELSE THE CODE WILL NOT FUNCTION

# You'll also notice I'm dropping the first column of data in the perov array, since that's now accounted for in q_sub
#Now let's go ahead and plot again to see if our trimming worked
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q_sub,perov_sub[:,1],'g-', label='$MAPbIBr_2$') #plot subfield of data
plt.xlabel('Q [$\AA^{-1}$]',size=12) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=12)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper right hand corner
plt.title('Final')

q_linear = np.hstack((q_sub[0:10], q_sub[-11:-1])) #I'm taking the starting and ending values
perov_linear = np.hstack((perov_sub[0:10,0], perov_sub[-11:-1,0])) #We'll use these to fit a straight line
slope, intercept = np.polyfit(q_linear, perov_linear, 1) #Do linear fit
back = slope*q_sub+intercept #Create background array of the form Background = Ax+B


# In[155]:


#Let's begin by getting our data ready to analyze
perov_fit = perov_sub[:,1]-back #We'll begin by subtracting the background we calculated for this piece of data

#Now let's define a function we'll want to fit to - this is analagous to the "straight-line-model" from tutorial 03
#We'll call our function gaussian and it will calculate the expression described above
def gaussian(x, a, b, c): 
    return a*np.exp(-(x - b)**2/(2*c**2))

#We'll also give an initial guess for our fits based off of a visual interpretaion of our data
p0 = [220, 1.05, 0.1]

#Use scipy.optimize.curve_fit to fit our desired data
popt, pcov = scipy.optimize.curve_fit(gaussian, q_sub, perov_fit, p0)


# In[156]:


#To confirm our fits it's always nice to plot our model versus our data.
plt.figure(figsize=(8,6)) #make plot larger
plt.plot(q_sub,perov_fit,'g-', label='$MAPbIBr_2$ Initial') #plot subfield of data
plt.plot(q_sub,gaussian(q_sub, *popt),'b--', label='Model') #plot best fit
plt.xlabel('Q [$\AA^{-1}$]',size=12) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=12)#Define y-axis label
plt.legend(loc="upper right")#Put legend in upper left hand corner


# In[157]:


#Note that our initial 2-theta values are in degrees, not radians
#We also need to define the X-ray wavelength we're using
wave = 0.9763 #wavelength used at SSRL in Angstroms
theta = perov[:,0]/2*np.pi/180 #convert from 2-theta to theta in radians
q = 4*np.pi*np.sin(theta)/wave #convert to q


# In[158]:


#With confidence in our fit we can now go ahead and print/make note of/tabulate our parameters of interest. 
#Print peak intensity
print('Final Intensity:', popt[0])

#Caculate and pring d-spacing
d = 2*np.pi/popt[1] #Applying d = 2*pi/Q
print('d-Spacing: ', d) 

#Print lattice constant
miller = [1, 0, 0] #need to guess miller indices of peak
a = d/np.sqrt(miller[0]**2+miller[1]**2+miller[2]**2) #calculate a using a = d/sqrt(h^2+k^2+l^2) for a cubic lattice
print('Lattice Spacing:', a)


# In[ ]:





# In[ ]:





# In[ ]:




