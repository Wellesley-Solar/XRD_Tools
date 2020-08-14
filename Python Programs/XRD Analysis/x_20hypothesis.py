#this program is currently set up to run after converting to q, trimming data, 
# and isolating a given peak
# should take a guess composition and see if that is fitting with the xrd structure
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize 
from scipy.optimize import curve_fit
import pandas as pd
from xrdfunctions import csv_to_np, two_to_q, trim_data, num_files, back_subtract, pvoigt, gaussian, three_gaussians, q_to_a
from chemistry import predict_q, predict_lattice, predict_chem

#%% Import Data
perov_import = csv_to_np('One_Sun_100_wedge.csv')
#%% Covert to Q
q = perov_import[:,1]
perov = perov_import[:,2:] #seperate XRD intensities from rest of data 
miller = [1,0,0]
#%% hypothesis testing
def limit_fit(q,data,initial,hypothesis,index):
    intensity = 100
    fwhm = .005
    flex = .001
    p0 = [intensity, predict_q(hypothesis,index), fwhm,
    intensity, predict_q(initial,index), fwhm,
    intensity, predict_q(0.9,index), fwhm]
    upper_limit = [intensity*10, predict_q(hypothesis,index)+flex, fwhm*20,
    intensity*10, predict_q(initial,index)+flex, fwhm*20,
    intensity*10, predict_q(1,index),fwhm*20]
    lower_limit = [0, predict_q(hypothesis,index)-flex, 0,
    0, predict_q(initial,index)-flex, 0,
    0, predict_q(initial,index),0]
    popt,pcov = curve_fit(three_gaussians, q, data, p0, bounds=(lower_limit, upper_limit), maxfev=8000)
    return popt, np.sqrt(np.diag(pcov))
#%% hypothesis screening testing
fits=[]
error = []
chem = []
hyp = []
q_avg = []
for x in range(200):
    predict = limit_fit(q,perov[:,1],0.5,x*.001,miller)
    fits.append(predict[0])
    error.append(sum(predict[1]))
    temp = predict[0]
    chem.append((temp[0]*predict_chem(temp[1],miller)+temp[3]*predict_chem(temp[4],miller)+temp[6]*predict_chem(temp[7],miller))/(temp[0]+temp[3]+temp[6])) 
    q_avg.append((temp[0]*temp[1]+temp[3]*temp[4]+temp[6]*temp[7])/(temp[0]+temp[3]+temp[6]))
    hyp.append(x*.001*100)   

plt.plot(hyp,chem)
#%%seeing if a single hypothesis works for all compositions
totalerror = []
hyp = []
for percent in range(8):
    guess = percent*.01+.04
    comp = [.33, .5, .67, .75]
    fits = []
    error = []
    q_avg = []
    q_init = []
    i_phase = []
    br_phase = []
    for x in range(num_files(perov)):
        predict = limit_fit(q,perov[:,x],comp[x],guess,miller)
        fits.append(predict[0])
        error.append(sum(predict[1]))
        temp = predict[0]
        q_avg.append((temp[0]*temp[1]+temp[3]*temp[4]+temp[6]*temp[7])/(temp[0]+temp[3]+temp[6]))
        q_init.append(predict_q(comp[x],miller))
        i_phase.append(temp[1])
        br_phase.append(temp[7])
    totalerror.append(sum(error))
    hyp.append(guess)

#%%
best = 0.04 #hyp[totalerror.index(min(totalerror))]
fig,a =  plt.subplots(2,2)
#x = 0.33
predict = limit_fit(q,perov[:,0],comp[0],best,miller)
fits = predict[0]
a[0][0].plot(q,perov[:,0],'b-', label='Data') #plot subfield of data
a[0][0].plot(q,three_gaussians(q, *fits),'c--', label='Model') #plot best fit
a[0][0].plot(q, gaussian(q, fits[0],fits[1],fits[2]), 'b-.', label='Peak 1')
a[0][0].plot(q, gaussian(q, fits[3],fits[4],fits[5]), 'b--', label='Peak 2')
a[0][0].plot(q, gaussian(q, fits[6],fits[7],fits[8]), 'b:', label='Peak 3')
#x=0.5
predict = limit_fit(q,perov[:,1],comp[1],best,miller)
fits = predict[0]
a[0][1].plot(q,perov[:,1],'b-', label='Data') #plot subfield of data
a[0][1].plot(q,three_gaussians(q, *fits),'c--', label='Model') #plot best fit
a[0][1].plot(q, gaussian(q, fits[0],fits[1],fits[2]), 'b-.', label='Peak 1')
a[0][1].plot(q, gaussian(q, fits[3],fits[4],fits[5]), 'b--', label='Peak 2')
a[0][1].plot(q, gaussian(q, fits[6],fits[7],fits[8]), 'b:', label='Peak 3')
#x=0.67
predict = limit_fit(q,perov[:,2],comp[2],best,miller)
fits = predict[0]
a[1][0].plot(q,perov[:,2],'b-', label='Data') #plot subfield of data
a[1][0].plot(q,three_gaussians(q, *fits),'c--', label='Model') #plot best fit
a[1][0].plot(q, gaussian(q, fits[0],fits[1],fits[2]), 'b-.', label='Peak 1')
a[1][0].plot(q, gaussian(q, fits[3],fits[4],fits[5]), 'b--', label='Peak 2')
a[1][0].plot(q, gaussian(q, fits[6],fits[7],fits[8]), 'b:', label='Peak 3')
#x=0.75
predict = limit_fit(q,perov[:,3],comp[3],best,miller)
fits = predict[0]
a[1][1].plot(q,perov[:,3],'b-', label='Data') #plot subfield of data
a[1][1].plot(q,three_gaussians(q, *fits),'c--', label='Model') #plot best fit
a[1][1].plot(q, gaussian(q, fits[0],fits[1],fits[2]), 'b-.', label='Peak 1')
a[1][1].plot(q, gaussian(q, fits[3],fits[4],fits[5]), 'b--', label='Peak 2')
a[1][1].plot(q, gaussian(q, fits[6],fits[7],fits[8]), 'b:', label='Peak 3')
plt.show()    

#%% individual guess
predict = limit_fit(q,perov[:,0],.33,0.04,miller)
fits = predict[0]
plt.plot(q,perov[:,0],'b-', label='Data') #plot subfield of data
plt.plot(q,three_gaussians(q, *fits),'c--', label='Model') #plot best fit
plt.plot(q, gaussian(q, fits[0],fits[1],fits[2]), 'b-.', label='Peak 1')
plt.plot(q, gaussian(q, fits[3],fits[4],fits[5]), 'b--', label='Peak 2')
plt.plot(q, gaussian(q, fits[6],fits[7],fits[8]), 'b:', label='Peak 3')
plt.xlabel('Q [$\AA^{-1}$]',size=14) #Define x-axis label
plt.ylabel('Intensity [a.u.]',size=14)#Define y-axis label
