import numpy as np
import csv
import os
import decimal
import matplotlib.pyplot as plt
from matplotlib import pyplot as mp
import pandas as pd
import math
import scipy
from scipy.optimize import curve_fit
from math import pi,sin
from heapq import nsmallest


def fwhm(y_values_q, x_values):
    y_values, q_l, q_r,q,inten = [], [], [],[],[]

    # To make 'y_values_temp', a numpy array, into a python list
    for x in range(0,len(y_values_q)):
        y_values.append(y_values_q[x])
    peak_height = max(y_values)
    print(peak_height)
    half_peak_height = max(y_values)/2
    print(half_peak_height)
    # Splitting the y_values data into before and after x_value at peak height
    y_l_q = y_values[0:y_values.index(peak_height)]
    #print(y_l_temp)
    y_r_q = y_values[y_values.index(peak_height):len(y_values)]
    #print(y_r_temp)
    # Finds 1st closest value to half_peak_height in y_l and y_r
    y_l = nsmallest(1, y_l_q, key=lambda x: abs(x-half_peak_height))
    print (y_l)
    y_r = nsmallest(1, y_r_q, key=lambda x: abs(x-half_peak_height))
    print (y_r)
    inten.append(y_l)
    inten.append(y_r)
    #print(inten)
    # Gets x_value pairs for y_l and y_r
    print(y_values)
    print(x_values)
    print (x_values[y_values.index(y_l)])
    q_l.append(x_values[y_values.index(y_l)])#[0]
    q_r.append(x_values[y_values.index(y_r)])#[0]+ len(y_l) -1
    fwhm_n = ((abs(q_l[0] - q_r[0])))
    q.append(q_l)
    q.append(q_r)
    #print (temp_l[0])
    print(q_l)
    print(q_r)
    #print ( fwhm_n)
    plt.plot(q,inten)
    return(float(fwhm_n))
#%%
#Defining all functions that will be used in program
def readcsv(filename):
    data = pd.read_csv(filename)
    return(np.array(data))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
 





