# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:01:59 2022

@author: Usasi.CHOWDHURY
"""


import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import math
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import curve_fit
from FetchData import accesspath, astropy_ndobj, createarray

def peaks(Y2,m):
    #for i in range(m):
    i=13
    peaks, _ = find_peaks(-Y2[:,i],height=None, prominence=200, distance=10000)
    widths = peak_widths(-Y2[:,i], peaks, rel_height=200)
    return peaks,_,widths
def plotmulti(xarr,yarr,m):
    i=13#for i in range(m):
    cm = plt.get_cmap('gist_rainbow')
    plt.plot(xarr[:,i], yarr[:,i],  label='KID %.d'%(i+1), color = cm(i*20),linewidth=2.0)
    plt.xlabel("Frequency(GHz)",fontsize=30)
    plt.ylabel("Frequenccy Shift(Hz)",fontsize=30)
    #plt.xlim(xmin,xmax)
    #plt.ylim(-60000)
    plt.tick_params(labelsize=25, length=6, width=2, grid_alpha=0.5)
    plt.legend(loc='upper left',fontsize='large')
def Gauss(x, A, mu, sigma,c):
    

    """
    The response function (Thesis: Flanigan_2018, Pg.-51)
    This function calculates the S21 regarding the input values

    Parameters
    ----------
    x : string
        input data(in GHz)
    A : amplitude
    mu : center
    sigma : sigma #FWHM = 2*sigma*sqrt(2*ln2) = 2.3548*sigma
    Returns
    -------
    y or S21 linear

    """
    y = c + (A/(sigma*math.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2*sigma**2))
    # y = (A/(np.pi))*(sigma/((x-mu)**2+(sigma**2)))
    return y
num = 234 #input("Run No: ")
filename = 'tt_2022_04_07_12h29m52' #input("Enter the FITS filename: ")
accesspath(num,filename)
data = astropy_ndobj(filename)
m = 15 #int(input("No of resonances:"))
xmax = 110 #int(input("The maximum value of in the x axis: "))
xmin = 80 #int(input("The minimum value of in the x axis: "))

xnew_arr,ynew_arr = createarray(data,m)
peaks,_,widths = peaks(ynew_arr,m)
sigma_arr = []
for i in range (m):
    x_kid1 = xnew_arr[:,i]
    y_kid1 = ynew_arr[:,i]
    plt.plot(x_kid1,y_kid1)
    min_in= np.argmin(y_kid1) #minimum index
    plt.plot(x_kid1[min_in],y_kid1[min_in],'o')
    # chosenPoints = np.where((x_kid1>96) & (x_kid1<97.5))
    
    x_1 = x_kid1[min_in-600:min_in+600]
    y_1 = y_kid1[min_in-600:min_in+600]
    # width = np.mean(widths[0])/100
    # xdata = xnew_arr[(peaks[2]-width):(peaks[2]+width),13]
    # ydata = ynew_arr[peaks[2]-width:peaks[2]+width]
    parameters, covariance = curve_fit(Gauss,x_1,y_1, p0 = [-14000,x_kid1[min_in],0.5,-1000])
    A = parameters[0]
    mu = parameters[1]
    sigma = parameters[2]
    c = parameters[3]
    fit_y = Gauss(x_1,A,mu,sigma,c)
    sigma_arr = np.insert(sigma_arr,i,sigma, axis=0)
    FWHM = sigma_arr * 2.3548 *1e3
    
    # for i in range(len(peaks)):
    #     plt.plot(xnew_arr[peaks[i],13],ynew_arr[peaks[i],13],'o')
    plt.plot(x_1,y_1)
    # plt.plot(x_1, ydata, '-', label='data')
   
    plt.plot(x_1, fit_y, '-', label='fit_KID_%s'%i)
    plt.legend()
KID= np.array(['NAN','NAN','NAN',105.360,86.33,103.520,87.680,89.100,101.920,98.3,100.02,93.590,95.090,96.750,90.920,92.0])
plt.figure()   
plt.plot(KID,'ro') 
plt.plot(FWHM,'ro')           
plt.xlabel("Frequency(GHz)",fontsize=30)
plt.ylabel("S21 Response",fontsize=30)
plt.legend(loc='upper left',fontsize='large')
  