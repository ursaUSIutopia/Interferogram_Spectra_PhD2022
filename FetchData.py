# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:03:24 2022

@author: Usasi.CHOWDHURY
"""

def accesspath(num,filename):
    """
    This function changes the filepath according to the RUN number and filename

    Parameters
    ----------
    num : string
         Run number
    filename : string
        FITS file name without the extension

    Returns
    -------
    None.

    """
    import os
    os.chdir(r"C:/Users/Usasi.CHOWDHURY/Desktop/Python_Program/Acquisition/NIKA0/Run_%s/"%num)
    return

def astropy_ndobj(filename):
    """
    
    This function changes FITS to NumPy object array.

    Parameters
    ----------
    filename : string
        Name of FITS file without extension

    Returns
    -------
    data : NumPy object array
        FITS data in accessible format

    """
    
    from astropy.io import fits
    hdul = fits.open("%s.fits"%filename) #fits data open
    hdul.info()
    data = hdul[1].data 
    return data

def createarray(data,m):
    import numpy as np
    xnew_arr=np.arange((np.size(data)*m), dtype=np.float64).reshape(np.size(data),m)
    ynew_arr=np.arange((np.size(data)*m), dtype=np.float64).reshape(np.size(data),m)
    for i in range(m):
        xnew_arr[:,i]= data['X_0%.2d'%i] / 1000 #[Make a array for x data] 
        ynew_arr[:,i]= data['Y_0%.2d'%i] #[Make a array for y data]
    return xnew_arr, ynew_arr