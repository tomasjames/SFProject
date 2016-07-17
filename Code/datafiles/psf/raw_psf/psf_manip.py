################################################################################
############################ 4th Year Project Code #############################
##############################  PSF Manipulation  ##############################
################################################################################

############################# Import statements ################################
# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Import time tools
import time

# Import units
from radmc3dPy.natconst import *

# Import numpy
import numpy as np

# Import csv
import csv

# Import matplotlib
from matplotlib.pyplot import *

# Import astropy
from astropy.io import fits

# Import scipy
import scipy.interpolate
import scipy.ndimage

# Import glob for file scraping
import glob

import random

from convolve.convolution import *

################################ Read in PSFs ################################

# Initialise a list to store the read in data
data, filenames = [], []

# Initiate counter to help store the passbands
count = -1

# Loop through all FITS files
for file in glob.glob('*.fits'):
    # Open them
    with open(file) as f:
        count += 1

        # Read in using loadtxt and then append the data to a list for use later
        contents = fits.open(file)

        # Append the name of the file to the position (count+len(data)) followed by
        # the file contents
        data.insert(count+len(data), file)
        data.append(contents)
        filenames.append(file)

# Assign to variables and put fluxes into an array (data plotted later)
# +1 term avoids the filename that is placed in the array for identification
blue = data[data.index(filenames[0])+1]
green = data[data.index(filenames[1])+1]
red = data[data.index(filenames[2])+1]
plw = data[data.index(filenames[3])+1]
pmw = data[data.index(filenames[4])+1]
psw = data[data.index(filenames[5])+1]

########################## Determine PSF sizes and resize ######################

blue_size, blue_pix = blue[0].header['NAXIS1'], blue[0].header['CDELT1']*(-1)
green_size, green_pix = green[0].header['NAXIS1'], green[0].header['CDELT1']*(-1)
red_size, red_pix = red[0].header['NAXIS1'], red[0].header['CDELT1']*(-1)
plw_size, plw_pix = plw[0].header['NAXIS1'], plw[0].header['CDELT1']
pmw_size, pmw_pix = pmw[0].header['NAXIS1'], pmw[0].header['CDELT1']
psw_size, psw_pix = psw[0].header['NAXIS1'], psw[0].header['CDELT1']

# Resize PSFs to be the same size as the smallest PSF
green = congrid(green[0].data, (blue_size,blue_size))
red = congrid(red[0].data, (blue_size,blue_size))
plw = congrid(plw[0].data, (blue_size,blue_size))
pmw = congrid(pmw[0].data, (blue_size,blue_size))
psw = congrid(psw[0].data, (blue_size,blue_size))

# Determine new data quantities
new_blue_pix = blue_pix
new_green_pix = green_pix*(blue_size/np.float64(green_size))
new_red_pix = red_pix*(blue_size/np.float64(red_size))
new_plw_pix = plw_pix*(blue_size/np.float64(plw_size))
new_pmw_pix = pmw_pix*(blue_size/np.float64(pmw_size))
new_psw_pix = psw_pix*(blue_size/np.float64(psw_size))

# Collect resized data into list
data_store = [[blue[0].data], [green], [red], [plw], [pmw], [psw]]
size_store = [new_blue_pix, new_green_pix, new_red_pix, new_plw_pix, new_pmw_pix, new_psw_pix]

################################## Save to FITS file ############################

# Change the folder to the parent directory
os.chdir('../')

# Loop through all of the filenames
for i in range(0,len(filenames)):

    # Do the write
    fname = str(filenames[i])
    head = fits.PrimaryHDU(data_store[i])
    head.writeto(fname)

    # Remove erroneous NAXIS3 keyword card if it exists
    statement = 'NAXIS3' in head.header
    if statement == True:
        del head.header['NAXIS3']

    # Write some header information
    # Write relevant information to the FITS header
    fits.setval(fname, 'CTYPE1', value='RA---TAN')
    fits.setval(fname, 'CTYPE2', value='DEC--TAN')
    fits.setval(fname, 'CDELT1', value=(-1)*size_store[i])
    fits.setval(fname, 'CDELT2', value=size_store[i])
    fits.setval(fname, 'CRPIX1', value=np.round(len(data_store[i][0])/2))
    fits.setval(fname, 'CRPIX2', value=np.round(len(data_store[i][0])/2))
    fits.setval(fname, 'CRVAL1', value=0.0)
    fits.setval(fname, 'CRVAL2', value=0.0)
