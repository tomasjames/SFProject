################################################################################
############################ 4th Year Project Code #############################
##############################  Functions to run  ##############################
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


######################### Function to determine number of beams ###########################
def beams_pp(filename):

    '''
    Determines the number of beams per pixel for any given filter.

    Simulations would normally expect 1 beam per pixel, so the number of beams
    aligns with the number of pixels. Increasing the number of beams per pixel provides
    a higher resoltion.

    Keywords
    d: the distance to the source (in pc)

    ins: the instrument under consideration. Should be declared as a string
        Options:
            SPIRE
            PACS

    filt: the filter that accompanies the instrument declared in ins keyword
        Options - SPIRE:
            psw (250 micron)
            pmw (350 micron)
            plw (500 micron)

        Options - PACS:
            70 (70 micron)
            100 (100 micron)
            160 (160 micron)

    extent: the image width (in AU)
    '''

    '''
    # Use a basic if statement to determine the correct FWHM for the filter
    if ins == 'SPIRE':
        if filt == 'psw':
            fwhm = 18.
            wav = 250.
        elif filt == 'pmw':
            fwhm == 25.
            wav = 350.
        elif filt == 'plw':
            fwhm == 37.
            wav = 500.

    elif ins == 'PACS':
        if filt == '70':
            fwhm = 5.2
            wav = 70.
        elif filt == '100':
            fwhm = 7.7
            wav = 100.
        elif filt == '160':
            fwhm = 12.
            wav = 160.
    '''
    # Determine the current working directory
    direc = os.getcwd()

    # Read in the data
    data = fits.open(str(direc)+str('/')+str(filename)+str('.fits'))

    # Extract header
    hdu = data[0].header

    # Assign pixel size and distance to variables
    theta = hdu['PIXSIZE']
    d = hdu['DISTANCE']

    # Repeat for image width
    image_width = hdu['IMGWIDTH']

    # And for pixel width
    pix_width = hdu['PIXWIDTH']

    # Determine the physical resolution of a pixel
    res_pixel = ((d/pc)*theta)/206256

    # Determine the number of beams
    b = ((pix_width)/pc)/res_pixel

    return b

def beams_req(d, img_width, npix):

    '''
    Determines the number of beams for any given filter.

    Simulations would normally expect 1 beam per pixel, so the number of beams
    aligns with the number of pixels. Increasing the number of beams per pixel provides
    a higher resolution.

    Keywords
    d: the distance to the source (in pc)

    filt: the filter that accompanies the instrument declared in ins keyword
        Options - SPIRE:
            psw (250 micron)
            pmw (350 micron)
            plw (500 micron)

        Options - PACS:
            blue (70 micron)
            green (100 micron)
            red (160 micron)

    extent: the image width (in AU)
    '''



    res_pixel = (pix_width/pc)/b

    b = (pix_width/pc)/(((d/pc)*theta)/206256)
