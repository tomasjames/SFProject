################################################################################
######################### PSF Quantity Determination ###########################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import numpy
import numpy as np

# Import matplotlib
from matplotlib.pyplot import *

# Import operating system manipulation tools
import os
import glob

# Import astropy
from astropy.io import fits

# Import radmc3dPy
from radmc3dPy.natconst import *

############################## Declare function #################################

def beams(ins, filter, extent):

    '''
    Determines the number of beams for any given filter.

    Simulations would normally expect 1 beam per pixel, so the number of beams
    aligns with the number of pixels. Increasing the number of beams per pixel provides
    a higher resoltion.

    Keywords
    ins: the instrument under consideration. Should be declared as a string
        Options:
            SPIRE
            PACS

    filter: the filter that accompanies the instrument declared in ins keyword
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

    # Use a basic if statement to determine the correct FWHM for the filter
    if ins == 'SPIRE':
        if filter == 'psw':
            fwhm = 18.
            wav = 250.
        elif filter == 'pmw':
            fwhm == 25.
            wav = 350.
        elif filter == 'plw':
            fwhm == 37.
            wav = 500.

    elif ins == 'PACS':
        if filter == '70':
            fwhm = 5.2
            wav = 70.
        elif filter == '100':
            fwhm = 7.7
            wav = 100.
        elif filter == '160':
            fwhm = 12.
            wav = 160.

    # Determine the physical resolution of the image
    res = (wav*fwhm)/206256

    # Determine the number of beams
    b = ((extent*au)/pc)/res

    return b
