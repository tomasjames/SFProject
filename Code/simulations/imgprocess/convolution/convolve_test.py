################################################################################
############################ 4th Year Project Code #############################
########################### Image Convolution Tests ############################
################################################################################

############################# Import statements ################################

# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Import units
from radmc3dPy.natconst import *

# Import numpy
import numpy as np

# Import matplotlib
from matplotlib.pyplot import *

# Import astropy
from astropy.io import fits

############################## Define functions ################################

# Function taken from http://scipy-cookbook.readthedocs.org/items/Rebinning.html
def rebin( a, newshape ):
        '''Rebin an array to a new shape.
        '''
        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = np.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        return a[tuple(indices)]

############################# Read in FITS data ################################

# Read in the PSF
psf_raw = fits.open('./theoretical_spire_beam_model_psw_V0_2.fits')

# Extract the data
psf = psf_raw[0].data
psf_header = psf_raw[0].header

# Read in the PSF
data_raw = fits.open('./HorseHead.fits')

# Extract the data
data = data_raw[0].data
data_header = data_raw[0].header

# Pull out data dimensions
data_x = data_header['NAXIS1'] # Number of pixels in X axis
data_y = data_header['NAXIS2'] # Number of pixels in Y axis

psf_x = psf_header['NAXIS1'] # Number of pixels in X axis
psf_y = psf_header['NAXIS2'] # Number of pixels in Y axis

if data_x != psf_x and data_y != psf_y:
    print 'Dimensions of the data image do not match the dimensions of the PSF image. The data image will now be rebinned to reduce its dimensions to that of the PSF image.\n'

# Plot the data for comparison purposes
imshow(psf,origin='lower')
savefig('psf.png',dpi=300)
close()

imshow(data,origin='lower')
savefig('HorseHead.png',dpi=300)
close()

############################### Rebin the data #################################

# Use the above rebin function to change the dimensions of the array
data_rebin = rebin(data,[psf_x,psf_y])

# Plot the resulting data to assess accuracy of rebin
imshow(data_rebin,origin='lower')
savefig('HorseHead_rebin.png',dpi=300)
close()

# Determine dimensions of the new image
data_rebin_x = len(data_rebin)
data_rebin_y = len(data_rebin[0])

print 'The image has now been rebinned to have dimensions', data_rebin_x, 'x', data_rebin_y, '. The new image has also been saved to the working directory.'
