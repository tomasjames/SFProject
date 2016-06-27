################################################################################
########################### Image Convolution Tests ############################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import numpy
import numpy as np

# Import matplotlib
from matplotlib.pyplot import *

# Import astropy
from astropy.io import fits
from astropy.convolution import convolve

# Import operating system manipulation tools
import os

# Import radmc3dPy tools
from radmc3dPy.natconst import *

############################## Define functions ################################

# Function taken from http://scipy-cookbook.readthedocs.org/items/Rebinning.html
def rebin(a, newshape):
        '''
        Rebin an array a to a new shape newshape.
        '''
        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = np.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        return a[tuple(indices)]

################################# Read data ####################################

# Read in the FITS file containing the simulation results
sim_result = fits.open('../../data/psw/background_15K/psw.fits')

# Split into image data and header data
data = sim_result[0].data
hdu = sim_result[0].header

# Pull out the total image width
image_width = hdu['IMGWIDTH']

# Pull out the pixel width
data_theta = hdu['PIXSIZE']

# Determine the size of one pixel
data_size = (250*data_theta)/206256

############################# Read in FITS data ################################

# Read in the PSF
psf_raw = fits.open('psf/theoretical_spire_beam_model_psw_V0_2.fits')

# Extract the data
psf = psf_raw[0].data
psf_header = psf_raw[0].header

# Assign cards to variables
psf_theta = psf_header['CDELT1'] # Pixel width in degrees

# Determine the size of one pixel
psf_size_pixel = (150*psf_theta*3600)/206256 # 150 is the distance in parsecs to the source
psf_size = psf_size_pixel*len(psf[0]) # Total PSF spans psf_size in parsecs

# Determine the total number of 'beams' required
beams = (image_width/au)/psf_size # Divide the required width by the size of the PSF

# Plot the data for comparison purposes
subplot2grid((10,8), (0,0), colspan=8,rowspan=8)
title(r'Raw SPIRE 250$\mu m$ PSF')
imshow(psf,origin='lower')
colorbar()
tight_layout
subplot2grid((10,8), (8,0), colspan=8,rowspan=2)
plot(np.linspace(0,len(psf),len(psf)),psf[len(psf)/2])
xlim(min(np.linspace(0,len(psf),len(psf))), max(np.linspace(0,len(psf),len(psf))))
ylabel('Pixel Value')
xlabel('Pixel')
tight_layout()
savefig('psf.pdf',dpi=300)
close()

imshow(data,origin='lower')
title(r'RADMC-3D Intensity Data for SPIRE 250$\mu m$')
colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
savefig('data.pdf',dpi=300)
close()

############################### Rebin the data #################################

# First rebin the initial images to have the correct number of beams (save computation time)
data_rebin = rebin(a=data,newshape=[beams*11,beams*11])

# Determine whether the pixel widths of the PSF and data differ
if data_size == psf_size:
    print 'The pixel sizes match!\n'
    psf_rebin = psf

else:
    print 'The pixel sizes do not match! Rebinning the PSF so that the pixel sizes match.\n'

    # Use the above rebin function to change the dimensions of the array
    psf_rebin = rebin(a=psf,newshape=[11,11])

# Plot the resulting data to assess accuracy of rebin
subplot2grid((10,8), (0,0), colspan=8,rowspan=8)
title(r'Rebinned SPIRE 250$\mu m$ PSF')
imshow(psf_rebin,origin='lower')
colorbar()
tight_layout
subplot2grid((10,8), (8,0), colspan=8,rowspan=2)
plot(np.linspace(0,len(psf_rebin[len(psf_rebin)/2]),len(psf_rebin[len(psf_rebin)/2])),psf_rebin[len(psf_rebin)/2])
xlim(min(np.linspace(0,len(psf_rebin),len(psf_rebin))), max(np.linspace(0,len(psf_rebin),len(psf_rebin))))
ylabel('Pixel Value')
xlabel('Pixel')
tight_layout()
savefig('psf_rebin.pdf',dpi=300)
close()

# Determine dimensions of the new image
psf_rebin_x = len(psf_rebin)
psf_rebin_y = len(psf_rebin[0])

print 'The PSF image has now been rebinned to have dimensions', psf_rebin_x, 'x', psf_rebin_y, '. The new image has also been saved to the working directory.\n'

############################### Convolve the data ##############################

print 'Currently convolving...'

# Use the convolution to convolve the rebinned data with the PSF
conv = convolve(data_rebin, psf_rebin, boundary='extend')

print 'Convolution complete!\nSaving the image...\n'

# Plot the resulting convolution
imshow(conv,origin='lower')
title(r'Convolved (and rebinned) Intensity Data for SPIRE 250$\mu m$')
colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
savefig('convolved.pdf',dpi=300)
close()

print 'The convolution has been performed and the convolved image saved to the working directory.\n'

############################### Check image power ##############################

# Begin with default, non altered image
raw_power = np.sum(np.sum(data_rebin))

# Repeat for convolved image
convolved_power = np.sum(np.sum(conv))

# Print a statement to the terminal to show overall power
print 'The initial power in the image is ', raw_power
print 'The power in the convolved image is ', convolved_power
