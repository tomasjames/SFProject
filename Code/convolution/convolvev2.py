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
from astropy.convolution import convolve_fft

# Import operating system manipulation tools
import os

# Import radmc3dPy tools
from radmc3dPy.natconst import *

# Import rebin
from rebin import rebin, congrid

################################# Read data ####################################

# Read in the FITS file containing the simulation results
sim_result = fits.open('../data/psw/background_15K/psw.fits')

# Split into image data and header data
data = sim_result[0].data
hdu = sim_result[0].header

# Pull out the total image width
image_width = hdu['IMGWIDTH'] # In cm

# Pull out the pixel width
data_theta = hdu['PIXSIZE'] # Pixel size in "

# Determine the size of one pixel
data_size = (250*data_theta)/206256 # Extent of one pixel

############################# Read in FITS data ################################

# Read in the PSF
psf_raw = fits.open('psf/theoretical_spire_beam_model_psw_V0_2.fits')

# Extract the data
psf = psf_raw[0].data
psf_header = psf_raw[0].header

# Assign cards to variables
psf_theta = psf_header['CDELT1'] # Pixel width in "
image_width = hdu['IMGWIDTH'] # Image width in cm

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
title(str('RADMC-3D Intensity Data for SPIRE ')+str(hdu['EFF'])+str('$ \mu m$'))
colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
savefig('data.pdf',dpi=300)
close()

############################### Rebin the data #################################

# Take the data pixel size and determine how much larger/smaller it is than the PSF
size_factor = data_theta/psf_theta

# Convolve requires that the array be odd in dimension
# Check to see if the size_factor is an odd number
if (np.round(size_factor) % 2 == 0):
    print 'The dimensions of the PSF are even; altering to make the dimensions odd...\n'
    size_factor = np.round(size_factor)
else:
    print 'The dimensions of the PSF are odd; proceeding as normal...\n'

# Determine whether the pixel widths of the PSF and data differ
if data_theta == psf_theta:
    print 'The pixel sizes match!\n'
    psf_rebin = psf
else:
    print 'The pixel sizes do not match! Rebinning the PSF so that the pixel sizes match.\n'

    # Use the above rebin function to change the dimensions of the array
    psf_rebin = congrid(psf,(psf_header['NAXIS1']/size_factor,psf_header['NAXIS2']/size_factor))

    new_psf_theta = psf_theta*size_factor

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

# Print statements for debug
print 'The raw PSF "/pixel:', psf_theta, '\n'
print 'The raw data "/pixel:', data_theta, '\n'
print 'The rebinned PSF "/pixel:', new_psf_theta, '\n'

print 'The PSF image has now been rebinned to have dimensions', psf_rebin_x, 'x', psf_rebin_y, '. The new image has also been saved to the working directory.\n'

############################### Convolve the data ##############################

print 'Currently convolving...'

# Use the convolution to convolve the rebinned data with the PSF
#conv = convolve(data, psf_rebin, boundary='extend')
conv = convolve_fft(data, psf_rebin)

print 'Convolution complete!\nSaving the image...\n'

# Plot the resulting convolution
#imshow(conv,origin='lower',vmin=0,vmax=1)
imshow(conv,origin='lower')
colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
title(str('Convolved Data for SPIRE ')+str(hdu['EFF'])+str(r'$ \mu m$'))
savefig('convolved.pdf',dpi=300)
close()

print 'The convolution has been performed and the convolved image saved to the working directory.\n'

############################### Check image power ##############################

# Begin with default, non altered image
raw_power = np.sum(np.sum(data))

# Repeat for convolved image
convolved_power = np.sum(np.sum(conv))

# Print a statement to the terminal to show overall power
print 'The initial power in the image is ', raw_power
print 'The power in the convolved image is ', convolved_power
