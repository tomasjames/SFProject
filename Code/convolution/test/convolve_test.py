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

# Import radmc3dPy
from radmc3dPy.natconst import *

# Import rebin
from rebin import rebin, congrid

########################## Produce basic test data #############################

# Produce a dummy array to house the image
data_raw = np.zeros([500,500])

# Instantiate rudimentary counter
count_x, count_y = -1, -1

# Loop through the image and produce a smaller square of ones within the image of zeros
# Loop through x pixels
for i in range(0,len(data_raw[0])):
    # Count the x counter by 1 and reset the y counter
    count_x += 1
    count_y = -1

    # Loop through y pixels
    for j in range(0,len(data_raw[0])):
        # Increase counter by 1
        count_y += 1

        # Determine whether pixel is within pre defined bounds
        if i >= (len(data_raw)/2 - 10) and i <= (len(data_raw)/2 + 10) and j >= (len(data_raw)/2 - 10) and j <= (len(data_raw)/2 + 10):
            data_raw[count_x,count_y] = 1

# Check to see whether file has already been written
# If so, delete the file (Python cannot overwrite a binary file) and rewrite
if os.path.isfile('test.fits') == True:
    os.remove('test.fits')

    # Write the test data to a FITS file for later convolution
    hdu = fits.PrimaryHDU(data_raw)
    hdu.writeto('test.fits')

else:
    # Write the test data to a FITS file for later convolution
    hdu = fits.PrimaryHDU(data_raw)
    hdu.writeto('test.fits')

# Distance to the source
d = 300*pc

# The width of the source
imag_width = 15000*au

# The width of one pixel
pix_width = (imag_width)/(len(data_raw[0]))

# The angular size of one pixel
theta_rad = pix_width/d # In radians
#data_theta = theta_rad*(360/2*np.pi)*3600 # Convert to degrees and then to arcseconds
data_theta = theta_rad*3600 # Convert to arcseconds

# The effective wavelength that the data is plotted at
eff = '250'

# Determine the size of one pixel
data_size = (250*data_theta)/206256

fits.setval('test.fits', 'DISTANCE', value=d)
fits.setval('test.fits', 'IMGWIDTH', value=imag_width)
fits.setval('test.fits', 'PIXWIDTH', value=pix_width)
fits.setval('test.fits', 'PIXSIZE', value=data_theta)
fits.setval('test.fits', 'EFF', value=eff)

############################# Read in FITS data ################################

# Read in the PSF
psf_raw = fits.open('../psf/psf_0250.fits')
test_raw = fits.open('test.fits')

# Extract the data
psf = psf_raw[0].data
psf_header = psf_raw[0].header

data_raw = test_raw[0].data
hdu = test_raw[0].header

# Assign cards to variables
psf_theta = psf_header['CDELT2']*3600 # Pixel width in arcseconds
image_width = hdu['IMGWIDTH'] # Image width in cm

# Determine the size of one pixel
psf_size_pixel = ((d/pc)*psf_theta)/206256 # d is the distance in parsecs to the source and 206256 is the number of arseconds in a radian
#psf_size = psf_size_pixel*len(psf[0]) # Total PSF spans psf_size in parsecs
psf_size = ((d/pc)*16)/(206256)

# Determine the total number of 'beams' required
beams = np.round((image_width/pc)/psf_size) # Divide required width by the size of the PSF

# Plot the data for comparison purposes
subplot2grid((10,8), (0,0), colspan=8,rowspan=8)
title(r'Raw SPIRE 250$\mu m$ PSF')
#imshow(psf,origin='lower',vmin=0,vmax=1)
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

imshow(data_raw,origin='lower')
colorbar()
title(str('Raw Data for SPIRE ')+str(hdu['EFF'])+str(r'$ \mu m$'))
savefig('data.pdf',dpi=300)
close()

############################### Rebin the data #################################

# Define the new data dimension
new_dim = 101

# Assign data_raw to data for consistency
data = data_raw

# First rebin the initial images to have the correct number of beams (save computation time)
data_rebin = congrid(data,(beams*new_dim,beams*new_dim))

# Determine whether the pixel widths of the PSF and data differ
if data_size == psf_size:
    print 'The pixel sizes match!\n'
    psf_rebin = psf

else:
    print 'The pixel sizes do not match! Rebinning the PSF so that the pixel sizes match.\n'

    # Use the above rebin function to change the dimensions of the array
    psf_rebin = congrid(psf,(new_dim,new_dim))

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
conv = convolve(data_raw, psf_rebin, boundary='extend')

print 'Convolution complete!\nSaving the image...\n'

# Plot the resulting convolution
#imshow(conv,origin='lower',vmin=0,vmax=1)
imshow(conv,origin='lower')
colorbar()
title(str('Convolved Data for SPIRE ')+str(hdu['EFF'])+str(r'$ \mu m$'))
savefig('convolved.pdf',dpi=300)
close()

print 'The convolution has been performed and the convolved image saved to the working directory.\n'

############################### Check image power ##############################

# Begin with default, non altered image
raw_power = np.sum(np.sum(data_raw))

# Repeat for convolved image
convolved_power = np.sum(np.sum(conv))

# Print a statement to the terminal to show overall power
print 'The initial power in the image is ', raw_power
print 'The power in the convolved image is ', convolved_power
