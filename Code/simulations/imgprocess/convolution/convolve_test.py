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
from astropy.convolution import convolve

# Import 3d plotting tools
from mpl_toolkits.mplot3d import Axes3D

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

############################# Read in FITS data ################################

# Read in the PSF
psf_raw = fits.open('./theoretical_spire_beam_model_psw_V0_2.fits')
#psf_raw = fits.open('./Corrected_PSF_SPIRE_250.fits')

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
    print 'Dimensions of the data image do not match the dimensions of the PSF image. The PSF image will now be rebinned so that its dimensions match the data image.\n'

# Plot the data for comparison purposes
subplot2grid((12,12), (0,2), colspan=10,rowspan=10)
imshow(psf,origin='lower')
colorbar()
tight_layout
subplot2grid((12,12), (10,0), colspan=12,rowspan=2)
plot(np.linspace(0,len(psf),len(psf)),psf[len(psf)/2])
xlim(min(np.linspace(0,len(psf),len(psf))), max(np.linspace(0,len(psf),len(psf))))
ylabel('Pixel Value')
xlabel('Pixel')
tight_layout()
savefig('psf.png',dpi=300)
close()

imshow(data,origin='lower')
colorbar()
savefig('HorseHead.png',dpi=300)
close()

############################### Rebin the data #################################

# Use the above rebin function to change the dimensions of the array
#psf_rebin = rebin(a=psf,newshape=[data_x,data_y])
psf_rebin = rebin(a=psf,newshape=[11,11])

# Plot the resulting data to assess accuracy of rebin
subplot2grid((12,12), (0,2), colspan=10,rowspan=10)
imshow(psf_rebin,origin='lower')
colorbar()
tight_layout
subplot2grid((12,12), (10,0), colspan=12,rowspan=2)
plot(np.linspace(0,len(psf_rebin[len(psf_rebin)/2]),len(psf_rebin[len(psf_rebin)/2])),psf_rebin[len(psf_rebin)/2])
xlim(min(np.linspace(0,len(psf_rebin),len(psf_rebin))), max(np.linspace(0,len(psf_rebin),len(psf_rebin))))
ylabel('Pixel Value')
xlabel('Pixel')
tight_layout()
savefig('psf_rebin.png',dpi=300)
close()

# Determine dimensions of the new image
psf_rebin_x = len(psf_rebin)
psf_rebin_y = len(psf_rebin[0])

print 'The PSF image has now been rebinned to have dimensions', psf_rebin_x, 'x', psf_rebin_y, '. The new image has also been saved to the working directory.\n'

########################### Determine pixel widths #############################

# The FWHM for the PSW band
fwhm = 18. # FWHM in arcseconds

# Determine the radius of the source
d = 150. # Distance to the source in parsecs
r = (fwhm/3600)*d # Calculates the radius based upon the small angle approximation

print 'The radius of the psf, assuming it subtends an angle of', fwhm, '\", is', r, 'pc.'

# Define data storage lists
data_store_index = []

# Determine the range of data which constitutes the PSF
for i in range(0,len(psf_rebin)):
    # Checks to see if the datapoint in the middle row is greater than 10 times the mean value of the first 100 points
    if psf_rebin[len(psf_rebin)/2][i] > np.mean(psf[len(psf_rebin)/2][0:100]*10):
        data_store_index.append(i)

psf_data_store = psf_rebin[len(psf_rebin)/2][data_store_index]

# Determine the number of pixels from edge to centre (the centre being the median of the array)
psf_r = np.round(np.median(data_store_index)) - psf_data_store[0]

print 'The number of pixels in the radius of the PSF is', psf_r, 'pixels.'

# Divide this by the radius to get the distance of 1 pixel
pix = r/psf_r

print 'This corresponds to a pixel width of', pix, 'pc.\n'

'''
answer = input('Plot 3D histogram of the PSF?\n')

if answer == 'Yes':
    print 'Entering the 3D histogram routine to visualise the PSF...'

    # Create a figure for plotting the data as a 3D histogram.
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create an X-Y mesh of the same dimension as the 2D data. You can
    # think of this as the floor of the plot.

    x_data, y_data = np.meshgrid( np.arange(psf.shape[1]),
                                  np.arange(psf.shape[0]) )
    #
    # Flatten out the arrays so that they may be passed to "ax.bar3d".
    # Basically, ax.bar3d expects three one-dimensional arrays:
    # x_data, y_data, z_data. The following call boils down to picking
    # one entry from each array and plotting a bar to from
    # (x_data[i], y_data[i], 0) to (x_data[i], y_data[i], z_data[i]).
    #
    x_data = x_data.flatten()
    y_data = y_data.flatten()
    z_data = psf.flatten()
    ax.bar3d( x_data,
              y_data,
              np.zeros(len(z_data)),
              50, 50, z_data, alpha=0.4, color='g')

    # Finally, display the plot.
    savefig('psf_hist.png',dpi=300)
    close()
'''

############################### Convolve the data ##############################

print 'Currently convolving...'

# Use the convolution to convolve the rebinned data with the PSF
conv = convolve(data, psf_rebin, boundary='extend')

print 'Convolution complete!\nSaving the image...\n'

# Plot the resulting convolution
imshow(conv,origin='lower')
colorbar()
savefig('convolved.png',dpi=300)
close()

print 'The convolution has been performed and the convolved image saved to the working directory.\n'
