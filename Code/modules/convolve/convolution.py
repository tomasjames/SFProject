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
from astropy.convolution import convolve, convolve_fft

# Import scipy
import scipy.interpolate
import scipy.ndimage

# Import glob for file scraping
import glob

import random

########################## Function to rebin the PSF ##########################

# Function taken from http://scipy-cookbook.readthedocs.io/items/Rebinning.html Example 3
def congrid(a, newdims, method='nearest', centre=True, minusone=True):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method, bounds_error=False, fill_value=0 )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method, bounds_error=False, fill_value=0  )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None


######################### Function to convolve data and PSF ###########################

def conv(filter):

    '''
    Takes data (in the form of a FITS file) and convolves with a PSF.

    Keywords
    filter: the filter that the data considered is imaged in
        Options - SPIRE:
            psw (250 micron)
            pmw (350 micron)
            plw (500 micron)

        Options - PACS:
            70 (70 micron)
            100 (100 micron)
            160 (160 micron)

    psf: the filename of the point spread function of the filter considered
        Example:
            SPIRE 250 micron PSF has filename 'theoretical_spire_beam_model_psw_V0_2'
        This file *must* be present in folder 'psf'

    '''

    ################################# Read data ####################################

    # Determine the current working directory
    direc = os.getcwd()

    # Read in the FITS file containing the simulation results
    sim_result = fits.open(str(direc)+str('/')+str(filter)+str('.fits'))

    # Split into image data and header data
    data = sim_result[0].data
    hdu = sim_result[0].header

    # Pull out the total image width
    image_width = hdu['IMGWIDTH'] # In cm

    # Pull out the pixel width
    data_theta = hdu['PIXSIZE'] # Pixel size in "/pixel

    # Determine the size of one pixel
    data_size = (250*data_theta)/206256 # Extent of one pixel

    ############################# Read in FITS data ################################

    # Read in the PSF
    # Given that the PSF has different file names depending on the filter under consideration, use this if statement to determine
    # which is the correct one
    if filter == 'psw':
        psf_raw = fits.open('/Users/tomasjames/Documents/University/Project/ZiggyStarDust/Code/datafiles/psf/ken/psf_0250.fits')

    elif filter == 'pmw':
        psf_raw = fits.open('/Users/tomasjames/Documents/University/Project/ZiggyStarDust/Code/datafiles/psf/theoretical_spire_beam_model_pmw_V0_2.fits')

    elif filter == 'plw':
        psf_raw = fits.open('/Users/tomasjames/Documents/University/Project/ZiggyStarDust/Code/datafiles/psf/theoretical_spire_beam_model_plw_V0_2.fits')

    # Extract the data
    psf = psf_raw[0].data
    psf_header = psf_raw[0].header

    # Assign cards to variables
    psf_theta_deg = psf_header['CDELT2'] # Pixel width in deg
    psf_theta = psf_theta_deg*3600.
    #image_width = hdu['IMGWIDTH'] # Image width in cm

    # Plot the data for comparison purposes
    subplot2grid((10,8), (0,0), colspan=8,rowspan=8)
    title(r'Raw SPIRE 250$\mu m$ PSF')
    matshow(psf,origin='lower')
    colorbar()
    xlabel('Pixel')
    ylabel('Pixel')
    '''
    #tight_layout
    #grid(color='g',linestyle='-',linewidth=1)
    subplot2grid((10,8), (8,0), colspan=8,rowspan=2)
    plot(np.linspace(0,len(psf),len(psf)),psf[len(psf)/2])
    xlim(min(np.linspace(0,len(psf),len(psf))), max(np.linspace(0,len(psf),len(psf))))
    ylabel('Pixel Value')
    xlabel('Pixel')
    #tight_layout
    '''
    savefig('psf.pdf',dpi=300)
    close()

    matshow(data,origin='lower')
    title(str('RADMC-3D Intensity Data for SPIRE ')+str(hdu['EFF'])+str('$ \mu m$'))
    #colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
    colorbar(label=r'$I$ [MJy/ster]')
    savefig('data.pdf',dpi=300)
    close()

    ############################### Rebin the PSF #################################

    # Take the data pixel size and determine how much larger/smaller it is than the PSF
    size_factor = data_theta/psf_theta

    # Convolve requires that the array be odd in dimension
    # Check to see if the size_factor is an odd number
    # This ensures the PSF is centered in the frame
    if (np.round(size_factor) % 2 == 0):
        print 'The dimensions of the PSF are even; altering to make the dimensions odd...\n'

    else:
        print 'The dimensions of the PSF are odd; proceeding as normal...\n'

    # Determine whether the pixel widths of the PSF and data differ
    if size_factor == 1:
        print 'The pixel sizes match!\n'
        psf_rebin = psf

    elif size_factor < 1:
        print 'The pixel sizes do not match, and the PSF has larger pixel size. Necessary to increase the size of the RADMC-3D image.\n'

        # This will resize the original data to be 1/2 the original size using the value closest to the original (method='neighbour')
        data = congrid(data, (len(data[0])*(1./2),len(data[:,0])*(1./2)),centre=True)
        imshow(data,origin='lower')
        title(r'Rebinned Data')
        #colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
        colorbar(label=r'$I$ [MJy/ster]')
        savefig('data_rebin.png')
        close()

        # Determine the new angular size of each pixel
        data_theta = data_theta/(1./2)

        # Adjust the size factor accordingly
        size_factor = data_theta/psf_theta

        # Use the above rebin function to change the dimensions of the PSF to
        psf_rebin_unnorm = congrid(psf, ((psf_header['NAXIS1']/size_factor)+1,(psf_header['NAXIS2']/size_factor)+1),centre=True)

        # Renormalise the PSF back to unity (or close to it)
        psf_rebin = psf_rebin_unnorm*(np.sum(psf)/np.sum(psf_rebin_unnorm))

    elif size_factor > 1:

        # Use the above rebin function to change the dimensions of the PSF to
        psf_rebin_unnorm = congrid(psf, ((psf_header['NAXIS1']/size_factor),(psf_header['NAXIS2']/size_factor)),centre=True)

        # Renormalise the PSF back to unity (or close to it)
        psf_rebin = psf_rebin_unnorm*(np.sum(psf)/np.sum(psf_rebin_unnorm))

    # Determine the new pixel angular size
    new_psf_theta = psf_theta*size_factor

    # Save rebinned PSF to fits file
    fname = str('psf_rebin.fits')
    head = fits.PrimaryHDU(psf_rebin)
    head.writeto(fname)

    # Assign info to header
    fits.setval(fname, 'DISTANCE', value=hdu['DISTANCE'])
    fits.setval(fname, 'PIXSIZE', value=new_psf_theta)
    fits.setval(fname, 'EFF', value=hdu['EFF'])

    # Plot the resulting data to assess accuracy of rebin
    subplot2grid((10,8), (0,0), colspan=6,rowspan=6)
    title(r'Rebinned SPIRE 250$\mu m$ PSF')
    matshow(psf_rebin,origin='lower')
    colorbar()
    xlabel('Pixel')
    ylabel('Pixel')
    '''
    # major ticks every 20, minor ticks every 5
    major_ticks = np.arange(0, len(psf_rebin[0]), 10)
    minor_ticks = np.arange(0, len(psf_rebin[0]), 1)

    xticks(major_ticks)
    xticks(minor_ticks)
    yticks(major_ticks)
    yticks(minor_ticks)

    grid(which='minor',linestyle='-',alpha=0.2)
    grid(which='major',linestyle='-',alpha=0.2)

    subplot2grid((10,8), (8,0), colspan=8,rowspan=2)
    plot(np.linspace(0,len(psf_rebin[len(psf_rebin)/2]),len(psf_rebin[len(psf_rebin)/2])),psf_rebin[len(psf_rebin)/2])
    xlim(min(np.linspace(0,len(psf_rebin),len(psf_rebin))), max(np.linspace(0,len(psf_rebin),len(psf_rebin))))
    ylabel('Pixel Value')
    xlabel('Pixel')
    #tight_layout
    '''
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
    matshow(conv,origin='lower')
    #colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
    colorbar(label=r'$I$ [MJy/ster]')
    title(str('Convolved Data for SPIRE ')+str(hdu['EFF'])+str(r'$ \mu m$'))
    savefig('convolved.pdf',dpi=300)
    close()


    # Save to a FITS file
    # Check to see whether file has already been written
    # If so, delete the file (Python cannot overwrite a binary file)
    if os.path.isfile(str(filter)+str('_convolved.fits')) == True:
        os.remove(str(filter)+str('_convolved.fits'))

    # Do the write
    fname = str(filter)+str('_convolved.fits')
    head = fits.PrimaryHDU(conv)
    head.writeto(fname)

    # Write some header information
    # Write relevant information to the FITS header
    fits.setval(fname, 'DISTANCE', value=hdu['DISTANCE'])
    fits.setval(fname, 'IMGWIDTH', value=hdu['IMGWIDTH'])
    fits.setval(fname, 'PIXWIDTH', value=hdu['PIXWIDTH'])
    fits.setval(fname, 'PIXSIZE', value=data_theta)
    fits.setval(fname, 'EFF', value=hdu['EFF'])

    print 'The convolution has been performed and the convolved image saved to the working directory.\n'

    ############################### Check image power ##############################

    # Begin with default, non altered image
    raw_power = np.sum(np.sum(data))

    # Repeat for convolved image
    convolved_power = np.sum(np.sum(conv))

    # Print a statement to the terminal to show overall power
    print 'The initial power in the image is ', raw_power
    print 'The power in the convolved image is ', convolved_power