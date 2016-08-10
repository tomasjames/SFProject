################################################################################
############################ 4th Year Project Code #############################
##############################  Functions to run  ##############################
################################################################################

############################# Import statements ################################
# Import matplotlib
import matplotlib as mpl
mpl.use('Qt4Agg')
from matplotlib.pyplot import *

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

# Import astropy
from astropy.io import fits

# Import scipy
import scipy.interpolate
import scipy.ndimage

# Import glob for file scraping
import glob

import random

# Import the necessary modules
from inputfile.datafilegen import *
from convolve.convolution import *

# Import photutils
import photutils

########################## Function to run initial sim ##########################
def simulation(mode, filt, npix, sizeau, d, mass, cloud_density, outside_density, cloud_temp, outside_temp, kappa_0, lambda_0, B, amr, dust):

    '''
    Begins the RADMC-3D simulation based on the keyword arguments supplied.

    Keywords
    mode: whether the code should generate all input files or just perform a raytrace
        Options:
            'd': run the data file script first AND THEN run the raytrace
            'r': run the raytrace ONLY

    filt: the filter under consideration
        Options - SPIRE:
            psw (250 micron)
            pmw (350 micron)
            plw (500 micron)

        Options - PACS:
            70 (70 micron)
            100 (100 micron)
            160 (160 micron)

    npix: the number of pixels in any given dimensions. The total number of pixels will be npix**3
        Options:
            Float

    sizeau: the image width in astronomical units
        Options:
            Float

    d: the distance to the source/object in parsecs
    '''

    ###################### Check for correct default.inp files #####################

    print '########################################################################'
    print '########################## Default File Check ##########################'
    print '########################################################################'

    # Check to see whether radmc3d.inp exists. If so, leave it alone and if not write it to the working directory
    if os.path.isfile('radmc3d.inp') == True:
        print '\nradmc3d.inp already exists; no need to write/overwrite\n'
    else:
        # This line writes the radmc3d.inp file and then immediately closes it
        open('radmc3d.inp', 'w').close()
        print '\n--->radmc3d.inp did not already exists; a blank file called radmc3d.inp has been written to the working directory\n'

    ############################## Set up initial model ############################

    #data_choice = input('Run datafilegen.py to generate input files, or go straight to raytrace? (\'d\' for datafilegen.py, \'r\' for raytrace)\n')
    data_choice = mode
    if data_choice == 'd':
        # Call the data file generation generation script to write the necessary files to the working directory
        datafilegen(m=mass, cloud=cloud_density, outside=outside_density, width=sizeau, l=npix, cloud_temperature=cloud_temp, outside_temperature=outside_temp, nlam=10000, opaclaw='H', kappa_0=kappa_0, lambda_0=lambda_0, B=B, amr_gen=amr, dust_gen=dust, points=100)

    ################################## Run raytrace ################################

    # Interface with operating system to run the Monte-Carlo sim (and allowing the
    # code to use wavelength_micron.inp)
    print(str('Now invoking the command: radmc3d image loadlambda npix ')+str(npix)+str(' sizeau ')+str(sizeau))
    os.system(str('radmc3d image loadlambda npix ')+str(npix)+str(' sizeau ')+str(sizeau))

    ########################### Account for transmission ###########################

    # Initialise the data from the ray trace
    # This is stored in the class radmc3dImage
    imag = radmc3dPy.image.readImage('image.out')

    # Read in the specially created file holding the interpolated transmissions
    trans_data = np.loadtxt('transmission.txt')

    # Write a new image.out file called image_trans.out that will contain all of the new intensity information
    with open('image_trans.out', 'w') as f:
        #image_trans = open('image_trans.out', 'w')
        image_trans = csv.writer(f, delimiter=' ')

        # Begin writing of the file by writing the format number
        #image_trans.writerows(str(1)+str('\n'))
        image_trans.writerow([1])

        # Write the number of pixels in x and y dimensions
        #image_trans.write("%.f" % imag.nx)
        #image_trans.write(" %.f \n" % imag.ny)
        image_trans.writerow([imag.nx, imag.ny])

        # Write the number of wavelength points
        #image_trans.write(str('          ')+str(imag.nwav)+str('\n'))
        #image_trans.writerow([imag.nwav])
        image_trans.writerow([1])

        # Write the pixel sizes
        #image_trans.write("%.f" % imag.sizepix_x)
        #image_trans.write(" %.f \n" % imag.sizepix_y)
        image_trans.writerow([imag.sizepix_x, imag.sizepix_y])

        # Because the image is a composite of all wavelengths, only write the average of the wavelength points
        # Define effective filter wavelengths
        if filt == 'psw':
            cenwav = 247.12451
        elif filt == 'pmw':
            cenwav = 346.71804
        elif filt == 'plw':
            cenwav = 496.10677
        elif filt == 'blue':
            cenwav = 68.92474
        elif filt == 'green':
            cenwav = 97.90361
        elif filt == 'red':
            cenwav = 153.94392

        image_trans.writerow(["%.15f" % cenwav])
        # Writes a blank line to seperate the wavelength points from the intensity points
        image_trans.writerow([])

        # Begin computing transmission weighted sum to produce composite image
        # This is to be achieved by looping through each pixel in the image and multiplying by each transmission coefficient at that wavelength. This is then summed along all wavelengths to give a 'composite' image

        # Declare a list to store the variables
        store, store_all, trans_store = [], [], []
        summation = []

        for i in range(0, imag.nx):
            for j in range(0, imag.ny):
                for k in range(0, imag.nwav):
                    #if i == 0. and j == 0. and k == 0.:
                        #image_trans.write('  \n')
                    trans_store.append(np.float64(trans_data[k][1])) # Appends the transmission to a list
                    store.append(np.float64(imag.image[j][i][k]*trans_data[k][1])) # Determines the transmission weighting
                    store_all.append(np.float64(imag.image[j][i][k]*trans_data[k][1]))
                summation.append(np.float64(np.sum(store)/np.sum(trans_store))) # Reduces that weighting to one discrete point
                #image_trans.write(str(np.sum(store))+str('\n'))
                store, trans_store = [], []

        image_trans.writerows(zip(np.float64(summation)))

    # Close the file for memory purposes
    #image_trans.close()

    # Save store_all to another .txt file to SED
    with open('image_trans_raw.txt', 'w') as g:
        image_trans_raw = csv.writer(g, delimiter=' ')

        image_trans_raw.writerows(zip(np.float64(store_all)))

    ########################### Save data to FITS file #######################

    # Open the image_trans_raw.txt file into an array and convert to MJy/sr
    image_trans_raw_1d = np.loadtxt('image_trans.out',skiprows=6)/(10**(-23))

    # Reshape the array to be 2-dimensional
    image_trans_raw = np.reshape(image_trans_raw_1d, (int(np.sqrt(len(image_trans_raw_1d))), int(np.sqrt(len(image_trans_raw_1d)))))

    # Begin defining information for the header
    # Define distance to the source
    d = d*pc # Distance to source plotted in pc and converted to cm

    # Determine image width
    #amr = np.genfromtxt('amr_grid.inp', skip_header=6, filling_values=np.NaN) # Read in the spatial grid
    #imag_width = amr[0][-1] - amr[0][0] # Determine the image width in cm
    imag_width = sizeau*au
    pix_width = imag_width/(npix) # Determine the pixel width in cm

    # Use small angle approximation to determine the angle subtended by the pixel
    theta_rad = np.arctan(pix_width/d) # In radians
    theta = theta_rad*(360/(2*np.pi))*3600 # Convert to degrees and then to arcseconds

    # Check to see whether file has already been written
    # If so, delete the file (Python cannot overwrite a binary file)
    if os.path.isfile(str(filt)+str('.fits')) == True:
        os.remove(str(filt)+str('.fits'))

        # Save the image as a FITS file (declare the file name first)
        fname = str(filt)+str('.fits')
        hdu = fits.PrimaryHDU(image_trans_raw)
        hdu.writeto(fname)

    else:
        # Save the image as a FITS file (declare the file name first)
        fname = str(filt)+str('.fits')
        hdu = fits.PrimaryHDU(image_trans_raw)
        hdu.writeto(fname)

    # Write relevant information to the FITS header
    fits.setval(fname, 'DISTANCE', value=d)
    fits.setval(fname, 'IMGWIDTH', value=imag_width)
    fits.setval(fname, 'PIXWIDTH', value=pix_width)
    fits.setval(fname, 'PIXSIZE', value=theta)
    fits.setval(fname, 'EFF', value=cenwav)

    # Close the file
    #hdu.close()

    # Read in the FITS file containing the simulation results
    sim_result = fits.open(str(filt)+str('.fits'))

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
    if filt == 'psw':
	    #psf_raw = fits.open('../../../../datafiles/psf/ken/psf_0250.fits')
        psf_raw = fits.open('../../../../../datafiles/psf/theoretical_spire_beam_model_psw_V0_2.fits')

    elif filt == 'pmw':
        psf_raw = fits.open('../../../../../datafiles/psf/theoretical_spire_beam_model_pmw_V0_2.fits')

    elif filt == 'plw':
        psf_raw = fits.open('../../../../../datafiles/psf/theoretical_spire_beam_model_plw_V0_2.fits')

    elif filt == 'red':
        psf_raw = fits.open('../../../../../datafiles/psf/PSF_red_slope-1_small.fits')

    elif filt == 'blue':
        psf_raw = fits.open('../../../../../datafiles/psf/PSF_blue_slope-1_small.fits')

    elif filt == 'green':
        psf_raw = fits.open('../../../../../datafiles/psf/PSF_green_slope-1_small.fits')

    # Extract the data
    psf = psf_raw[0].data[0]
    psf_header = psf_raw[0].header

    # Ensure that the data is normalised
    if np.round(np.sum(psf)) == '1.':
        print 'PSF power is normalised!\n'
    else:
        print 'PSF power is not normalising; normalising now!\n'

        # Determine the total power
        tot_power = np.sum(psf)

        # Normalise the PSF
        psf = psf/tot_power

    # Assign cards to variables
    psf_theta_deg = psf_header['CDELT2'] # Pixel width in deg
    psf_theta = psf_theta_deg*3600.
    #image_width = hdu['IMGWIDTH'] # Image width in cm

    # Plot the data for comparison purposes
    #subplot2grid((10,8), (0,0), colspan=8,rowspan=8)
    title(r'Raw SPIRE 250$\mu m$ PSF')
    #matshow(psf,origin='lower')
    imshow(psf,interpolation='nearest',origin='lower')
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

    #matshow(data,origin='lower')
    imshow(data,interpolation='nearest',origin='lower')
    title(str('RADMC-3D Intensity Data for SPIRE ')+str(hdu['EFF'])+str('$ \mu m$'))
    #colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]')
    colorbar(label=r'$I$ [MJy/ster]')
    savefig('data.pdf',dpi=300)
    close()

    ############################### Rebin the PSF #################################

    # Take the data pixel size and determine how much larger/smaller it is than the PSF
    size_factor = data_theta/psf_theta

    # Determine data dimension
    dimen = np.round(psf_header['NAXIS1']/size_factor)

    # Convolve requires that the array be odd in dimension
    # Check to see if the size_factor is an odd number
    # This ensures the PSF is centered in the frame
    if (dimen % 2 == 0):
        print 'The dimensions of the PSF are even; altering to make the dimensions odd...\n'
        dimen = dimen+1
    else:
        print 'The dimensions of the PSF are odd; proceeding as normal...\n'

    # Determine whether the pixel widths of the PSF and data differ
    if size_factor == 1:
        print 'The pixel sizes match!\n'
        psf_rebin = psf

    elif size_factor < 1:
        print 'The pixel sizes do not match, and the PSF has larger pixel size. Necessary to increase the size of the RADMC-3D image.'
        print 'Note: this is not strictly accurate; rebinning designed to decrease dimensions rather than increase!\n'

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
        psf_rebin_unnorm = congrid(psf, (dimen,dimen), centre=True)

        # Renormalise the PSF back to unity (or close to it)
        psf_rebin = psf_rebin_unnorm*(np.sum(psf)/np.sum(psf_rebin_unnorm))

    elif size_factor > 1:

        # Use the above rebin function to change the dimensions of the PSF to
        psf_rebin_unnorm = congrid(psf, (dimen,dimen), centre=True)

        # Renormalise the PSF back to unity (or close to it)
        psf_rebin = psf_rebin_unnorm*(np.sum(psf)/np.sum(psf_rebin_unnorm))

    # Determine the new pixel angular size
    new_psf_theta = psf_theta*size_factor

    # Check to see whether file has already been written
    # If so, delete the file (Python cannot overwrite a binary file)
    if os.path.isfile('psf_rebin.fits') == True:
        os.remove('psf_rebin.fits')

        # Save rebinned PSF to fits file
        fname = str('psf_rebin.fits')
        head = fits.PrimaryHDU(psf_rebin)
        head.writeto(fname)

    else:
        # Save rebinned PSF to fits file
        fname = str('psf_rebin.fits')
        head = fits.PrimaryHDU(psf_rebin)
        head.writeto(fname)

    # Assign info to header
    fits.setval(fname, 'DISTANCE', value=hdu['DISTANCE'])
    fits.setval(fname, 'PIXSIZE', value=new_psf_theta)
    fits.setval(fname, 'EFF', value=hdu['EFF'])

    # Plot the resulting data to assess accuracy of rebin
    #subplot2grid((10,8), (0,0), colspan=6,rowspan=6)
    title(r'Rebinned SPIRE 250$\mu m$ PSF')
    #matshow(psf_rebin,origin='lower')
    imshow(psf_rebin,interpolation='nearest',origin='lower')
    colorbar()
    xlabel('Pixel')
    ylabel('Pixel')
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
    conv = convolve(data, psf_rebin, boundary='extend', normalize_kernel=True)
    #conv = convolve_fft(data, psf_rebin, boundary='extend')

    print 'Convolution complete!\nSaving the image...\n'

    # Plot the resulting convolution
    #matshow(conv,origin='lower')
    imshow(conv,interpolation=None,origin='lower')
    #colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]',use_gridspec=False)
    colorbar(label=r'$I$ [MJy/ster]')
    title(str('Convolved Data for SPIRE ')+str(hdu['EFF'])+str(r'$ \mu m$'))
    savefig('convolved.pdf',dpi=300)
    close()

    # Save to a FITS file
    # Check to see whether file has already been written
    # If so, delete the file (Python cannot overwrite a binary file)
    if os.path.isfile(str(filt)+str('_convolved.fits')) == True:
        os.remove(str(filt)+str('_convolved.fits'))

    # Do the write
    fname = str(filt)+str('_convolved.fits')
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

    # Write another FITS file (one that adheres to the FITS standard) for use with PPMAP
    # Check to see whether file has already been written
    # If so, delete the file (Python cannot overwrite a binary file)
    if os.path.isfile(str(filt)+str('_ppmap.fits')) == True:
        os.remove(str(filt)+str('_ppmap.fits'))

        # Save the image as a FITS file (declare the file name first)
        fname = str(filt)+str('_ppmap.fits')
        hdu = fits.PrimaryHDU(conv)
        hdu.writeto(fname)

    else:
        # Save the image as a FITS file (declare the file name first)
        fname = str(filt)+str('_ppmap.fits')
        hdu = fits.PrimaryHDU(conv)
        hdu.writeto(fname)

    # Write relevant information to the FITS header
    fits.setval(fname, 'CDELT1', value=(-1)*np.float64(new_psf_theta)/3600) # Degrees per pixel (CDELT1 has to be negative)
    fits.setval(fname, 'CDELT2', value=np.float64(new_psf_theta)/3600)
    fits.setval(fname, 'CRPIX1', value=len(conv[0])/2) # Pulls out reference pixel
    fits.setval(fname, 'CRPIX2', value=len(conv[:,0])/2) # Pulls out reference pixel
    fits.setval(fname, 'CRVAL1', value=159.4507) # Coordinate of that pixel
    fits.setval(fname, 'CRVAL2', value=-19.8196)
    fits.setval(fname, 'CTYPE1', value='GLON-CAR') # Changes the coordinate system to Galactic Longitude
    fits.setval(fname, 'CTYPE2', value='GLAT-CAR') # Changes the coordinate system to Galactic Latitude

    ############################### Check image power ##############################

    # Begin with default, non altered image
    raw_power = np.sum(np.sum(data))

    # Repeat for convolved image
    convolved_power = np.sum(np.sum(conv))

    # Print a statement to the terminal to show overall power
    print 'The initial power in the image is ', raw_power
    print 'The power in the convolved image is ', convolved_power

    ################################################################################
    ############################## Commonise image scales ##########################
    ################################################################################

    # Check whether the filter under consideration is within the list of used filters
    # If so, run the convolution routine
    if filt in ('psw','pmw','blue','green','red'):

        # Read in the kernel
        if filt == 'psw':
            kern = fits.open('../../../../../datafiles/psf/comconv/spire250_to_spire500.fits')
        elif filt == 'pmw':
            kern = fits.open('../../../../../datafiles/psf/comconv/spire350_to_spire500.fits')
        elif filt == 'blue':
            kern = fits.open('../../../../../datafiles/psf/comconv/pacs70_to_spire500.fits')
        elif filt == 'green':
            kern = fits.open('../../../../../datafiles/psf/comconv/pacs100_to_spire500.fits')
        elif filt == 'red':
            kern = fits.open('../../../../../datafiles/psf/comconv/pacs160_to_spire500.fits')

        # Extract data and header
        kern_vals = kern[0].data
        kern_hdu = kern[0].header

        # Plot the initial kernel
        #subplot2grid((10,8), (0,0), colspan=6,rowspan=6)
        title(str(r'Conversion Kernel'))
        #matshow(kern_vals,origin='lower')
        imshow(kern_vals,interpolation=None,origin='lower')
        colorbar()
        xlabel('Pixel')
        ylabel('Pixel')
        savefig('kern.pdf',dpi=300)
        close()

        # Assume that CD1_1 and CD2_2 give the pixel scales in RA and DEC respectively
        # This is in degrees per pixel so requires conversion to "
        kern_theta = kern_hdu['CD1_1']*3600

        # Read in the convolved data to determine new dimensions
        convolved_hdu = fits.open(str(filt)+str('_convolved.fits'))[0].header

        # Extract pixel size
        data_theta = convolved_hdu['PIXSIZE']

        # Determine the size factor
        size_factor = data_theta/kern_theta

        # Determine data dimension
        dimen = np.round(kern_hdu['NAXIS1']/size_factor)

        # Convolve requires that the array be odd in dimension
        # Check to see if the size_factor is an odd number
        # This ensures the PSF is centered in the frame
        if (dimen % 2 == 0):
            print 'The dimensions of the PSF are even; altering to make the dimensions odd...\n'
            dimen = dimen+1

        else:
            print 'The dimensions of the PSF are odd; proceeding as normal...\n'

        # Use the above rebin function to change the dimensions of the PSF to
        kern_rebin_unnorm = congrid(kern_vals, (dimen,dimen), centre=True)

        # Renormalise the PSF back to unity (or close to it)
        kern_rebin = kern_rebin_unnorm*(np.sum(kern_vals)/np.sum(kern_rebin_unnorm))

        # Plot the resulting data to assess accuracy of rebin
        #subplot2grid((10,8), (0,0), colspan=6,rowspan=6)
        title(str(r'Rebinned Conversion Kernel'))
        #matshow(kern_rebin,origin='lower')
        imshow(kern_rebin,interpolation='nearest',origin='lower')
        colorbar()
        xlabel('Pixel')
        ylabel('Pixel')
        savefig('kern_rebin.pdf',dpi=300)
        close()

        # Determine the new pixel angular size
        new_kern_theta = kern_theta*size_factor

        # Determine dimensions of the new image
        kern_rebin_x = len(kern_rebin)
        kern_rebin_y = len(kern_rebin[0])

        # Print statements for debug
        print 'The raw PSF "/pixel:', kern_theta, '\n'
        print 'The raw data "/pixel:', data_theta, '\n'
        print 'The rebinned PSF "/pixel:', new_kern_theta, '\n'

        print 'The PSF image has now been rebinned to have dimensions', kern_rebin_x, 'x', kern_rebin_y, '. The new image has also been saved to the working directory.\n'

        # Check to see whether file has already been written
        # If so, delete the file (Python cannot overwrite a binary file)
        if os.path.isfile('kern_rebin.fits') == True:
            os.remove('kern_rebin.fits')

            # Save rebinned PSF to fits file
            fname = str('kern_rebin.fits')
            head = fits.PrimaryHDU(kern_rebin)
            head.writeto(fname)

        else:
            # Save rebinned PSF to fits file
            fname = str('kern_rebin.fits')
            head = fits.PrimaryHDU(kern_rebin)
            head.writeto(fname)

        # Save the common resolution image to file
        print 'Currently convolving to common resolution...'

        # Perform the convolution
        common_conv = convolve(conv, kern_rebin, boundary='extend', normalize_kernel=True)

        print 'Convolution complete!\nSaving the image...\n'

    # Otherwise, just leave the data as is
    else:
        print 'Common convolution has not been performed as the filter being considered is the PLW filter - i.e.'
        print 'the data is already at the correct resolution.\n'

        # Instantiate a 'fake' variable
        # This is because the PLW has not undergone common convolution, but further scripts will look for that filename
        # Therefore assign the variable (and save the file) in the same way to make further scripts more uniform
        common_conv = conv

    # Reopen a file to extract the necessary metadata for images
    hdu = fits.open(str(filt)+str('_convolved.fits'))[0].header

    # Plot the resulting convolution
    #matshow(common_conv,origin='lower')
    imshow(common_conv,interpolation=None,origin='lower')
    #colorbar(label=r'$I_{\nu}$ [erg/s/cm/cm/Hz/ster]',use_gridspec=False)
    colorbar(label=r'$I$ [MJy/ster]')
    title(str(r'Convolved Data to SPIRE 500$\mu m$ for SPIRE ')+str(hdu['EFF'])+str(r'$ \mu m$'))
    savefig('common_convolved.pdf',dpi=300)
    close()


    # Save to a FITS file
    # Check to see whether file has already been written
    # If so, delete the file (Python cannot overwrite a binary file)
    if os.path.isfile(str(filt)+str('_common_convolved.fits')) == True:
        os.remove(str(filt)+str('_common_convolved.fits'))

    # Do the write
    fname = str(filt)+str('_common_convolved.fits')
    head = fits.PrimaryHDU(common_conv)
    head.writeto(fname)

    # Write some header information
    # Write relevant information to the FITS header
    fits.setval(fname, 'DISTANCE', value=hdu['DISTANCE'])
    fits.setval(fname, 'IMGWIDTH', value=hdu['IMGWIDTH'])
    fits.setval(fname, 'PIXWIDTH', value=hdu['PIXWIDTH'])
    fits.setval(fname, 'PIXSIZE', value=data_theta)
    fits.setval(fname, 'EFF', value=hdu['EFF'])

    print 'The convolved image saved to the working directory.\n'

    ############################### Check image power ##############################

    # Begin with default, non altered image
    raw_power = np.sum(np.sum(conv))

    # Repeat for convolved image
    convolved_power = np.sum(np.sum(common_conv))

    # Print a statement to the terminal to show overall power
    print 'The initial power in the image is ', raw_power
    print 'The power in the convolved image is ', convolved_power

    return common_conv

    '''
    ################################################################################
    ################################# Generate the SED #############################
    ################################################################################

    ################################## Open the image ##############################

    # Read in the convolved and transmission weighted image
    img = fits.open(str(filt)+str('_common_convolved.fits'))
    img_data = img[0].data
    img_header = img[0].header

    ######################### Determine frequency of the image #####################

    # Determine parameters
    wav = imag.wav

    # Loop through possible filters to find effective wavelength
    if filt == 'psw':
        w_cen = 247.124514
    elif filt == 'pmw':
        w_cen = 346.71804
    elif filt == 'plw':
        w_cen = 496.10677
    elif filt == 'blue':
        w_cen = 68.92474
    elif filt == 'green':
        w_cen = 97.90361
    elif filt == 'red':
        w_cen = 153.94392

    # Convert to frequency
    v = cc/(wav*10**-4)
    v_cen = cc/(w_cen*10**-4)

    ############################ Sample flux at each pixel #########################

    # Check to see if folder exists before entering file write
    if os.path.isdir(str('../../../curvefitting/')+str(sim_name)) == True:
        print 'The folder', sim_name, 'already exists. Moving on to file write...'
    else:
        print 'The folder', sim_name, 'does not already exist. Will now create it...'
        os.makedirs(str('../../../curvefitting/')+str(sim_name))

    # Check to see if file already exists
    if os.path.isfile(str('../../../curvefitting/')+str(sim_name)+str(filt)+str('_average_data.txt')) == True:
        print 'Data storage file already exists; opening now\n'
        save_data = open(str('../../../curvefitting/')+str(sim_name)+str('/')+str(filt)+str('_average_data.txt'), 'a+')
        data_store = csv.writer(save_data, delimiter=' ')

    else:
        print 'Data storage file does not already exist; writing now\n'
        save_data = open(str('../../../curvefitting/')+str(sim_name)+str('/')+str(filt)+str('_average_data.txt'), 'w+')
        data_store = csv.writer(save_data, delimiter=' ')

    # Invoke a counter
    count = 0

    # Loop through the pixels in the image
    for i in range(0,len(img_data[0])):
        for j in range(0,len(img_data[:,0])):

            # Add to the counter
            count =+ 1

            # Pull out the pixel value (img_data[i,j]) and write to a row
            data_store.writerow([count, v_cen, img_data[i,j]])

    # Close the file
    save_data.close()
    '''

########################## Function to generate the SED ##########################
def sedGeneration(filt, sim_name, kappa_0, lambda_0, B, withPSF=True):

    '''
    Generates one point of an SED on a pixel by pixel basis for the RADMC-3D output returned in the function sim().
    The function will write the SED points to a text file, <filter>_average_data.txt.

    For a full SED, sedGeneration() should be called for all filters imaged in, building a bank of <filter>_average_data.txt files.

    Keywords
    filt: the filter under consideration
        Options - SPIRE:
            psw (250 micron)
            pmw (350 micron)
            plw (500 micron)

        Options - PACS:
            70 (70 micron)
            100 (100 micron)
            160 (160 micron)


    sim_name: a name to be given to the simulation in order to allow it to write to organisational folders of this name
        Options:
            'cloud'
            'sphdata'
            'herschel_snaps'

    limit: the pixel value required in which to mask the data. This should be the lowest pixel value of the source encountered.
        Options:
            Float

    B: the dust emissitivity index
        Options:
            '1.8'
            '2.0'
            '2.2'

    withPSF: boolean argument that will generate an SED for the PSF convolved data OR the original data
        Options:
            True (for PSF)
            False (without PSF)
    '''

    ################################## Open the image ##############################

    if withPSF == True:
        # Read in the convolved and transmission weighted image
        img = fits.open(('{}_common_convolved.fits').format(filt))
    else:
        # Read in the convolved and transmission weighted image
        img = fits.open(('{}.fits').format(filt))


    img_data = img[0].data
    img_header = img[0].header

    ######################### Determine frequency of the image #####################

    # Loop through possible filters to find effective wavelength
    if filt == 'psw':
        w_cen = 247.124514
    elif filt == 'pmw':
        w_cen = 346.71804
    elif filt == 'plw':
        w_cen = 496.10677
    elif filt == 'blue':
        w_cen = 68.92474
    elif filt == 'green':
        w_cen = 97.90361
    elif filt == 'red':
        w_cen = 153.94392

    # Convert to frequency
    v_cen = cc/(w_cen*10**-4)

    ############################ Sample flux at each pixel #########################

    '''
    # Check to see if folder exists before entering file write
    if os.path.isdir(str('../../../../curvefitting/B=')+str(B)+str('/')+str(sim_name)) == True:
        print 'The folder', sim_name, 'already exists. Moving on to file write...'
    else:
        print 'The folder', sim_name, 'does not already exist. Will now create it...'
        os.makedirs(str('../../../../curvefitting/')+str(B)+str('/')+str(sim_name))
    '''

    # Check to see if file already exists
    if os.path.isfile(str('../../../../curvefitting/')+str(sim_name)+str('/B=')+str(B)+str('/')+str(filt)+str('_average_data.txt')) == True:
        print 'Data storage file already exists; opening now\n'
        save_data = open(str('../../../../curvefitting/')+str(sim_name)+str('/B=')+str(B)+str('/')+str(filt)+str('_average_data.txt'), 'a+')
        data_store = csv.writer(save_data, delimiter=' ')

    else:
        print 'Data storage file does not already exist; writing now\n'
        save_data = open(str('../../../../curvefitting/')+str(sim_name)+str('/B=')+str(B)+str('/')+str(filt)+str('_average_data.txt'), 'w+')
        data_store = csv.writer(save_data, delimiter=' ')

    # Assign the flux calibration error based on the filter by checking whether the filter is PACS (first if) of SPIRE (else)
    if filt > ['blue', 'green', 'red']:
        sigma_percent = 0.15
    else:
        sigma_percent = 0.10

    # Invoke a counter
    count = 0

    # Loop through the pixels in the image
    for i in range(0,len(img_data[0])):
        for j in range(0,len(img_data[:,0])):

            # Add to the counter
            count =+ 1

            # Pull out the pixel value (img_data[i,j]) and write to a row
            data_store.writerow([count, v_cen, img_data[i,j], sigma_percent*img_data[i,j]])

    # Close the file
    save_data.close()
