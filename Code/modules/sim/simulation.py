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

########################## Function to run initial sim ##########################
def simulation(filt, npix, sizeau, d):

    '''
    Begins the RADMC-3D simulation based on the keyword arguments supplied.

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

    data_choice = input('Run datafilegen.py to generate input files, or go straight to raytrace? (\'d\' for datafilegen.py, \'r\' for raytrace)\n')
    if data_choice == 'd':
        # Call the data file generation generation script to write the necessary files to the working directory
        execfile('/Users/tomasjames/Documents/University/Project/ZiggyStarDust/Code/datafiles/vanilla/datafilegen.py')

    ################################## Run raytrace ################################

    # Interface with operating system to run the Monte-Carlo sim (and allowing the
    # code to use wavelength_micron.inp)
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
        #avwav = np.mean(imag.wav)

        # Define effective filter wavelengths
        if filt == 'psw':
            cenwav = 247.12451
        elif filt == 'pmw':
            cenwav = 346.71804
        elif filt == 'plw':
            cenwav = 496.10677

        image_trans.writerow(["%.15f" % cenwav])
        # Writes a blank line to seperate the wavelength points from the intensity points
        image_trans.writerow([])

        # Begin computing transmission weighted sum
        # This is to be achieved by looping through each pixel in the image and multiplying by each transmission coefficient at that wavelength. This is then summed along all wavelengths to give a 'composite' image

        # Declare a list to store the variables
        store, store_all = [], []
        summation = []

        for i in range(0, imag.nx):
            for j in range(0, imag.ny):
                for k in range(0, imag.nwav):
                    #if i == 0. and j == 0. and k == 0.:
                        #image_trans.write('  \n')
                    store.append(np.float64(imag.image[i][j][k]*trans_data[k][1]))
                    store_all.append(np.float64(imag.image[i][j][k]*trans_data[k][1]))
                summation.append(np.float64(np.sum(store)))
                #image_trans.write(str(np.sum(store))+str('\n'))
                store = []

        image_trans.writerows(zip(np.float64(summation)))

    # Close the file for memory purposes
    #image_trans.close()

    # Save store_all to another .txt file to SED
    with open('image_trans_raw.txt', 'w') as g:
        image_trans_raw = csv.writer(g, delimiter=' ')

        image_trans_raw.writerows(zip(np.float64(store_all)))

    ########################### Save data to FITS file #######################

    # Open the image_trans_raw.txt file into an array
    image_trans_raw_1d = np.loadtxt('image_trans.out',skiprows=6)

    # Reshape the array to be 2-dimensional
    image_trans_raw = np.reshape(image_trans_raw_1d, (np.sqrt(len(image_trans_raw_1d)), np.sqrt(len(image_trans_raw_1d))))

    # Begin defining information for the header
    # Define distance to the source
    d = d*pc # Distance to source plotted in pc and converted to cm

    # Determine image width
    amr = np.loadtxt('amr_grid.inp', skiprows=6) # Read in the spatial grid
    imag_width = amr[0][-1] - amr[0][0] # Determine the image width in cm
    pix_width = imag_width/((len(amr[0]))) # Determine the pixel width in cm

    # Use small angle approximation to determine the angle subtended by the pixel
    theta_rad = pix_width/d # In degrees
    #theta = theta_rad*(360/2*np.pi)*3600 # Convert to degrees and then to arcseconds
    theta = theta_rad*3600 # Convert to degrees and then to arcseconds

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

    ########################## Plot the resulting data #######################

    # Query whether it should be plotted locally
    p = input('Shall the image be plotted?\n')

    if p == 'yes':
        # Initialise the image
        imag_trans = radmc3dPy.image.readImage('image_trans.out', binary=False)

        # Plot the image in a matplotlib figure (ifreq is the index of the lambdarange to plot)
        radmc3dPy.image.plotImage(imag_trans, arcsec=False, au=True, dpc=d, log=False, bunit='inu')

    print '\n######################################################################'
    print 'Please run the command \"viewimage\" in the Terminal at this point to'
    print 'start the GUI image viewer'
    print '########################################################################'



########################## Function to generate the SED ##########################
def sedGeneration(filt, sim_name):

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
            String
    '''

    ############################### Read in image data ##############################

    # Read in the transmission weighted image devised earlier
    image_trans_raw = np.loadtxt('image_trans_raw.txt')

    # Read in the original image file
    imag = radmc3dPy.image.readImage('image.out')
    #imag = np.reshape(imag.image, (1,len(imag.image)*len(imag.image[0])))

    # Resize the transmission weighted data
    image_trans_raw = np.reshape(image_trans_raw, (len(imag.image),len(imag.image[0]),len(imag.image[0][0])))

    print 'Data dimensions dictate that there are', imag.nx*imag.ny, 'rows before a new wavelength data entry begins.\n This means that there are', imag.nwav, 'wavelength fluxes at the central pixel.\n'

    ########################## Read in the transmission data #######################

    # Reads the transmission data
    trans_data = np.loadtxt('transmission.txt')

    # Assign to variables
    trans_wav = trans_data[:,0]
    trans_v = cc/(trans_wav*10**-4)
    trans = trans_data[:,1]

    ######################### Determine frequency of the image #####################

    # Determine parameters
    wav = imag.wav
    w_cen = 247.124514
    #flux = middle_flux_trans

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
    if os.path.isfile(str('../../../curvefitting/')+str(sim_name)+str('/psw_average_data.txt')) == True:
        print 'Data storage file already exists; opening now\n'
        save_data = open(str('../../../curvefitting/')+str(sim_name)+str('/')+str(filt)+str('_average_data.txt'), 'a+')
        data_store = csv.writer(save_data, delimiter=' ')

    else:
        print 'Data storage file does not already exist; writing now\n'
        save_data = open(str('../../../curvefitting/')+str(sim_name)+str('/')+str(filt)+str('_average_data.txt'), 'w+')
        data_store = csv.writer(save_data, delimiter=' ')

    # Sum integer to track the loop
    count = 0

    # Instantiate lists to store values
    flux, flux_trans, flux_trans_index = [], [], []

    # Loop through the array to find the middle pixel and log its flux
    for x in range(0,imag.nx):
        for y in range(0,imag.ny):
            for l in range(0,imag.nwav):

                # Assign values to the pixel numbers to track position in loop
                xpix = x
                ypix = y
                zpix = l

                raw_flux = imag.image[x][y][l]
                raw_flux_trans = image_trans_raw[x][y][l]

                flux.append(raw_flux)
                flux_trans.append(raw_flux_trans)

            # Determine the flux 'seen' by SPIRE
            weighted_flux_mean = np.sum(flux_trans)/np.sum(trans)
            weighted_flux_std = np.std(flux_trans)

            # Invoke a counter
            count += 1

            # Append th determined information to the file
            data_store.writerow([count, xpix, ypix, zpix, v_cen, weighted_flux_mean, weighted_flux_std])

            # Reset the lists to 0
            flux, flux_trans, flux_trans_index = [], [], []

    if sum == imag.nx*imag.ny*imag.nwav:
        print 'The loop has been executed over all of the elements.\n'

    # Close the file
    save_data.close()
