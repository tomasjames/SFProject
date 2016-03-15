################################################################################
############################ 4th Year Project Code #############################
############################## Initial Simulation ##############################
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

from matplotlib.pyplot import *

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
    print '\n--->radmc3d.inp did not alreadt exists; a blank file called radmc3d.inp has been written to the working directory\n'

############################## Set up initial model ############################
# Call the data file generation generation script to write the necessary files to the working directory
execfile('../../../../datafiles/vanilla/datafilegen.py')

############################## Set up initial model ############################
# Writes the default parameter file for the 2d sphere model
radmc3dPy.analyze.writeDefaultParfile('3d_cloud')

# Setup the dust module with the ascii input files
radmc3dPy.setup.problemSetupDust('3d_cloud', binary=False, nx=128, ny=128, nz=128, xbound=[-15000*au,15000*au], ybound=[-15000*au,15000*au], zbound=[-15000*au,15000*au], nphot=2000000.)

########################### Run Monte-Carlo simulation ########################

# Interface with operating system to run the Monte-Carlo sim (and allowing the
# code to use wavelength_micron.inp)
os.system('radmc3d image loadlambda')

############################# Plot the resulting data ##########################

# Define wavelength ranges of spire to plot (PSW, PMW and PLW)
#spire = [[196.5351,298.1259],[277.3117,423.4707],[386.6218,679.3126]]
plw_ext = [391.4346,690.8139]

# Plot image for first SPIRE wavelength band (PSW)
#radmc3dPy.image.makeImage(npix=150000, sizeau=20000, incl=90., lambdarange=plw_ext, nlam=60)

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
    avwav = np.mean(imag.wav)
    image_trans.writerow(["%.15f" % avwav])
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

# Close the file for memory purposes
#image_trans.close()

########################## Plot the resulting data #######################

# Initialise the image
imag_trans = radmc3dPy.image.readImage('image_trans.out', binary=False)

# Plot the image in a matplotlib figure (ifreq is the index of the lambdarange to plot)
radmc3dPy.image.plotImage(imag_trans, arcsec=False, au=True, dpc=150., log=False, bunit='inu')


print '\n######################################################################'
print 'Please run the command \"viewimage\" in the Terminal at this point to'
print 'start the GUI image viewer'
print '########################################################################'
