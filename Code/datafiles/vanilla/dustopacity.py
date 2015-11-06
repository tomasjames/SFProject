################################################################################
############################ 4th Year Project Code #############################
########################## Dust Opacity Determination ##########################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import Numpy
import numpy as np

########################## Define dust opacity function ########################

def dustopacity(nlam, kappa_0, w_ref, B, opaclaw='H', w_0=0.1, w_nlam=10000):

    '''
    Uses the dust opacity power law defined in Hildebrand (1983) to determine
    the dust opacity at a range of wavelengths from w_0 (0.1um by default) to
    w_nlam (1000um by default) given some reference wavelength w_ref and dust
    spectral index B. Recommended that kappa_0 = 1 by default.

    The function will generate the range of wavelengths from 0.1um to 1000um
    with nlam evenly spaced intervals. The function will then write the required
    dustkappa_silicate.inp file to the working directory. The format of this
    file follows the standard specified in the RADMC-3D manual.
    '''

    ############################ Define dust opacity ###########################

    # This is the frequency dependent dust opacity defined in Hildebrand (1983)
    # Parameterising the original developed by Hildebrand gives the below. This
    # was found in The Physics and Chemistry of the ISM (AGGM Tielens) in
    # wavelength form
    #                       opac = kappa_0*(v/v_0)**B

    ########################### Create necessary files #########################

    # Creates the dust opacity file (dustopac.inp) and opens it in write mode
    f = open('dustkappa_silicate.inp', 'w')

    ######################## Read 'wavelength_micron.inp' ######################
    '''
    # Read the file into an array
    w = np.loadtxt('wavelength_micron.inp')

    # Take the first entry (the number of wavelength points)
    nlam = w[0]

    # Take the second and last entries (i.e. the beginning and end points for
    # the wavelength range) to use later
    lambda_init = w[1]
    lambda_fin = w[-1]
    '''

    # The start and end points of the wavelength range in microns
    lambda_init = w_0
    lambda_fin = w_nlam

    ################# Evaluate opacities over wavelength range #################

    # Create array of wavelengths from range given and reshape
    w = np.vstack(np.linspace(lambda_init, lambda_fin, nlam))

    if opaclaw == 'H':
        # Evaluate opacities
        opacity = np.vstack(kappa_0*(w/w_ref)**B)

    # Concantenate arrays to create 2D array of all data in format ready to be
    # written to the .inp file
    data = np.concatenate((w,opacity), axis=1)

    ########################## Write explanatory notes #########################

    f.write('# This files contains all of the dust opacities for the power law\n')

    if opaclaw == 'H':
        f.write('# given in Hildebrand (1983): kappa_abs = kappa_0(v/v_0)**B\n')

    ########################## Write the iformat integer #######################

    # Writes the 'iformat' integer. This should always be 1 for lambda and
    # kappa_abs columns to be written
    f.write('1               # Format number of this file\n')

    ########################## Write the nlam integer ##########################

    # Writes the nlam integer (i.e. how many wavelength points are in the file)
    f.write(str(nlam) + '               # Nr of wavelength points in the file\n')

    ################################ Write data ################################

    for row in data:
        for column in row:
            f.write('%14.8f' % column)
        f.write('\n')
    f.close()

    print 'dustkappa_silicate.inp written to the working directory'
