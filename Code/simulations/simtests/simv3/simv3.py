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

from radmc3dPy.natconst import *

############################## Set up initial model ############################
# Writes the default parameter file for the 2d sphere model
radmc3dPy.analyze.writeDefaultParfile('3d_cloud')

# Setup the dust module with the ascii input files
radmc3dPy.setup.problemSetupDust('3d_cloud', binary=False, tstar='1.0*ts', rstar='1.0*rs', nx=128, ny=128, nz=128, xbound=[-10000*au,10000*au], ybound=[-10000*au,10000*au], zbound=[-10000*au,10000*au], nphot=2000000.)

########################### Run Monte-Carlo simulation ########################

# Interface with operating system to run the Monte-Carlo sim (and allowing the
# code to use wavelength_micron.inp)
os.system('radmc3d image')

############################# Plot the resulting data ##########################

# Define wavelength ranges of spire to plot (PSW, PMW and PLW)
spire = [[196.5351,298.1259],[277.3117,423.4707],[386.6218,679.3126]]

# Plot image for first SPIRE wavelength band (PSW)
radmc3dPy.image.makeImage(npix=100000, sizeau=20000, incl=90., lambdarange=[196.5351,298.1259], nlam=60)

# Initialise the image
imag = radmc3dPy.image.readImage()

# Plot the image in a matplotlib figure (ifreq is the index of the lambdarange to plot)
radmc3dPy.image.plotImage(imag, arcsec=False, au=False, dpc=150., log=False, bunit='inu')

print '\n######################################################################'
print 'Please run the command \"viewimage\" in the Terminal at this point to'
print 'start the GUI image viewer'
print '########################################################################'

############################# Misc. print statements ###########################
