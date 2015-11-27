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
radmc3dPy.setup.problemSetupDust('3d_cloud', binary=False, tstar='1.0*ts', rstar='1.0*rs',
                                    nx=128, ny=128, nz=128, xbound=[-10000*au,10000*au],
                                    ybound=[-10000*au,10000*au], zbound=[-10000*au,10000*au],
                                    nphot=2000000.)

################################# Pause script #################################
# This pauses the script until the user presses enter. This is in place so that
# user is able to edit the .inp files to modify the simulation paramaters
raw_input(
'\n######################################################################################\n'
'Script has been paused, and all default .inp files have been written'
' in the working directory. If amendments are required, please make'
' them now and press enter to continue:\n')

############################ Run Monte-Carlo simulation ########################

# Interface with operating system to run the Monte-Carlo sim (and allowing the
# code to use wavelength_micron.inp)
os.system('radmc3d image')

############################# Plot the resulting data ##########################
# Generate a canvas to plot over
#radmc3dPy.image.makeImage(npix=1000, sizeau=45000., incl=90., pointau=[15, 15, 15], lambdarange=[5.,100.], nlam=96)

radmc3dPy.image.makeImage(npix=10000, sizeau=20000., wav=850., incl=90.)

# Initialise the image
imag = radmc3dPy.image.readImage()

# Plot the image in a matplotlib figure (ifreq is the index of the lambdarange to plot)
radmc3dPy.image.plotImage(imag, arcsec=False, au=False, dpc=150., log=False, bunit='inu')

print '\n######################################################################'
print 'Please run the command \"viewimage\" in the Terminal at this point to'
print 'start the GUI image viewer'
print '########################################################################'

############################# Misc. print statements ###########################
