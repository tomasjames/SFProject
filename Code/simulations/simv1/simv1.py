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

############################## Set up initial model ############################
# Writes the default parameter file for the 2d sphere model
radmc3dPy.analyze.writeDefaultParfile('spher2d_1')

# Setup the dust module with the ascii input files
radmc3dPy.setup.problemSetupDust('spher2d_1', binary=False,
                                    tstar='0.003*ts', crd_sys='sph', nx=10,
                                    ny=10, nz=10, nphot=2000000.)

'''
# Copy the dust opacity and data files from the datafiles directory
print '\nCopying datafiles from datafiles directory...\n'
os.system('cp -v ../../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../../datafiles/molecule_co.inp .')
print '\nCopy complete'
'''

################################# Pause script #################################
# This pauses the script until the user presses enter. This is in place so that
# user is able to edit the .inp files to modify the simulation paramaters
raw_input(
'\n######################################################################################\n'
'Script has been paused, and all default .inp files have been written'
' in the working directory. If amendments are required, please make'
' them now and press enter to continue:\n')

############################ Run Monte-Carlo simulation ########################
# Start a timer
start = time.time()

# Interface with operating system to run the Monte-Carlo sim
os.system('radmc3d mctherm')

# End a timer
end = time.time()

############################# Plot the resulting data ##########################
# Generate a canvas to plot over
radmc3dPy.image.makeImage(npix=1000, sizeau=200., wav=10000., incl=45.)

# Initialise the image
imag = radmc3dPy.image.readImage()

# Plot the image in a matplotlib figure
radmc3dPy.image.plotImage(imag, arcsec=True, dpc=150., log=True, maxlog=5)

print '\n######################################################################'
print 'Please run the command \"viewimage\" in the Terminal at this point to'
print 'start the GUI image viewer'
print '########################################################################'

############################# Misc. print statements ###########################
print '\nThe thermal Monte-Carlo simulation took', end-start, 's\n'
