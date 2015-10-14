################################################################################
############################ 4th Year Project Code #############################
############################## Initial Simulation ##############################
################################################################################

############################# Import statements ################################
# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

############################## Set up initial model ############################
# Writes the default parameter file for the 2d sphere model
radmc3dPy.analyze.writeDefaultParfile('spher2d_1')

# Setup the dust module with the ascii input files
radmc3dPy.setup.problemSetupDust('spher2d_1', binary=False, tstar='0.003*ts', nphot=20000.)

# Copy the dust opacity and data files from the datafiles directory
print '\nCopying datafiles from datafiles directory...\n'
os.system('cp -v ../../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../../datafiles/molecule_co.inp .')
print '\nCopy complete'

################################# Pause script #################################
# This pauses the script until the user presses enter. This is in place so that
# user is able to edit the .inp files to modify the simulation paramaters
raw_input(
'\n######################################################################################\n'
'Script has been paused, and all default .inp files have been written'
' in the working directory. If amendments are required, please make'
' them now and press enter to continue:\n')
