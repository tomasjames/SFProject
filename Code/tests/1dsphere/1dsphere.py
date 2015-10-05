################################################################################
############################## RADMC-3D Tests ##################################
############################# 1D Sphere Setup ##################################
################################################################################

############################# Import statements ################################

# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Write a parameter file containing all of the default parameters
radmc3dPy.analyze.writeDefaultParfile('spher1d_1')

# Setup the dust module with the ascii input files
radmc3dPy.setup.problemSetupDust('spher1d_1', binary=False)

# Copy the dust opacity and data files from the datafiles directory
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../datafiles/molecule_co.inp .')
