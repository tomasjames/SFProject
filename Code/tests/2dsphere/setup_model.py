################################################################################
############################## RADMC-3D Tests ##################################
############################# 2D Sphere Setup ##################################
################################################################################

############################# Import statements ################################
# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Import matplotlib
import matplotlib.pyplot as plt

# Write a parameter file containing all of the default parameters
radmc3dPy.analyze.writeDefaultParfile('spher2d_1')

# Setup the dust module with the ascii input files
radmc3dPy.setup.problemSetupDust('spher2d_1', binary=False)

# Copy the dust opacity and data files from the datafiles directory
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../datafiles/molecule_co.inp .')

# Run thermal Monte-Carlo simulation
os.system('radmc3d mctherm')

# Plot image
radmc3dPy.image.makeImage(npix=400, sizeau=400, wav=800., incl=45., posang=43.)

imag = radmc3dPy.image.readImage()
radmc3dPy.image.plotImage(imag, arcsec=True, dpc=140., log=True, maxlog=5)
