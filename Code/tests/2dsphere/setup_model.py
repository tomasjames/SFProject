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
radmc3dPy.setup.problemSetupDust('spher2d_1', binary=False, tstar='0.003*ts',
                                    nx=50, ny=50, nphot=20000.)

# Copy the dust opacity and data files from the datafiles directory
os.system('cp -v ../datafiles/dustkappa_silicate.inp .')
os.system('cp -v ../datafiles/molecule_co.inp .')

# Run thermal Monte-Carlo simulation
os.system('radmc3d mctherm')

# Plot image
radmc3dPy.image.makeImage(npix=400, sizeau=400, wav=10000., incl=90., posang=43.)
imag_1 = radmc3dPy.image.readImage()
radmc3dPy.image.plotImage(imag_1, arcsec=True, dpc=140., log=True, maxlog=5)

radmc3dPy.image.makeImage(npix=400, sizeau=400, wav=10000., incl=0., posang=43.)
imag_2 = radmc3dPy.image.readImage()
radmc3dPy.image.plotImage(imag_2, arcsec=True, dpc=140., log=True, maxlog=5)
