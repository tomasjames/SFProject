################################################################################
############################ 4th Year Project Code #############################
########################## Custom Data File Generation #########################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import Numpy
import numpy as np

print '\n########################################################################'
print 'This script will request some information before writing multiple required'
print 'custom files to the working directory for use with RADMC-3D. These should'
print 'be copied and pasted to the directory in which RADMC-3D\'s thermal Monte'
print 'Carlo simulation or raytrace is run.'
print '########################################################################\n'

################################################################################
############################### Set up amr_grid.inp ############################
################################################################################

print '########################################################################'
print '                             amr_grid.inp                                 '
print '########################################################################\n'

######################## Write file to working directory #######################

# Writes a file called amr_grid.inp to the working directory and allows it to be
# written to
amr = open('amr_grid.inp', 'w')

############################ Define axes to work with ##########################

# Define the axis length (i.e. the number of pixels in each axis)
l = input('How many pixels should I assign to a dimension? (128 recommended) \n')

# Define the width of each pixel
dpix = input('And what width would you like each pixel to be? (This should be in AU) \n')

# Create an array of pixels in each dimension
x,y,z = np.linspace(0,l*dpix,l+1), np.linspace(0,l*dpix,l+1), np.linspace(0,l*dpix,l+1)

################# Write integers for RADMC-3D to learn about grid ##############

# Write the iformat number
amr.write('1         # The same as most other files\n')

# Write the grid style
amr.write('0         # 0 = regular\n')

# Write coordinate system
amr.write('1         # < 100 = Cartesian\n')

# Write gridinfo integer
amr.write('0         # If 1, abundant information is written regarding the grid\n')

# Write integers to activate dimensions
amr.write('1 1 1 \n')

################# Loop through all grid values and write to file ###############

# Write the number of base cells in each dimension
amr.write(str(l)+ ' ' + str(l) + ' ' +str(l) + '\n')

# Loop through each pixel and write to the file the edges of each pixel
for i in x:
    amr.write(str(i) + str(' '))

    if i == x[-1]:
        amr.write(str(i) + str('\n'))

for j in y:
    amr.write(str(j) + str(' '))

    if j == y[-1]:
        amr.write(str(j) + str('\n'))

for k in z:
    amr.write(str(k) + str(' '))

    if k == z[-1]:
        amr.write(str(k) + str('\n'))

print '\'amr_grid.inp\' has been written to the working directory\n'

################################################################################
############################# Set up dust_density.inp ##########################
################################################################################

print '\n########################################################################'
print '                           dust_density.inp                               '
print '########################################################################\n'

density = open('dust_density.inp', 'w')

# Information from the writing of amr_grid.inp is needed. Firstly, find the
# centre of the 3D grid
centre = [x[l/2], y[l/2], z[l/2]]

# Query for radius of the cloud
r = input('What radius cloud should be modelled? (This should be in AU)\n')

# Determine how many of the pixels lies within the cloud
r_pix = r/dpix

# Determine which points lie within the cloud



################################################################################
######################### Set up dustkappa_silicate.inp ########################
################################################################################

print '\n########################################################################'
print '                         dustkappa_silicate.inp                           '
print '########################################################################\n'

# This is the frequency dependent dust opacity defined in Hildebrand (1983)
# Parameterising the original developed by Hildebrand gives the below. This
# was found in The Physics and Chemistry of the ISM (AGGM Tielens) in
# wavelength form
#                       opac = kappa_0*(v/v_0)**B

############################# Create necessary files ###########################

# Creates the dust opacity file (dustopac.inp) and opens it in write mode
silicate = open('dustkappa_silicate.inp', 'w')

########################### Instantiate constant terms #########################

# The start and end points of the wavelength range in microns
lambda_init = 0.1 # 0.1 um
lambda_fin = 10000 # 10000 um

# Query for reference wavelength
w_ref = input('Which reference wavelength should I use to evaluate the opacities?\n')

# Ask for reference opacity
kappa_0 = input('What reference intensity should be used? (1 is a safe default)\n')

# Ask for the number of wavelength points
nlam = input('How many wavelength points do you want me to use in construction?\n')

##################### Evaluate opacities over wavelength range #################

# Create array of wavelengths from range given and reshape
w = np.vstack(np.linspace(lambda_init, lambda_fin, nlam))

# Ask for opacity law
opaclaw = input('Which opacity law should I use? Answer with: \'H\' (Hildebrand) \n')

# Ask for spectral index
B = input('What dust spectral index should be used? (1.7 is recommended)\n')

if opaclaw == 'H':
    # Evaluate opacities
    opacity = np.vstack(kappa_0*(w/w_ref)**B)

# Concantenate arrays to create 2D array of all data in format ready to be
# written to the .inp file
data = np.concatenate((w,opacity), axis=1)

########################## Write explanatory notes #########################

silicate.write('# This files contains all of the dust opacities for the power law\n')

if opaclaw == 'H':
    silicate.write('# given in Hildebrand (1983): kappa_abs = kappa_0(v/v_0)**B\n')

########################## Write the iformat integer #######################

# Writes the 'iformat' integer. This should always be 1 for lambda and
# kappa_abs columns to be written
silicate.write('1               # Format number of this file\n')

########################## Write the nlam integer ##########################

# Writes the nlam integer (i.e. how many wavelength points are in the file)
silicate.write(str(nlam) + '               # Nr of wavelength points in the file\n')

################################ Write data ################################

for row in data:
    for column in row:
        silicate.write('%14.8f' % column)
    silicate.write('\n')
silicate.close()

print '\'dustkappa_silicate.inp\' written to the working directory\n'
