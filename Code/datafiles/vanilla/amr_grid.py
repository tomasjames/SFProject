################################################################################
############################ 4th Year Project Code #############################
########################## Custom AMR Grid Write File ##########################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import Numpy
import numpy as np

######################## Write file to working directory #######################

# Writes a file called amr_grid.inp to the working directory and allows it to be
# written to
f = open('amr_grid.inp', 'w')

print '########################################################################'
print 'This script will request some information before writing the returned'
print 'values to the file \'amr_grid.inp\' for RADMC-3D to deal with. At the'
print 'moment this code will only handle 3-dimensional grids.'
print '########################################################################\n'

############################ Define axes to work with ##########################

# Define the axis length (i.e. the number of pixels in each axis)
l = input('How many pixels should I assign to a dimension? \n')

# Define the width of each pixel
dpix = input('And what width would you like each pixel to be? \n')

# Create an array of pixels in each dimension
x,y,z = np.linspace(0,l*dpix,l+1), np.linspace(0,l*dpix,l+1), np.linspace(0,l*dpix,l+1)

###################### Write basic comments to explain file ####################

f.write('# This is a custom amr_grid.inp file for use with RADMC-3D. It is 3D.\n')

################# Write integers for RADMC-3D to learn about grid ##############

# Write the iformat number
f.write('1         # The same as most other files\n')

# Write the grid style
f.write('0         # 0 = regular\n')

# Write coordinate system
f.write('1         # < 100 = Cartesian\n')

# Write gridinfo integer
f.write('0         # If 1, abundant information is written regarding the grid\n')

# Write integers to activate dimensions
f.write('1 1 1 \n')

################# Loop through all grid values and write to file ###############

# Write the number of base cells in each dimension
f.write(str(l)+ ' ' + str(l) + ' ' +str(l) + '\n')

# Loop through each pixel and write to the file the edges of each pixel
for i in x:
    f.write(str(i) + str(' '))

    if i == x[-1]:
        f.write(str(i) + str('\n'))

for j in y:
    f.write(str(j) + str(' '))

    if j == y[-1]:
        f.write(str(j) + str('\n'))

for k in z:
    f.write(str(k) + str(' '))

    if k == z[-1]:
        f.write(str(k) + str('\n'))
