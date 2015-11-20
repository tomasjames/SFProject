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
############################ Preliminary Information ###########################
################################################################################

################################ Define constants ##############################

# Define solar mass (in g) and solar radius (in cm)
m_sun = 1.988e33
r_sun = 6.955e10

# Define proton mass (in g/cm^3)
m_p = 1.67e-24

# Define mean molecular weight
u = 2.33

# cm to AU conversation factor
au_factor = 6.685e-14

######################## Ask user for intended quantities ######################

# Ask for the cloud mass and convert to grams
mass = input('What cloud mass should be used? (Answer should be in Solar Masses)\n')
mass = mass*(m_sun)
print 'The mass of the cloud considered is', mass, 'g.\n'

# Ask for the cloud number density and convert to g/cm^3
cloud_density = input('What number density should be assigned to the cloud?\n')
cloud_density = cloud_density*u*m_p
print 'The cloud density is', cloud_density, 'g/cm^3.\n'

# Repeat for outside density
outside_density = input('What number density should be assigned outside of the cloud?\n')
outside_density = outside_density*u*m_p
print 'The outside density is', outside_density, 'g/cm^3.\n'

# Calculate radius of the cloud in cm
r = ((3./(4*np.pi))*(mass/cloud_density))**(1./3.)
print 'The radius of the cloud is', r, 'cm (or', r*au_factor, ' AU).\n'

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

# Define the image width
width = input('What image width would you like to use? The radius of your cloud is ' +str(r*au_factor) +' AU - image width should be appropriately larger than this (This should be in AU) \n')

# Define the axis length (i.e. the number of pixels in each axis)
l = input('How many pixels should I assign to a dimension? (128 recommended) \n')

# Define the width of each pixel in AU
dpix = width/l
print 'This corresponds to a pixel width of', dpix,'AU.'

# Create an array of pixels in each dimension
x,y,z = np.linspace(0,width,l+1), np.linspace(0,width,l+1), np.linspace(0,width,l+1)
#x,y,z = np.arange(0,(l*dpix+1)), np.arange(0,(l*dpix+1)), np.arange(0,(l*dpix+1))

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
    if i == x[-1]:
        amr.write(str(i) + str('\n'))
    else:
        amr.write(str(i) + str(' '))

for j in y:
    if j == y[-1]:
        amr.write(str(j) + str('\n'))

    else:
        amr.write(str(j) + str(' '))

for k in z:
    if k == z[-1]:
        amr.write(str(k) + str('\n'))

    else:
        amr.write(str(k) + str(' '))

print '\'amr_grid.inp\' has been written to the working directory\n'

################################################################################
############# Set up dust_density.inp and dust_temperature.inp #################
################################################################################

print '\n########################################################################'
print '                dust_density.inp and dust_temperature.dat                 '
print '########################################################################\n'

# Writes the files to the working directory
density = open('dust_density.inp', 'w')
temperature = open('dust_temperature.dat', 'w')

# Writes the format number
density.write('1         # The same as most other files\n')
temperature.write('1         # The same as most other files\n')

# Write the number of cells
density.write(str(l**3)+'         # Number of cells\n')
temperature.write(str(l**3)+'         # Number of cells\n')

# Write the number of dust species
density.write('1         # Number of dust species\n')
temperature.write('1         # Number of dust species\n')

# Firstly, find the centre of the 3D grid
centre = [x[l/2], y[l/2], z[l/2]]

# Determine how many of the pixels lies within the cloud (convert the radius to
# AU and then divide by the width of one pixel in AU)
r_pix = (r*au_factor)/dpix

# Ask for the temperature to be assigned to the cloud
print 'Dust density in the cloud was supplied previously. It is taken to be'
print cloud_density, 'g/cm^3.\n'
cloud_temperature = input('What temperature should be assigned to the cloud?\n')

# Ask for the temperature outside
print 'Dust density outside of the cloud was supplied previously. It is taken to be'
print outside_density, 'g/cm^3.\n'
outside_temperature = input('What temperature should be assigned outside of the cloud?\n')

# Define empty lists to store densities
density_cube = np.zeros([l,l,l])

x_cube = np.zeros([l,l,l])
y_cube = np.zeros([l,l,l])
z_cube = np.zeros([l,l,l])

# Determine where centre of cloud lies in model space
for n in range(0,len(z)-1):
    for m in range(0,len(y)-1):
        for l in range(0,len(x)-1):
            if np.sqrt((x[l]-centre[0])**2 + (y[m]-centre[1])**2 + (z[n]-centre[2])**2) <= r*au_factor:
                if n == (len(z)-1):
                    density.write(str(cloud_density))
                    temperature.write(str(cloud_temperature))
                else:
                    density.write(str(cloud_density)+'\n')
                    temperature.write(str(cloud_temperature)+'\n')
                '''
                density_cube[l,m,n] = cloud_density
                x_cube[l,m,n] = x[l]
                y_cube[l,m,n] = y[m]
                z_cube[l,m,n] = z[n]
                '''
            else:
                if n == (len(z)-1):
                    density.write(str(outside_density))
                    temperature.write(str(outside_temperature))
                else:
                    density.write(str(outside_density)+'\n')
                    temperature.write(str(outside_temperature)+'\n')
                '''
                density_cube[l,m,n] = outside_density
                x_cube[l,m,n] = x[l]
                y_cube[l,m,n] = y[m]
                z_cube[l,m,n] = z[n]
                '''

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
