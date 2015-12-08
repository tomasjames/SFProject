################################################################################
############################ 4th Year Project Code #############################
########################## Custom Data File Generation #########################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import Numpy
import numpy as np

# Import constants from radmc3dPy
from radmc3dPy.natconst import *

# Import glob for dredging folders for files
import glob

# Import matplotlib to plot
import matplotlib.pyplot as plt

# Import interpolation modules
from scipy.interpolate import interp1d

print '\n########################################################################'
print 'This script will request some information before writing multiple required'
print 'custom files to the working directory for use with RADMC-3D. These should'
print 'be copied and pasted to the directory in which RADMC-3D\'s thermal Monte'
print 'Carlo simulation or raytrace is run.'
print '########################################################################\n'

################################################################################
############################ Preliminary Information ###########################
################################################################################


######################## Ask user for intended quantities ######################

# Define dust-to-gas ratio
d2g = 0.01

# Ask for the cloud mass and convert to grams
m = input('What cloud mass should be used? (Answer should be in Solar Masses)\n')
dust_mass = d2g**m*ms
print 'The mass of the cloud considered is', m*ms, 'g. The dust mass is ', dust_mass, 'g.\n'

# Ask for the cloud number density and convert to g/cm^3
cloud = input('What number density should be assigned to the cloud?\n')
cloud_density = cloud*muh2*mp*d2g
print 'The cloud density is', cloud_density, 'g/cm^3.\n'

# Repeat for outside density
outside = np.float64(input('What number density should be assigned outside of the cloud?\n'))
outside_density = outside*muh2*mp*d2g
print 'The outside density is', outside_density, 'g/cm^3.\n'

# Calculate radius of the cloud in cm
r = ((3./4)*(1./np.pi)*(dust_mass/cloud_density))**(1./3)
print 'The radius of the cloud is', r, 'cm (or', r/au, ' AU).\n'

#print 'The dust mass in the cloud is ', ((4./3.)*(np.pi)*(r**3)*(cloud_density))/ms, 'MSun.'

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
width = input('What image width would you like to use? The radius of your cloud is ' +str(r/au) +' AU - image width should be appropriately larger than this (This should be in AU) \n')
print 'This is a width of ', width*au, 'cm.\n'

# Define the axis length (i.e. the number of pixels in each axis)
l = input('How many pixels should I assign to a dimension? (128 recommended) \n')

# Define the width of each pixel in AU
dpix = (width*au)/l
print 'This corresponds to a pixel width of', dpix,'cm.'

# Create an array of pixels in each dimension
x,y,z = np.linspace(-(width*au)/2,(width*au)/2,l+1), np.linspace(-(width*au)/2,(width*au)/2,l+1), np.linspace(-(width*au)/2,(width*au)/2,l+1)
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
        amr.write(str(np.float64(i)) + str('\n'))
    else:
        amr.write(str(np.float64(i)) + str(' '))

for j in y:
    if j == y[-1]:
        amr.write(str(np.float64(j)) + str('\n'))

    else:
        amr.write(str(np.float64(j)) + str(' '))

for k in z:
    if k == z[-1]:
        amr.write(str(np.float64(k)) + str('\n'))

    else:
        amr.write(str(np.float64(k)) + str(' '))

print '\'amr_grid.inp\' has been written to the working directory\n'

amr.close()

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
#centre = [x[l/2], y[l/2], z[l/2]]
centre = [0,0,0]

# Determine how many of the pixels lies within the cloud (convert the radius to
# cm and then divide by the width of one pixel in cm)
r_pix = (r*au)/dpix

# Ask for the temperature to be assigned to the cloud
print 'Dust density in the cloud was supplied previously. It is taken to be'
print cloud_density, 'g/cm^3.\n'
cloud_temperature = np.float64(input('What temperature should be assigned to the cloud?\n'))

# Ask for the temperature outside
print 'Dust density outside of the cloud was supplied previously. It is taken to be'
print outside_density, 'g/cm^3.\n'
outside_temperature = np.float64(input('What temperature should be assigned outside of the cloud?\n'))

# Define empty lists to store densities
density_cube = np.zeros([l,l,l])

x_cube = np.zeros([l,l,l])
y_cube = np.zeros([l,l,l])
z_cube = np.zeros([l,l,l])

# Determine where centre of cloud lies in model space
for n in range(0,len(z)-1):
    for m in range(0,len(y)-1):
        for l in range(0,len(x)-1):
            if (x[l]-centre[0])**2 + (y[m]-centre[1])**2 + (z[n]-centre[2])**2 <= r**2:
                if n == (len(z)-1):
                    density.write(str(cloud_density))
                    temperature.write(str(cloud_temperature))
                else:
                    density.write(str(cloud_density)+'\n')
                    temperature.write(str(cloud_temperature)+'\n')

                density_cube[l,m,n] = cloud_density
                x_cube[l,m,n] = x[l]
                y_cube[l,m,n] = y[m]
                z_cube[l,m,n] = z[n]

            else:
                if n == (len(z)-1):
                    density.write(str(outside_density))
                    temperature.write(str(outside_temperature))
                else:
                    density.write(str(outside_density)+'\n')
                    temperature.write(str(outside_temperature)+'\n')

                density_cube[l,m,n] = outside_density
                x_cube[l,m,n] = x[l]
                y_cube[l,m,n] = y[m]
                z_cube[l,m,n] = z[n]

density.close()
temperature.close()

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
lambda_fin = 10000. # 10000 um

# Query for reference wavelength
w_ref = input('Which reference wavelength should I use to evaluate the opacities?\n')

# Ask for reference opacity
kappa_0 = input('What reference intensity should be used? (1 is a safe default)\n')

# Ask for the number of wavelength points
nlam = input('How many wavelength points do you want me to use in construction?\n')

##################### Evaluate opacities over wavelength range #################

# Create array of wavelengths from range given and reshape
w = np.linspace(lambda_init, lambda_fin, nlam)

# Ask for opacity law
opaclaw = input('Which opacity law should I use? Answer with: \'H\' (Hildebrand) \n')

# Ask for spectral index
B = input('What dust spectral index should be used? (1.7 is recommended)\n')

if opaclaw == 'H':
    # Evaluate opacities
    opacity = kappa_0*(w/w_ref)**B

# Concantenate arrays to create 2D array of all data in format ready to be
# written to the .inp file
#data = np.concatenate((w,opacity), axis=1)

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

for o in range(0,len(w)):
    if i == len(w):
        silicate.write(str(w[o])+str('    ')+str(opacity[o]))
    else:
        silicate.write(str(w[o])+str('    ')+str(opacity[o])+str('\n'))

print '\'dustkappa_silicate.inp\' written to the working directory\n'

silicate.close()

################################################################################
######################### Set up wavelength_micron.inp #########################
################################################################################

print '\n########################################################################'
print '       wavelength_micron.inp and camera_wavelength_micron.inp       '
print '########################################################################\n'

# Write the file
wavelength_micron = open('wavelength_micron.inp', 'w')

# Scan the directory to find file containing .txt
# This works by matching any filename that has PSW before the file extension
files = glob.glob('./*Herschel*')

# Define empty list to store file contents
spire, spire_wav, spire_trans = [], [], []

# Determine whether
if files == ['./Herschel_SPIRE.PSW_ext.dat.txt']:
    with open('Herschel_SPIRE.PSW_ext.dat.txt') as f:
        spire_psw = f.readlines()
        print 'I am reading from \'Herschel_SPIRE.PSW_ext.dat.txt\' and have assigned its output to psw_band.\n'

        # Because the content is not read in using proper formatting this splits each row into columns that allow the first column to be wavelength points for the passband and the second column to be the transmission coefficint at that wavelength
        for row in spire_psw:
            spire.append(row.split())

elif files == ['./Herschel_SPIRE.PMW_ext.dat.txt']:
    with open('Herschel_SPIRE.PMW_ext.dat.txt') as f:
        spire_pmw = f.readlines()
        print 'I am reading from \'Herschel_SPIRE.PMW_ext.dat.txt\' and have assigned its output to pmw_band.\n'

        # Because the content is not read in using proper formatting this splits each row into columns that allow the first column to be wavelength points for the passband and the second column to be the transmission coefficint at that wavelength
        for row in spire_pmw:
            spire.append(row.split())

elif files == ['./Herschel_SPIRE.PLW_ext.dat.txt']:
    with open('Herschel_SPIRE.PLW_ext.dat.txt') as f:
        spire_plw = f.readlines()
        print 'I am reading from \'Herschel_SPIRE.PLW_ext.dat.txt\' and have assigned its output to plw_band.\n'

        # Because the content is not read in using proper formatting this splits each row into columns that allow the first column to be wavelength points for the passband and the second column to be the transmission coefficint at that wavelength
        for row in spire_plw:
            spire.append(row.split())

for a in range(0,len(spire)):
    spire_wav.append(np.float64(spire[a][0]))
    spire_trans.append(np.float64(spire[a][1]))

print 'The correct format data has now been assigned to spire with the wavelength component being in spire_wav and the transmission component being in spire_trans.\n'

# Determine the lower and upper bounds of the passband
lower = np.float64(spire[0][0])
higher = np.float64(spire[-1][0])

# Determine the number of wavelength points in the passband
nwav = len(spire)

# Save the number of wavelength points
wavelength_micron.write(str(nwav)+str('\n'))

# Instantiate the number of wavelength points`
#wavs = np.linspace(np.log10(1e-1), np.log10(1e3), nwav)
wavs = np.linspace(1e-1, 1e3, nwav)

# Takes each element of the list and writes it to the wavelength_micron file
for p in wavs:
    if p == wavs[-1]:
        #wavelength_micron.write(str(10**(p)))
        wavelength_micron.write(str(p))
    else:
        #wavelength_micron.write(str(10**(p)) + str('\n'))
        wavelength_micron.write(str(p)+str('\n'))

# Close the file
wavelength_micron.close()

# Write the file
camera = open('camera_wavelength_micron.inp', 'w')

print 'There will be', nwav, 'points between', lower, 'um and', higher, 'um.\n'
print ''

# This determines whether the spacing is linear or logarithmic
if spire_wav[2] - spire_wav[1] != spire_wav[1] - spire_wav[0]:
    print 'Spacing in wavelength is not constant; interpolation is required.\n'

    # Generate a function that best fits the data given in the download
    interp = interp1d(spire_wav, spire_trans, kind='cubic')

    # Generate a linear array of points to feed into it
    linear_wav = np.linspace(lower,higher,nwav)

    # Feed in and determine the new, linearly spaced transmission points
    interp_trans = interp(linear_wav)

    plt.plot(spire_wav, spire_trans, 'r+', label=str(files))
    plt.plot(linear_wav, interp_trans, 'b--', label='Interpolated Points')
    plt.xlabel('$\lambda$ ($\AA$)')
    plt.ylabel('Transmission')
    plt.legend(loc='best')
    plt.savefig('interpolated.png', bbox_inches='tight')
    plt.close()

    print 'Interpolation complete. Plot of the original function and the interpolated function have been saved to the working directory for accuracy inspection.\n'
    print 'Data will now be written to \'camera_wavelength_micron.inp\'\n'

    # Save the number of wavelength points
    camera.write(str(nwav)+str('\n'))

    # Writes wavelength points
    for q in range(0,len(linear_wav)):
        if q == len(linear_wav):
            camera.write(str(linear_wav[q]*10**-4))
        else:
            camera.write(str(linear_wav[q]*10**-4) + str('\n'))

else:
    # Save the number of wavelength points
    camera.write(str(nwav)+str('\n'))

    # Writes wavelength points
    for q in range(0,len(spire)):
        if q == len(spire):
            camera.write(str(spire[q][0]*10**-4))
        else:
            camera.write(str(spire[q][0]*10**-4) + str('\n'))

camera.close()

'''
################################################################################
############################ Set up external_source.inp ########################
################################################################################

print '\n########################################################################'
print '                           external_source.inp                            '
print '########################################################################\n'

# Write the file
external = open('external_source.inp', 'w')

# Write the format number
external.write('2        # This is the format number (1=Hertz i.e. frequency)\n')

# Write number of wavelength points
external.write(str(nwav)+str('\n'))

isrf = []

# Loop through the wavelengths as defined earlier for wavelength_micron.inp and save to external_source.inp as well as evaluate the intensity (using the Planck function)
for r in wavs:
    external.write(str(10**(r)) + str('\n'))

    # Split Planck function into seperate quotients to make error tracking easier
    a = 2.0*hh*cc**2
    #b = hh*cc/(((10**r)*10**-4)*kk*outside_temperature)
    #isrf.append(a/((((10**r)*10**-4)**5)*(np.exp(b) - 1.0)))
    b = hh*cc/(((10**r))*kk*outside_temperature)
    isrf.append(a/((((10**r))**5)*(np.exp(b) - 1.0)))

# Write the intensities to the end of external_source.inp
for s in isrf:
    if s == isrf[-1]:
        external.write(str(s))
    else:
        external.write(str(s)+str('\n'))

external.close()

################################################################################
############################### Set up stars.inp ###############################
################################################################################

print '\n########################################################################'
print '                               stars.inp                                  '
print '########################################################################\n'

# Determine whether the user would like stars
userStar = input('Do you want to put any stars in the simulation?')

if userStar == 'Yes':

    # Write the file
    stellar = open('stars.inp', 'w')

    # Write the format number
    stellar.write('2        # Format number: should always be 2')

    # Write the number of stars (assume 1 for now)
    stellar.write('1        # Number of stars in the field')

    # Write the radius of the star along with its mass and coordinates (assume Sun)
    stellar.write('1    1    0    0    0')

else:
    print 'No stars will be considered (and stars.inp has not been written)\n'
'''
