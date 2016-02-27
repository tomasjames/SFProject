################################################################################
############################ 4th Year Project Code #############################
###################### SED Plotting and Manipulation Script ####################
################################################################################

############################# Import statements ################################
# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Import time tools
import time

# Import units
from radmc3dPy.natconst import *

# Import numpy
import numpy as np

# Import csv
import csv

from matplotlib.pyplot import *

############################### Generate the SED ###############################

# Read in the transmission weighted image devised earlier
image_trans_raw = np.loadtxt('image_trans_raw.txt')

# Read in the original image file
imag = radmc3dPy.image.readImage('image.out')
#imag = np.reshape(imag.image, (1,len(imag.image)*len(imag.image[0])))

# Resize the transmission weighted data
image_trans_raw = np.reshape(image_trans_raw, (len(imag.image),len(imag.image[0]),len(imag.image[0][0])))

print 'Data dimensions dictate that there are', imag.nx*imag.ny, 'rows before a new wavelength data entry begins.\n This means that there are', imag.nwav, 'wavelength fluxes at the central pixel.\n'
'''
# Reshape the raw transmission based image based on the dimensions in imag ie zeroth row is every x pixel value at y=0, first row is every x pixel at y=1 and so on)
image = np.reshape(imag.image, (imag.nwav,imag.nx*imag.ny))
image_trans_raw = np.reshape(image_trans_raw, (imag.nwav,imag.nx*imag.ny))

print 'Image data has been restructured such that there are now', imag.nwav, 'rows with', imag.nx*imag.ny, 'columns in each row.\n'

# Loop through the reshaped array to find the value of the middle pixel at every wavelength
middle_flux, middle_flux_trans = [], []

# Because the middle pixel is (imag.nx-1)/2, loop through every y value and find the value of the middle pixel and append. Should end up with ~100 values, 1 for each wavelength
for i in range(0,imag.nwav):
    if image[i][np.round((imag.nx*imag.ny)/2)+imag.nx/2] == 0:
        middle_flux.append(np.NaN)
    else:
        middle_flux.append(image[i][np.round((imag.nx*imag.ny)/2)+imag.nx/2])
    if image_trans_raw[i][np.round((imag.nx*imag.ny)/2)+imag.nx/2] != 0:
        middle_flux_trans.append(image_trans_raw[i][np.round((imag.nx*imag.ny)/2)+imag.nx/2])
    else:
        middle_flux_trans.append(np.NaN)
'''

sum = 0
middle_flux, middle_flux_trans, middle_flux_trans_index  = [], [], []

for l in range(0,imag.nwav):
    for y in range(0,imag.ny):
        for x in range(0,imag.nx):
            sum += 1
            if x == imag.nx/2 and y == imag.ny/2:
                middle_flux_trans.append(image_trans_raw[x][y][l])
                middle_flux.append(imag.image[x][y][l])
                middle_flux_trans_index.append(sum)

if sum == len(image_trans_raw):
    print 'The loop has been executed over all of the elements.\n'

################# Read the spectrum.out file and assign variables ##############

# Determine parameters
wav = imag.wav
w_cen = 496.106774
flux = middle_flux_trans

v = cc/(wav*10**-4)
v_cen = cc/(w_cen*10**-4)

########################## Read in the transmission data #######################

# Reads the transmission data
trans_data = np.loadtxt('transmission.txt')

# Assign to variables
trans_wav = trans_data[:,0]
trans_v = cc/(trans_wav*10**-4)
trans = trans_data[:,1]

# Determine the flux 'seen' by SPIRE
weighted_flux_mean = np.sum(middle_flux_trans)/np.sum(trans)
weighted_flux_std = np.std(middle_flux_trans)

'''
# Sanity check: check to see if frequencies correspond
if trans_v.all() == v.all():
    print 'We are a GO to proceed!\n'
else:
    print 'ABORT ABORT ABORT\n'

# Weight the fluxes by multplying by the transmission coefficient at each point
weighted_flux = trans*flux

# Determine the flux 'seen' by SPIRE
#weighted_flux_mean = np.mean(weighted_flux)
weighted_flux_std = np.std(weighted_flux)

weighted_flux_mean = np.sum(weighted_flux)/np.sum(trans)
'''
################################# Save the data ################################

# Check to see if file already exists
if os.path.isfile('../../../curvefitting/SPIRE_average_data.txt') == True:
    print 'Data storage file already exists; opening now\n'
    save_data = open('../../../curvefitting/SPIRE_average_data.txt', 'a+')
    data_store = csv.writer(save_data, delimiter=' ')

else:
    print 'Data storage file does not already exist; writing now\n'
    save_data = open('../../../curvefitting/SPIRE_average_data.txt', 'w+')
    data_store = csv.writer(save_data, delimiter=' ')

# Append to the file the passband
data_store.writerow([v_cen,np.float64(weighted_flux_mean),np.float64(weighted_flux_std)])

# Close the file
save_data.close()

################################# Plotting routine #############################

subplot2grid((12,12), (0,0), colspan=6,rowspan=12)
plot(v,middle_flux,'r--',label='SED from image.out')
plot(v,middle_flux_trans,'b--',label='SED from image_trans.out')
errorbar(v_cen,weighted_flux_mean,yerr=weighted_flux_std,fmt='go',label='Flux \'seen\' by SPIRE')
xlabel(r'$\nu(Hz)$')
ylabel('Intensity $(erg/s/cm^{2}/Hz/sr)$')
legend(loc='best',prop={'size':9})
subplot2grid((12,12), (0,6), colspan=6,rowspan=12)
plot(v,trans,'g--',label='Transmission Curve')
xlabel(r'$\nu(Hz)$')
ylabel('Transmission')
suptitle('Spectral Energy Distribution for PLW Band Synthetic Data\n')
legend(loc='best',prop={'size':9})
tight_layout()
savefig('spectrum_plw_unweighted_v2.png', dpi=300, bbox_inches='tight')
close()

'''
figure(2)
plot(trans_v,weighted_flux,'g--')
errorbar(v_cen,weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PLW Band Synthetic Data as seen by SPIRE\n')
savefig('spectrum_plw_weighted_v2.png', dpi=300, bbox_inches='tight')
close()

figure(3)
errorbar(v_cen,weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PLW Band Synthetic Data as seen by SPIRE\n')
savefig('spectrum_plw_seen_v2.png', dpi=300, bbox_inches='tight')
close()
'''
