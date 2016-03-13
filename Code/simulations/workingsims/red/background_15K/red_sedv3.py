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

############################### Read in image data ##############################

# Read in the transmission weighted image devised earlier
image_trans_raw = np.loadtxt('image_trans_raw.txt')

# Read in the original image file
imag = radmc3dPy.image.readImage('image.out')
#imag = np.reshape(imag.image, (1,len(imag.image)*len(imag.image[0])))

# Resize the transmission weighted data
image_trans_raw = np.reshape(image_trans_raw, (len(imag.image),len(imag.image[0]),len(imag.image[0][0])))

print 'Data dimensions dictate that there are', imag.nx*imag.ny, 'rows before a new wavelength data entry begins.\n This means that there are', imag.nwav, 'wavelength fluxes at the central pixel.\n'

########################## Read in the transmission data #######################

# Reads the transmission data
trans_data = np.loadtxt('transmission.txt')

# Assign to variables
trans_wav = trans_data[:,0]
trans_v = cc/(trans_wav*10**-4)
trans = trans_data[:,1]

######################### Determine frequency of the image #####################

# Determine parameters
wav = imag.wav
w_cen = 166.069819
#flux = middle_flux_trans

v = cc/(wav*10**-4)
v_cen = cc/(w_cen*10**-4)

############################ Sample flux at each pixel #########################

# Check to see if file already exists
if os.path.isfile('../../../curvefitting/red_average_data.txt') == True:
    print 'Data storage file already exists; opening now\n'
    save_data = open('../../../curvefitting/red_average_data.txt', 'a+')
    data_store = csv.writer(save_data, delimiter=' ')

else:
    print 'Data storage file does not already exist; writing now\n'
    save_data = open('../../../curvefitting/red_average_data.txt', 'w+')
    data_store = csv.writer(save_data, delimiter=' ')

# Sum integer to track the loop
sum = 0

# Instantiate lists to store values
flux, flux_trans, flux_trans_index = [], [], []

# Loop through the array to find the middle pixel and log its flux
for x in range(0,imag.nx):
    for y in range(0,imag.ny):
        for l in range(0,imag.nwav):
            sum += 1

            # Assign values to the pixel numbers to track position in loop
            xpix = x
            ypix = y
            zpix = l

            raw_flux = imag.image[x][y][l]
            raw_flux_trans = image_trans_raw[x][y][l]

            flux.append(raw_flux)
            flux_trans.append(raw_flux_trans)

        # Determine the flux 'seen' by SPIRE
        weighted_flux_mean = np.sum(flux_trans)/np.sum(trans)
        weighted_flux_std = np.std(flux_trans)

        # Append th determined information to the file
        data_store.writerow([sum, xpix, ypix, zpix, v_cen, weighted_flux_mean, weighted_flux_std])

        # Reset the lists to 0
        flux, flux_trans, flux_trans_index = [], [], []

if sum == imag.nx*imag.ny*imag.nwav:
    print 'The loop has been executed over all of the elements.\n'

# Close the file
save_data.close()

################################# Plotting routine #############################

'''
subplot2grid((12,12), (0,0), colspan=6,rowspan=12)
plot(v,middle_flux,'r--',label='SED from image.out')
plot(v,middle_flux_trans,'b--',label='SED from image_trans.out')
errorbar(v_cen,weighted_flux_mean,yerr=weighted_flux_std,fmt='go',label='Flux \'seen\' by PACS')
xlabel(r'$\nu(Hz)$')
ylabel('Intensity $(erg/s/cm^{2}/Hz/sr)$')
legend(loc='best',prop={'size':9})
subplot2grid((12,12), (0,6), colspan=6,rowspan=12)
plot(v,trans,'g--',label='Transmission Curve')
xlabel(r'$\nu(Hz)$')
ylabel('Transmission')
suptitle('Spectral Energy Distribution for Red Band Synthetic Data\n')
legend(loc='best',prop={'size':9})
tight_layout()
savefig('spectrum_red_unweighted_v2.png', dpi=300, bbox_inches='tight')
close()

figure(2)
plot(trans_v,weighted_flux,'g--')
errorbar(v_cen,weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data as seen by SPIRE\n')
savefig('spectrum_psw_weighted_v2.png', dpi=300, bbox_inches='tight')
close()

figure(3)
errorbar(v_cen,weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data as seen by SPIRE\n')
savefig('spectrum_psw_seen_v2.png', dpi=300, bbox_inches='tight')
close()
'''
