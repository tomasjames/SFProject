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

###################### Run the command to generate the SED #####################

# Plot SED (results are saved to spectrum.out)
os.system('radmc3d sed loadlambda')

################# Read the spectrum.out file and assign variables ##############

# Read in the spectrum
spectrum = np.loadtxt('spectrum.out', skiprows=2)

# Plot the spectrum
wav = spectrum[:,0]
flux = spectrum[:,1]

v = cc/(wav*10**-4)

########################## Read in the transmission data #######################

# Reads the transmission data
trans_data = np.loadtxt('transmission.txt')

# Assign to variables
trans_wav = trans_data[:,0]
trans_v = cc/(trans_wav*10**-4)
trans = trans_data[:,1]

# Sanity check: check to see if frequencies correspond
if trans_v.all() == v.all():
    print 'We are a GO to proceed!\n'
else:
    print 'ABORT ABORT ABORT\n'

# Weight the fluxes by multplying by the transmission coefficient at each point
weighted_flux = trans*flux

# Determine the flux 'seen' by SPIRE
weighted_flux_mean = np.mean(weighted_flux)
weighted_flux_std = np.std(weighted_flux)

################################# Save the data ################################

# Check to see if file already exists
if os.path.isfile('../../../stats/SPIRE_average_data.txt') == True:
    print 'Data storage file already exists; opening now\n'
    save_data = open('../../../stats/SPIRE_average_data.txt', 'a+')
    data_store = csv.writer(save_data, delimiter=' ')

else:
    print 'Data storage file does not already exist; writing now\n'
    save_data = open('../../../stats/SPIRE_average_data.txt', 'w+')
    data_store = csv.writer(save_data, delimiter=' ')

# Append to the file the passband
data_store.writerow([np.mean(v),np.float64(weighted_flux_mean),np.float64(weighted_flux_std)])

# Close the file
save_data.close()

################################# Plotting routine #############################

figure(1)
plot(v,flux,'b--')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data\n')
savefig('spectrum_psw_unweighted.png', dpi=300, bbox_inches='tight')
close()

figure(2)
plot(trans_v,weighted_flux,'g--')
errorbar(np.mean(v),weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data as seen by SPIRE\n')
savefig('spectrum_psw_weighted.png', dpi=300, bbox_inches='tight')
close()

figure(3)
errorbar(np.mean(v),weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
xlabel(r'$\nu(Hz)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data as seen by SPIRE\n')
savefig('spectrum_psw_seen.png', dpi=300, bbox_inches='tight')
close()
