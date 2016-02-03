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

################################# Plotting routine #############################

################################# Plotting routine #############################

figure(1)
plot(wav, flux, 'b--')
xlabel('Wavelength $(\mu m)$')
ylabel('Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data\n')
savefig('spectrum_psw.png', dpi=300, bbox_inches='tight')
close()

figure(2)
plot(wav, flux/max(flux), 'g--')
xlabel('Wavelength $(\mu m)$')
ylabel('Normalised Flux')
title('\nSpectral Energy Distribution for PSW Band Synthetic Data\n')
savefig('spectrum_norm_psw.png', dpi=300, bbox_inches='tight')
close()
