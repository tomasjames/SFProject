################################################################################
############################ 4th Year Project Code #############################
################################# Plot Maps ####################################
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

# Read in both data types
inp_data = np.loadtxt('../../curvefitting/cloud/datafeed.txt')
chi_data = np.loadtxt('../../curvefitting/cloud/chi.txt')

# Split data types into plottable quantities
inp_N = inp_data[:,2]
inp_T = inp_data[:,3]

chi_N = chi_data[:,1]
chi_T = chi_data[:,2]

################################# Plot image data ################################

# Reshape the data such that x and y pixels correspond
N_pix = chi_N.reshape(np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))
T_pix = chi_T.reshape(np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))

# Plot the 2 images side by side for comparison
figure(1)
imshow(N_pix,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('A Map of the $\chi^{2}$ Recovered $N$\n')
savefig('map_N.png', dpi=300)
close()

figure(2)
imshow(T_pix,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('A Map of the $\chi^{2}$ Recovered $T$\n')
tight_layout()
savefig('map_T.png', dpi=300)
close()
