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
inp_data = np.loadtxt('../../curvefitting/sphdata/datafeed.txt')
chi_data = np.loadtxt('../../curvefitting/sphdata/chi.txt')

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
subplot2grid((6,12), (0,0), colspan=6,rowspan=6)
imshow(N_pix,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('$N$')
subplot2grid((6,12), (0,6), colspan=6,rowspan=6)
imshow(T_pix,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('$T$')
suptitle('Maps of the $\chi^{2}$ Recovered $N$ and $T$\n')
tight_layout()
savefig('map.png', dpi=300)
close()
