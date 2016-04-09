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

from astropy.io import fits

############################### Read in image data ##############################

# Read in both data types
inp_data = np.loadtxt('../../curvefitting/sphdata/datafeed.txt')
chi_data = np.loadtxt('../../curvefitting/sphdata/chi.txt',skiprows=1)

# Split data types into plottable quantities
inp_N = inp_data[:,1]
inp_T = inp_data[:,2]

chi_N = chi_data[:,1]
chi_T = chi_data[:,2]

################################# Plot image data ################################

# Reshape the data such that x and y pixels correspond
N_chi_inp = np.reshape(chi_N, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))
T_chi_inp = np.reshape(chi_T, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))

N_data_inp = np.reshape(inp_N, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))
T_data_inp = np.reshape(inp_T, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))

# Plot the 2 images side by side for comparison
figure(1)
imshow(N_chi_inp,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('A Map of the $\chi^{2}$ Recovered $N$\n')
savefig('map_N_chi.png', dpi=300)
close()

figure(2)
imshow(T_chi_inp,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('A Map of the $\chi^{2}$ Recovered $T$\n')
tight_layout()
savefig('map_T_chi.png', dpi=300)
close()

figure(3)
imshow(N_data_inp,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('A Map of the $N$ from the User Defined Values\n')
tight_layout()
savefig('map_N_data.png', dpi=300)
close()

figure(4)
imshow(T_data_inp,origin='lower')
colorbar()
xlabel('X')
ylabel('Y')
title('A Map of the $T$ from the User Defined Values\n')
tight_layout()
savefig('map_T_data.png', dpi=300)
close()

################################# Save image data ################################

# Collate all of the data into one array
combined = [N_chi_inp, T_chi_inp, N_data_inp, T_data_inp]

# Define the names of the individual data sets
combined_names = ['N_chi_inp', 'T_chi_inp', 'N_data_inp', 'T_data_inp']

for i in range(0,len(combined)):
    # Define a new header file
    hdul = fits.HDUList()

    # Append to a primary header
    hdul.append(fits.PrimaryHDU())

    # Append the data
    hdul.append(fits.ImageHDU(data=combined[i]))

    # Write the data to a fits file with name of the data array
    hdul.writeto(str(combined_names[i])+str('.fits'))
