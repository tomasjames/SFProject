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
inp_data = np.loadtxt('../../curvefitting/filament/datafeed.txt')
chi_data = np.loadtxt('../../curvefitting/filament/chi.txt',skiprows=1)

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

# Plot the data along with a PDF
figure()
subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
imshow(np.log10(N_chi_inp),origin='lower')
colorbar(label='$log_{10}N\/(g\/cm^{-3})$')
xlabel('X')
ylabel('Y')
title('A Map of the $\chi^{2}$ Recovered $N$\n')

subplot2grid((6,6), (5,0), colspan=4,rowspan=1)
hist(np.log10(N_chi_inp),bins=50,log=True,normed=True)
title('PDF of $\chi^{2}$ Recovered $N$')
xlabel('$log_{10}N\/(g\/cm^{-3})$')
ylabel('Frequency')

savefig('map_N_chi.png', dpi=300)
close()

# Plot the data along with a PDF
figure()
subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
imshow(T_chi_inp,origin='lower')
colorbar(label='$T\/(K)$')
xlabel('X')
ylabel('Y')
title('A Map of the $\chi^{2}$ Recovered $T$\n')

subplot2grid((6,6), (5,0), colspan=4,rowspan=1)
hist(T_chi_inp,bins=50,log=True,normed=True)
title('PDF of $\chi^{2}$ Recovered $T$')
xlabel('$T\/(K)$')
ylabel('Frequency')

savefig('map_T_chi.png', dpi=300)
close()

# Plot the data along with a PDF
figure()
N_data_inp[N_data_inp == 0] = np.nan
subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
imshow(np.log10(N_data_inp),origin='lower')
colorbar(label='$log_{10}N\/(g\/cm^{-3})$')
xlabel('X')
ylabel('Y')
title('A Map of the Data Input $N$\n')

subplot2grid((6,6), (5,0), colspan=4,rowspan=1)
hist(np.log10(N_data_inp),bins=50,log=True,normed=True)
title('PDF of Data Input $N$')
xlabel('$log_{10}N\/(g\/cm^{-3})$')
ylabel('Frequency')

savefig('map_N_data.png', dpi=300)
close()

# Plot the data along with a PDF
figure()
T_data_inp[T_data_inp == 0] = np.nan
subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
imshow(T_data_inp,origin='lower')
colorbar(label='$T\/(K)$')
xlabel('X')
ylabel('Y')
title('A Map of the Data Input $T$\n')

subplot2grid((6,6), (5,0), colspan=4,rowspan=1)
hist(T_data_inp,bins=50,log=True,normed=True)
title('PDF of Data Input $T$')
xlabel('$T\/(K)$')
ylabel('Frequency')

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
