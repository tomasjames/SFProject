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
chi_data = np.loadtxt('../../curvefitting/sphdata/chi_fine.txt',skiprows=1)

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
subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
imshow(np.log10(N_chi_inp),origin='lower',vmin=np.log10(min(chi_N)),vmax=np.log10(max(inp_N)))
colorbar(label='$log_{10}N\/(g\/cm^{-2})$')
xlabel('X (pixels)')
ylabel('Y (pixels)')
title('A Map of the $\chi^{2}$ Recovered $N$\n')

subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
#hist(np.log10(N_chi_inp),bins=5,normed=True)
h, b = np.histogram(np.log10(N_chi_inp), bins=len(N_chi_inp), density=True)
#bar(b[:-1],h,width=(max(h)-min(h))/len(N_chi_inp))
bar(b[:-1],h,width=0)
xlim(np.log10(min(chi_N)),np.log10(max(inp_N)))
title('PDF of $\chi^{2}$ Recovered $N$')
xlabel('$log_{10}N\/(g\/cm^{-2})$')
ylabel('Normalised \n Frequency')

savefig('map_N_chi.png', dpi=300)
close()

# Plot the data along with a PDF
figure()
subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
imshow(T_chi_inp,origin='lower',vmin=min(inp_T),vmax=max(chi_T))
colorbar(label='$T\/(K)$')
xlabel('X (pixels)')
ylabel('Y (pixels)')
title('A Map of the $\chi^{2}$ Recovered $T$\n')

subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
#hist(T_chi_inp,bins=5,normed=True)
h, b = np.histogram(T_chi_inp, bins=len(T_chi_inp), density=True)
#bar(b[:-1],h,width=(max(h)-min(h))/len(T_chi_inp))
bar(b[:-1],h,width=0)
xlim(min(inp_T),max(inp_T))
title('PDF of $\chi^{2}$ Recovered $T$')
xlabel('$T\/(K)$')
ylabel('Normalised \n Frequency')

savefig('map_T_chi.png', dpi=300)
close()

# Plot the data along with a PDF
figure()
N_data_inp[N_data_inp == 0] = np.nan
subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
imshow(np.log10(N_data_inp),origin='lower',vmin=np.log10(min(chi_N)),vmax=np.log10(max(inp_N)))
colorbar(label='$log_{10}N\/(g\/cm^{-2})$')
xlabel('X (pixels)')
ylabel('Y (pixels)')
title('A Map of the Data Input $N$\n')

subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
#hist(np.log10(N_data_inp),bins=5,normed=True)
h, b = np.histogram(np.log10(N_data_inp), bins=len(N_data_inp), density=True)
#bar(b[:-1],h,width=(max(h)-min(h))/len(N_data_inp))
bar(b[:-1],h,width=0)
xlim(np.log10(min(chi_N)),np.log10(max(inp_N)))
title('PDF of Data Input $N$')
xlabel('$log_{10}N\/(g\/cm^{-2})$')
ylabel('Normalised \n Frequency')

savefig('map_N_data.png', dpi=300)
close()

# Plot the data along with a PDF
figure()
T_data_inp[T_data_inp == 0] = np.nan
subplot2grid((6,7), (0,0), colspan=4,rowspan=4)
imshow(T_data_inp,origin='lower',vmin=min(inp_T),vmax=max(chi_T))
colorbar(label='$T\/(K)$')
xlabel('X (pixels)')
ylabel('Y (pixels)')
title('A Map of the Data Input $T$\n')

subplot2grid((6,7), (5,0), colspan=4,rowspan=2)
#hist(T_data_inp,bins=5,normed=True)
h, b = np.histogram(T_data_inp, bins=len(T_data_inp), density=True)
#bar(b[:-1],h,width=(max(h)-min(h))/len(T_data_inp))
bar(b[:-1],h,width=0)
xlim(min(inp_T),max(inp_T))
title('PDF of Data Input $T$')
xlabel('$T\/(K)$')
ylabel('Normalised \n Frequency')

savefig('map_T_data.png', dpi=300)
close()

##################### Determine line of sight T variations #######################

# This is essentially the standard deviated weighted mean of the temperatures
sigma_T_inp = (((inp_T - chi_T)**2)*inp_N)/(inp_N)
sigma_T_chi = (((inp_T - chi_T)**2)*chi_N)/(chi_N)

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

################################# Plot contours ################################

figure(1)
plot(chi_N, chi_T, 'g.', label=r'$\chi^{2}$ Recovered Values')
xlabel(r'$N_{\chi^{2}}$')
ylabel(r'$T_{\chi^{2}}$')
title(r'A Plot Showing the $\chi^{2}$ Recovered Values of $T$ and $N$')
legend(loc='best')
savefig('contours.png')
close()

######################## Plot further graphs for report #########################

# Plot the recovered temperatures against the input temperatures
figure(2)
plot(chi_T, inp_T, 'b.', markersize=2)
xlabel(r'$T_{\chi^{2}}\/(K)$')
ylabel(r'$T_{N}\/(K)$')
xlim(min(inp_T),max(chi_T))
ylim(min(inp_T),max(chi_T))
title(r'A Graph Comparing the $\chi^{2}$ Recovered $T$ to the data input $T$')
savefig('T.png', dpi=300)
close()

figure(3)
subplot(2,1,1)
plot(np.log10(inp_N), sigma_T_inp, 'g.', markersize=2)
xlabel(r'$N_{input}\/(g\/cm^{-2})$')
ylabel(r'$\sigma^{2}_{T_{N}}$')
title('A Graph Comparing the Line of Sight \n Temperature Variations to the Input Column Density')

subplot(2,1,2)
plot(inp_T, sigma_T_inp, 'g.', markersize=2)
xlabel(r'$T_{input}\/(K)$')
ylabel(r'$\sigma^{2}_{T_{N}}$')
title('A Graph Comparing the Line of Sight \n Temperature Variations to the Input Temperature')
tight_layout()
savefig('sigma_T_inp.png', dpi=300)
close()

figure(4)
subplot(2,1,1)
plot(np.log10(chi_N), sigma_T_chi, 'g.', markersize=2)
xlabel(r'$N_{\chi^{2}}\/(g\/cm^{-2})$')
ylabel(r'$\sigma^{2}_{T_{\chi^{2}}}$')
title('A Graph Comparing the Line of Sight \n Temperature Variations to the $\chi^{2}$ Column Density')

subplot(2,1,2)
plot(chi_T, sigma_T_chi, 'g.', markersize=2)
xlabel(r'$T_{\chi^{2}}\/(K)$')
ylabel(r'$\sigma^{2}_{T_{\chi^{2}}}$')
title('A Graph Comparing the Line of Sight \n Temperature Variations to the $\chi^{2}$ Temperature')
tight_layout()
savefig('sigma_T_chi.png', dpi=300)
close()

figure(5)
plot(np.log10(inp_N), (inp_T/chi_T), 'b.', markersize=2)
xlabel(r'$N_{input}\/(g\/cm^{-2})$')
ylabel(r'$\frac{T_{input}}{T_{\chi^{2}}}$')
title('A Graph Comparing the Input Column Density to the Temperature Ratio')
savefig('T_ratio_inp.png', dpi=300)
close()

figure(6)
plot(np.log10(chi_N), (inp_T/chi_T), 'b.', markersize=2)
xlabel(r'$N_{\chi^{2}}\/(g\/cm^{-2})$')
ylabel(r'$\frac{T_{input}}{T_{\chi^{2}}}$')
title(r'A Graph Comparing the $\chi^{2}$ Column Density to the Temperature Ratio')
savefig('T_ratio_chi.png', dpi=300)
close()
