################################################################################
############################ 4th Year Project Code #############################
############################## Plot Dendrograms ################################
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

# Import plotting tools
from matplotlib.pyplot import *
from astropy.io import fits
from astrodendro import Dendrogram

# Import file mining
import glob

# Import random
import random

############################### Read in image data ##############################

# Read in both data types
inp_data = np.loadtxt('../../curvefitting/sphdata/datafeed.txt')
chi_data = np.loadtxt('../../curvefitting/sphdata/chi_fine.txt',skiprows=1)

# Split data types into plottable quantities
inp_N = inp_data[:,1]
inp_T = inp_data[:,2]

chi_N, chi_N_error = chi_data[:,1], chi_data[:,-2]
chi_T, chi_T_error = chi_data[:,2], chi_data[:,-1]

# Reshape the data such that x and y pixels correspond
N_chi_inp, N_chierror_inp = np.reshape(chi_N, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))), np.reshape(chi_N_error, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))
T_chi_inp, T_chierror_inp = np.reshape(chi_T, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data)))), np.reshape(chi_T_error, (np.sqrt(len(chi_data)),np.sqrt(len(chi_data))))

N_data_inp = np.reshape(inp_N, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))
T_data_inp = np.reshape(inp_T, (np.sqrt(len(inp_data)),np.sqrt(len(inp_data))))

######################## Analyse data with dendrograms ############################

# Define naive tracking lists
count, data, filenames = 0, [], []

# Define data storage
M_chi_list, M_data_list = [], []
gpe_chi, gpe_data = [], []
ratio_chi, ratio_data = [], []

# Import the imag to allow image dimensions to be used
imag = radmc3dPy.image.readImage('../../data/blue/dust_project/image.out')

# Define the number of pixels in one leaf to be regarded a detection
condition = 4.

# Loop through all fits files
for file in glob.glob('*.fits'):
    count += 1

    # Read in using astropy and then append the data to a list for use later
    contents = fits.open(str(file))[1]

    # Define a dendrogram instance
    d = Dendrogram.compute(contents.data, verbose=True)

    # Let the dendrogram become an interactive plot
    #p = d.plotter()

    if 'N_chi' in file:
        N_chi_indices, N_chi_vals = [], []

        N_chi_data = contents.data

        # Save a file for storage of determined values
        N_chi_mass = open('N_chi_mass.txt', 'w')
        cs = csv.writer(N_chi_mass, delimiter=' ')

        # Save a line to allow better understanding of the columns
        cs_towrite = ['Index', 'U', 'GPE', 'Mass/ms', 'Ratio']
        cs.writerow(cs_towrite)

        # Open the temperature file and extract the numerical data
        T_chi = fits.open(str('T_chi_inp.fits'))[1].data

        # Add the data to a figure
        fig = figure()

        # Let the dendrogram become an interactive plot
        p = d.plotter()

        # Define subplots to plot with
        ax1 = subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
        ax2 = subplot2grid((6,6), (5,0), colspan=4,rowspan=1)

        ax1.imshow(contents.data, origin='lower', interpolation='nearest')

        p.plot_contour(ax1, color='black')

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Get the mask of the region highlighted
            mask = struct.get_mask()

            # Create FITS HDU objects
            mask_hdu = fits.PrimaryHDU(mask.astype('short'), contents.header)

            # Determine the type of structure
            if struct.is_leaf:

                # Determine whether the number of pixels in the leaf is likely to be a core
                if np.float64(len(struct.indices(subtree=True)[0])) >= np.float64(condition) and np.sum(N_chi_inp[struct.indices(subtree=True)]/N_chierror_inp[struct.indices(subtree=True)]) >= np.float64(2.):
                #if np.sum(N_chi_inp[struct.indices(subtree=True)]/N_chierror_inp[struct.indices(subtree=True)]) >= np.float64(3.):
                    N = struct.values(subtree=True)
                    ind_N = struct.indices(subtree=True)

                    # Extract the indices of the core and its value
                    N_chi_indices.append(ind_N)
                    N_chi_vals.append(N)

                    # Determine the mass of the gas (assuming each core is a sphere as the area is pi*r**2) at these points
                    M = np.sum(N*(muh2*mp)*(np.pi*(imag.sizepix_x/2)**2))
                    M_chi_list.append(M)

                    # Extract the chi-squared recovered temperatures at these points
                    T_weight = T_chi[ind_N]
                    T = np.sum(T_weight*N)/np.sum(N)

                    # Determine the gpe (GMm/r). Neglect the second M term (the term coming from infinity) as this cancels when determining the virial ratio later
                    gpe = (-gg)*M/((imag.sizepix_x/2))
                    gpe_chi.append(gpe)

                    # Determine the internal energy per unit mass
                    u = (3./2)*((kk*T)/(muh2*mp))

                    # Virial theorem: 2K+U=0, -U/2K=1. U is the internal energy of the overall mass, i.e. U=u*M.
                    ratio = (u)/(-gpe)
                    ratio_chi.append(ratio)

                    if ratio < 1:
                        p.plot_contour(ax1, structure=struct, lw=2, colors='red')
                    elif ratio > 1:
                        p.plot_contour(ax1, structure=struct, lw=2, colors='green')

                    cs_towrite = [struct.indices(subtree=True), gpe, u, M/ms, ratio]

                    # Write to file
                    cs.writerow(cs_towrite)

        # Plot the entire tree in black
        p.plot_tree(ax2, color='black')

        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$N\/(g\/cm^{-3})$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        # Save the image and then close to save memory and allow the next image to take its place
        fig.savefig(str(file)+str('.jpg'), dpi=300)
        close()

        N_chi_mass.close()

    elif 'N_data' in file:

        N_data_indices, N_data_vals = [], []

        N_data_data = contents.data

        # Save a file for storage of determined values
        N_data_mass = open('N_data_mass.txt', 'w')
        cs = csv.writer(N_data_mass, delimiter=' ')

        # Save a line to allow better understanding of the columns
        cs_towrite = ['Index', 'U', 'GPE', 'Mass/ms', 'Ratio']
        cs.writerow(cs_towrite)

        # Open the temperature file and extract the numerical data
        T_data = fits.open(str('T_data_inp.fits'))[1].data

        # Add the data to a figure
        fig = figure()

        # Let the dendrogram become an interactive plot
        p = d.plotter()

        # Define subplots to plot with
        ax1 = subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
        ax2 = subplot2grid((6,6), (5,0), colspan=4,rowspan=1)

        ax1.imshow(contents.data, origin='lower', interpolation='nearest')

        p.plot_contour(ax1, color='black')

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Get the mask of the region highlighted
            mask = struct.get_mask()

            # Create FITS HDU objects
            mask_hdu = fits.PrimaryHDU(mask.astype('short'), contents.header)

            # Determine the type of structure
            if struct.is_leaf:

                # Determine whether the number of pixels in the leaf is likely to be a core
                if np.float64(len(struct.indices(subtree=True)[0])) >= np.float64(condition):
                    N = struct.values(subtree=True)
                    ind_N = struct.indices(subtree=True)

                    # Extract the indices of the core and its value
                    N_data_indices.append(ind_N)
                    N_data_vals.append(N)

                    # Determine the mass of the gas (assuming each core is a sphere as the area is pi*r**2) at these points
                    M = np.sum(N*(muh2*mp)*(np.pi*(imag.sizepix_x/2)**2))
                    M_data_list.append(M)

                    # Extract the chi-squared recovered temperatures at these points
                    T_weight = T_data[ind_N]
                    T = np.sum(T_weight*N)/np.sum(N)

                    # Determine the gpe
                    gpe = (-gg)*M/((imag.sizepix_x/2))
                    gpe_data.append(gpe)

                    # Determine the internal energy per unit mass
                    u = (3./2)*((kk*T)/(muh2*mp))

                    # Virial theorem: 2K+U=0, -U/2K=1. U is the internal energy of the overall mass, i.e. U=u*M.
                    ratio = (u)/(-gpe)
                    ratio_data.append(ratio)

                    if ratio < 1:
                        p.plot_contour(ax1, structure=struct, lw=2, colors='red')
                    elif ratio > 1:
                        p.plot_contour(ax1, structure=struct, lw=2, colors='green')

                    cs_towrite = [struct.indices(subtree=True), gpe, u, M/ms, ratio]

                    # Write to file
                    cs.writerow(cs_towrite)

        # Plot the entire tree in black
        p.plot_tree(ax2, color='black')

        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$N\/(g\/cm^{-3})$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        # Save the image and then close to save memory and allow the next image to take its place
        fig.savefig(str(file)+str('.jpg'), dpi=300)
        close()

        N_data_mass.close()

    else:
        # Let the dendrogram become an interactive plot
        p = d.plotter()

        # Add the data to a figure
        fig = figure()

        # Define subplots to plot with
        ax1 = subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
        ax2 = subplot2grid((6,6), (5,0), colspan=4,rowspan=1)
        #ax13 = subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
        #ax14 = subplot2grid((3, 3), (1, 2), rowspan=2)

        ax1.imshow(contents.data, origin='lower', interpolation='nearest')
        '''
        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Get the mask of the region highlighted
            mask = struct.get_mask()

            # Create FITS HDU objects
            mask_hdu = fits.PrimaryHDU(mask.astype('short'), contents.header)

            ax1.imshow(contents.data, origin='lower', interpolation='nearest')

            # Show contour for ``min_value``
            p.plot_contour(ax1, color='black')

            # Determine the type of structure
            if struct.is_leaf:
                #leaves.append(struct)

                # Extract the indices of the core and its value
                indices = struct.indices(subtree=True)
                vals = struct.values(subtree=True)

                # Highlight branches
                #p.plot_contour(ax1, structure=struct, lw=3, colors="#%06x" % random.randint(0, 0xFFFFFF))
                p.plot_contour(ax1, structure=struct, lw=3)
        '''
        # Plot the entire tree in black
        p.plot_tree(ax2, color='black')

        # Add labels to the plot
        if 'T_chi' in file:
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_title('A Map Showing the Structure\n of $T$ for the $\chi^{2}$ Recovered Values')

            ax2.set_xlabel('Structure')
            ax2.set_ylabel('$T\/(K)$')
            ax2.set_title('A Dendrogram Showing the Structure\n of $T$ for the $\chi^{2}$ Recovered Values')

        elif 'T_data' in file:
            ax1.set_xlabel('X')
            ax1.set_ylabel('Y')
            ax1.set_title('A Map Showing the Structure\n of $T$ for the Data Input Values')

            ax2.set_xlabel('Structure')
            ax2.set_ylabel('$T\/(K)$')
            ax2.set_title('A Dendrogram Showing the Structure\n of $T$ for the Data Input Values')

        #ax1.imshow(contents, origin='lower', interpolation='nearest', cmap=cm.Blues, vmax1=4.0)

        # Plot black contours on to the data
        # Show contour for ``min_value``
        #p.plot_contour(ax1, color='black')

        # Save the image and then close to save memory and allow the next image to take its place
        fig.savefig(str(file)+str('.jpg'), dpi=300)
        close()
