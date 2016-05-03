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

########################### Determine expected column density ######################

# Define the density and temperature expected of a core
rho_0 = (10**4)*mp*muh2
T_0 = 15

# Determine the Jeans mass and radius and use these to define the column density
#M_j_SI = (((5*k*T_0)/(G*muh2*mh))**(3./2)*(3./(4*np.pi*rho_0))**(1./2))*1000
M_j = (((5*kk*T_0)/(gg*muh2*mp))**(3./2))*(3./(4*np.pi*rho_0))**(1./2)
#R_j = (np.sqrt((15*k*(T_0))/(4*np.pi*G*muh2*mh*rho_0)))*100
R_j = ((15*kk*(T_0))/(4*np.pi*gg*muh2*mp*rho_0))**(1./2)
N_j = (M_j/(muh2*mp*(R_j)**2))

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
    print 'File is...', file

    count += 1

    # Read in using astropy and then append the data to a list for use later
    contents = fits.open(str(file))[1]

    if 'N_chi_inp' in file:

        # Define a dendrogram instance
        d = Dendrogram.compute(contents.data, min_value=N_j/100, min_npix=4, verbose=True)

        N_chi_indices, N_chi_vals = [], []

        N_chi_data = contents.data

        # Save a file for storage of determined values
        N_chi_mass = open('N_chi_mass.txt', 'w')
        cs = csv.writer(N_chi_mass, delimiter=' ')

        # Save a line to allow better understanding of the columns
        cs_towrite = ['First x Pixel', 'First y Pixel', 'U', 'GPE', 'Mass/ms', 'Ratio']
        cs.writerow(cs_towrite)

        # Open the temperature file and extract the numerical data
        T_data = fits.open(str('T_chi_inp.fits'))[1].data

        # Add the data to a figure
        fig = figure()

        # Let the dendrogram become an interactive plot
        p = d.plotter()

        # Define subplots to plot with
        ax1 = subplot2grid((6,6), (0,0), colspan=4,rowspan=4)
        ax2 = subplot2grid((6,6), (5,0), colspan=4,rowspan=1)

        ax1.imshow(contents.data, origin='lower', interpolation='nearest')

        #p.plot_contour(ax1, color='black')

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Get the mask of the region highlighted
            mask = struct.get_mask()

            # Create FITS HDU objects
            mask_hdu = fits.PrimaryHDU(mask.astype('short'), contents.header)

            # Determine the type of structure
            if struct.is_leaf:

                # Extract the numerical values at the leaf
                N = struct.values(subtree=True)
                ind_N = struct.indices(subtree=True)

                # Determine the number density of the leaf
                num_pix = len(struct.indices(subtree=True)[0]) # Determine the number of pixels in the leaf
                tot_area = num_pix*(imag.sizepix_x*imag.sizepix_y) # Determine the total area (area of one pixel multiplied by the number of pixels)
                R_core = np.sqrt(tot_area/np.pi) # If those pixels in the leaf are spherical then they would have an area piR**2 = area, so R = (area/pi)**0.5
                M = np.mean(N*100)*(muh2*mp)*(np.pi*(R_core)**2) # Determine the mass

                # Pull temperature at the given indices
                T_weight = T_data[ind_N]
                T = np.sum(T_weight*N*100)/np.sum(N*100)

                # Extract the indices of the core and its value
                N_chi_indices.append(ind_N)
                N_chi_vals.append(N*100)

                # Append the earlier determined mass to the list
                M_chi_list.append(M)

                # Determine the gpe (GMm/r). Neglect the second M term (the term coming from infinity) as this cancels when determining the virial ratio later
                gpe = (3./5)*(-gg)*M/((R_core))
                gpe_chi.append(gpe)

                # Determine the internal energy per unit mass
                u = (3./2)*((kk*T)/(muh2*mp))

                # Virial theorem: 2K+U=0, -U/2K=1. U is the internal energy of the overall mass, i.e. U=u*M.
                ratio = (u)/(-gpe)
                ratio_chi.append(ratio)

                if ratio > 1:
                    p.plot_contour(ax1, structure=struct, lw=2, colors='red')
                    bound_store = 'Unbound'
                elif ratio <= 1:
                    p.plot_contour(ax1, structure=struct, lw=2, colors='green')
                    bound_store = 'Bound'

                    cs_towrite = [ind_N[0][0], ind_N[1][0], gpe, u, M/ms, ratio]

                    # Write to file
                    cs.writerow(cs_towrite)

        # Plot the entire tree in black
        p.plot_tree(ax2, color='black')

        print 'Evaluated leaves...'

        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$N_{dust}\/(cm^{-2})$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        print 'Saving...'

        # Save the image and then close to save memory and allow the next image to take its place
        fig.savefig(str(file)+str('.jpg'), dpi=300)
        close()

        print 'Saved!'

        N_chi_mass.close()

        print 'Closed...'

    if 'N_data_inp' in file:

        # Define a dendrogram instance
        d = Dendrogram.compute(contents.data, min_value=N_j/100, min_delta=3, min_npix=4, verbose=True)

        N_data_indices, N_data_vals = [], []

        N_data_data = contents.data

        # Save a file for storage of determined values
        N_data_mass = open('N_data_mass.txt', 'w')
        cs = csv.writer(N_data_mass, delimiter=' ')

        # Save a line to allow better understanding of the columns
        cs_towrite = ['First x Pixel', 'First y Pixel','U', 'GPE', 'Mass/ms', 'Ratio']
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

        #p.plot_contour(ax1, color='black')

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Get the mask of the region highlighted
            mask = struct.get_mask()

            # Create FITS HDU objects
            mask_hdu = fits.PrimaryHDU(mask.astype('short'), contents.header)

            # Determine the type of structure
            if struct.is_leaf:

                # Extract the numerical values at the leaf
                N = struct.values(subtree=True)
                ind_N = struct.indices(subtree=True)

                # Determine the number density of the leaf
                num_pix = len(struct.indices(subtree=True)[0]) # Determine the number of pixels in the leaf
                tot_area = num_pix*(imag.sizepix_x*imag.sizepix_y) # Determine the total area (area of one pixel multiplied by the number of pixels)
                R_core = np.sqrt(tot_area/np.pi) # If those pixels in the leaf are spherical then they would have an area piR**2 = area, so R = (area/pi)**0.5
                M = np.mean(N*100)*(muh2*mp)*(np.pi*(R_core)**2) # Determine the mass

                # Pull temperature at the given indices
                T_weight = T_data[ind_N]
                T = np.sum(T_weight*N*100)/np.sum(N*100)

                # Extract the indices of the core and its value
                N_data_indices.append(ind_N)
                N_data_vals.append(N*100)

                # Append the earlier determined mass to the list
                M_data_list.append(M)

                # Determine the gpe (GMm/r). Neglect the second M term (the term coming from infinity) as this cancels when determining the virial ratio later
                gpe = (3./5)*(-gg)*M/((R_core))
                gpe_chi.append(gpe)

                # Determine the internal energy per unit mass
                u = (3./2)*((kk*T)/(muh2*mp))

                # Virial theorem: 2K+U=0, -U/2K=1. U is the internal energy of the overall mass, i.e. U=u*M.
                ratio = (u)/(-gpe)
                ratio_chi.append(ratio)

                if ratio > 1:
                    p.plot_contour(ax1, structure=struct, lw=2, colors='red')
                    bound_store = 'Unbound'
                elif ratio <= 1:
                    p.plot_contour(ax1, structure=struct, lw=2, colors='green')
                    bound_store = 'Bound'

                    cs_towrite = [ind_N[0][0], ind_N[1][0], gpe, u, M/ms, ratio]

                    # Write to file
                    cs.writerow(cs_towrite)

        # Plot the entire tree in black
        p.plot_tree(ax2, color='black')

        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $N$ for the Data Derived Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$N_{dust}\/(cm^{-2})$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the Data Derived Values')

        # Save the image and then close to save memory and allow the next image to take its place
        fig.savefig(str(file)+str('.jpg'), dpi=300)
        close()

        N_data_mass.close()

    else:
        # Define a dendrogram instance
        d = Dendrogram.compute(contents.data, verbose=True)

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

######################## Plot recovered masses ##########################

# Read in the data
chi = np.loadtxt('N_chi_mass.txt',skiprows=1)
data = np.loadtxt('N_data_mass.txt',skiprows=1)

psw_raw = np.loadtxt('../../data/psw/dust_project/image_trans.out',skiprows=6)
eff = 247.12451

psw = np.reshape(psw_raw, (np.sqrt(len(psw_raw)),np.sqrt(len(psw_raw))))

# Pull out the pixel locations
xpix_chi, ypix_chi = chi[:,0], chi[:,1]
xpix_data, ypix_data = data[:,0], data[:,1]

# Pull out the masses
mass_chi = chi[:,4]
mass_data = data[:,4]

figure(1)
imshow(psw,origin='left')
colorbar()
plot(xpix_chi, ypix_chi, 'c.', label=r'$\chi^{2}$')
plot(xpix_data, ypix_data, 'm.', label=r'$Data$')
xlim(0,200)
ylim(0,200)
xlabel('X (Pixels)')
ylabel('Y (Pixels)')

savefig('masses.png',dpi=300)
close()
