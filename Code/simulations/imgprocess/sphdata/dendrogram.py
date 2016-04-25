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

############################# Read in the data ################################

count, data, filenames = 0, [], []

# Loop through all fits files
for file in glob.glob('*.fits'):
    count += 1

    # Read in using astropy and then append the data to a list for use later
    contents = fits.open(str(file))[1]

    '''
    # Append the name of the file to the position (count+len(data)) followed by
    # the file contents
    data.insert(count+len(data), file)
    data.append(contents)
    filenames.append(file)
    '''

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
            leaves.append(struct)

            # Extract the indices of the core and its value
            indices = struct.indices(subtree=True)
            vals = struct.values(subtree=True)

            # Highlight branches
            p.plot_contour(ax1, structure=struct, lw=3, colors="#%06x" % random.randint(0, 0xFFFFFF))

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
    elif 'N_chi' in file:
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$N\/(g\/cm^{-3})$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the $\chi^{2}$ Recovered Values')
    elif 'T_data' in file:
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $T$ for the Data Input Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$T\/(K)$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $T$ for the Data Input Values')
    elif 'N_data' in file:
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_title('A Map Showing the Structure\n of $N$ for the Data Input Values')

        ax2.set_xlabel('Structure')
        ax2.set_ylabel('$N\/(g\/cm^{-3})$')
        ax2.set_title('A Dendrogram Showing the Structure\n of $N$ for the Data Input Values')

    #ax1.imshow(contents, origin='lower', interpolation='nearest', cmap=cm.Blues, vmax1=4.0)

    # Plot black contours on to the data
    # Show contour for ``min_value``
    #p.plot_contour(ax1, color='black')

    # Save the image and then close to save memory and allow the next image to take its place
    fig.savefig(str(file)+str('.jpg'), dpi=300)
    close()
