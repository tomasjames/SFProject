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

############################# Read in the data ################################

count, data, filenames = 0, [], []

# Loop through all fits files
for file in glob.glob('*.fits'):
    count += 1

    # Read in using astropy and then append the data to a list for use later
    contents = fits.getdata(str(file))

    '''
    # Append the name of the file to the position (count+len(data)) followed by
    # the file contents
    data.insert(count+len(data), file)
    data.append(contents)
    filenames.append(file)
    '''

    # Define a dendrogram instance
    d = Dendrogram.compute(contents, verbose=True)

    # Let the dendrogram become an interactive plot
    p = d.plotter()

    # Add the data to a figure
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)

    # Plot the entire tree in black
    p.plot_tree(ax, color='black')

    # Add labels to the plot
    if 'T_chi' in file:
        ax.set_xlabel('Structure')
        ax.set_ylabel('$T\/(K)$')
        title('A Dendrogram Showing the Structure of $T$ for the $\chi^{2}$ Recovered Values')
    elif 'N_chi' in file:
        ax.set_xlabel('Structure')
        ax.set_ylabel('$N\/(g\/cm^{-3})$')
        title('A Dendrogram Showing the Structure of $N$ for the $\chi^{2}$ Recovered Values')
    elif 'T_data' in file:
        ax.set_xlabel('Structure')
        ax.set_ylabel('$T\/(K)$')
        title('A Dendrogram Showing the Structure of $T$ for the Data Input Values')
    elif 'N_chi' in file:
        ax.set_xlabel('Structure')
        ax.set_ylabel('$N\/(g\/cm^{-3})$')
        title('A Dendrogram Showing the Structure of $N$ for the Data Input Values')

    #ax.imshow(contents, origin='lower', interpolation='nearest', cmap=cm.Blues, vmax=4.0)

    # Plot black contours on to the data
    # Show contour for ``min_value``
    #p.plot_contour(ax, color='black')

    # Save the image and then close to save memory and allow the next image to take its place
    fig.savefig(str(file)+str('.jpg'), dpi=300)
    close()
