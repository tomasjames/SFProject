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

    # Define a dendrogram instance
    d = Dendrogram.compute(contents.data, verbose=True)

    # Let the dendrogram become an interactive plot
    p = d.plotter()

    if 'T_chi' in file:
        T_chi_indices, T_chi_vals = [], []

        T_chi_data = contents.data

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Determine the type of structure
            if struct.is_leaf:

                # Extract the indices of the core and its value
                T_chi_indices.append(struct.indices(subtree=True))
                T_chi_vals.append(struct.values(subtree=True))

    elif 'N_chi' in file:
        N_chi_indices, N_chi_vals = [], []

        N_chi_data = contents.data

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Determine the type of structure
            if struct.is_leaf:

                # Extract the indices of the core and its value
                N_chi_indices.append(struct.indices(subtree=True))
                N_chi_vals.append(struct.values(subtree=True))

    elif 'T_data' in file:
        T_data_indices, T_data_vals = [], []

        T_data_data = contents.data

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Determine the type of structure
            if struct.is_leaf:

                # Extract the indices of the core and its value
                T_data_indices.append(struct.indices(subtree=True))
                T_data_vals.append(struct.values(subtree=True))

    elif 'N_data' in file:
        N_data_indices, N_data_vals = [], []

        N_data_data = contents.data

        # Find the structure within the data
        for i in range(0,len(d)):
            struct = d[i]

            # Determine the type of structure
            if struct.is_leaf:

                # Extract the indices of the core and its value
                N_data_indices.append(struct.indices(subtree=True))
                N_data_vals.append(struct.values(subtree=True))

####################### Determine core properties #########################

#dataset = ['T_chi_indices', 'T_chi_vals', 'N_chi_indices', 'N_chi_vals', 'T_data_indices', 'T_data_vals', 'N_data_indices', 'N_data_vals']
dataset = ['T_chi_indices', 'N_chi_indices', 'T_data_indices', 'N_data_indices']

imag = radmc3dPy.image.readImage('../../data/blue/dust_project/image.out')

T_chi_T, N_chi_T, T_chi_N, N_chi_N, T_data_T, N_data_T, T_data_N, N_data_N = [],[],[],[],[],[],[],[]
u_T_chi_T, M_chi_T, gpe_chi_T, ratio_store = [], [], [], []

# Loop through each dataset
for i in dataset:

    # Determine which dataset is under consideration
    if dataset[0] == i:

        # Take the values at each indices and store
        for j in range(0,len(T_chi_indices)):

            # Take the values of the temperature at the pre defined points
            T = T_chi_data[T_chi_indices[j]]
            T_chi_T.append(T)
            # Determine the internal energy per unit mass at these points
            u = (3./2)*((kk*T)/(muh2*mp))
            u_T_chi_T.append(u)

            # Take the values of the column density at the same points
            N = N_chi_data[T_chi_indices[j]]
            N_chi_T.append(N)
            # Determine the mass of the gas (assuming each core is a sphere as the area is pi*r**2) at these points
            M = N*(muh2*mp)*(np.pi*(imag.sizepix_x/2)**2)
            M_chi_T.append(M)
            # Determine the gpe
            gpe = (-gg)*M/((imag.sizepix_x/2))

            # GPE is GMm/r however because of the Virial Theorem of 2K+U=0, -U/2K = 1. U is the internal energy of the overall mass
            # i.e. U = u*M. Dividing the 2 cancels the mass.
            gpe_chi_T.append(gpe)
            ratio = (u)/(-gpe)
            ratio_store.append(ratio)

    elif dataset[1] == i:
        for j in range(0,len(N_chi_indices)):
            T_chi_N.append(T_chi_data[N_chi_indices[j]])
            N_chi_N.append(N_chi_data[N_chi_indices[j]])
    elif dataset[2] == i:
        for j in range(0,len(T_data_indices)):
            T_data_T.append(T_data_data[T_data_indices[j]])
            N_data_T.append(N_data_data[T_data_indices[j]])
    elif dataset[3] == i:
        for j in range(0,len(N_data_indices)):
            T_data_N.append(T_data_data[N_data_indices[j]])
            N_data_N.append(N_data_data[N_data_indices[j]])
