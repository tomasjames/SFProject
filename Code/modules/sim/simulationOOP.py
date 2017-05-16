import urllib2  # the lib that handles the url stuff

import numpy as np

import radmc3dPy # Import RADMC-3D

import csv # Import csv

class Simulation:

    def __init__(self, cenwav, trans_data, imag):
        self.cenwav = cenwav
        self.trans_data = trans_data
        self.imag = imag

    def transmission(self, trans_data, imag):
        '''
        Method to take observational data and degrade it using filter transmission curve
        
        Arguments:
            trans_data: 2D array containing transmission coefficient and wavelength/frequency at each coefficient
            imag: radmc3dPy instance of the raytraced imaged, i.e. radmc3dPy.image.readImage()
        '''

        with open('image_trans.out', 'w') as f:
            image_trans = csv.writer(f, delimiter=' ')

            image_trans.writerow([1]) # Begin writing of the file by writing the format number
            image_trans.writerow([imag.nx, imag.ny]) # Write the number of pixels in x and y dimensions
            image_trans.writerow([1]) # Write the number of wavelength points
            image_trans.writerow([imag.sizepix_x, imag.sizepix_y]) # Write the pixel sizes
            image_trans.writerow(["%.15f" % self.cenwav]) # # Write the filter's central wavelength
            image_trans.writerow([]) # Writes a blank line to seperate the wavelength points from the intensity points

            # Begin computing transmission weighted sum to produce composite image
            # This is to be achieved by looping through each pixel in the image and multiplying by each transmission coefficient at that wavelength. This is then summed along all wavelengths to give a 'composite' image

            # Declare a list to store the variables
            store, store_all, trans_store = [], [], []
            summation = []

            for i in range(0, imag.nx):
                for j in range(0, imag.ny):
                    for k in range(0, imag.nwav):
                        trans_store.append(np.float64(trans_data[k][1])) # Appends the transmission to a list
                        store.append(np.float64(imag.image[j][i][k]*trans_data[k][1])) # Determines the transmission weighting
                        store_all.append(np.float64(imag.image[j][i][k]*trans_data[k][1]))
                    summation.append(np.float64(np.sum(store)/np.sum(trans_store))) # Reduces that weighting to one 
                    store, trans_store = [], []

            image_trans.writerows(zip(np.float64(summation)))

            # Save store_all to another .txt file to SED
            with open('image_trans_raw.txt', 'w') as g:
                image_trans_raw = csv.writer(g, delimiter=' ')

                image_trans_raw.writerows(zip(np.float64(store_all)))
