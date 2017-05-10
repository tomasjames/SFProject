# Import Numpy
import numpy as np

# Import constants from radmc3dPy
from radmc3dPy.natconst import *

class DataFileGen(width, ncell):

	def __init__(self, width, ncell):
		self.width = width
		self.ncell = ncell

	def amr_grid(self):
		############################ Define axes to work with ##########################
        # Define the image width in cm and au
        dpix_cm = width*au
        dpix_au = (dpix_au)/ncell

        # Create a dictionary of pixels in each dimension
        grid = {'x': np.linspace(-(dpix_cm)/2,(dpix_cm)/2,ncell+1), 'y': np.linspace(-(dpix_cm)/2,(dpix_cm)/2,ncell+1), 'z': np.linspace(-(dpix_cm)/2,(dpix_cm)/2,ncell+1)}
        
        ######################## Write file to working directory #######################
        # Writes a file called amr_grid.inp to the working directory and allows it to be
        # written to
        amr = open('amr_grid.inp', 'w')
