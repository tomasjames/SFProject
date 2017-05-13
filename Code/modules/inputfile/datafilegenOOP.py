# Import Numpy
import numpy as np

# Import constants from radmc3dPy
from radmc3dPy.natconst import *


class DataFileGen:

    def __init__(self, width, ncell):
        self.width = width
        self.ncell = ncell

        self.image_width_cm = self.width*au # Converts width in au to width in cm
        self.cell_width_cm = self.image_width_cm/self.ncell # Determines the width of a cell in cm

    def dimensions(self): 
        '''
        Returns a dictionary of pixels in x, y and z dimensions
        '''
        return {
            'x': np.linspace(-(self.image_width_cm)/2,(self.image_width_cm)/2,self.ncell+1), 
            'y': np.linspace(-(self.image_width_cm)/2,(self.image_width_cm)/2,self.ncell+1), 
            'z': np.linspace(-(self.image_width_cm)/2,(self.image_width_cm)/2,self.ncell+1)
        }

    def write_amr_grid(self):

        amr = open('amr_grid.inp', 'w') # Writes file to the current directory
        amr.write('1 # The same as most other files\n') # Write the iformat number
        amr.write('0 # 0 = regular\n') # Write the grid style
        amr.write('1 # < 100 = Cartesian\n') # Write coordinate system
        amr.write('0 # If 1, abundant information is written regarding the grid\n') # Write gridinfo integer
        amr.write('1 1 1 \n') # Write integers to activate dimensions
        amr.write(str(self.ncell)+ ' ' + str(self.ncell) + ' ' +str(self.ncell) + '\n') # Write the number of base cells in each dimension

        for dimension in self.dimensions():
            print dimension
            for i in self.dimensions()[dimension]:
                if i == self.dimensions()[dimension][-1]:
                    print 'In the last pixel'
                    amr.write(str(np.float64(i)) + str('\n'))
                else:
                    amr.write(str(np.float64(i)) + str(' '))

        amr.close()
