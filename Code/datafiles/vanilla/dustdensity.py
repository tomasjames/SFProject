################################################################################
############################ 4th Year Project Code #############################
########################## Dust Denstiry Assignment ############################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import Numpy
import numpy as np

############## Read 'amr_grid.inp' to get grid information #####################

# Read amr_grid.inp and skip first
file = open('amr_grid.inp', 'r')
f = file.read()
