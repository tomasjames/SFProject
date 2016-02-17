################################################################################
############################## RADMC-3D Tests ##################################
############################## Gas Model Setup #################################
################################################################################

############################# Import statements ################################
import radmc3dPy
import os

########################### Create parameter file ##############################
radmc3dPy.analyze.writeDefaultParfile('ppdisk')

############################## Setup the model #################################
radmc3dPy.setup.problemSetupDust('ppdisk')

############################# Run thermal MC sim ###############################
os.system('radmc3d mctherm')

############################ Set up the gas model ##############################
radmc3dPy.setup.problemSetupGas('ppdisk')

########################### Calculate channel map ##############################
os.system('radmc3d image npix 400 sizeau 200 incl 45. phi 0. posang 43. iline 3 vkms 1.0')
