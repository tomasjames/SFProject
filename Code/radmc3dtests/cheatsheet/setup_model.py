################################################################################
############################## RADMC-3D Tests ##################################
############################# Cheatsheet Setup #################################
################################################################################

############################# Import statements ################################
import radmc3dPy
import os

########################### Create parameter file ##############################
radmc3dPy.analyze.writeDefaultParfile('ppdisk')

############################## Setup the model #################################
radmc3dPy.setup.problemSetupDust('ppdisk', mdisk='0.01*ms')

############################# Run thermal MC sim ###############################
os.system('radmc3d mctherm')

######################## Generate image from MC sim ############################
radmc3dPy.image.makeImage(npix=400, sizeau=200, wav=800., incl=45, posang=43.)

############################# Read image and plot ##############################
imag = radmc3dPy.image.readImage()
radmc3dPy.image.plotImage(imag, arcsec=True, dpc=140., log=True, maxlog=5)
