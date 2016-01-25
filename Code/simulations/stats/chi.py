################################################################################
############################ 4th Year Project Code #############################
############################# Chi Squared Routine ##############################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import radmc3dPy
from radmc3dPy.natconst import *

# Import numpy
import numpy as np

# Import csv
import csv

# Import plotting tools
from matplotlib.pyplot import *

from scipy.stats import chisquare

############################### Define functions ###############################

# Chi squared
def chi(O,E,sigma):
    '''
    Returns simply the minimised chi squared value for a given dataset.

    O is the flux from the actual data.
    Sigma is the square root of the variance of P (assume 10% of P)
    E is the data that is expected
    '''
    return sum((((O-E)/sigma)**2))

# Modified blackbody
def mbb(N,opac,wav,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity per unit solid angle, i.e. the flux observed.
    '''
    return (N*opac*((2*(hh)*(cc)**2)/wav**5)*(1/(np.exp(hh*cc/(wav*kk*T))-1)))

############################ Read in the spectral data #########################

# Read in the spectrum
spectrum = np.loadtxt('spectrum.out', skiprows=2)

wav = spectrum[:,0]
flux = spectrum[:,1]/max(spectrum[:,1])

######################### Determine (again) the opacity ########################

band = input('Which passband is to be considered? Answer with: \'PSW\', \'PMW\', \'PLW\'\n')

if band == 'PSW':
    # The start and end points of the wavelength range in microns
    lambda_init = 199.4540
    lambda_fin = 298.5657

    # Query for reference wavelength
    w_ref = 249.00983
    v_ref = cc/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'microns which equates to a frequency of', v_ref, 'Hz.\n'

# Ask for reference opacity
#kappa_0 = input('What reference intensity should be used? (1 is a safe default)\n')
kappa_0 = 0.009
print ('\nThe k_0 opacity is taken to be 0.009 cm^2/g as per Ossenkopf and Henning, 1994.\n')

# Ask for the number of wavelength points
#nlam = input('How many wavelength points do you want me to use in construction?\n')
nlam = len(spectrum)
print 'There are ', nlam, 'points in the opacity spectrum.\n'

##################### Evaluate opacities over wavelength range #################

# Create array of wavelengths from range given and reshape
w = np.linspace(lambda_init, lambda_fin, nlam)
v = np.linspace(cc/lambda_init, cc/lambda_fin, nlam)

# Ask for opacity law
opaclaw = input('Which opacity law should I use? Answer with: \'H\' (Hildebrand) \n')

# Ask for spectral index
B = input('What dust spectral index should be used? (1.7 is recommended)\n')

if opaclaw == 'H':
    # Evaluate opacities
    #opacity = kappa_0*(v/v_ref)**B
    opacity = kappa_0*(w_ref/w)**B

###################### Determine the modified black body #######################

# Generate list of column densities
N = np.linspace(1e20,1e24,1e4)

# Loop through each column density and determine modified black body curve
mod = []

for i in range(0,len(N)):
    mod.append(mbb(N[i],opacity,wav,T=10))

# Plot the figure and save for reference
figure(1)
plot(wav,mod)
xlabel('Wavelength ($\mu m$)')
ylabel('Flux')
title(str('\nThe Modified Black Body Curve for N=')+str(N[0])+str('$cm^{-1}$ to ')+str(N[-1])+str('$cm^{-1}$\n'))
savefig(str(band)+str('_modbb.png'),dpi=300)
close()

######################### Apply Chi Squared Routine ############################

# The data in mod is the expected data whilst the data read in (spectrum) is the real data
chivals = []

# Loop through each value of the expected data (i.e. each value of the column density) and determine the chi squared value for each value
for j in range(0,len(N)):
    chivals.append(chi(flux,mod[i]/max(mod[i]),sigma=(0.1*mod[i]/max(mod[i]))))

# Plot the figure and save for reference
figure(2)
errorbar(N,chivals)
xlabel('N $(cm^{-1})$')
ylabel('$\chi^2$')
title(str('\nThe $\chi^2$ Distribution for N=')+str(N[0])+str('$cm^{-1}$ to ')+str(N[-1])+str('$cm^{-1}$\n'))
savefig(str(band)+str('_chisquared_N.png'),dpi=300)
close()
