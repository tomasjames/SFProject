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
    Sigma is the square root of the variance of O (assume 10% of O)
    E is the data that is expected
    '''
    return sum((((O-E)/sigma)**2))

# Modified blackbody
def mbb(N,opac,wav,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity per unit solid angle, i.e. the flux observed.
    '''
    a = (2*(6.626e-34)*(3e8)**2)/(wav**5)
    b = ((6.6262e-34)*(3e8))/(wav*(1.38e-23)*T)
    #return N*opac*(a*(1./(np.exp(b)-1)))
    return opac*N*(a*(1./(np.exp(b)-1)))

############################ Read in the spectral data #########################

# Read in the spectrum
spectrum = np.loadtxt('spectrum.out', skiprows=2)

wav = spectrum[:,0]
flux = spectrum[:,1]

######################### Determine (again) the opacity ########################

band = input('Which passband is to be considered? Answer with: \'PSW\', \'PMW\', \'PLW\'\n')

if band == 'PSW':
    # The start and end points of the wavelength range in microns
    lambda_init = 199.4540e-6
    lambda_fin = 298.5657e-6

    # Query for reference wavelength
    w_ref = 1.30e9
    v_ref = cc*10**-2/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'microns which equates to a frequency of', v_ref, 'Hz.\n'

if band == 'PMW':
    # The start and end points of the wavelength range in microns
    lambda_init = 281.6949
    lambda_fin = 424.7548

    # Query for reference wavelength
    w_ref = 1.30e9
    v_ref = cc*10**-2/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'microns which equates to a frequency of', v_ref, 'Hz.\n'

# Ask for reference opacity
#kappa_0 = input('What reference intensity should be used? (1 is a safe default)\n')
kappa_0 = 3.09e-1
print ('\nThe k_0 opacity is taken to be 3.09e-1 cm^2/g as per Ossenkopf and Henning, 1994.\n')

# Ask for the number of wavelength points
#nlam = input('How many wavelength points do you want me to use in construction?\n')
nlam = len(spectrum)
print 'There are ', nlam, 'points in the opacity spectrum.\n'

##################### Evaluate opacities over wavelength range #################

# Create array of wavelengths from range given and reshape
w = np.linspace(lambda_init, lambda_fin, nlam)

# Ask for opacity law
opaclaw = input('Which opacity law should I use? Answer with: \'H\' (Hildebrand) \n')

# Ask for spectral index
B = input('What dust spectral index should be used? (1.7 is recommended)\n')

if opaclaw == 'H':
    # Evaluate opacities
    opacity = kappa_0*(w_ref/w)**B

###################### Determine the modified black body #######################

# Generate list of column densities
#N = np.linspace(1e21,9e22,10)
N = [1e21]

# Loop through each column density and determine modified black body curve
mod = []

# Plot the figure and save for reference
figure(1)

for i in range(0,len(N)):
    blackbody = mbb(N[i],opacity,w,T=20)
    mod.append(blackbody)

    plot(np.log10(wav),np.log10(blackbody),label=str('N=')+str(N[i]))

xlabel('Wavelength ($\mu m$)')
ylabel('Flux')
legend(loc='best')
title(str('\nThe Modified Black Body Curve for N=')+str(N[0])+str('$cm^{-2}$ to ')+str(N[-1])+str('$cm^{-2}$\n'))
savefig(str(band)+str('_modbb.png'),dpi=300)
close()

######################### Apply Chi Squared Routine ############################

# The data in mod is the expected data whilst the data read in (spectrum) is the real data
chivals = []

# Loop through each value of the expected data (i.e. each value of the column density) and determine the chi squared value for each value
for j in range(0,len(N)):
    chivals.append(chi(flux,mod[i],sigma=(0.1*mod[i])))

# Plot the figure and save for reference
figure(2)
plot(N,chivals)
xlabel('N $(cm^{-1})$')
ylabel('$\chi^2$')
title(str('\nThe $\chi^2$ Distribution for N=')+str(N[0])+str('$cm^{-1}$ to ')+str(N[-1])+str('$cm^{-1}$\n'))
savefig(str(band)+str('_chisquared_N.png'),dpi=300)
close()
