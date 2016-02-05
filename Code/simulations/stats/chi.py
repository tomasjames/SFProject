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
def mbb(sigma,N,opac,v,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity per unit solid angle, i.e. the flux observed.
    '''
    #a = (2*(hh)*(cc)**2)/(wav**5)
    #b = ((hh)*(cc))/(wav*(kk)*T)
    a = (2*(hh)*(v**3))/(cc**2)
    b = (hh*v)/(kk*T)
    return sigma*N*opac*a*(1./(np.exp(b)-1))

############################ Read in the spectral data #########################

# Read in the spectrum
spectrum = np.loadtxt('spectrum.out', skiprows=2)

wav = spectrum[:,0]*10**-4
flux = spectrum[:,1]

v = cc/wav

######################### Determine (again) the opacity ########################

band = input('Which passband is to be considered? Answer with: \'PSW\', \'PMW\', \'PLW\'\n')

if band == 'PSW':
    # The start and end points of the wavelength range in microns
    lambda_init = 199.4540e-4
    lambda_fin = 298.5657e-4
    v_init = cc/lambda_init
    v_fin = cc/lambda_fin

    # Query for reference wavelength
    w_ref = 1.30e-5
    v_ref = cc/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

    # Solid angle of the beam
    sigma_arc = 465 # 465 square arcseconds
    sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

if band == 'PMW':
    # The start and end points of the wavelength range in microns
    lambda_init = 281.6949e-4
    lambda_fin = 424.7548e-4
    v_init = cc/lambda_init
    v_fin = cc/lambda_fin

    # Query for reference wavelength
    w_ref = 1.30e-4
    v_ref = cc/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

    # Solid angle of the beam
    sigma_arc = 822 # 465 square arcseconds
    sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

if band == 'PLW':
    # The start and end points of the wavelength range in microns
    lambda_init = 391.4346e-4
    lambda_fin = 690.8139e-4
    v_init = cc/lambda_init
    v_fin = cc/lambda_fin

    # Query for reference wavelength
    w_ref = 1.30e-4
    v_ref = cc/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

    # Solid angle of the beam
    sigma_arc = 1768 # 465 square arcseconds
    sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

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
v = np.linspace(v_init, v_fin, nlam)

# Ask for opacity law
opaclaw = input('Which opacity law should I use? Answer with: \'H\' (Hildebrand) \n')

# Ask for spectral index
B = input('What dust spectral index should be used? (1.7 is recommended)\n')

if opaclaw == 'H':
    # Evaluate opacities
    #opacity = kappa_0*(w_ref/w)**B
    opacity = kappa_0*(v/v_ref)**B

#################### Determine the expected column density #####################
'''
# Define dust-to-gas ratio
d2g = 0.01

# Ask for the cloud mass and convert to grams
m = input('What cloud mass should be used? (Answer should be in Solar Masses)\n')
dust_mass = d2g*m*ms
print 'The mass of the cloud considered is', m*ms, 'g. The dust mass is ', dust_mass, 'g.\n'

# Ask for the cloud number density and convert to g/cm^3
cloud = input('What number density should be assigned to the cloud?\n')
cloud_density = cloud*muh2*mp*d2g
print 'The cloud density is', cloud_density, 'g/cm^3.\n'

T = input('What temperature should be assigned to the cloud?\n')
print 'The cloud temperature is', T, 'K.\n'

# Determine the Jeans' mass
M_j = (((5*kk*T)/(gg*muh2*mp))**(3./2))*(3./(4*np.pi*cloud_density))**(1./2)
print 'The Jeans\' mass was determined to be', M_j, 'g.\n'

# Determine the Jeans' radius
R = np.sqrt((15*kk*T)/(4*np.pi*cloud_density*gg*muh2*mp))
print 'The Jeans\' length was determined to be', M_j, 'cm.\n'

# Calculate the column density
col = M_j/(np.pi*(R/2)**2)
print 'The column density was determined to be', col, 'g/cm^2.'
'''

N = np.linspace(1e24,9e24,9)

###################### Determine the modified black body #######################

# Loop through each column density and determine modified black body curve
mod = []

# Plot the figure and save for reference
figure(1)

for i in range(0,len(N)):
    blackbody = mbb(sigma,N[i],opacity,v,T=10)
    mod.append(blackbody)

    plot(v,blackbody,label=str('N=')+str(N[i]))
    #loglog(v,blackbody,label=str('N=')+str(N[i]))

xlabel(r'$\nu(Hz)$')
ylabel('Intensity $(erg/cm^{2}/s/Hz/ster)$')
legend(loc='best')
title(str('\nThe Modified Black Body Curve for N=')+str(N[0])+str('$cm^{-2}$ to ')+str(N[-1])+str('$cm^{-2}$\n'))
#xlim(min(v),max(v))
#ylim(min(blackbody),max(blackbody))
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
