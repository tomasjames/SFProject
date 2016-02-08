################################################################################
############################ 4th Year Project Code #############################
########################### Chi Squared Routine v2 #############################
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
def mbb(N,opac,v,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity per unit solid angle, i.e. the flux observed.
    '''
    #a = (2*(hh)*(cc)**2)/(wav**5)
    #b = ((hh)*(cc))/(wav*(kk)*T)
    a = (2*(hh)*(v**3))/(cc**2)
    b = (hh*v)/(kk*T)
    return N*opac*a*(1./(np.exp(b)-1))

################### Read the average data determined earlier ###################

# Reads the data
average_data = np.loadtxt('SPIRE_average_data.txt', delimiter=" ")

# Sort the data
average_data_sort = average_data[np.argsort(average_data[:,0])]

# Assign to variables (data plotted later)
plw,pmw,psw = average_data_sort[0], average_data_sort[1], average_data_sort[2]

######################### Determine (again) the opacity ########################

# Query for reference wavelength
w_ref = 1.30e-5
v_ref = cc/w_ref
print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

# Ask for reference opacity
kappa_0 = 3.09e-1
print ('\nThe k_0 opacity is taken to be 3.09e-1 cm^2/g as per Ossenkopf and Henning, 1994.\n')

# Define the start and end points of the spectrum
lambda_init = 199.4540e-4
lambda_fin = 690.8139e-4
v_init = cc/lambda_init
v_fin = cc/lambda_fin
nlam = 300

# Create array of wavelengths from range given and reshape
w = np.linspace(lambda_init, lambda_fin, nlam)
v = np.linspace(v_init, v_fin, nlam)

# Solid angle of the beam
sigma_arc = 1768 # 465 square arcseconds
sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

# Ask for spectral index
#B = input('What dust spectral index should be used? (1.7 is recommended)\n')
B = 2.0

# Evaluate opacities
#opacity = kappa_0*(w_ref/w)**B
opacity = kappa_0*(v/v_ref)**B

'''
band = input('Which passband is to be considered? Answer with: \'PSW\', \'PMW\', \'PLW\'\n')

if band == 'PSW':
    # The start and end points of the wavelength range in microns
    lambda_init = 199.4540e-4
    lambda_fin = 298.5657e-4
    lambda_eff = 247.12451e-4
    v_init = cc/lambda_init
    v_fin = cc/lambda_fin
    v_eff = cc/lambda_eff

    # Solid angle of the beam
    sigma_arc = 465 # 465 square arcseconds
    sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

if band == 'PMW':
    # The start and end points of the wavelength range in microns
    lambda_init = 281.6949e-4
    lambda_fin = 424.7548e-4
    lambda_eff = 346.71804e-4
    v_init = cc/lambda_init
    v_fin = cc/lambda_fin
    v_eff = cc/lambda_eff

    # Solid angle of the beam
    sigma_arc = 822 # 465 square arcseconds
    sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

if band == 'PLW':
    # The start and end points of the wavelength range in microns
    lambda_init = 391.4346e-4
    lambda_fin = 690.8139e-4
    lambda_eff = 496.10677e-4
    v_init = cc/lambda_init
    v_fin = cc/lambda_fin
    v_eff = cc/lambda_eff

    # Solid angle of the beam
    sigma_arc = 1768 # 465 square arcseconds
    sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

# Ask for the number of wavelength points
nlam = 100
print 'There are ', nlam, 'points in the opacity spectrum.\n'

#################### Determine the expected column density #####################

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

# Determine the Jeans' length
R = np.sqrt((15*kk*T)/(4*np.pi*cloud_density*gg*muh2*mp))
print 'The Jeans\' length was determined to be', M_j, 'cm.\n'

# Calculate the column density
col = M_j/(np.pi*(R/2)**2)
print 'The column density was determined to be', col, 'g/cm^2.'
'''
#N = np.linspace(1e24,9e24,9)
N = [125]

###################### Determine the modified black body #######################

# Loop through each column density and determine modified black body curve
mod = []

# Plot the figure and save for reference
figure(1)

for i in range(0,len(N)):
    blackbody = mbb(N[i],opacity,v,T=10)
    mod.append(blackbody)

    loglog(v,blackbody,label=str('N=')+str(N[i]))
    #loglog(v,blackbody,label=str('N=')+str(N[i]))

# Plot the data
loglog(psw[0],psw[1],'bx',label='PSW')
loglog(pmw[0],pmw[1],'gx',label='PMW')
loglog(plw[0],plw[1],'rx',label='PLW')
xlabel(r'$\nu (Hz)$')
ylabel(r'Flux $(erg/cm^{2}/s/Hz/ster)$')
legend(loc='best')
savefig('SPIRE_averages.png',dpi=300)
close()

'''
xlabel(r'$\nu(Hz)$')
ylabel('Intensity $(erg/cm^{2}/s/Hz/ster)$')
legend(loc='best')
title(str('\nThe Modified Black Body Curve for N=')+str(N[0])+str('$cm^{-2}$ to ')+str(N[-1])+str('$cm^{-2}$\n'))
#xlim(min(v),max(v))
#ylim(min(blackbody),max(blackbody))
savefig(str(band)+str('_modbb.png'),dpi=300)
close()


####################### Perform weighted average routine #######################

# Multiply each intensity by the transmission at that wavelength/flux
weighted_flux = blackbody*trans
weighted_v = v*trans

# Average the above array to determine the actual flux that Herschel 'sees'
weighted_flux_mean = np.mean(weighted_flux)
weighted_flux_std = np.std(weighted_flux)
weighted_v_mean = np.mean(weighted_v)
weighted_v_std = np.std(weighted_v)

figure(2)
errorbar(v_eff,weighted_flux_mean,yerr=weighted_flux_std,fmt='x')
#errorbar(v,weighted_flux,'g--',yerr=np.std(weighted_flux))
plot(v,weighted_flux,'g--')
xlabel(r'$\nu (Hz)$')
ylabel(r'Flux $(erg/cm^{2}/s/Hz/ster)$')
savefig(str(band)+str('_average.png'),dpi=300)
close()

######################### Apply Chi Squared Routine ############################

# The data in mod is the expected data whilst the data read in (spectrum) is the real data
chivals = []

# Loop through each value of the expected data (i.e. each value of the column density) and determine the chi squared value for each value
for j in range(0,len(N)):
    chivals.append(chi(flux,mod[i],sigma=(0.1*mod[i])))

# Plot the figure and save for reference
figure(3)
plot(N,chivals)
xlabel('N $(cm^{-1})$')
ylabel('$\chi^2$')
title(str('\nThe $\chi^2$ Distribution for N=')+str(N[0])+str('$cm^{-1}$ to ')+str(N[-1])+str('$cm^{-1}$\n'))
savefig(str(band)+str('_chisquared_N.png'),dpi=300)
close()
'''
