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
def mbb(N,sigma,opac,v,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity per unit solid angle, i.e. the flux observed.
    '''
    #a = (2*(hh)*(cc)**2)/(wav**5)
    #b = ((hh)*(cc))/(wav*(kk)*T)
    a = (2*(hh)*(v**3))/(cc**2)
    b = (hh*v)/(kk*T)
    return N*sigma*opac*a*(1./(np.exp(b)-1))

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
#sigma_arc = 1768 # 465 square arcseconds
#sigma = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

# Ask for spectral index
#B = input('What dust spectral index should be used? (1.7 is recommended)\n')
B = 2.0

# Evaluate opacities
#opacity = kappa_0*(w_ref/w)**B
opacity = kappa_0*(v/v_ref)**B

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
M_j = ((5./2)*(np.sqrt(15./8*np.pi))*(kk/gg)**(3./2))*(dust_mass**-2)*((T**3)/(cloud))**1./2
print 'The Jeans\' mass was determined to be', M_j, 'g.\n'

# Determine the Jeans' length
R_j = np.sqrt((15./8*np.pi)*((kk*T)/(gg*(dust_mass**2)*cloud)))
print 'The Jeans\' length was determined to be', R_j, 'cm.\n'

# Calculate the column density
col = (M_j/(np.pi*(R_j**2)))
print 'The column density was determined to be', col, 'g/cm^2.'

#N = np.linspace(1e24,9e24,9)
N = np.linspace(0,2*col,200)
r = 7171 # Radius of the cloud in AU
d = 1 # Distance to source in parsecs
sigma = (np.pi*(r*1.496e+11)**2)/(d*3.086e+16)**2
#N = np.linspace(1,400,3990)

###################### Determine the modified black body #######################

# Loop through each column density and determine modified black body curve
mod = []

# Plot the figure and save for reference
figure(1)

for i in range(0,len(N)):
    blackbody = mbb(N[i],sigma,opacity,v,T=10)
    mod.append(blackbody)

    #loglog(v,blackbody,label=str('N=')+str(N[i]))

# Plot the data
loglog(psw[0],psw[1],'bx',label='PSW')
loglog(pmw[0],pmw[1],'gx',label='PMW')
loglog(plw[0],plw[1],'rx',label='PLW')
xlabel(r'$\nu \/(Hz)$')
ylabel(r'Intensity $(erg/cm^{2}/s/Hz/ster)$')

######################### Apply Chi Squared Routine ############################

fit_storage = []

# The data in mod is the expected data whilst the data read in (spectrum) is the real data
chivals = []

for i in range(0,len(N)):

    # Find the element in the theoretical SED closest to the band points
    # This will return (index,value)
    to_fit_psw = min(enumerate(mod[i]), key=lambda x: abs(x[1]-psw[1]))
    to_fit_pmw = min(enumerate(mod[i]), key=lambda y: abs(y[1]-pmw[1]))
    to_fit_plw = min(enumerate(mod[i]), key=lambda z: abs(z[1]-plw[1]))

    print to_fit_psw, to_fit_pmw, to_fit_plw

    # Pull out the frequency at each of those points
    v_psw = v[to_fit_psw[0]]
    v_pmw = v[to_fit_pmw[0]]
    v_plw = v[to_fit_plw[0]]

    to_fit = np.array([to_fit_psw[1],to_fit_pmw[1],to_fit_plw[1]])
    ps = np.array([psw[1],pmw[1],plw[1]])

    # Loop through each value of the expected data (i.e. each value of the column density) and determine the chi squared value for each value
    chivals.append(chi(to_fit,ps,sigma=(0.1*ps)))

# Plot the figure and save for reference
figure(2)
plot(N,chivals)
grid(True,which='both')
xlabel('N $(cm^{-1})$')
ylabel('$\chi^2$')
title(str('\nThe $\chi^2$ Distribution for N=')+str(N[0])+str('$cm^{-1}$ to ')+str(N[-1])+str('$cm^{-1}$\n'))
savefig('chisquared_N.png',dpi=300)
close()

# Determine the best fit blackbody
chi_min_index = chivals.index(min(chivals))
chi_min_blackbody = mod[chi_min_index]

# Plot the data overlayed to the points
figure(1)
loglog(v,chi_min_blackbody,label=str('$\chi^2$')+str(' Minimum Column Density=')+str(N[chi_min_index]))
grid(True,which='both')
legend(loc='best')
savefig('SPIRE_averages.png',dpi=300)
close()
