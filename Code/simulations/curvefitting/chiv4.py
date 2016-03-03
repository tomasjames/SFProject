################################################################################
############################ 4th Year Project Code #############################
########################### Chi Squared Routine v4 #############################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

# Import radmc3dPy
import radmc3dPy
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
def mbb(N,dust_mass,opac,v,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity of the source over the beam.
    '''
    #a = (2*(hh)*(cc)**2)/(wav**5)
    #b = ((hh)*(cc))/(wav*(kk)*T)
    a = (2*(hh)*(v**3))/(cc**2)
    b = (hh*v)/(kk*T)
    #return N*dust_mass*opac*(a*(1./(np.exp(b)-1)))/100
    return (N*muh2*mp*opac*(a*(1./(np.exp(b)-1))))/100

################### Read the average data determined earlier ###################

# Reads the data
average_data = np.loadtxt('SPIRE_average_data.txt', delimiter=" ")

# Sort the data
average_data_sort = average_data[np.argsort(average_data[:,0])]

# Assign to variables and put fluxes into an array (data plotted later)
plw,pmw,psw = average_data_sort[0], average_data_sort[1], average_data_sort[2]
ps = np.array([psw[1],pmw[1],plw[1]])

# Read the initial radmc3dPy output to get image dimensions and info
imag = radmc3dPy.image.readImage('../workingsims/plw/background_15K/image.out')

######################### Determine (again) the opacity ########################

# Query for reference wavelength
w_ref = 250e-4
v_ref = cc/w_ref
print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

# Ask for reference opacity
kappa_0 = 4.0
print '\nThe k_0 opacity is taken to be', kappa_0, 'cm^2/g as per Ossenkopf and Henning, 1994. According to http://arxiv.org/pdf/1302.5699v1.pdf the value of B=2.08 fits the power law.\n'

# Define the start and end points of the spectrum
lambda_init = 199.4540e-4
lambda_fin = 690.8139e-4
v_init = cc/lambda_init
v_fin = cc/lambda_fin
nlam = 300

# Create array of wavelengths from range given and reshape
w = np.linspace(lambda_init, lambda_fin, nlam)
v = np.linspace(v_init, v_fin, nlam)

# Find the index of the frequency that most closely matches the frequency of the band (i.e. find the difference between the two and find the index at which this is the minimum)
psw_index = min(enumerate(v), key=lambda x: abs(x[1]-psw[0]))
pmw_index = min(enumerate(v), key=lambda y: abs(y[1]-pmw[0]))
plw_index = min(enumerate(v), key=lambda z: abs(z[1]-plw[0]))

# Solid angle of the beam
sigma_arc = 1768 # 465 square arcseconds
sigma_beam = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

# Ask for spectral index
#B = input('What dust spectral index should be used? (1.7 is recommended)\n')
B = 2.0

# Evaluate opacities
#opacity = kappa_0*(w_ref/w)**B
opacity = kappa_0*(v/v_ref)**B

#################### Determine the expected column density #####################

print 'Entering the column density determination script\n'

# Dust to gas ratio
d2g = 0.01

# Determine the mass of the cloud
m = input('What cloud mass should be used? (Answer should be in Solar Masses)\n')
cloud_dust_mass = d2g*m*ms # 0.01 is the dust to gas ratio
print 'The mass of the cloud considered is', m*ms, 'g. The dust mass is ', cloud_dust_mass, 'g.\n'

# Determine the mass of one dust grain
print 'Presume gas is Hydrogen, therefore dust mass is 1/100 the mass of the gas.\n'
cloud = input('What number density should be assigned to the cloud?\n')
cloud_density = cloud*muh2*mp*d2g
dust_mass = muh2*mp*d2g
print 'The mass of 1 dust grain is', dust_mass, 'g.\n'

# Calculate radius of the cloud in cm
r = ((3./4)*(1./np.pi)*(cloud_dust_mass/cloud_density))**(1./3)
print 'The radius of the cloud is', r, 'cm (or', r/au, ' AU).\n'
V = (4./3)*(np.pi)*(r**3)

# Determine the number of dust grains in the cloud by taking the mass of the cloud and dividing it by the mass of 1 dust grain
N_d = (cloud_dust_mass/dust_mass)/((2*r)/(imag.sizepix_x))**2
#N_d = N_d_cloud/(V/(imag.sizepix_x/2)**3)

# Ask for distance to source
d = input('At what distance is the source? This answer should be in parsecs.\n')
D = np.float64(d)*pc # Convert to cm

# Determine solid angle of the source and pixel
sigma_source = (np.pi*r**2)/D**2
sigma_pix = (imag.sizepix_x*imag.sizepix_y)/D**2

# From Ward-Thompson and Whitworth, column density is the number of dust grains per unit area
col = N_d/((D**2)*(sigma_pix))

print 'This produces a column density of ', col, '/cm**2.'

N = np.linspace(10**21, 10**23, 1000)
T = np.linspace(9,11,200)
#N = np.linspace(np.round(col,-22),2*np.round(col,-22),10000)
#A = (imag.sizepix_x*imag.sizepix_y) # Area of one pixel in cm

###################### Determine the modified black body #######################

# Loop through each column density and determine modified black body curve
mod,chivals,N_index,T_index = [],[],[],[]

# Define lists to store band fluxes
to_fit_psw_list, to_fit_pmw_list, to_fit_plw_list = [], [], []

# Loop through both N and T to determine the blackbody
for i in range(0,len(N)):
    for j in range(0,len(T)):
        blackbody = mbb(N[i],dust_mass,opacity,v,T=T[j])
        mod.append(blackbody)
        N_index.append(N[i])
        T_index.append(T[j])

# Loop through all values of the modified blackbody
for k in range(0,len(mod)):

    # Find the flux at the index determined earlier
    to_fit_psw = mod[k][psw_index[0]]
    to_fit_pmw = mod[k][pmw_index[0]]
    to_fit_plw = mod[k][plw_index[0]]

    # Append these fluxes to the lists
    to_fit_psw_list.append(to_fit_psw)
    to_fit_pmw_list.append(to_fit_pmw)
    to_fit_plw_list.append(to_fit_plw)

    # Put the fitting fluxes into an array
    to_fit = np.array([to_fit_psw,to_fit_pmw,to_fit_plw])

    # Append the chi squared value
    chivals.append(chi(to_fit,ps,sigma=[psw[2],pmw[2],plw[2]]))

# Determine the chi squared minimum
chi_min_index = chivals.index(min(chivals))
chi_min_blackbody = mod[chi_min_index]

# Plot the data
figure(1)
errorbar(psw[0],psw[1],yerr=psw[2],fmt='bo',label='PSW')
errorbar(pmw[0],pmw[1],yerr=pmw[2],fmt='go',label='PMW')
errorbar(plw[0],plw[1],yerr=plw[2],fmt='ro',label='PLW')
plot(v,chi_min_blackbody,label=str('$\chi^2$')+str(' Minimum:\n $N$=')+str(N_index[chi_min_index])+str('$cm^{-2}$ \n')+str(' $T$=')+str(np.float(T_index[chi_min_index]))+str('$K$'))
xlabel(r'$\nu \/(Hz)$')
ylabel(r'Intensity $(erg/cm^{2}/s/Hz/ster)$')
xscale("log", nonposx='clip')
yscale("log", nonposx='clip')
grid(True,which='both')
legend(loc='best')
title('The $\chi^{2}$ Minimised Best Fit SED for PSW, PMW and PLW Bands\n')
savefig('SPIRE_averages_v4.png',dpi=300)
close()

'''
# Plot the figure and save for reference
figure(2)
plot(N,chivals)
grid(True,which='both')
xlabel('$N\,(cm^{-2})$')
ylabel('$\chi^2$')
title(str('\nThe $\chi^2$ Distribution for $N$=')+str(N[0])+str('$cm^{-2}$ to ')+str(N[-1])+str('$cm^{-2}$\n'))
savefig('chisquared_N_v3.png',dpi=300)
close()
'''
