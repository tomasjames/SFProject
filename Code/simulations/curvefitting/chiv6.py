################################################################################
############################ 4th Year Project Code #############################
########################### Chi Squared Routine v5 #############################
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

# Import glob
import glob

############################### Define functions ###############################

# Chi squared
def chi(O,E,sigma):
    '''
    Returns simply the minimised chi squared value for a given dataset.

    O is the flux from the actual data.
    Sigma is the square root of the variance of O
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
    return (N*muh2*mp*opac*(a*(1./(np.exp(b)-1))))/100

################### Read the average data determined earlier ###################

# Initialise a list to store the read in data
data, filenames = [], []

# Initiate counter to help store the passbands
count = -1

# Loop through all files that contain average data
for file in glob.glob('*_average_data.txt'):
    # Open them
    with open(file) as f:
        count += 1

        # Read in using loadtxt and then append the data to a list for use later
        contents = np.loadtxt(file, delimiter=" ")

        # Append the name of the file to the position (count+len(data)) followed by
        # the file contents
        data.insert(count+len(data), file)
        data.append(contents)
        filenames.append(file)

# Assign to variables and put fluxes into an array (data plotted later)
blue = data[data.index(filenames[0])+1]
green = data[data.index(filenames[1])+1]
plw = data[data.index(filenames[2])+1]
pmw = data[data.index(filenames[3])+1]
psw = data[data.index(filenames[4])+1]
red = data[data.index(filenames[5])+1]

# Extract meaningful data from input files
v_data = np.array([psw[:,4],pmw[:,4],plw[:,4],blue[:,4],green[:,4],red[:,4]])
flux = np.array([psw[:,5],pmw[:,5],plw[:,5],blue[:,5],green[:,5],red[:,5]])
flux_error = np.array([psw[:,6],pmw[:,6],plw[:,6],blue[:,6],green[:,6],red[:,6]])

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
lambda_init = 55.670637e-4
lambda_fin = 690.8139e-4
v_init = cc/lambda_init
v_fin = cc/lambda_fin
nlam = 600

# Create array of wavelengths from range given and reshape
w = np.linspace(lambda_init, lambda_fin, nlam)
v = np.linspace(v_init, v_fin, nlam)

# Find the index of the frequency that most closely matches the frequency of the band (i.e. find the difference between the two and find the index at which this is the minimum)
psw_index = min(enumerate(v), key=lambda x: abs(x[1]-psw[0,4]))
pmw_index = min(enumerate(v), key=lambda y: abs(y[1]-pmw[0,4]))
plw_index = min(enumerate(v), key=lambda z: abs(z[1]-plw[0,4]))
blue_index = min(enumerate(v), key=lambda a: abs(a[1]-blue[0,4]))
green_index = min(enumerate(v), key=lambda b: abs(b[1]-green[0,4]))
red_index = min(enumerate(v), key=lambda c: abs(c[1]-red[0,4]))

# Solid angle of the beam
sigma_arc = 1768 # 465 square arcseconds
sigma_beam = sigma_arc*(1/4.24e10) # 465 square arcseconds in steradians

# Ask for spectral index
B = 2.0

# Evaluate opacities
opacity = kappa_0*(v/v_ref)**B

#################### Determine the expected column density #####################

print 'Entering the column density determination script\n'

# Dust to gas ratio
d2g = 0.01

# Read in the dust density information
dust_density = np.loadtxt('.././workingsims/blue/background_15K/dust_density.inp', skiprows=3)
dust_temperature = np.loadtxt('.././workingsims/blue/background_15K/dust_temperature.dat', skiprows=3)

# Create file to store the values
datafeed_store = open('datafeed.txt', 'w')
df = csv.writer(datafeed_store, delimiter=' ')

# Because of the format of the dust_density file, create a list of the indices that correspond to the pixel values. The file is based on a xpix**3 providing xpix is the number of pixels in all 3 dimensions
npix = np.round(len(dust_density)**(1./3))
xpix, ypix, zpix = np.arange(0,npix), np.arange(0,npix), np.arange(0,npix)

# Create list to store values of dust density summed along every x,y coordinate
dust_density_line, T_line, col_full, T_full = [], [], [], []

# Ask for distance to source
#d = input('At what distance is the source? This answer should be in parsecs.\n')
d = 150
D = np.float64(d)*pc # Convert to cm

# The mass of 1 dust grain is simply 1./100th the mass of Hydrogen
dust_mass = muh2*mp*d2g

# Determine solid angle of the pixel
sigma_pix = (imag.sizepix_x*imag.sizepix_y)/D**2

# Loop over the pixels in dust_density
for x in xpix:
    for y in ypix:
        # Reset the dust storage list
        dust_cumulative, T_cumulative = [], []
        for z in zpix:
            # Append the dust storage list with the value of the every dust density along the z axis
            dust_cumulative.append(dust_density[x+y+z])
            T_cumulative.append(dust_temperature[x+y+z])

        # Store the sum
        dust_density_line.append(sum(dust_cumulative))
        T_line.append(np.mean(T_cumulative))

        # The dust density is dust_density_line and so therefore the dust mass in one pixel along the line of sight is dust_density_line*volume
        dust_mass_pixel = (sum(dust_cumulative))*(imag.sizepix_x*imag.sizepix_y*(npix*imag.sizepix_y))

        # Determine the number of dust grains in the pixel
        N_d = (dust_mass_pixel/dust_mass)

        # From Ward-Thompson and Whitworth, column density is the number of dust grains per unit area
        col = N_d/((D**2)*(sigma_pix))

        # Assign all of the writable items to a variable for easier write
        df_towrite = [x, y, col, np.mean(T_cumulative)]

        # Save to a file
        df.writerow(df_towrite)

        col_full.append(col)
        T_full.append(T_line)

N = np.linspace(np.log10(min(col_full)/10), np.log10(max(col_full)*10), 40)
T = np.linspace(8,20,40)

datafeed_store.close()

print 'Column densities have been evaluated and have been saved to the file datafeed_store.txt,\n'

#N = np.linspace(np.round(col,-22),2*np.round(col,-22),10000)
#A = (imag.sizepix_x*imag.sizepix_y) # Area of one pixel in cm

###################### Determine the modified black body #######################

print 'Now determining the modified blackbody curves.\n'

chi_min_index_all, chi_min_blackbody_all = [], []

# Save this to a file for storage
chi_store = open('chi.txt', 'w')
cs = csv.writer(chi_store, delimiter=' ')

for g in range(0,len(flux)):
    for h in range(0,len(flux[g])):

        # Loop through each column density and determine modified black body curve
        mod,chivals,N_index,T_index = [],[],[],[]

        # Loop through both N and T to determine the blackbody
        for i in range(0,len(N)):
            for j in range(0,len(T)):
                blackbody = mbb(10**N[i],dust_mass,opacity,v,T=T[j])

                # Append the value of the given blackbody to a list for storage
                mod.append(blackbody)

                # Append the given values of N and T
                N_index.append(10**N[i])
                T_index.append(T[j])

                #print 'N=',10**N[i],'and T=',T[i]

                # Define lists to store band fluxes
                to_fit_psw_list, to_fit_pmw_list, to_fit_plw_list, to_fit_blue_list, to_fit_green_list, to_fit_red_list = [], [], [], [], [], []

                # Find the flux at the index determined earlier
                to_fit_psw = blackbody[psw_index[0]]
                to_fit_pmw = blackbody[pmw_index[0]]
                to_fit_plw = blackbody[plw_index[0]]
                to_fit_blue = blackbody[blue_index[0]]
                to_fit_green = blackbody[green_index[0]]
                to_fit_red = blackbody[red_index[0]]

                # Put the fitting fluxes into an array
                to_fit = np.array([to_fit_plw, to_fit_pmw, to_fit_psw, to_fit_red, to_fit_green, to_fit_blue])
                #to_fit = np.array([to_fit_psw,to_fit_pmw,to_fit_plw,to_fit_blue])

                # Takes the 6 data points for each pixel and puts them on a list to allow easier assignment
                points = np.array([flux[0][h],flux[1][h],flux[2][h],flux[3][h],flux[4][h],flux[5][h]])
                points_error = np.array([flux_error[0][h],flux_error[1][h],flux_error[2][h],flux_error[3][h],flux_error[4][h],flux_error[5][h]])

                # Append the chi squared value
                chivals.append(chi(to_fit,points,sigma=points_error))
                #chivals.append(chi(to_fit,ps,sigma=[psw[2],pmw[2],plw[2],blue[2]]))

                #print str('Found the chi-squared landscape. Moving to the next values...\n')

        # Determine the chi squared minimum
        chi_min_index = chivals.index(min(chivals))
        chi_min_blackbody = mod[chi_min_index]

        # Append this value to a list to store the values
        chi_min_index_all.append(chi_min_index)
        chi_min_blackbody_all.append(chi_min_blackbody)

        cs_towrite = [g, h, N_index[chi_min_index], T_index[chi_min_index], min(chivals)]

        # Write to file
        cs.writerow(cs_towrite)
        #print 'Writing row to datafile...\n'

chi_store.close()

'''
# Plot the data
figure(1)
errorbar(psw[0],psw[1],yerr=psw[2],fmt='co',label='SPIRE: PSW')
errorbar(pmw[0],pmw[1],yerr=pmw[2],fmt='yo',label='SPIRE: PMW')
errorbar(plw[0],plw[1],yerr=plw[2],fmt='mo',label='SPIRE: PLW')
errorbar(blue[0],blue[1],yerr=blue[2],fmt='bo',label='PACS: Blue')
errorbar(green[0],green[1],yerr=green[2],fmt='go',label='PACS: Green')
errorbar(red[0],red[1],yerr=red[2],fmt='ro',label='PACS: Red')
plot(v,chi_min_blackbody,label=str('$\chi^2$')+str(' Minimum:\n $N$=')+str(N_index[chi_min_index])+str('$cm^{-2}$ \n')+str(' $T$=')+str(np.float(T_index[chi_min_index]))+str('$K$'))
xlabel(r'$\nu \/(Hz)$')
ylabel(r'Intensity $(erg/cm^{2}/s/Hz/ster)$')
xscale("log", nonposx='clip')
yscale("log", nonposx='clip')
grid(True,which='both')
legend(loc='best')
title('The $\chi^{2}$ Minimised Best Fit SED for PACS and SPIRE Bands\n')
savefig('SPIRE_averages_v4.png',dpi=300)
close()

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
