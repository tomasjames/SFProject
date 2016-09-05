################################################################################
############################ 4th Year Project Code #############################
##############################  Functions to run  ##############################
################################################################################

############################# Import statements ################################
# Import RADMC-3D
import radmc3dPy

# Import OS manipulation tools
import os

# Import units
from radmc3dPy.natconst import *

# Import numpy
import numpy as np

# Import csv
import csv

# Import matplotlib
import matplotlib as mpl
mpl.use('Qt4Agg')
from matplotlib.pyplot import *

# Import astropy
from astropy.io import fits

# Import glob for file scraping
import glob

import random

############################ Define Chi squared test ############################

# Chi squared
def chi(O,E,sigma):

    '''
    Returns simply the minimised chi squared value for a given dataset.

    O is the flux from the actual data.
    Sigma is the square root of the variance of O
    E is the data that is expected
    '''

    return sum((((O-E)/sigma)**2))

######################### Define a modified blackbody ############################

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
    return (N*muh2*mp*opac*(a*(1./(np.exp(b)-1))))/(10**(-23))

def dataDerive(data_type, region, kappa_0, lambda_0, B, d, npix, sizeau):

    '''
    Function to take input quantities for RADMC-3D and derive from them 'true' values of the dust column density and dust temperature

    Keywords:
    data_type: the type of data being used
        Options:
            'radmc'
            'arepo'
            'herschel_snaps'
    '''

    #################### Determine the expected column density #####################

    B_val = 'B={}'.format(B)

    print 'Entering the column density determination script\n'

    # Dust to gas ratio
    d2g = 0.01

    if data_type == 'herschel_snaps':

        # Define the first folder
        first_folder = 'herschel_snaps_psf'
        second_folder = region

        dust_density = np.loadtxt(('/export/home/c1158976/Code/simulations/{}/{}/{}/psw/dust_density.inp').format(first_folder,second_folder,B_val), skiprows=3)
        dust_temperature = np.loadtxt(('/export/home/c1158976/Code/simulations/{}/{}/{}/psw/dust_temperature.dat').format(first_folder,second_folder,B_val), skiprows=3)
        imag = fits.open(('/export/home/c1158976/Code/simulations/{}/{}/{}/psw/psw_common_convolved.fits').format(first_folder,second_folder,B_val))

    elif data_type == 'radmc':
        dust_density = np.loadtxt(('../../../workingsims_psf/{}/psw/background_15K/dust_density.inp').format(B_val), skiprows=3)
        dust_temperature = np.loadtxt(('../../../workingsims_psf/{}/psw/background_15K/dust_temperature.dat').format(B_val), skiprows=3)
        imag = fits.open(('../../../workingsims_psf/{}/psw/background_15K/psw_common_convolved.fits').format(B_val))

    elif data_type == 'arepo':
        dust_density = np.loadtxt(('../../../data_psf/{}/psw/dust_project/dust_density.inp').format(B_val), skiprows=3)
        dust_temperature = np.loadtxt(('../../../data_psf/{}/psw/dust_project/dust_temperature.dat').format(B_val), skiprows=3)
        imag = fits.open(('../../../data_psf/{}/psw/dust_project/psw_common_convolved.fits').format(B_val))

    '''
    if data_type == 'radmc':
        # Read in the dust density information
        #dust_density = np.loadtxt(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/dust_density.inp'), skiprows=3)
        #dust_temperature = np.loadtxt(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/dust_temperature.dat'), skiprows=3)
        #imag = radmc3dPy.image.readImage(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/image.out'))
        #imag = fits.open(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/')+str(filter)+str('_common_convolved.fits'))
        dust_density = np.loadtxt(('../../../workingsims_psf/{}/psw/background_15K/dust_density.inp').format(B_val), skiprows=3)
        dust_temperature = np.loadtxt(('../../../workingsims_psf/{}/psw/background_15K/dust_temperature.dat').format(B_val), skiprows=3)
        imag = fits.open(('../../../workingsims_psf/{}/psw/background_15K/psw_common_convolved.fits').format(B_val))

    elif data_type == 'arepo':
        #dust_density = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_density.inp'), skiprows=3)
        #dust_temperature = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_temperature.dat'), skiprows=3)
        #imag = radmc3dPy.image.readImage(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/image.out'))
        #imag = fits.open(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/')+str(filter)+str('_common_convolved.fits'))
        dust_density = np.loadtxt(('../../../data_psf/{}/psw/dust_project/dust_density.inp').format(B_val), skiprows=3)
        dust_temperature = np.loadtxt(('../../../data_psf/{}/psw/dust_project/dust_temperature.dat').format(B_val), skiprows=3)
        imag = fits.open(('../../../data_psf/{}/psw/dust_project/psw_common_convolved.fits').format(B_val))

    elif data_type == 'herschel_snaps':
        #dust_density = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_density.inp'), skiprows=3)
        #dust_temperature = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_temperature.dat'), skiprows=3)
        #imag = radmc3dPy.image.readImage(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/image.out'))
        #imag = fits.open(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/')+str(filter)+str('_common_convolved.fits'))
        dust_density = np.loadtxt(('../../../../herschel_snaps_psf/x1_comp_region1/{}/psw/dust_density.inp').format(B_val), skiprows=3)
        dust_temperature = np.loadtxt(('../../../../herschel_snaps_psf/x1_comp_region1/{}/psw/dust_temperature.dat').format(B_val), skiprows=3)
        imag = fits.open(('../../../../herschel_snaps_psf/x1_comp_region1/{}/psw/psw_common_convolved.fits').format(B_val))
    '''

    # Create file to store the values
    datafeed_store = open('datafeed.txt', 'w')
    df = csv.writer(datafeed_store, delimiter=' ')

    # Determine data dimensions
    nx = len(imag[0].data[0])
    ny = len(imag[0].data[:,0])

    # Because of the format of the dust_density file, create a list of the indices that correspond to the pixel values. The file is based on a xpix**3 providing xpix is the number of pixels in all 3 dimensions
    xpix, ypix, zpix = np.arange(0,nx), np.arange(0,ny), np.arange(0,(len(dust_density)/(nx*ny)))

    # Create list to store values of dust density summed along every x,y coordinate
    dust_density_line, T_line, col_full, T_full = [], [], [], []

    # Ask for distance to source
    #d = input('At what distance is the source? This answer should be in parsecs.\n')
    D = np.float64(d)*pc # Convert to cm

    # The mass of 1 dust grain is simply 1./100th the mass of Hydrogen
    dust_mass = muh2*mp*d2g

    # Instantiate a counter and a number of lists to store values (and for diagnostics)
    count = 0
    #store_density, store_temp = [], []

    # Loop over the 2D image square cube
    for i in range(0,len(xpix)*len(ypix)):

        # Add to the counter
        count += 1

        # Finds the locations of every each pixel's Z values
        loc = np.arange(i,len(zpix)*(len(xpix)*len(ypix)),len(xpix)*len(ypix))

        # Sums these values (i.e. sums along the line of sight) and stores the locations
        #store_density.append(sum(dust_density[loc]))

        # The dust density is dust_density along the line of sight and so therefore the dust mass in one pixel along the line of sight is dust_density_line*volume
        #dust_mass_pixel = (sum(dust_density[loc]))*(imag.sizepix_x*imag.sizepix_y)*(len(dust_density[loc])*imag.sizepix_y)
        #dust_mass_cell = sum(((dust_density[loc]))*(imag.sizepix_x*imag.sizepix_y*imag.sizepix_y))

        # Account for the column-weighted temperature
        col_T = np.sum((dust_temperature[loc]*dust_density[loc]))/(np.sum(dust_density[loc]))

        # Repeat similar procedure for the temperature
        #store_temp.append(col_T)

        # Determine the number of dust grains in the pixel
        #N_d = (dust_mass_pixel/(dust_mass*len(dust_density[loc])))
        #N_d = dust_mass_pixel/dust_mass

        # From Ward-Thompson and Whitworth, column density is the number of dust grains per unit area
        #col = N_d/((D**2)*(sigma_pix))
        pixwidth = (sizeau*au)/npix
        col = np.sum((dust_density[loc]*pixwidth)/(muh2*mp))
        #col = ((dust_mass_cell/(muh2*mp*imag.sizepix_x*imag.sizepix_y)))

        # Assign all of the writable items to a variable for easier write
        df_towrite = [i, col, col_T]

        # Save to a file
        df.writerow(df_towrite)

        col_full.append(col)
        T_full.append(col_T)

    datafeed_store.close()

    print 'Column densities have been evaluated and have been saved to the file datafeed.txt,\n'

    return col_full, T_full

########################## Define wrapper to apply to data ########################

def chiTest(data_type, region, output_name, N, T, kappa_0, lambda_0, B, d):

    '''
    A function to run the Chi-squared minimisation analysis over the SEDs generated.

    Keywords
    data_type: whether the data being considered is RADMC-3D data, or Arepo data
        Options:
            'radmc': spherical, isothermal cloud SEDs
            'arepo': the arepo simulation SEDs

    output_name: the name of the output file to store the values
        Options:
            'chi_coarse.txt' or 'chi_fine.txt'

    N: an array of specific length containing all of the values of dust column density to be used in the sim
        Values MUST be logarithmic, i.e. 1x10**17 would be 17
        Options:
            N = np.array([17,18,19,20])

    T: an array of specific length containing all of the values of dust temperature to be used in the sim
        Values NEED NOT be logarithmic
        Options:
            T = np.array([10,11,12])
    '''

    ################### Read the average data determined earlier ###################

    # Instantiate the B val
    B_val = 'B={}'.format(B)

    '''
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

            # Append the file contents
            data.append(contents)
    '''

    # Assign to variables and put fluxes into an array (data plotted later)
    # +1 term avoids the filename that is placed in the array for identification
    blue = np.loadtxt('blue_average_data.txt', delimiter=" ")
    green = np.loadtxt('green_average_data.txt', delimiter=" ")
    plw = np.loadtxt('plw_average_data.txt', delimiter=" ")
    pmw = np.loadtxt('pmw_average_data.txt', delimiter=" ")
    psw = np.loadtxt('psw_average_data.txt', delimiter=" ")
    red = np.loadtxt('red_average_data.txt', delimiter=" ")

    # Extract meaningful data from input files
    #v_data = np.array([psw[:,1],pmw[:,1],plw[:,1],blue[:,1],green[:,1],red[:,1]])
    flux = np.array([psw[:,2],pmw[:,2],plw[:,2],blue[:,2],green[:,2],red[:,2]])
    flux_error = np.array([psw[:,3],pmw[:,3],plw[:,3],blue[:,3],green[:,3],red[:,3]])

    # Pull out the discrete frequencies used (band[:,1] pulls out all of the frequency data)
    v = np.unique(np.array([psw[:,1],pmw[:,1],plw[:,1],blue[:,1],green[:,1],red[:,1]]))
    #v = np.unique(np.array([psw[:,1],pmw[:,1],plw[:,1]]))

    # Read the initial radmc3dPy output to get image dimensions and info
    #imag = radmc3dPy.image.readImage('../../data/blue/dust_project/image.out')

    ######################### Determine (again) the opacity ########################

    # Query for reference wavelength
    #w_ref = 250e-4
    w_ref = lambda_0*(1e-4)
    v_ref = cc/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

    # Ask for reference opacity
    #kappa_0 = 4.0
    print '\nThe k_0 opacity is taken to be', kappa_0, 'cm^2/g as per Ossenkopf and Henning, 1994. According to http://arxiv.org/pdf/1302.5699v1.pdf the value of B=2.08 fits the power law.\n'

    # Evaluate opacities
    #opacity = kappa_0*(v/v_ref)**B
    opacity = kappa_0*(v/v_ref)**B

    '''
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
    psw_index = min(enumerate(v), key=lambda x: abs(x[1]-psw[0][1]))
    pmw_index = min(enumerate(v), key=lambda y: abs(y[1]-pmw[0,1]))
    plw_index = min(enumerate(v), key=lambda z: abs(z[1]-plw[0,1]))
    blue_index = min(enumerate(v), key=lambda a: abs(a[1]-blue[0,1]))
    green_index = min(enumerate(v), key=lambda b: abs(b[1]-green[0,1]))
    red_index = min(enumerate(v), key=lambda c: abs(c[1]-red[0,1]))
    '''

    ############################# Read in necessary values ##########################

    # Define the dust to gas ratio
    d2g = 0.01

    if data_type == 'herschel_snaps':

        # Define the first folder
        first_folder = 'herschel_snaps_psf'
        second_folder = region

        imag = fits.open(('/export/home/c1158976/Code/simulations/{}/{}/B=2.0/psw/psw_common_convolved.fits').format(first_folder,second_folder))

    elif data_type == 'radmc':
        imag = fits.open(('../../../workingsims_psf/B=2.0/psw/background_15K/psw_common_convolved.fits'))

    elif data_type == 'arepo':
        imag = fits.open(('../../../data_psf/B=2.0/psw/dust_project/psw_common_convolved.fits'))

    # Determine data dimensions
    nx = len(imag[0].data[0])
    ny = len(imag[0].data[:,0])

    # Ask for distance to source
    #d = input('At what distance is the source? This answer should be in parsecs.\n')
    #d = 300
    D = np.float64(d)*pc # Convert to cm

    # The mass of 1 dust grain is simply 1./100th the mass of Hydrogen
    dust_mass = muh2*mp*d2g

    ###################### Determine the modified black body #######################

    print 'Now determining the modified blackbody curves.\n'

    # Create 2 2D arrays of the data to track the progress of the loop
    T_mesh, N_mesh = np.meshgrid(T,N)

    chi_min_index_all, chi_min_blackbody_all = [], []

    # Save this to a file for storage
    chi_store = open(str(output_name), 'w')
    cs = csv.writer(chi_store, delimiter=' ')

    # Save a line to allow better understanding of the columns
    cs_towrite = ['Index', 'Column Density', 'Temperature', 'Minimised Chi-Squared', 'One-Sigma N', 'One-Sigma T', 'Two-Sigma N', 'Two-Sigma T', 'Three-Sigma N', 'Three-Sigma T']
    cs.writerow(cs_towrite)

    # Loop through all pixels in the image
    for h in range(0,nx*ny):

        # Loop through each column density and determine modified black body curve
        mod,chivals,N_index,T_index = [],[],[],[]

        # Takes the 6 data points for each pixel and puts them on a list to allow easier assignment
        #points = np.array([flux[2][h],flux[1][h],flux[0][h]])
        #points_error = np.array([flux_error[2][h],flux_error[1][h],flux_error[0][h]])
        points = np.array([flux[2][h],flux[1][h],flux[0][h],flux[5][h],flux[4][h],flux[3][h]])
        points_error = np.array([flux_error[2][h],flux_error[1][h],flux_error[0][h],flux_error[5][h],flux_error[4][h],flux_error[3][h]])

        # Loop through both N and T to determine the blackbody
        for i in range(0,len(N)):
            for j in range(0,len(T)):
                blackbody = mbb(10**N[i],dust_mass,opacity,v,T=T[j])

                # Append the value of the given blackbody to a list for storage
                mod.append(blackbody)

                # Append the given values of N and T
                N_index.append(N[i])
                T_index.append(T[j])

                # Append the chi squared value
                #chisquared = chi(points,blackbody,sigma=points_error)
                chisquared = chi(points,blackbody,sigma=points_error)

                # This accounts for any infinite values of the Chi-squared coefficient
                if chisquared == np.inf:
                    chivals.append(np.nan)
                else:
                    chivals.append(chisquared)

        # Determine the chi squared minimum
        chi_min_index = chivals.index(min(chivals))
        chi_min_blackbody = mod[chi_min_index]

        ############## Error analysis begins here #############

        # Reshape the chivals to fit the format of N_mesh and T_mesh
        chivals_mesh = np.reshape(chivals, (len(T),len(N)))

        # Determine error in the value by finding location of points that lie within 3 sigma
        # Start by finding the location in the mesh of the minimum chi-squared value
        val = zip(*np.where(chivals_mesh == np.amin(chivals_mesh)))[0]

        # According to Data Analysis notes 2D distribution has the following values of delta-chi-squared
        delta_1 = 2.3 #(one_sigma)
        delta_2 = 6.17 #(two_sigma)
        delta_3 = 11.8 #(three_sigma)

        # Find the locations of the points within the one sigma bound
        one_sigma = zip(*np.where((chivals_mesh >= chivals_mesh[val]-delta_1/2) & (chivals_mesh <= chivals_mesh[val]+delta_1/2)))

        T_temp, N_temp = [], []
        T_error, N_error = [], []

        # Pull the values of N and T at these points
        for a in one_sigma:
            T_temp.append(T_mesh[a])
            N_temp.append(N_mesh[a])

        T_error.append(max(T_temp) - min(T_temp))
        N_error.append(max(N_temp) - min(N_temp))

        # Append this value to a list to store the values
        chi_min_index_all.append(chi_min_index)
        chi_min_blackbody_all.append(chi_min_blackbody)

        '''
        print '============== N ================'
        print 'Coarse N:', N_index[chi_min_index]
        print 'One sigma N:', min(N_temp), ' to ', max(N_temp)

        print '============== T ================'
        print 'Coarse T:', T_index[chi_min_index]
        print 'One sigma T:', min(T_temp), ' to ', max(T_temp)
        '''

        # From here, repeat the Chi-squared analyses over the one sigma range by first defining the new range to look over
        #N_refined = np.linspace(N_index[chi_min_index]-(max(N_temp) -min(N_temp)),N_index[chi_min_index]+(max(N_temp) - min(N_temp)),100)
        #T_refined = np.linspace(T_index[chi_min_index]-(max(T_temp) - min(T_temp)),T_index[chi_min_index]+(max(T_temp) - min(T_temp)),100)

        N_refined = np.linspace(min(N_temp),max(N_temp),100)
        T_refined = np.linspace(min(T_temp),max(T_temp),100)

        # Loop through each column density and determine modified black body curve
        mod_refined,chivals_refined,N_index_refined,T_index_refined = [],[],[],[]

        # Loop through both N and T to determine the blackbody
        for i in range(0,len(N_refined)):
            for j in range(0,len(T_refined)):
                blackbody_refined = mbb(10**N_refined[i],dust_mass,opacity,v,T=T_refined[j])

                # Append the value of the given blackbody to a list for storage
                mod_refined.append(blackbody_refined)

                # Append the given values of N and T
                N_index_refined.append(10**N_refined[i])
                T_index_refined.append(T_refined[j])

                # Append the chi squared value (using the same points and points_error determined earlier)
                #chisquared_refined = chi(points,blackbody_refined,sigma=points_error)
                chisquared_refined = chi(points,blackbody_refined,sigma=points_error)

                # This accounts for any infinite values of the Chi-squared coefficient
                if chisquared_refined == np.inf:
                    chivals_refined.append(np.nan)
                else:
                    chivals_refined.append(chisquared_refined)

        # Determine the chi squared minimum
        chi_min_index_refined = chivals_refined.index(min(chivals_refined))
        chi_min_blackbody_refined = mod_refined[chi_min_index_refined]

        ############## Error analysis begins here #############

        # Mesh N and T for error determination
        T_mesh_refined, N_mesh_refined = np.meshgrid(T_refined,N_refined)

        # Reshape the chivals to fit the format of N_mesh and T_mesh
        chivals_mesh_refined = np.reshape(chivals_refined, (len(T_refined),len(N_refined)))

        # Determine error in the value by finding location of points that lie within 3 sigma
        # Start by finding the location in the mesh of the minimum chi-squared value
        val_refine = zip(*np.where(chivals_mesh_refined == np.amin(chivals_mesh_refined)))[0]

        # Find the locations of the points within these bounds
        one_sigma_refined = zip(*np.where((chivals_mesh_refined >= chivals_mesh_refined[val_refine]-delta_1) & (chivals_mesh_refined <= chivals_mesh_refined[val_refine]+delta_1)))
        two_sigma_refined = zip(*np.where((chivals_mesh_refined >= chivals_mesh_refined[val_refine]-delta_2) & (chivals_mesh_refined <= chivals_mesh_refined[val_refine]+delta_2)))
        three_sigma_refined = zip(*np.where((chivals_mesh_refined >= chivals_mesh_refined[val_refine]-delta_3) & (chivals_mesh_refined <= chivals_mesh_refined[val_refine]+delta_3)))

        T_temp_refined, N_temp_refined = [], []
        T_error_refined, N_error_refined = [], []

        # Pull the values of N and T at these points
        for b in [one_sigma_refined, two_sigma_refined, three_sigma_refined]:
            for a in b:
                T_temp_refined.append(T_mesh_refined[a])
                N_temp_refined.append(10**N_mesh_refined[a])

            T_error_refined.append(max(T_temp_refined) - min(T_temp_refined))
            N_error_refined.append(max(N_temp_refined) - min(N_temp_refined))

        # Append this value to a list to store the values
        #chi_min_index_all.append(chi_min_index)
        #chi_min_blackbody_all.append(chi_min_blackbody)

        cs_towrite = [h, N_index_refined[chi_min_index_refined], T_index_refined[chi_min_index_refined], min(chivals_refined), N_error_refined[0], T_error_refined[0], N_error_refined[1], T_error_refined[1], N_error_refined[2], T_error_refined[2]]

        # Write to file
        cs.writerow(cs_towrite)
        #print 'Writing row to datafile...\n'

        if h == 0:
        # Plot the data
            figure(1)
            errorbar(psw[0][1],psw[0][2],yerr=psw[0][-1],fmt='co',label='SPIRE: PSW')
            errorbar(pmw[0][1],pmw[0][2],yerr=pmw[0][-1],fmt='yo',label='SPIRE: PMW')
            errorbar(plw[0][1],plw[0][2],yerr=plw[0][-1],fmt='mo',label='SPIRE: PLW')
            errorbar(blue[0][1],blue[0][2],yerr=blue[0][-1],fmt='bx',label='PACS: Blue')
            errorbar(green[0][1],green[0][2],yerr=green[0][-1],fmt='gx',label='PACS: Green')
            errorbar(red[0][1],red[0][2],yerr=red[0][-1],fmt='rx',label='PACS: Red')
            plot(v,chi_min_blackbody,label=str('$\chi^2$')+str(' Minimum:\n $N$=')+str(N_index[chi_min_index])+str('+/-')+str(N_error[0])+str('$cm^{-2}$ \n')+str(' $T$=')+str(np.float(T_index[chi_min_index]))+str('+/-')+str(T_error[0])+str('$\/K$'))
            xlabel(r'$\nu \/(Hz)$')
            #ylabel(r'Intensity $(erg/cm^{2}/s/Hz/ster)$')
            ylabel(r'Intensity $[MJy/ster]$')
            xscale("log", nonposx='clip')
            yscale("log", nonposx='clip')
            grid(True,which='both')
            legend(loc='best',prop={'size':8})
            title(str('The $\chi^{2}$ Minimised Best Fit SED for PACS and SPIRE Bands\n at the ')+str(h)+str('th Pixel\n'))
            savefig(str('averages_coarse_')+str(h)+str('.png'),dpi=300)
            close()

        # Also plot probability contour
            figure(2)
            #plot(10**N_mesh,T_mesh, 'b.')
            plot(N_index[chi_min_index], T_index[chi_min_index], 'ro', label=str(r'$\chi^{2}$ Minimum=')+str(min(chivals)))
            CS = contour(10**N_mesh,T_mesh,chivals_mesh)
            clabel(CS)
            xlabel(r'$N\/(cm^{-2})$')
            ylabel(r'$T\/(K)$')
            title(r'$\chi^{2} Contours$')
            legend(loc='best',prop={'size':8})
            savefig(str('contours_')+str(h)+str('.png'),dpi=300)
            close()

        # And the chi-squared landscape
            figure(3)
            subplot(2,1,1)
            plot(T_index, chivals, 'b.')
            plot(T_index[chi_min_index], chivals[chi_min_index], 'ro')
            xlabel(r'$T\/(K)$')
            ylabel(r'$\chi^{2}$')

            subplot(2,1,2)
            plot(N_index, chivals, 'b.')
            plot(N_index[chi_min_index], chivals[chi_min_index], 'ro')
            xlabel(r'$N\/(cm^{-2})$')
            ylabel(r'$\chi^{2}$')

            savefig(str('landscape_')+str(h)+str('.png'),dpi=300)
            close()

    chi_store.close()
