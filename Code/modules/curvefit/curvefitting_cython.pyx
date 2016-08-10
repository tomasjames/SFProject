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
cimport numpy as np

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
def chi(double O, double E, double sigma):

    '''
    Returns simply the minimised chi squared value for a given dataset.

    O is the flux from the actual data.
    Sigma is the square root of the variance of O
    E is the data that is expected
    '''

    return sum((((O-E)/sigma)**2))

######################### Define a modified blackbody ############################

# Modified blackbody
def mbb(double N, double dust_mass, double opac,double v, double T):

    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity of the source over the beam.
    '''

    # Define variables in Cython/C syntax
    cdef double a
    cdef double b

    a = (2*(hh)*(v**3))/(cc**2)
    b = (hh*v)/(kk*T)

    return (N*muh2*mp*opac*(a*(1./(np.exp(b)-1))))/(10**(-23))

def dataDerive(data_type, double kappa_0, double lambda_0, double B):

    '''
    Function to take input quantities for RADMC-3D and derive from them 'true' values of the dust column density and dust temperature

    Keywords:
    data_type: the type of data being used
        Options:
            'radmc' or 'arepo'
    '''

    ########################## Declare variables for Cython ########################
    '''
    cdef double d2g
    cdef double d
    cdef double D
    cdef double dust_mass
    cdef int count
    '''
    #################### Determine the expected column density #####################

    B_val = str('B=')+str(B)

    print 'Entering the column density determination script\n'

    # Dust to gas ratio
    d2g = 0.01

    if data_type == 'radmc':
        # Read in the dust density information
        dust_density = np.loadtxt(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/dust_density.inp'), skiprows=3)
        dust_temperature = np.loadtxt(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/dust_temperature.dat'), skiprows=3)
        #imag = radmc3dPy.image.readImage(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/image.out'))
        imag = fits.open(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/')+str(filter)+str('_common_convolved.fits'))

    elif data_type == 'arepo':
        dust_density = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_density.inp'), skiprows=3)
        dust_temperature = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_temperature.dat'), skiprows=3)
        #imag = radmc3dPy.image.readImage(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/image.out'))
        imag = fits.open(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/')+str(filter)+str('_common_convolved.fits'))

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
    d = 300
    D = np.float64(d)*pc # Convert to cm

    # The mass of 1 dust grain is simply 1./100th the mass of Hydrogen
    dust_mass = muh2*mp*d2g

    # Instantiate a counter and a number of lists to store values (and for diagnostics)
    count = 0
    store_density, store_loc, store_temp = [], [], []

    # Loop over the 2D image square cube
    for i in range(0,len(xpix)*len(ypix)):

        # Add to the counter
        count += 1

        # Finds the locations of every each pixel's Z values
        loc = np.arange(i,len(zpix)*(len(xpix)*len(ypix)),len(xpix)*len(ypix))

        # Sums these values (i.e. sums along the line of sight) and stores the locations
        store_density.append(sum(dust_density[loc]))
        store_loc.append(loc)

        # The dust density is dust_density along the line of sight and so therefore the dust mass in one pixel along the line of sight is dust_density_line*volume
        #dust_mass_pixel = (sum(dust_density[loc]))*(imag.sizepix_x*imag.sizepix_y)*(len(dust_density[loc])*imag.sizepix_y)
        #dust_mass_cell = sum(((dust_density[loc]))*(imag.sizepix_x*imag.sizepix_y*imag.sizepix_y))

        # Account for the column-weighted temperature
        col_T = np.sum((dust_temperature[loc]*dust_density[loc]))/(np.sum(dust_density[loc]))

        # Repeat similar procedure for the temperature
        store_temp.append(col_T)

        # Determine the number of dust grains in the pixel
        #N_d = (dust_mass_pixel/(dust_mass*len(dust_density[loc])))
        #N_d = dust_mass_pixel/dust_mass

        # From Ward-Thompson and Whitworth, column density is the number of dust grains per unit area
        #col = N_d/((D**2)*(sigma_pix))
        col = np.sum((dust_density[loc]*imag.sizepix_y)/(muh2*mp))
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

def chiTest(data_type, output_name, double N, double T, double kappa_0, double lambda_0, double B):

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
    B_val = 'B = %s' % (B)

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
    # +1 term avoids the filename that is placed in the array for identification
    blue = data[data.index(filenames[0])+1]
    green = data[data.index(filenames[1])+1]
    plw = data[data.index(filenames[2])+1]
    pmw = data[data.index(filenames[3])+1]
    psw = data[data.index(filenames[4])+1]
    red = data[data.index(filenames[5])+1]

    # Extract meaningful data from input files
    v_data = np.array([psw[:,1],pmw[:,1],plw[:,1],blue[:,1],green[:,1],red[:,1]])
    flux = np.array([psw[:,2],pmw[:,2],plw[:,2],blue[:,2],green[:,2],red[:,2]])
    flux_error = np.array([psw[:,3],pmw[:,3],plw[:,3],blue[:,3],green[:,3],red[:,3]])

    # Read the initial radmc3dPy output to get image dimensions and info
    #imag = radmc3dPy.image.readImage('../../data/blue/dust_project/image.out')

    ######################### Determine (again) the opacity ########################

    # Query for reference wavelength
    #w_ref = 250e-4
    w_ref = lambda_0
    v_ref = cc/w_ref
    print 'The reference wavelength is taken to be', w_ref, 'cm which equates to a frequency of', v_ref, 'Hz.\n'

    # Ask for reference opacity
    #kappa_0 = 4.0
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
    psw_index = min(enumerate(v), key=lambda x: abs(x[1]-psw[0][1]))
    pmw_index = min(enumerate(v), key=lambda y: abs(y[1]-pmw[0,1]))
    plw_index = min(enumerate(v), key=lambda z: abs(z[1]-plw[0,1]))
    blue_index = min(enumerate(v), key=lambda a: abs(a[1]-blue[0,1]))
    green_index = min(enumerate(v), key=lambda b: abs(b[1]-green[0,1]))
    red_index = min(enumerate(v), key=lambda c: abs(c[1]-red[0,1]))

    # Evaluate opacities
    opacity = kappa_0*(v/v_ref)**B

    ############################# Read in necessary values ##########################

    # Dust to gas ratio
    d2g = 0.01

    if data_type == 'radmc':
      # Read in the dust density information
      #dust_density = np.loadtxt(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/dust_density.inp'), skiprows=3)
      #dust_temperature = np.loadtxt(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/dust_temperature.dat'), skiprows=3)
      #imag = radmc3dPy.image.readImage(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/image.out'))
      #imag = fits.open(('../../../workingsims_psf/')+str(B_val)+str('/psw/background_15K/')+str(filter)+str('_common_convolved.fits'))
      dust_density = np.loadtxt(('../../../workingsims_psf/ %s /psw/background_15K/dust_density.inp' % (B_val)), skiprows=3)
      dust_temperature = np.loadtxt(('../../../workingsims_psf/ %s /psw/background_15K/dust_temperature.dat' % (B_val)), skiprows=3)
      imag = fits.open(('../../../workingsims_psf/ %s /psw/background_15K/ %s _common_convolved.fits' % (B_val, filter)))

    elif data_type == 'arepo':
      dust_density = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_density.inp'), skiprows=3)
      dust_temperature = np.loadtxt(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/dust_temperature.dat'), skiprows=3)
      #imag = radmc3dPy.image.readImage(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/image.out'))
      imag = fits.open(('../../../data_psf/')+str(B_val)+str('/psw/dust_project/')+str(filter)+str('_common_convolved.fits'))

    # Create file to store the values
    datafeed_store = open('datafeed.txt', 'w')
    df = csv.writer(datafeed_store, delimiter=' ')

    # Determine data dimensions
    nx = len(imag[0].data[0])
    ny = len(imag[0].data[:,0])

    # Because of the format of the dust_density file, create a list of the indices that correspond to the pixel values. The file is based on a xpix**3 providing xpix is the number of pixels in all 3 dimensions
    xpix, ypix, zpix = np.arange(0,nx), np.arange(0,ny), np.arange(0,(len(dust_density)/(nx*ny)))

    # Ask for distance to source
    #d = input('At what distance is the source? This answer should be in parsecs.\n')
    d = 300
    D = np.float64(d)*pc # Convert to cm

    # The mass of 1 dust grain is simply 1./100th the mass of Hydrogen
    dust_mass = muh2*mp*d2g

    ###################### Determine the modified black body #######################

    print 'Now determining the modified blackbody curves.\n'

    chi_min_index_all, chi_min_blackbody_all = [], []

    # Save this to a file for storage
    chi_store = open(str(output_name), 'w')
    cs = csv.writer(chi_store, delimiter=' ')

    # Save a line to allow better understanding of the columns
    cs_towrite = ['Index', 'Column Density', 'Temperature', 'Minimised Chi-Squared', 'One-Sigma N', 'One-Sigma T', 'Two-Sigma N', 'Two-Sigma T', 'Three-Sigma N', 'Three-Sigma T']
    cs.writerow(cs_towrite)

    for h in range(0,imag.nx*imag.ny):

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
                #to_fit = np.array([to_fit_plw, to_fit_pmw, to_fit_psw])

                # Takes the 6 data points for each pixel and puts them on a list to allow easier assignment
                #points = np.array([flux[2][h],flux[1][h],flux[0][h]])
                #points_error = np.array([flux_error[2][h],flux_error[1][h],flux_error[0][h]])
                points = np.array([flux[2][h],flux[1][h],flux[0][h],flux[5][h],flux[4][h],flux[3][h]])
                points_error = np.array([flux_error[2][h],flux_error[1][h],flux_error[0][h],flux_error[5][h],flux_error[4][h],flux_error[3][h]])

                # Append the chi squared value
                chisquared = chi(to_fit,points,sigma=points_error)
                #chisquared = chi(to_fit,points)

                if chisquared == np.inf:
                    chivals.append(np.nan)
                else:
                    chivals.append(chisquared)
                #chivals.append(chi(to_fit,ps,sigma=[psw[2],pmw[2],plw[2],blue[2]]))

                #print str('Found the chi-squared landscape. Moving to the next values...\n')

        if h == imag.nx*imag.ny/10:
            print '\r[=>         ] 10%'
        elif h == 2*(imag.nx*imag.ny)/10:
            print '\r[==>        ] 20%'
        elif h == 3*(imag.nx*imag.ny)/10:
            print '\r[===>       ] 30%'
        elif h == 4*(imag.nx*imag.ny)/10:
            print '\r[====>      ] 40%'
        elif h == 5*(imag.nx*imag.ny)/10:
            print '\r[=====>     ] 50%'
        elif h == 6*(imag.nx*imag.ny)/10:
            print '\r[======>    ] 60%'
        elif h == 7*(imag.nx*imag.ny)/10:
            print '\r[=======>   ] 70%'
        elif h == 8*(imag.nx*imag.ny)/10:
            print '\r[========>  ] 80%'
        elif h == 9*(imag.nx*imag.ny)/10:
            print '\r[=========> ] 90%'
        elif h == (imag.nx*imag.ny)-1:
            print '\r[==========>] 100%'

        # Determine the chi squared minimum
        chi_min_index = chivals.index(min(chivals))
        chi_min_blackbody = mod[chi_min_index]

        ############## Error analysis begins here #############

        # Create 2 2D arrays of the data to track the progress of the loop
        T_mesh, N_mesh = np.meshgrid(T,N)

        # Reshape the chivals to fit the format of N_mesh and T_mesh
        chivals_mesh = np.reshape(chivals, (len(T),len(N)))

        # Determine error in the value by finding location of points that lie within 3 sigma
        # Start by finding the location in the mesh of the minimum chi-squared value
        val = zip(*np.where(chivals_mesh == np.amin(chivals_mesh)))

        # According to Data Analysis notes 2D distribution has the following values of delta-chi-squared
        delta_1 = 2.3 #(one_sigma)
        delta_2 = 6.17 #(two_sigma)
        delta_3 = 11.8 #(three_sigma)

        # Find the locations of the points within these bounds
        one_sigma = zip(*np.where((chivals_mesh >= chivals_mesh[val[0]]-delta_1) & (chivals_mesh <= chivals_mesh[val[0]]+delta_1)))
        two_sigma = zip(*np.where((chivals_mesh >= chivals_mesh[val[0]]-delta_2) & (chivals_mesh <= chivals_mesh[val[0]]+delta_2)))
        three_sigma = zip(*np.where((chivals_mesh >= chivals_mesh[val[0]]-delta_3) & (chivals_mesh <= chivals_mesh[val[0]]+delta_3)))

        T_temp, N_temp = [], []
        T_error, N_error = [], []

        # Pull the values of N and T at these points
        for b in [one_sigma, two_sigma, three_sigma]:
            for a in b:
                T_temp.append(T_mesh[a])
                N_temp.append(10**N_mesh[a])

            T_error.append(max(T_temp) - min(T_temp))
            N_error.append(max(N_temp) - min(N_temp))

        # Append this value to a list to store the values
        chi_min_index_all.append(chi_min_index)
        chi_min_blackbody_all.append(chi_min_blackbody)

        cs_towrite = [h, N_index[chi_min_index], T_index[chi_min_index], min(chivals), N_error[0], T_error[0], N_error[1], T_error[1], N_error[2], T_error[2]]

        # Write to file
        cs.writerow(cs_towrite)
        #print 'Writing row to datafile...\n'

        # Determine location to save images
        if 'coarse' in str(output_name):
            save_path = 'imgdump/coarse/'
        else:
            save_path = 'imgdump/fine/'

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
            savefig(str(save_path)+str('averages_coarse_')+str(h)+str('.png'),dpi=300)
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
            savefig(str(save_path)+str('contours_')+str(h)+str('.png'),dpi=300)
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

            savefig(str(save_path)+str('landscape_')+str(h)+str('.png'),dpi=300)
            close()

    chi_store.close()
