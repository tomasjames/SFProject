################################################################################
############################ 4th Year Project Code #############################
##############################  Functions to run  ##############################
################################################################################

############################# Import statements ################################

# Import simulation run script
from sim.simulation import *

# Import Chi-squared bits
from curvefit.curvefitting import *

# Import image making scripts
from imgprocess.img import *

# This should be declared in all other imported files, but redeclare here just in case
import matplotlib as mpl
mpl.use('Qt4Agg')
from matplotlib.pyplot import *

print '\n######################################################################'
print '############################ handler.py ##############################'
print '######################################################################\n'

print '#######################################################################'
print 'This is the handler script for the dust emission project. This script '
print 'will execute all files required to run the simulation and analysis the '
print 'data with minimal user input, unless when required.'
print 'These steps will run the simulation and generate the pixel by pixel SEDs '
print 'for each band. It will then run the pixel-by-pixel Chi-Squared routine '
print 'before then visualising that data with maps of each data type. Finally, '
print 'the code will generate dendrograms and contour plots for the data.'
print '######################################################################\n'

############################# Run simulations ################################

# Define variables needed for function calls
mode = 'd'
npix = 128
sizeau = 15000
d = 300.
mass = 1.
cloud_density = 1e5
outside_density = 1e2
cloud_temp = 10.
outside_temp = 15.
amr = True
dust = True

# Define filters to run over
filt = np.array(['blue', 'green', 'plw', 'pmw', 'psw', 'red'])

print '#######################################################################'
#data_type = input('Are you analysing ideal simulation (\'sim\') data, filament simulation (\'filament\') or SPH/Arepo (\'sph\') data?\n')
#B = input('What value of B are you using? Select from: B=1.8, B=2.0, B=2.2\n')
#kappa_0 = input('What value of kappa_0 should be used? B=2.0 has kappa_0=4\n')
print '######################################################################\n'

data_type = 'sph'
#kappa_0 = 0.042
#lambda_0 = 300e-4 # Reference wavelength should be in cm
kappa_0 = 4.
lambda_0 = 250e-4

# Different values of dust emissivity index
B = input('What value of B are you using? Select from: B=1.8, B=2.0, B=2.2\n')

if data_type == 'sim':

    # Loop through all values of B
    for B in beta:

        # Fix the string so that it coincides with folder names
        B_val = str('B=')+str(B)

        # Run through the files
        for suffix in filt:

            # Print the intentions of the script to keep track of code's location
            print(str('Changing the directory to:')+str(' simulations/workingsims_psf/')+str(B_val)+str('/')+str(suffix)+str('/background_15K/\n'))

            # Change the directory to the directory housing the band data
            os.chdir('simulations/workingsims_psf/'+str(B_val)+str('/')+str(suffix)+str('/background_15K/'))

            # Print to tell simulation is being run
            print(str('Now running the simulation for ')+str(suffix)+str('\n'))

            # Run the simulation itself
            sim_data = simulation(mode=mode ,filt=str(suffix), npix=npix, sizeau=sizeau, d=d, mass=mass, cloud_density=cloud_density, outside_density=outside_density, cloud_temp=cloud_temp, outside_temp=outside_temp, kappa_0=kappa_0, lambda_0=lambda_0, B=B, amr=amr, dust=dust, sim_name='cloud')

            # Print to tell simulation is being run
            print(str('Now determining the SED for')+str(suffix)+str('\n'))

            # Determine the limit to filter for
            min_val = np.amin(sim_data)
            max_val = np.amax(sim_data)
            std_val = np.std(sim_data)

            # Determine the SED
            sedGeneration(filt=str(suffix), sim_name='cloud', kappa_0=kappa_0, lambda_0=lambda_0, B=B, withPSF=True)

            # Change the directory back to the folder containing this file
            os.chdir('../../../../../')

            # Code will now go back to the top of the loop and execute over the next band

        ############################# Run the Chi-squared #############################

        print '#######################################################################'
        print 'All simulations have now been run and their SEDs computed. The code '
        print 'will now move on to the Chi-squared analysis.'
        print '######################################################################\n'

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/curvefitting/cloud/')+str(B_val)+str('/\n'))

        # Once all bands have been executed over, need to run the Chi-Squared routine
        # Change directory to the directory that houses the Chi-Squared routine
        os.chdir(('simulations/curvefitting/cloud/')+str(B_val))

        # Another print statement
        print(str('Executing the datafeed run.\n'))

        dataDerive(data_type='radmc', kappa_0=kappa_0, lambda_0=lambda_0, B=B)

        # Another print statement
        print(str('Executing the coarse Chi-squared run. This may take some time\n'))

        # Determine quantities to loop through for N and T
        N = np.linspace(17,24,50)
        T = np.linspace(6,16,50)

        # Execute the chi-squared routine
        chiTest(data_type='radmc', output_name='chi_coarse.txt', N=N, T=T, kappa_0=kappa_0, lambda_0=lambda_0, B=B)

        # Another print statement
        print(str('Executing the fine Chi-squared run. This may take some time\n'))

        # Determine quantities to loop through for N and T
        N_data = np.log10(np.loadtxt('chi_coarse.txt',skiprows=1)[:,1])
        T_data = np.loadtxt('chi_coarse.txt',skiprows=1)[:,2]

        N = np.linspace(min(N_data),max(N_data),100)
        T = np.linspace(min(T_data),max(T_data),100)

        # Execute the chi-squared routine
        chiTest(data_type='radmc', output_name='chi_fine.txt', N=N, T=T, kappa_0=kappa_0, lambda_0=lambda_0, B=B)

        # Change the directory back to the folder containing this file
        os.chdir('../../../../../')

        ################### Run the mapping and dendrogram suites ######################

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/imgprocess/cloud_psf/')+str(B_val)+str('/\n'))

        # Go to the image processing folder
        os.chdir(('simulations/imgprocess/cloud_psf/')+str(B_val)+str('/'))

        # Another print statement
        print(str('Executing mapping script\n'))

        # Execute the script
        mapMaker(data_type='radmc')

        # Another print statement
        print(str('Executing dendrogram script\n'))

        # Execute the script
        dendrogram()

        print '#######################################################################'
        print 'handler.py has now finished running all necessary files. The results '
        print 'can be found in simulations/imgprocess/cloud/.'
        print '######################################################################\n'

elif data_type == 'sph':

    # Loop through all values of B
    for B in beta:

        # Fix the string so that it coincides with folder names
        B_val = str('B=')+str(B)

        # Run through the files
        for suffix in filt:
            '''
            # Print the intentions of the script to keep track of code's location
            print(str('Changing the directory to:')+str(' simulations/data_psf/')+str(B_val)+str('/')+str(suffix)+str('/dust_project/\n'))

            # Change the directory to the directory housing the band data
            os.chdir('simulations/data_psf/'+str(B_val)+str('/')+str(suffix)+str('/dust_project/'))

            # Print to tell simulation is being run
            print(str('Now running the simulation for ')+str(B_val)+str('/')+str(suffix)+str('\n'))

            # Run the simulation itself
            sim_data = simulation(mode='r', filt=str(suffix), npix=200, sizeau=133692, d=d, mass=mass, cloud_density=cloud_density, outside_density=outside_density, cloud_temp=cloud_temp, outside_temp=outside_temp, kappa_0=kappa_0, lambda_0=lambda_0, B=B, amr=False, dust=False, sim_name='sphdata')
            '''
            # Print to tell simulation is being run
            print(str('Now determining the SED for')+str(B_val)+str('/')+str(suffix)+str('\n'))

            # Determine the SED
            sedGeneration(filt=str(suffix), sim_name='sphdata', kappa_0=kappa_0, lambda_0=lambda_0, B=B, withPSF=True)

            # Change the directory back to the folder containing this file
            os.chdir('../../../../../')

            # Code will now go back to the top of the loop and execute over the next band

        ############################# Run the Chi-squared #############################

        print '#######################################################################'
        print 'All simulations have now been run and their SEDs computed. The code '
        print 'will now move on to the Chi-squared analysis.'
        print '######################################################################\n'

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/curvefitting/sphdata/')+str(B_val)+str('/\n'))

	    os.chdir(('simulations/curvefitting/sphdata/')+str(B_val)+str('/'))

        # Another print statement
        print(str('Executing the datafeed run.\n'))

        dataDerive(data_type='arepo', kappa_0=kappa_0, lambda_0=lambda_0, B=B)

        # Another print statement
        print(str('Executing the coarse Chi-squared run. This may take some time\n'))

        # Determine quantities to loop through for N and T
        N = np.linspace(17,24,50)
        T = np.linspace(6,16,50)

        # Execute the chi-squared routine
        chiTest(data_type='arepo', output_name='chi_coarse.txt', N=N, T=T, kappa_0=kappa_0, lambda_0=lambda_0, B=B)
        '''
        # Another print statement
        print(str('Executing the fine Chi-squared run. This may take some time\n'))

        # Determine quantities to loop through for N and T
        N_data = np.log10(np.loadtxt('chi_coarse.txt',skiprows=1)[:,1])
        T_data = np.loadtxt('chi_coarse.txt',skiprows=1)[:,2]

        N = np.linspace(min(N_data),max(N_data),100)
        T = np.linspace(min(T_data),max(T_data),100)

        # Execute the chi-squared routine
        chiTest(data_type='arepo', output_name='chi_fine.txt', N=N, T=T, kappa_0=kappa_0, lambda_0=lambda_0, B=B)
        '''
        # Change the directory back to the folder containing this file
        os.chdir('../../../../../')

        ################### Run the mapping and dendrogram suites ######################

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/imgprocess/sphdata_psf')+str(B_val)+str('\n'))

        # Go to the image processing folder
        os.chdir(('simulations/imgprocess/sphdata_psf/')+str(B_val))

        # Another print statement
        print(str('Executing mapping script\n'))

        # Execute the script
        mapMaker(data_type='arepo')

        # Another print statement
        print(str('Executing dendrogram script\n'))

        # Execute the script
        dendrogram()

        print '#######################################################################'
        print 'handler.py has now finished running all necessary files. The results '
        print 'can be found in simulations/imgprocess/sphdata/.'
        print '######################################################################\n'
