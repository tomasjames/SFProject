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
data_type = input('Are you analysing ideal simulation (\'sim\') data, filament simulation (\'filament\') or SPH/Arepo (\'sph\') data?\n')
print '######################################################################\n'

if data_type == 'sim':

    # Run through the files
    for suffix in filt:

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/workingsims_psf/')+str(suffix)+str('/background_15K/\n'))

        # Change the directory to the directory housing the band data
        os.chdir('simulations/workingsims_psf/'+str(suffix)+str('/background_15K/'))

        # Print to tell simulation is being run
        print(str('Now running the simulation for ')+str(suffix)+str('\n'))

        # Run the simulation itself
        simulation(mode=mode ,filt=str(suffix), npix=npix, sizeau=sizeau, d=d, mass=mass, cloud_density=cloud_density, outside_density=outside_density, cloud_temp=cloud_temp, outside_temp=outside_temp, amr=amr, dust=dust, sim_name='cloud')

        # Print to tell simulation is being run
        #print(str('Now determining the SED for')+str(suffix)+str('\n'))

        # Determine the SED
        #sedGeneration(filt=str(suffix), sim_name='cloud')

        # Change the directory back to the folder containing this file
        os.chdir('../../../../')

        # Code will now go back to the top of the loop and execute over the next band

    ############################# Run the Chi-squared #############################

    print '#######################################################################'
    print 'All simulations have now been run and their SEDs computed. The code '
    print 'will now move on to the Chi-squared analysis.'
    print '######################################################################\n'

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/curvefitting/cloud/\n'))

    # Once all bands have been executed over, need to run the Chi-Squared routine
    # Change directory to the directory that houses the Chi-Squared routine
    os.chdir('simulations/curvefitting/cloud/')

    # Another print statement
    print(str('Executing the Chi-squared run. This may take some time\n'))

    # Execute the chi-squared routine
    chiTest(data_type='radmc')

    # Change the directory back to the folder containing this file
    os.chdir('../../../')

    ################### Run the mapping and dendrogram suites ######################

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/imgprocess/cloud_psf/\n'))

    # Go to the image processing folder
    os.chdir('simulations/imgprocess/cloud_psf/')

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

elif data_type == 'filament':
    '''
    # Run through the files
    for suffix in band:

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/filament/')+str(suffix)+str('/core/\n'))

        # Change the directory to the directory housing the band data
        os.chdir('simulations/filament/'+str(suffix)+str('/core/'))

        # Another print statement
        print(str('Executing:')+str(' sim_')+str(suffix)+str('.py\n'))

        # Execute the simulation
        execfile(str('sim_')+str(suffix)+str('.py'))

        # Another print statement
        print(str('Executing:')+str(suffix)+str('_sedv3.py\n'))

        # Run the SED generation script
        execfile(str(suffix)+str('_sedv3')+str('.py'))

        # Change the directory back to the folder containing this file
        os.chdir('../../../../')

        # Code will now go back to the top of the loop and execute over the next band
    '''
    ############################# Run the Chi-squared #############################

    print '#######################################################################'
    print 'All simulations have now been run and their SEDs computed. The code '
    print 'will now move on to the Chi-squared analysis.'
    print '######################################################################\n'

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/curvefitting/filament/\n'))

    # Once all bands have been executed over, need to run the Chi-Squared routine
    # Change directory to the directory that houses the Chi-Squared routine
    os.chdir('simulations/curvefitting/filament/')

    # Another print statement
    print(str('Executing:')+str('chiv6.py\n'))

    # Execute the chi-squared routine
    execfile(str('chiv6.py'))

    # Change the directory back to the folder containing this file
    os.chdir('../../../')

    ################### Run the mapping and dendrogram suites ######################

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/imgprocess/filament/\n'))

    # Go to the image processing folder
    os.chdir('simulations/imgprocess/filament/')

    # Another print statement
    print(str('Executing:')+str(' map.py'))

    # Execute the script
    execfile(str('map.py'))

    # Another print statement
    print(str('Executing:')+str(' dendrogram.py'))

    # Execute the script
    execfile(str('dendrogram.py'))

    print '#######################################################################'
    print 'handler.py has now finished running all necessary files. The results '
    print 'can be found in simulations/imgprocess/filament/.'
    print '######################################################################\n'

elif data_type == 'sph':

    # Run through the files
    for suffix in filt:

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/data_psf/')+str(suffix)+str('/dust_project/\n'))

        # Change the directory to the directory housing the band data
        os.chdir('simulations/data_psf/'+str(suffix)+str('/dust_project/'))

        # Print to tell simulation is being run
        print(str('Now running the simulation for ')+str(suffix)+str('\n'))

        # Run the simulation itself
        simulation(mode=mode, filt=str(suffix), npix=200, sizeau=133692, d=d, mass=mass, cloud_density=cloud_density, outside_density=outside_density, cloud_temp=cloud_temp, outside_temp=outside_temp, amr=False, dust=False, sim_name='sphdata')

        # Print to tell simulation is being run
        print(str('Now determining the SED for')+str(suffix)+str('\n'))

        # Determine the SED
        sedGeneration(filt=str(suffix), sim_name='sphdata')

        # Change the directory back to the folder containing this file
        os.chdir('../../../../')

        # Code will now go back to the top of the loop and execute over the next band

    ############################# Run the Chi-squared #############################
    '''
    print '#######################################################################'
    print 'All simulations have now been run and their SEDs computed. The code '
    print 'will now move on to the Chi-squared analysis.'
    print '######################################################################\n'

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/curvefitting/sphdata/\n'))

    # Once all bands have been executed over, need to run the Chi-Squared routine
    # Change directory to the directory that houses the Chi-Squared routine
    os.chdir('simulations/curvefitting/sphdata/')

    # Another print statement
    print(str('Executing the Chi-squared run. This may take some time\n'))

    # Execute the chi-squared routine
    chiTest(data_type='arepo')

    # Change the directory back to the folder containing this file
    os.chdir('../../../')

    ################### Run the mapping and dendrogram suites ######################

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/imgprocess/sphdata_psf/\n'))

    # Go to the image processing folder
    os.chdir('simulations/imgprocess/sphdata_psf/')

    # Another print statement
    print(str('Executing mapping script\n'))

    # Execute the script
    mapMaker(data_type='arepo')

    # Another print statement
    print(str('Executing dendrogram script\n'))

    # Execute the script
    dendrogram()
    '''
    print '#######################################################################'
    print 'handler.py has now finished running all necessary files. The results '
    print 'can be found in simulations/imgprocess/sphdata/.'
    print '######################################################################\n'
