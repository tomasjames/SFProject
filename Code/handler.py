################################################################################
############################ 4th Year Project Code #############################
############################### Script Handler #################################
################################################################################

############################# Import statements ################################

import os

############################# Run the simulation ###############################

print '\n######################################################################'
print '############################ handler.py ##############################'
print '######################################################################\n'

print '#######################################################################'
print 'This is the handler script for the dust emission project. This script '
print 'will execute all files required to run the simulation and analysis the '
print 'data with minimal user input, unless when required. The steps run run '
print 'are as follows: '
print 'run: sim_band.py'
print 'run: band_sedv3.py'
print 'run: chiv6.py'
print 'run: map.py'
print 'run: dendrogram.py'
print 'These steps will run the simulation and generate the pixel by pixel SEDS '
print 'for each band. It will then run the pixel-by-pixel Chi-Squared routine '
print 'before then visualising that data with maps of each data type. Finally, '
print 'the code will generate dendrograms and contour plots for the data.'
print '######################################################################\n'

# Store a list of the bands to loop through
band = ['blue', 'green', 'psw', 'pmw', 'plw', 'red']

print '#######################################################################'
data_type = input('Are you analysing ideal simulation (\'sim\') data, filament simulation (\'filament\') or SPH/Arepo (\'sph\') data?\n')
print '######################################################################\n'

if data_type == 'sim':
    '''
    # Run through the files
    for suffix in band:

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/workingsims/')+str(suffix)+str('/background_15K/\n'))

        # Change the directory to the directory housing the band data
        os.chdir('simulations/workingsims/'+str(suffix)+str('/background_15K/'))

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
    print(str('Changing the directory to:')+str(' simulations/curvefitting/cloud/\n'))

    # Once all bands have been executed over, need to run the Chi-Squared routine
    # Change directory to the directory that houses the Chi-Squared routine
    os.chdir('simulations/curvefitting/cloud/')

    # Another print statement
    print(str('Executing:')+str('chiv6.py\n'))

    # Execute the chi-squared routine
    execfile(str('chiv6.py'))

    # Change the directory back to the folder containing this file
    os.chdir('../../../')

    ################### Run the mapping and dendrogram suites ######################

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/imgprocess/cloud/\n'))

    # Go to the image processing folder
    os.chdir('simulations/imgprocess/cloud/')

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
    '''
    # Run through the files
    for suffix in band:

        # Print the intentions of the script to keep track of code's location
        print(str('Changing the directory to:')+str(' simulations/data/')+str(suffix)+str('/dust_project/\n'))

        # Change the directory to the directory housing the band data
        os.chdir('simulations/data/'+str(suffix)+str('/dust_project/'))

        # Another print statement
        print(str('Executing:')+str(' sim_')+str(suffix)+str('.py\n'))

        # Execute the simulation
        execfile(str('sim_')+str(suffix)+str('.py'))

        # Another print statement
        print(str('sim_')+str(suffix)+str('.py\n'))

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
    print(str('Changing the directory to:')+str(' simulations/curvefitting/sphdata/\n'))

    # Once all bands have been executed over, need to run the Chi-Squared routine
    # Change directory to the directory that houses the Chi-Squared routine
    os.chdir('simulations/curvefitting/sphdata/')

    # Another print statement
    print(str('Executing:')+str('chiv6.py\n'))

    # Execute the chi-squared routine
    execfile(str('chiv6.py'))

    # Change the directory back to the folder containing this file
    os.chdir('../../../')

    ################### Run the mapping and dendrogram suites ######################

    # Print the intentions of the script to keep track of code's location
    print(str('Changing the directory to:')+str(' simulations/imgprocess/sphdata/\n'))

    # Go to the image processing folder
    os.chdir('simulations/imgprocess/sphdata/')

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
    print 'can be found in simulations/imgprocess/sphdata/.'
    print '######################################################################\n'
