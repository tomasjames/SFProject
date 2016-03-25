################################################################################
############################ 4th Year Project Code #############################
############################### Script Handler #################################
################################################################################

############################# Import statements ################################

# Store a list of the bands to loop through
band = ['blue', 'green', 'psw', 'pmw', 'plw', 'red']

# Run through the files
for i in band:
    if i == 'blue':
        os.system('cd simulations/workingsims/')+str(i)+str('/background_15K/)
        execfile(str('sim_')+str(i)+str('.py'))
    else:
        os.system('cd ../../simulations/workingsims/')+str(i)+str('/background_15K/)
        execfile(str('sim_')+str(i)+str('.py'))
