################################################################################
############################ 4th Year Project Code #############################
########################## Dust Opacity Determination ##########################
################################################################################

############################# Import statements ################################

# Import OS manipulation tools
import os

########################## Define dust opacity function ########################

def dustopacity(kappa_0, v, v_0, B):

    '''
    Uses the dust opacity power law defined in Shetty et al (2009) to determine
    the dust opacity at a range of frequencies v, given some initial frequency
    v_0 and dust spectral index B.
    '''

    # This is the frequency dependent dust opacity defined in Shetty et al (2009)
    opac = kappa_0*(v/v_0)**B
