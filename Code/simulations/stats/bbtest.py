################################################################################
############################ 4th Year Project Code #############################
############################# Chi Squared Routine ##############################
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

############################### Define functions ###############################

# Modified blackbody
def mbb(N,opac,wav,T):
    '''
    Function that returns the modified black body law that best describes cold dust emission as stated in Kelly et. al. (2012).

    This returns the intensity per unit solid angle, i.e. the flux density.
    '''
    a = (2*(hh)*(cc)**2)/(wav**5)
    b = ((hh)*(cc))/(wav*(kk)*T)
    #return N*opac*(a*(1./(np.exp(b)-1)))
    return N*opac*a*(1./(np.exp(b)-1))

################################## Spoof data ##################################

# Define wavelengths
wav = np.linspace(60e-4,1400e-4,1701)

# Assign spectral index
B = [0,2]

# Assign column density
N = 1e24

# Define reference quantities
w_ref = 1.303e-3
kappa_0 = 3.09e-1

opacity = []
fluxdensity = []

# Loop through a number of different spextral indices
for i in range(0,2):
    # Evaluate opacities
    opacity.append(kappa_0*(w_ref/wav)**B[i])

    plot(wav,opacity[i],'b--')
    xlabel('Wavelength $(m)$')
    ylabel('Dust Opacity $(g/cm^{3})$')
    savefig(str('dust_opacity_B=')+str(B[i])+str('.png'))
    close()

    # Evaluate fluxes
    fluxdensity.append(mbb(N,opacity[i],wav,T=20))

    loglog(wav,fluxdensity[i],'b--')
    xlabel('Wavelength $(cm)$')
    ylabel('Flux Density')
    savefig(str('flux_density_B=')+str(B[i])+str('.png'))
    close()

loglog(wav,fluxdensity[0],'b--',label=str('B=')+str(B[0]))
loglog(wav,fluxdensity[1],'g--',label=str('B=')+str(B[1]))
grid(True,which="majorminor",ls="-",color='0.65')
xlabel('Wavelength $(cm)$')
ylabel('Flux Density')
legend(loc='best')
savefig('flux_density_combined.png')
close()
