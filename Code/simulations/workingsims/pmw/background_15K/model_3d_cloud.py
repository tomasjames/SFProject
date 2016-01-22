"""
This is a radmc3dPy model template
A radmc3dPy model file can contain any / all of the functions below

    getDefaultParams()
    getModelDesc()
    getDustDensity()
    getDustTemperature()
    getGasAbundance()
    getGasDensity()
    getGasTemperature()
    getVelocity()
    getVTurb()

The description of the individual functions can be found in the docstrings below the function name.
If a model does not provide a variable or the variable should be calculated by RADMC-3D
(e.g. dust temperature) the corresponding function (e.g. get_dust_temperature) should be removed from
or commented out in the model file.

NOTE: When using this template it is strongly advised to renme the template model (to e.g. model_mydisk.py)
as the get_model_names() function in the setup module removes the name 'template' from the list of available
models.

"""

try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'

from radmc3dPy.natconst import *



# ============================================================================================================================
#
# ============================================================================================================================
def getModelDesc():
    """
    Function to provide a brief description of the model
    """

    return "Example model template"


# ============================================================================================================================
#
# ============================================================================================================================
def getDefaultParams():
    """
    Function to provide default parameter values

    OUTPUT:
    -------

    Returns a list whose elements are also lists with three elements:
    1) parameter name, 2) parameter value, 3) parameter description
    All three elements should be strings. The string of the parameter
    value will be directly written out to the parameter file if requested,
    and the value of the string expression will be evaluated and be put
    to radmc3dData.ppar. The third element contains the description of the
    parameter which will be written in the comment field of the line when
    a parameter file is written.
    """

    defpar = {}

    # Declare information about the cloud to determine radius for later assignment
    #cloud_mass = ms
    #cloud_density = 100000*muh2*mp

    defpar = [['mstar', '1.0*ms', 'Mass of the star(s)'],
    ['pstar', '[0., 0., 0.]', 'Position of the star(s) (cartesian coordinates)'],
    ['rstar', '1.0*rs', 'Radius of the star(s)'],
    ['tstar', '1.0*ts', 'Effective temperature of the star(s)'],
    ['crd_sys', "'car'", 'Coordinate system used (car/sph)'],
    ['nx', '128', 'Number of grid points in the first dimension'],
    ['ny', '128', 'Number of grid points in the second dimension'],
    ['nz', '128', 'Number of grid points in the third dimension'],
    ['xbound', '[0*au, 30000.0*au]', 'Boundaries for the x-grid'],
    ['ybound', '[0*au, 30000.0*au]', 'Boundaries for the y-grid'],
    ['zbound', '[0*au, 30000.0*au]', 'Boundaries for the z-grid'],
    ['nw', '[150]', 'Number of points in the wavelength grid'],
    ['wbound', '[0.1, 1e5]', 'Boundaries for the wavelength grid'],
    ['dustkappa_ext', "['silicate']", 'Dust opacity file name extension'],
    ['nphot', '100000', 'Number of photons in the thermal Monte Carlo simulation'],
    ['scattering_mode_max', '0', '0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering'],
    ['prho', '-2.0', ' Power exponent of the radial density distribution'],
    ['hpr', '0.1', ' Vertical scale height'],
    ['rho0', '1e-16*10000', 'Central density']]
    '''
    ['cloud_density', '100000*muh2*mp', 'Density within the cloud (g/cm^3)'],
    ['outside_density', 100*muh2*mp, 'Density outside of the cloud (g/cm^3)'],
    ['cloud_mass', ms, 'Mass of the cloud (g)'],
    ['cloud_radius', ((3./(4*np.pi))*(cloud_mass/cloud_density))**(1./3.), 'Cloud radius (cm)']]
    '''
    return defpar

'''
# ============================================================================================================================
#
# ============================================================================================================================
def getDustTemperature(grid=None, ppar=None):
    """
    Function to calcualte/set the dust temperature

    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model

    OUTPUT:
    -------
        returns the dust temperature in K
    """

    tdust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) + 1.0
    return tdust
# ============================================================================================================================
#
# ============================================================================================================================
def getGasAbundance(grid=None, ppar=None, ispec=''):
    """
    Function to create the conversion factor from volume density to number density of molecule ispec.
    The number density of a molecule is rhogas * abun

    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model
        ispec - The name of the gas species whose abundance should be calculated

    OUTPUT:
    -------
        returns the abundance as a Numpy array
    """

    gasabun = -1
    if ppar['gasspec_mol_name'].__contains__(ispec):
        ind = ppar['gasspec_mol_name'].index(ispec)
        gasabun[:,:,:] = ppar['gasspec_mol_abun'][ind]

    elif ppar['gasspec_colpart_name'].__contains__(ispec):
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
        ind = ppar['gasspec_colpart_name'].index(ispec)
        gasabun[:,:,:] = ppar['gasspec_colpart_abun'][ind]
    else:
        print 'ERROR'
        print ' The abundance of "'+ispec+'" is not specified in the parameter file'

    return gasabun

# ============================================================================================================================
#
# ============================================================================================================================
def getGasDensity(grid=None, ppar=None):
    """
    Function to create the total gas density distribution

    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model

    OUTPUT:
    -------
        returns the volume density in g/cm^3
    """

    rhogas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + 1e-20

    return rhogas

# ============================================================================================================================
#
# ============================================================================================================================
def getDustDensity(grid=None, ppar=None):
    """
    Function to create the dust density distribution

    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model

    OUTPUT:
    -------
        returns the volume density in g/cm^3
    """

    rhogas  = get_gas_density(grid=grid, ppar=ppar)
    rhodust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64)
    rhodust[:,:,:,0] = rhogas * ppar['dusttogas']
    return rhodust

# ============================================================================================================================
#
# ============================================================================================================================
def getVTurb(grid=None, ppar=None):
    """
    Function to create the turbulent velocity field

    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model

    OUTPUT:
    -------
        returns the turbulent velocity in cm/s
    """

    vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['gasspec_vturb']
    return vturb

# ============================================================================================================================
#
# ============================================================================================================================
def getVelocity(grid=None, ppar=None):
    """
    Function to create the turbulent velocity field

    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model

    OUTPUT:
    -------
        returns the turbulent velocity in cm/s
    """

    vel = np.zeros([grid.nx, grid.ny, grid.nz, 3], dtype=np.float64)
    return vel
'''
