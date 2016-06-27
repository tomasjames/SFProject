"How good is dust emission as a tracer of structure in star-forming molecular clouds?"
================================================================================

<p align="center">
  <img src="https://github.com/tomasjames/ZiggyStarDust/blob/master/Code/logo.JPG?raw=true" alt="Ziggy Stardust"/>
</p>

Project Overview
----------------

This project probes the structure of star forming molecular clouds by using dust emission as a tracer. The ‘formal’ assessment section of the project was completed on 24th May 2016 with a presentation (found in this repo) to peers on the work completed. 

**The project has also been extended beyond this assessment period to encompass 8 weeks of the Summer. This work involves contributions from Dr. Paul Clark, Prof. Anthony Whitworth and Dr. Ken Marsh.

This repository contains all code (written exclusively in Python) pertaining to the project. At present, the only specialist package used is RADMC-3D, a raytracing radiative transfer code (http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/).

Project Aims
------------

The (non-exhaustive) aims of the Summer extension are:

- better simulate real data by:
    - convolving the data with a PSF for the filters considered
    - simulating some form of photon shot noise
- run PPMAP code over FITS files to see how this simulation differs from the SED fitting routine

The primary (simplified) aims of this project were:

- to simulate a molecular cloud using the Arepo simulation containing x amount of pre-stellar cores (i.e. generate synthetic data)
- recover the values of dust column density, N and dust temperature, T for each pixel in the resulting data
- also use initial simulation input parameters to derive values of N and T for comparison purposes
- to produce dendrograms (https://dendrograms.readthedocs.org/en/latest/) of the resulting simulation recovered N to probe the structure of the simulated cloud
    - repeat this procedure for initial input derived N     
- Extract the core locations and masses in both dendrograms and compare populations (i.e. their locations, masses and number of cores detected)

Project Progress
----------------

27/06/2016: One week of work has culminated in a working PSF convolution script. At this point, it convolves any RADMC-3D FITS output with the PSF for a given filer. Now moving on to determine number of beams required to provide better resolution of Arepo data, as well as refactoring code into functions to aid productivity.

09/05/16: Final commit prior to signing off on development. All code has now been written, including the supplementary dissertation. To run code, simply use run the `handler.py` file. This will ask for the simulation to be run: working simulations include the `sim` and `sph` responses (filament still non-functional). If the user wishes, an `inputs.txt` file can be directed to the code such the user inputs from the code will come directly from this file, e.g. to run the basic isothermal simulation one would run `python handler.py < inputs_sim.txt` - this would run the code from start to finish. 

20/03/16: The machinery to analyse data is now *all* in place and functioning. All SEDs produced are the correct form, and the chi-squared routine, that now successfully loops over each pixel, works correctly. The next step in producing further results is to:
- write a brief script to plot maps of the recovered column density and temperature
  - plot contour plots should time allow
- apply this machinery to SPH simulation data and recover values of N and T
- push forwards to analyse data dendrogramically
- write up!

08/02/16: Code has been reformatted since last major additions, and a workflow has now been devised. To properly run the code, the following commands must be executed in the following order:
- `cd` to a simulation directory, eg `psw/background_15K`
- run `python sim_psw.py` to run the main raytracing simulation. This will also run `datafilegen.py` and save the `.inp` files to the working directory to be read in immediately after the script has been completed and the ray trace initiated. Once ray trace has been completed, save output (an image) if required (not necessary, though makes a nice wall hanging picture)
- run `python psw_sed.py` to compute the object's Spectral Energy Distribution. This script will not plot in any active windows but will save all required plots to the working directory. It will also write the essential SED info to a `.txt` file in the `stats` directory for later use, appending each time a `band_sed.py` script is run.
- repeat the above process until all 3 bands (PSW, PMW, PLW) have been completed. For ease of use, code should be run in PSW, PMW and PLW order (though is not essential).
- `cd` to `../../../stats` and run `python chiv2.py` to plot the averaged SED information and to fit theoretical SED to it.

NOTE: the chi-squared routine is not yet complete, and shall we completed as the next task however the basic data generation, simulation and reduction routines have now been completed and work.

04/12/15: The model has since been adapted to plot the three passbands on SPIRE in order to begin generating the synthetic data required to compare with the real data from Herschel SPIRE. However, transmission has not been accounted for in these models to date - the .txt files containing the transmission curve do not have regular wavelength spacing, necessitating the use of interpolation to generate a linearly spaced wavelength-transmission parameter space.

Once this is complete plots of each wavelength are multiplied by their transmission coefficient and summed (weighted) over all possible wavelengths in the passband to give the final image as seen by SPIRE. At this point, an SED can be determined by invoking `radmc3d spectra` (check that this is the correct command).

20/11/15: Custom data files for RADMC-3D are now generated with a Python script that asks for user defined quantities. It is possible to compute the radius of the molecular cloud from its mass and density (the code requests number density and converts this to g/cm^3). The resulting radius is then used to help define the image width (i.e. the user is told what the cloud radius is so as to make appropriate adjustments to the image width). Furthermore this width is split into n user defined pixels and each pixel in each dimension looped over to assign density and temperature.

These files are then copied over to the directory that the image raytrace is performed. By invoking `radmc3d image` the raytrace is then performed (providing all the information required is written to the `problem-params.inp` file - this is usually performed in model setup, ie `radmc3dPy.analyze.writeDefaultParfile('model_name')`). Allowing the code to run produces an image as seen below:

![A 1 MSun cloud (dust mass=0.01 Msun) with a Sun like star in the centre](Code/simulations/workingsims/psw/MSun_psw.png "A 1 MSun cloud with dust mass=0.01 Msun - this model has no source of photons.")

IMPORTANT: ALL UNITS MUST BE CGS (centimetre-grams-second)

As of 11/10/15 minor tests have been performed of RADMC-3D in Python
to ascertain its viability for the project. Those tests setup various simple
problems and simulations, the most useful of which is the creation of a 2D (as well
as 1D) sphere of uniform temperature.

As of 14/10/15 more complex simulations have been set up ([#1](https://github.com/tomasjames/ZiggyStarDust/issues/1) handles this).

Project 'To-Do List'
-------------------
This list is not exhaustive of the remaining tasks in the project, and is updated as and when
new ideas are added. Furthermore, this is not representative of the project's progress. More
often than not, for larger implementations more detailed issues will be created in addition
to the notes found here. Where possible, those issues will be linked next to the entry.

- [x] Get RADMC-3D up and running on OS X 10.11 El Capitan
- [x] Run basic tests included in cheatsheet (i.e. set up gas and dust continuum models)
- [x] Generate sphere of uniform temperature:
  - [x] in 1D &,
  - [x] 2D
- [x] Begin collecting BiBTeX library of references for report
- [x] Generate a 3D sphere and place a source (e.g. star) behind to assess radiative transfer. Use:
  - [x] Different opacity law (see notes in notebook for papers)
  - [x] Standard 1/r**2 density profile (i.e. decreasing with radius)
  - [x] Standard molecular cloud temperature profile (cooler in the cloud relative to the exterior)
    (See [#1](https://github.com/tomasjames/ZiggyStarDust/issues/1) for updates)
  - [x] Investigate the `ValueError: zero-size array to reduction operation minimum which has no identity` error that RADMC-3D is throwing when trying to run the raytrace with all custom files in place.
  - ~~[x] Implement IRF~~
- [x] Expand the code to take SPIRE waveband inputs
  - [x] PSW, PMW and PLW bands need to be coded so that simulations can be run in each of their respective passbands
  - [x] Extract SED from the simulations by using `radmc3d sed loadlambda`
  - [x] Convolve these SEDs with the SPIRE transmission curves
    - [x] Weighted average these to generate the flux 'seen' by SPIRE
  - [x] Fit an SED to the resulting datapoints to extract the column density and B values
    - [x] Assess quality of the fit using chi-squared routine
  - [x] Plot probability contours to verify the banana shaped contour as seen in Shetty et. al.
