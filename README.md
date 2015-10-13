"How good is dust emission as a tracer of structure in star-forming molecular clouds?"
================================================================================

Project Overview
----------------

This project probes the structure of star forming regions in molecular clouds
by using dust emission as a tracer.

This repository contains all code (written exclusively in Python) pertaining to
the project. At present, the only specialist package used is RADMC-3D, a raytracing
radiative transfer code (http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/).

Project Aims
------------

The primary aims of this project are:

- to simulate a molecular cloud containing x amount of pre-stellar cores (i.e. generate synthetic data)
- to produce a dendogram (https://dendrograms.readthedocs.org/en/latest/) of the
  resulting simulation to probe the structure of the simulated cloud
- compare with other simulations

Project Progress
----------------

As of last commit (11/10/15) minor tests have been performed of RADMC-3D in Python
to ascertain its viability for the project. Those tests setup various simple
problems and simulations, the most useful of which is the creation of a 2D (as well
as 1D) sphere of uniform temperature.

Project 'To-Do List'
-------------------
This list is not exhaustive of the remaining tasks in the project, and is updated as and when
new ideas are added. Furthermore, this is not representative of the project's progress.

- ~~Get RADMC-3D up and running on OS X 10.11 El Capitan~~
- ~~Run basic tests included in cheatsheet (i.e. set up gas and dust continuum models)~~
- ~~Generate sphere of uniform temperature:~~
  - ~~in 1D &,~~
  - ~~2D~~
- Adapt 2D sphere model to account for gas AND dust
- Look into generated some form of rudimentary stucture in 2D sphere
