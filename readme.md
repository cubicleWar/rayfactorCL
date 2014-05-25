# RayFatorCL

RayFactorCL is a primitive based ray tracer for the calculation of radiative view factors.

RayFactorCL utilises OpenCL for executing all the grunt work with the Philox pseduo random number generator from the [Random 123](http://www.thesalmons.org/john/random123/releases/1.06/docs/) library lying at the heart of the Monte Carlo simulation.

The code was developed targeting the NVidia Fermi and Kepler architectures will run on all OpenCL 1.1 enabled devices.

Although this repo is an XCode project, it can be compiled and run on Windows (A visual studio project file is planed for the future).