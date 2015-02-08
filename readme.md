# RayFatorCL

RayFactorCL is a primitive based ray tracer for the calculation of radiative view factors. It was developed as part of the research for the PhD Thesis ["The use of primitives in the calculation of radiative view factors"](http://ses.library.usyd.edu.au/bitstream/2123/10275/1/walker_tj_thesis.pdf). This repo is intended primarily as an example of an OpenCL project and secondly to host the code for follow on research.

Although at its heart RayFactor is a ray tracer for arbitrary 3D geometries it has been somewhat tailored for further research it the radiative modelling of fibre drawing furnaces (in particular the bounding surface functionality).

RayFactorCL utilises OpenCL for executing all the grunt work with the Philox pseudo random number generator from the [Random 123](http://www.thesalmons.org/john/random123/releases/1.06/docs/) library lying at the heart of the Monte Carlo simulation.

The code was developed targeting the NVidia Fermi and Kepler architectures will run on all OpenCL 1.1 enabled devices.

Although this repo is an XCode project, it can be compiled and run on Windows (A visual studio project file is planed for the future).

## Notes on Bounding <s>Volume</s> Surface Support

Very rudimentary bounding surface support has been added to RayFactorCL, this is essentially a rough hack that improves the speed of simulations of geometries containing objects which could be represented by a single object but are discretised for various reasons i.e. a cylindrical furnace wall which is broken up into 180 cylinders for heat flux calculations.

However, bounding surface specification is not required and orginal format input files can be used.

## Roadmap

Future developments as required are as follows:

- Roll in triangle FEM and bounding surface primitive support. 