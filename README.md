# CMEPS
NUOPC based Community Mediator for Earth Prediction Systems

## Overview and resources

The Community Mediator for Earth Prediction Systems (CMEPS) is a
NUOPC-compliant Mediator component used for coupling Earth system
model components. It is currently being used in NCAR's Community
Earth System Model (CESM) and NOAA's UFS subseasonal-to-seasonal
coupled system application.

For documentation see

https://escomp.github.io/CMEPS/  

In order to build the package, the NCAR [ParallelIO package](https://github.com/NCAR/ParallelIO) must be installed and an environment variable PIO=${PIO_DIRECTORY} set. [PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) is optional.

To build stand-alone libraries, run:
```
cmake .
make
```
