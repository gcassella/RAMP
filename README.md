# RAMP: Raytracing Achieved via Massive Parallelisation

*RAMP* is a Monte Carlo raytracing package for the simulation of neutron instrumentation, written in Python and parallelized via the OpenCL API.

![LET Vis](/images/LET_readme_pic.png)

## Documentation

The documentation for RAMP is [hosted online on readthedocs](https://ramp-mcr.readthedocs.io/en/latest/index.html) and may also be found in the `doc` folder of this repository.

## Installation

If the user has their OpenCL SDK already available on the system, RAMP can simply be installed using

```
 $ python setup.py install
```

The package is also available via the Python package index

```
 $ pip install mcramp
```

If errors are encountered relating to "CL\cl.h" or "OpenCL.lib", specify the include directory for your OpenCL installation to pip

```
 $ pip install --global-option=build_ext --global-option="-I/path/to/CL_includes" --global-option="-L/path/to/CL_lib" pyopencl
 $ pip install mcramp
```

If the above is meaningless to you, more complete installation instructions can be found [in the docs](https://ramp-mcr.readthedocs.io/en/latest/user/installation.html).