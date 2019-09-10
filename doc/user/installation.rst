Installation
============

Prerequisites
-------------
Running RAMP requires a working Python 3 installation. If you are unfamiliar with \
Python, this can be easily obtianed from the `Python website <https://www.python.org/>`_ \
or as part of the `Anaconda package <https://www.anaconda.com/>`_  which comes \
preloaded with many useful libraries.

In addition to Python 3, RAMP is parallelized using the OpenCL API. To execute \
OpenCL kernels on your device, it is necessary to have the appropriate OpenCL \
runtime installed. For older devices, this can be something of a dark art, however \
most modern GPUs should provide an OpenCL runtime as part of their drivers. If you \
are experiencing errors related to "CL.h" when trying to run RAMP, it is likely that \
you do not have an OpenCL runtime properly configured.

Obtaining RAMP
--------------
RAMP is most easily obtained by cloning the public GitHub repository::

 $ git clone https://github.com/gcassella/RAMP.git

which can then be installed via::

 $ cd RAMP
 $ python setup.py install