Installation
============

Prerequisites
-------------
Running RAMP requires a working Python 3 installation. If you are unfamiliar with \
Python, this can be easily obtianed as part of the `Anaconda package <https://www.anaconda.com/>`_ \
which comes preloaded with many useful libraries.

In addition to Python 3, RAMP is parallelized using the OpenCL API. To install the \
Python port of the OpenCL bindings and compile the OpenCL kernels used by RAMP, it \
will be necessary to install the OpenCL SDK appropriate to the vendor of the hardware \
you intend to run RAMP on. At the time of writing, these SDKs can be obtained from

    - `Intel OpenCL SDK <https://software.intel.com/en-us/intel-opencl/>`_
    - `NVIDIA CUDA Toolkit <https://developer.nvidia.com/cuda-downloads/>`_
    - `AMD GPUOpen SDK <https://gpuopen.com/compute-product/opencl-sdk/>`_

Please note: **you only need to install the SDK corresponding to the manufacturer of \
the device you wish to execute OpenCL code on**.

Installing PyOpenCL
-------------------

As OpenCL implementations are platform dependant, there is no one simple path to \
installing PyOpenCL for your specific combination of platform and device. The following \
are some broad pointers from experience that make the process easier.

For NVIDIA device users: The version of the CUDA toolkit required to obtain the correct \
OpenCL drivers for your device can be found \
`here  <https://docs.nvidia.com/deploy/cuda-compatibility/index.html#binary-compatibility__table-toolkit-driver>`_.

Installing PyOpenCL via `pip` requires some manual handling of prerequisites, however \
this can be avoided by installing via `conda` instead by invoking::

 $ conda config --add channels conda-forge
 $ conda install pyopencl=2019.1.2

Obtaining RAMP
--------------

RAMP is most easily obtained from the Python Package Index via `pip`,

 $ pip install mcramp

RAMP can also obtained by cloning the public GitHub repository to get the most up to date, \
but potentially unstable, features::

 $ git clone https://github.com/gcassella/RAMP.git

The package can then be installed via::

 $ cd RAMP
 $ python setup.py install

To ensure the package has been installed correctly, try running one of the examples \
found in the 'examples' directory::

 $ cd examples
 $ python LET_RAMP.py

If this prompts you to select your OpenCL device and then spits out a handful of \
detector image plots, everything is installed correctly!

Jupyter notebooks
-----------------

Some of the instructional examples for RAMP are presented as `Jupyter <https://jupyter.org/>`_ \
notebooks. These are interactive Python environments that allow code to be presented, edited, \
and ran alongside relevant text and graphical outputs. The easiest way to get up and running \
with RAMP notebooks is via Google Colab, see :ref:`google-colab`, however if you wish to run the notebooks locally \
follow the instructions here.

To install Jupyter (assuming you followed the instructions above and have a functional \
Python installation) run::

 $ pip install jupyterlab

Jupyter can then be run via::

 $ cd %your_ramp_directory_here%/notebooks/%example_you_want_to_run%
 $ jupyter notebook

This will open a page in your default web browser from which you can select the \
example notebook and run the interactive Python therein.

An example installation
-----------------------

The following sequence of commands resulted in a correctly configured installation \
on Windows Subsystem for Linux (Ubuntu 20.04) using an NVIDIA 1060 GPU.

Install Miniconda::

 $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
 $ chmod +x ./Miniconda3-latest-Linux-x86_64.sh
 $ ./Miniconda3-latest-Linux-x86_64.sh
 
Update gcc, required for PyOpenCL installation::

 $ sudo apt update
 $ sudo apt install gcc

Installed CUDA toolkit, required version determined as suggested above::

 $ wget http://developer.download.nvidia.com/compute/cuda/11.0.2/local_installers/cuda_11.0.2_450.51.05_linux.run
 $ sudo sh cuda_11.0.2_450.51.05_linux.run

Installed PyOpenCL via conda::

 $ conda config --add channels conda-forge
 $ conda install pyopencl=2019.1.2

Installed RAMP via pip::

 $ pip install mcramp