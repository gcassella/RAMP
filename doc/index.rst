RAMP: Raytracing Achieved via Massive Parallelisation
=====================================================

**RAMP** is a Monte Carlo raytracing package for the simulation of neutron \
instrumentation, written in Python and parallelized via the OpenCL API.

-------------------

.. _google-colab:

Quick start
-----------

The quickest way to get hands on with RAMP is running Jupyter notebooks \
via Google Colab, a service that allows Python code to be executed in the \
cloud with access to enterprise grade GPUs. To access RAMP via Colab:

1. Go to https://colab.research.google.com/
2. Log into a Google Account 
3. Go to 'File -> Open Notebook'
4. Go to the 'GitHub' tab
5. Enter 'https://github.com/gcassella/RAMP' and click the search icon
6. Select the notebook you wish to run - 'intro_example.ipynb' is a good place to start
7. Go to 'Runtime -> Change runtime type'
8. Select 'GPU' from the 'Hardware accelerator' dropdown and click 'Save'
9. Run code cells by selecting them and clicking the play icon on the left or pressing Ctrl + Enter 

User guide
----------

.. toctree::
    :maxdepth: 2

    user/installation
    user/hello
    user/idf
    user/running

Kernel list
-----------

.. toctree::
   :maxdepth: 3

   kernels/geom_kernels
   kernels/scat_kernels
