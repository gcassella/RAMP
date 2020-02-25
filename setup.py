# setup.py
from distutils.core import setup
from setuptools import find_packages

# Obtain the relative path for all of the headers to include for OpenCL as they
# aren't in package directories

import os
root_dir = os.getcwd()
file_set = set()

for dir_, _, files in os.walk(root_dir):
      for file_name in files:
            if ".h" in file_name or ".hpp" in file_name:
                  rel_dir = os.path.relpath(dir_, root_dir)
                  rel_file = os.path.join(rel_dir, file_name)
                  file_set.add(rel_file)

# Run setup, OpenCL kernels are specified as additional data files

setup(name="MCRAMP",
      version="0.0.1a",
      packages=find_packages(),
      package_data={'': ['*.cl']},
      data_files=[('headers', list(file_set))],
      include_package_data=True,
      
      install_requires=[
            'pyopencl',
            'numpy',
            'matplotlib',
            'nbsphinx',
            'numpydoc',
            'pandoc',
            'mcpl'
      ]
)