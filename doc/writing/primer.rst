(very) Brief OpenCL Primer
==========================

A full account of the OpenCL kernel language is beyond the scope of this document, however some basics that are important for RAMP will be outlined here for reference. For a more complete account of the OpenCL API, `see here <https://www.khronos.org/opencl/>`_. I have found `this reference card <https://www.khronos.org/files/opencl-1-2-quick-reference-card.pdf>`_ exceptionally helpful, and have learned almost everything I needed to know to write RAMP so far using only my knowledge of C and this reference card, as the OpenCL kernel language closely mirrors C99 and the majority of syntax is identical.

Manipulating OpenCL objects from Python is accomplished via the PyOpenCL library. The documentation for PyOpenCL can be found `here <https://documen.tician.de/pyopencl/>`_. 

Kernel program declarations
---------------------------

While normal C99 functions can be declared in OpenCL kernels as usual, the programs that are invoked from the host device must be prepended with the ``__kernel`` keyword.

In addition to this, there are some important keywords to bear in mind when declaring program arguments. Buffers should be passed as ``__global`` pointers (global memory is not the only option in OpenCL, but it is the only option we will care about here). Kernel parameters should be passed as ``const``.

Buffers
-------

Data which is to be operated on in a parallel fashion (i.e. a kernel invocation for each entry in an array) is stored in a buffer. A good example of how buffers are instantiated from the host device can be found in the `Detector1D scattering kernel <>`_.

Vectors and swizzling
---------------------

A number of non-standard vector types are natively supported by OpenCL. The most relevant of these for RAMP are ``float3``, ``float4``, ``float8``, and ``float16``.

Standard vector operations e.g. ``+, -, *``, ``dot(v1, v2)``, ``cross(v1, v2)`` exist for these types.

Vector elements are accessed as ``v.sX`` where X is a number (zero-indexed) corresponding to the element to access. For example, if ``float3 v = (float3){ 0.0f, 1.0f, 0.2f }`` then ``v.s0 = 0.0f``, ``v.s1 = 1.0f``, and ``v.s2 = 0.2f``. Vector indices can also be 'swizzled', this means that sub-vectors can be sliced by arranging element indices in any order e.g. ``float2 w = v.s20`` yields ``w.s0 = 0.2f`` and ``w.s1 = 0.0f``.