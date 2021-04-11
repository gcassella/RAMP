Conventions: Scattering kernels
===============================

The building blocks of instruments in RAMP are OpenCL kernels. These mini programs are executed by the RAMP architecture to manipulate the neutron buffer in the order specified by the instrument definition file. Most useful applications of RAMP will require some custom components, and the process of creating user defined kernels that can be used in RAMP is outlined here.

Template
--------

The basic template for a scattering kernel ``mykernel.cl`` is

.. code-block:: C

    #include "consts.h"

    __kernel void myprog(__global float16* neutrons,
        __global float8* intersections, __global uint* iidx,
        uint const comp_idx) {

      uint global_addr        = get_global_id(0);
      float16 neutron         = neutrons[global_addr];
      float8 intersection = intersections[global_addr];
      uint this_iidx          = iidx[global_addr];

      /* Check we are scattering from the intersected component */
      if (!(this_iidx == comp_idx)) {
          return;
      }

      /* Check termination flag */
      if (NEUTRON_DIE > 0.f)  {
          return;
      }

      /* Perform scattering here */
    

      /* ----------------------- */

      /* Update global memory and reset intersection */
      iidx[global_addr] = 0;
      neutron.sc = comp_idx;

      neutrons[global_addr]      = neutron;
      intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                           0.0f, 0.0f, 0.0f, 100000.0f );

    }

The basic template for a scattering kernel Python class ``MyKernel.py`` is

.. code-block:: Python

    from mcramp.scat.sprim import SPrim
    import mcramp
    import numpy as np

    class MyKernel(SPrim):
        def __init__(self, idx=0, ctx=0, **kwargs):

            self.prg = mcramp.build_kernel('mykernel.cl', ctx)
            self.idx = np.uint32(idx)

        def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
            self.prg.myprog(queue, (N, ),
                            None,
                            neutron_buf,
                            intersection_buf,
                            iidx_buf,
                            self.idx)

With these files in the same location as the instrument definition file, this kernel can be included in the instrument using the Python class name ``MyKernel`` in the ``name`` field of a component's scattering kernel entry.

Example: Chopper kernel
-----------------------

The file ``<RAMP dir>/mcramp/scat/chopper.cl`` defines a simple disk chopper kernel.