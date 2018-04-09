import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class SPowder():
    def __init__(self, fn, idx, ctx):
        reflections = self._LoadLAZ(fn)

        self.nreflections = len(reflections)
        self.idx = idx
        mf = cl.mem_flags

        self.reflections_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=reflections)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'powder.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.powder_scatter(queue, (N,),
                                   None,
                                   neutron_buf,
                                   intersection_buf,
                                   iidx_buf,
                                   np.uint32(self.idx),
                                   self.reflections_opencl,
                                   np.uint32(self.nreflections)).wait()

    def _LoadLAZ(self, fn):
        reflections_temp = np.loadtxt(fn)

        reflections = np.empty((0,), dtype=clarr.vec.float3)

        sum_intensity = sum(reflections_temp[:,11])

        reflections = np.array([(ref[4], ref[5], ref[11]/sum_intensity, 0.) for ref in reflections_temp],
                               dtype=clarr.vec.float3 )

        return reflections
