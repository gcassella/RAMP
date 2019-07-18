from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class EMon(SPrim):
    def __init__(self, binning=(0, 0, 0), idx=0, ctx=None):
        self.binning     = binning
        self.idx         = idx

        self.num_bins    = np.ceil((binning[2] - binning[0])/binning[1]).astype(np.uint32)
        self.histo       = np.zeros((self.num_bins, ), dtype=np.float32)

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                    mf.WRITE_ONLY,
                                    self.histo.nbytes)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e_mon.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.histo_cl,
                          self.binning)

        cl.enqueue_copy(queue, self.histo, self.histo_cl, is_blocking=True).wait()

    @property
    def binning(self):
        return self._binning

    @binning.setter
    def binning(self, val):
        self._binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)