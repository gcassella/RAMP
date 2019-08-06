from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class SDelta(SPrim):
    def __init__(self, twotheta=0.0, phi=0.0, deltaE=0.0,
                 idx=0, ctx=None):

        self.twotheta = np.float32(np.radians(twotheta))
        self.phi = np.float32(np.radians(phi))
        self.deltaE = np.float32(deltaE)
        self.idx = np.uint32(idx)
    
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'delta.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.delta(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.twotheta,
                          self.phi,
                          self.deltaE).wait()