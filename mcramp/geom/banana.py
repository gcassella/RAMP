from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class GBanana(GPrim):
    def __init__(self, radius=0, height=0, mintheta=0, maxtheta=0, idx=0, ctx=None):
        self.radius     = np.float64(radius)
        self.height     = np.float64(height)
        self.mintheta   = np.float64(mintheta)
        self.maxtheta   = np.float64(maxtheta)
        self.idx        = idx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'banana.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_banana(queue, (N, ),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.radius,
                                  self.height,
                                  self.mintheta,
                                  self.maxtheta)