from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class GBanana(GPrim):
    def __init__(self, radius, position, height, mintheta, maxtheta, idx, ctx):
        self.radius     = np.float32(radius)
        self.height     = np.float32(height)
        self.mintheta   = np.float32(mintheta)
        self.maxtheta   = np.float32(maxtheta)
        self.position   = position
        self.idx        = idx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'banana.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_banana(queue, (N, ),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.position,
                                  self.radius,
                                  self.height,
                                  self.mintheta,
                                  self.maxtheta).wait()
