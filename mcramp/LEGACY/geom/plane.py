from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class GPlane(GPrim):
    def __init__(self, width=0, height=0, position=[0, 0, 0], idx=0, ctx=None):
        self.position   = position
        self.width      = np.float64(width)
        self.height     = np.float64(height)
        self.idx        = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plane.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_plane(queue, (N, ),
                                 None,
                                 neutron_buf,
                                 intersection_buf,
                                 iidx_buf,
                                 self.idx,
                                 self.position,
                                 self.width,
                                 self.height).wait()
