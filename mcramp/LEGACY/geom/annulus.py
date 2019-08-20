from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class GAnnulus(GPrim):
    def __init__(self, inner_radius=0, outer_radius=0, position=(0, 0, 0), height=0, idx=0, ctx=None):
        self.inner_radius     = np.float32(inner_radius)
        self.outer_radius     = np.float32(outer_radius)
        self.height     = np.float32(height)
        self.position   = position
        self.idx        = idx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'annulus.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_annulus(queue, (N, ),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.position,
                                  self.inner_radius,
                                  self.outer_radius,
                                  self.height).wait()

    def lines(self):
        lines = []

        fmt = 'b-'

        theta = np.linspace(0, 2*np.pi, num=128)

        # Top curve
        x = self.inner_radius*np.sin(theta) + self.position['x']
        y = self.height / 2 * np.ones(x.shape) + self.position['y']
        z = self.inner_radius*np.cos(theta) + self.position['z']
        lines.append((x, y, z, fmt))
        # Bottom curve
        y = -self.height / 2 * np.ones(x.shape) + self.position['y']
        lines.append((x, y, z, fmt))

        return lines
