from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class GBanana(GPrim):
    def __init__(self, radius=0, position=(0, 0, 0), height=0, mintheta=0, maxtheta=0, idx=0, ctx=None):
        self.radius     = np.float32(radius)
        self.height     = np.float32(height)
        self.mintheta   = np.float32(mintheta)
        self.maxtheta   = np.float32(maxtheta)
        self.position   = position
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
                                  self.position,
                                  self.radius,
                                  self.height,
                                  self.mintheta,
                                  self.maxtheta).wait()

    def lines(self):
        lines = []

        fmt = 'b-'

        theta = np.linspace(self.mintheta, self.maxtheta)

        # Top curve
        x = self.radius*np.sin(theta) + self.position['x']
        y = self.height / 2 * np.ones(x.shape) + self.position['y']
        z = self.radius*np.cos(theta) + self.position['z']
        lines.append((x, y, z, fmt))
        # Bottom curve
        y = -self.height / 2 * np.ones(x.shape) + self.position['y']
        lines.append((x, y, z, fmt))
        # Side line A
        y = np.linspace(-self.height / 2, self.height/2) + self.position['y']
        x = (self.radius*np.sin(self.mintheta) + self.position['x']) * np.ones(y.shape)
        z = (self.radius*np.cos(self.mintheta) + self.position['z']) * np.ones(z.shape)
        lines.append((x, y, z, fmt))
        # Side line B
        x = (self.radius*np.sin(self.maxtheta) + self.position['x']) * np.ones(y.shape)
        z = (self.radius*np.cos(self.maxtheta) + self.position['z']) * np.ones(z.shape)
        lines.append((x, y, z, fmt))

        return lines
