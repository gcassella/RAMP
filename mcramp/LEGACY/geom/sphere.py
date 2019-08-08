from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class GSphere(GPrim):
    def __init__(self, radius=0, position=(0, 0, 0), idx=0, ctx=None):
        self.radius     = np.float64(radius)
        self.position   = position
        self.idx        = idx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sphere.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_sphere(queue, (N,),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.position,
                                  self.radius).wait()

    def lines(self):
        lines = []

        fmt = 'r-'

        theta = np.linspace(0, 2*np.pi)
        x = self.radius*np.cos(theta) + self.position['x']
        y = self.radius*np.sin(theta) + self.position['y']
        z = self.position['z'] * np.ones(theta.shape)
        lines.append((x,y,z, fmt))

        x = self.radius*np.cos(theta) + self.position['x']
        y = self.position['y'] * np.ones(theta.shape)
        z = self.radius*np.sin(theta) + self.position['z']
        lines.append((x, y, z, fmt))

        x = self.position['x'] * np.ones(theta.shape)
        y = self.radius*np.cos(theta) + self.position['y']
        z = self.radius*np.sin(theta) + self.position['z']
        lines.append((x, y, z, fmt))

        return lines