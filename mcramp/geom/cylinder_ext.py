from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class GCylinderExt(GPrim):
    """
    Geometry kernel for cylinder geometry. Intersects with the exterior of the
    cylinder, i.e. first intersection must be at positive time for scattering to
    occur.

    ...

    Parameters
    ----------
    radius : float
        The radius of the cylinder
    height : float
        The height of the cylinder

    Notes
    -----
    Intersection 1 :
        Negative time intersection with the portion of the banana behind the\
        neutron trajectory.
    Intersection 2 :
        Positive time intersection with the portion of the banana ahead of the\
        neutron trajectory.

    Methods
    -------
    None
    """

    def __init__(self, radius=0, height=0, idx=0, ctx=None):
        self.radius     = np.float32(radius)
        self.height     = np.float32(height)
        self.idx        = idx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cylinder_ext.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect(queue, (N, ),
                           None,
                           neutron_buf,
                           intersection_buf,
                           iidx_buf,
                           np.uint32(self.idx),
                           self.radius,
                           self.height)

    def lines(self):
        angles = np.linspace(0, np.pi)
        h_2 = self.height / 2.0

        x = []
        y = []
        z = []

        for ang in angles:
            x.append(self.radius*np.sin(ang))
            y.append(h_2)
            z.append(self.radius*np.cos(ang))

        for ang in reversed(angles):
            x.append(self.radius*np.sin(ang))
            y.append(-h_2)
            z.append(self.radius*np.cos(ang))

        x.append(x[0])
        y.append(y[0])
        z.append(z[0])

        angles = np.linspace(0, -np.pi)

        for ang in angles:
            x.append(self.radius*np.sin(ang))
            y.append(h_2)
            z.append(self.radius*np.cos(ang))

        for ang in reversed(angles):
            x.append(self.radius*np.sin(ang))
            y.append(-h_2)
            z.append(self.radius*np.cos(ang))

        x.append(x[0])
        y.append(y[0])
        z.append(z[0])

        return [x, y, z]

        
