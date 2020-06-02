from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class GSphere(GPrim):
    """
    Geometry kernel for 'sphere' geometry. Intersects with the exterior of the
    sphere, i.e. first intersection time must be positive for scattering to
    occur.

    Parameters
    ----------
    radius : float
        The radius of the sphere

    Notes
    -----
    Intersection 1 :
        First point of intersection with the sphere geometry - 'entering' sphere.
    Intersection 2 :
        Second point of intersection with the sphere geometry - 'exiting' sphere.

    Methods
    -------
    None
    """

    def __init__(self, radius=0, idx=0, ctx=None):
        self.radius     = np.float32(radius)
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
                                  self.radius)

    def lines(self):
        x = []
        y = []
        z = []

        angles = np.linspace(0, 2*np.pi)
        for ang in angles:
            x.append(0)
            y.append(self.radius*np.cos(ang))
            z.append(self.radius*np.sin(ang))

        x.append(x[0])
        y.append(y[0])
        z.append(z[0])

        for ang in angles:
            x.append(self.radius*np.sin(ang))
            y.append(self.radius*np.cos(ang))
            z.append(0)

        x.append(x[0])
        y.append(y[0])
        z.append(z[0])

        for ang in np.linspace(0, np.pi / 2):
            x.append(self.radius*np.sin(ang))
            y.append(self.radius*np.cos(ang))
            z.append(0)

        for ang in angles:
            x.append(self.radius*np.cos(ang))
            y.append(0)
            z.append(self.radius*np.sin(ang))

        return [x, y, z]
