from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class GBanana(GPrim):
    """
    Geometry kernel for 'banana' geometry. Intersects with the interior of the
    banana, i.e. first intersection must be at negative time for scattering to
    occur.

    ...

    Parameters
    ----------
    radius : float
        The radius of the banana
    height : float
        The height of the banana
    mintheta : float
        The minimum valid intersection angle with the banana
    maxtheta : float
        The maximum valid intersection angle with the banana

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

    def __init__(self, radius=0, height=0, mintheta=0, maxtheta=0, idx=0, ctx=None):
        self.radius     = np.float32(radius)
        self.height     = np.float32(height)
        self.mintheta   = np.float32(mintheta)
        self.maxtheta   = np.float32(maxtheta)
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

    def lines(self):
        radial_intervals = 50
        mintheta_r = np.radians(self.mintheta)
        maxtheta_r = np.radians(self.maxtheta)
        a_step = (maxtheta_r - mintheta_r) / float(radial_intervals)
        x = []
        y = []
        z = []
        for i in reversed(range(radial_intervals)):
            y.append(self.height / 2.0)
            x.append(self.radius*np.sin(mintheta_r + i*a_step))
            z.append(self.radius*np.cos(mintheta_r + i*a_step))

        for i in range(radial_intervals):
            y.append(-self.height / 2.0)
            x.append(self.radius*np.sin(mintheta_r + i*a_step))
            z.append(self.radius*np.cos(mintheta_r + i*a_step))

        x.append(x[0])
        y.append(y[0])
        z.append(z[0])

        return [x, y, z]
