from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class GBox(GPrim):
    """
    Geometry kernel for 'box' geometry. Intersects with the exterior of the
    box, i.e. first intersection must be at positive time for scattering to occur.

    ...

    Parameters
    ----------
    width : float
        The width of the box
    height : float
        The height of the box
    depth : float
        The depth of the box

    Notes
    -----
    Intersection 1 :
        First point of intersection with the box geometry - 'entering' box.
    Intersection 2 :
        Second point of intersection with the box geometry - 'exiting' box.

    Methods
    -------
    None
    """

    def __init__(self, width=0, height=0, depth=0, idx=0, ctx=None):
        self.width = np.float32(width)
        self.height = np.float32(height)
        self.depth = np.float32(depth)
        self.idx        = idx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'box.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_box(queue, (N,),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.width,
                                  self.height,
                                  self.depth)