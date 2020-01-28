from .gprim import GPrim  #pylint: disable=E0401

import os
import pyopencl as cl
import numpy as np

class GArm(GPrim):
    """
    Geometry kernel for Arm component - does not store a real intersection.

    Parameters
    ----------
    None

    Notes
    -----
    None

    Methods
    -------
    None

    """

    def __init__(self, idx=0, ctx=None):
        self.idx = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'arm.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.arm(queue, (N, ),
                      None,
                      neutron_buf,
                      intersection_buf,
                      iidx_buf,
                      self.idx)