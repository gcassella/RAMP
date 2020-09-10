from .sprim import SPrim

import os
import pyopencl as cl
import numpy as np

class SArm(SPrim):
    """
    Scattering kernel for Arm component - does not alter the neutron state.

    Parameters
    ----------
    None

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, idx=0, ctx=None, **kwargs):
        self.idx = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'arm.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.arm(queue, (N, ),
                      None,
                      neutron_buf,
                      intersection_buf,
                      iidx_buf,
                      self.idx)