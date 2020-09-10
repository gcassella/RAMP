from .sprim import SPrim

import pyopencl as cl
import pyopencl.array as clarr
import numpy as np
import os

class SPolariser(SPrim):
    """
    Scattering kernel for ideal polariser component. Sets neutron
    polarisation to value of polarisation parameter.

    Parameters
    ----------
    polarisation : 3-tuple of floats
        Value to set neutron polarisation to
    

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, polarisation = [0, 0, 0], idx=0, ctx=None, filename=None, **kwargs):
        mf = cl.mem_flags
        self.idx = np.uint32(idx)

        self.polarisation = np.array((polarisation[0], polarisation[1], polarisation[2], 0.0), dtype=clarr.vec.float3)
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'polariser.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.polariser(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.polarisation)