from .sprim import SPrim

import pyopencl as cl
import pyopencl.array as clarr
import numpy as np
import os

class SBeamstop(SPrim):
    """
    Scattering kernel for the beamstop component. Absorbs neutrons.

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

    def __init__(self, idx=0, ctx=None, filename=None, **kwargs):
        self.idx = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'beamstop.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.beamstop(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx)