from .sprim import SPrim

import pyopencl as cl
import pyopencl.array as clarr
import numpy as np
import os

class SAnalyser(SPrim):
    """
    Scattering kernel for ideal analyser component. Adjusts neutron weight
    based on dot product between neutron polarisation and polarisation
    parameter, and sets neutron polarisation to the analysed direction.

    Parameters
    ----------
    polarisation : 3-tuple of floats
        Vector against which neutron polarisation should be analysed
    

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
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'analyser.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.analyser(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.polarisation)