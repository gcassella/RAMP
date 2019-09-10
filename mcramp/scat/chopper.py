from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class SChopper(SPrim):
    """
    Scattering kernel for Chopper component. Replicates the functionality of the
    DiskChopper component in McStas. Neutrons with velocities that are not
    permitted by the chopper are terminated.

    Parameters
    ----------
    slit_width : float
        Width of the chopper slits in meters
    radius : float
        Radius of the chopper disc in meters
    freq : float
        Angular frequency of the chopper in radians per second
    n_slits : float
        Number of chopper slits
    phase : float
        Initial phase of the chopper in radians
    jitter : float
        Jitter in the chopper phase in radians

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, slit_width=0.0, radius=0.0, freq=0.0,
                 n_slits = 0, phase = 0.0, jitter = 0.0, idx=0, ctx=0, **kwargs):

        self.slit_width = np.float32(slit_width)
        self.radius = np.float32(radius)
        self.freq = np.float32(freq)
        self.n_slits = np.uint32(n_slits)
        self.phase = np.float32(phase)
        self.jitter = np.float32(jitter)

        self.idx = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'chopper.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.chopper(queue, (N, ),
                                None,
                                neutron_buf,
                                intersection_buf,
                                iidx_buf,
                                self.idx,
                                self.slit_width,
                                self.radius,
                                self.freq,
                                self.n_slits,
                                self.phase,
                                self.jitter)