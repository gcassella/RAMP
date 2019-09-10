from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class SLinearCollimator(SPrim):
    """
    Scattering kernel for Linear Collimator component. Recreates the functionality
    of the Collimator_linear component in McStas. Neutrons with divergence
    exceeding that permitted by the collimator are terminated.

    Parameters
    ----------
    length : float
        Flight path length of the collimator
    divergence_H : float
        Maximum horizontal divergence accepted by the collimator
    divergence_V : float
        Maximum vertical divergence accepted by the collimator
    transmission : float
        Transmission coefficient of the collimator

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, length=0.0, divergence_H=0.0, divergence_V=0.0, transmission=1.0, idx=0, ctx=0,
                **kwargs):
        
        min2rad = lambda x: x * np.pi / (60.0 * 180.0)
        
        self.length = np.float32(length)
        self.slope_H = np.float32(np.tan(min2rad(divergence_H)))
        self.slope_V = np.float32(np.tan(min2rad(divergence_V)))
        self.transmission = np.float32(transmission)
        self.idx    = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'collimator_lin.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.collimator(queue, (N, ),
                                None,
                                neutron_buf,
                                intersection_buf,
                                iidx_buf,
                                self.idx,
                                self.length,
                                self.slope_H,
                                self.slope_V,
                                self.transmission)