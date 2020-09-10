from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class SGuide(SPrim):
    """
    Scattering kernel for tapered rectangular Guide. Recreates the functionality
    of the Guide component in McStas. The path of the neutron through the guide
    is numerically simulated and its weight adjusted according to the reflectivity
    function of the guide walls.

    Intersection is taken as the point at which the neutron enters the guide and
    the guide geometry is taken to lie centered along the z axis.

    Parameters
    ----------
    w1 : float
        Width of the guide entrance in meters
    h1 : float
        Height of the guide entrance in meters
    w2 : float
        Width of the guide exit in meters
    h2 : float
        Height of the guide exit in meters
    l : float
        Length of the guide in meters
    R0 : float
        Low-angle reflectivity of the guide
    Qc : float
        Critical scattering vector of the guide
    alpha : float
        Slope of the reflectivity
    m : float
        m-value of the guide coating
    W : float
        Width of the guide supermirror cutoff
    max_bounces : float
        Cutoff to prevent infinite scattering due to numerical error in the kernel

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, w1=0, h1=0, w2=0, h2=0, l=0, 
                 R0=0, Qc=0, alpha=0, m=1, W=0, idx=0, ctx=0, max_bounces=50,
                 **kwargs):

        if "CONFIG_USE_DOUBLE=1" in os.environ["PYOPENCL_BUILD_OPTIONS"]:
            self.w1     = np.float64(w1)
            self.h1     = np.float64(h1)
            self.w2     = np.float64(w2)
            self.h2     = np.float64(h2)
            self.l      = np.float64(l)
            self.R0     = np.float64(R0)
            self.Qc     = np.float64(Qc)
            self.alpha  = np.float64(alpha)
            self.m      = np.float64(m)
            self.W      = np.float64(W)
        else:
            self.w1     = np.float32(w1)
            self.h1     = np.float32(h1)
            self.w2     = np.float32(w2)
            self.h2     = np.float32(h2)
            self.l      = np.float32(l)
            self.R0     = np.float32(R0)
            self.Qc     = np.float32(Qc)
            self.alpha  = np.float32(alpha)
            self.m      = np.float32(m)
            self.W      = np.float32(W)
        self.idx    = np.uint32(idx)
        self.max_bounces = np.uint32(max_bounces)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'guide.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.guide_scatter(queue, (N, ),
                                None,
                                neutron_buf,
                                intersection_buf,
                                iidx_buf,
                                self.idx,
                                self.w1,
                                self.h1,
                                self.w2,
                                self.h2,
                                self.l,
                                self.R0,
                                self.Qc,
                                self.alpha,
                                self.m,
                                self.W,
                                self.max_bounces)

    def lines(self):
        w1_2 = self.w1 / 2
        h1_2 = self.h1 / 2
        w2_2 = self.w2 / 2
        h2_2 = self.h2 / 2
        l = self.l

        x_arr = [-w1_2, -w2_2, w2_2, w1_2, -w1_2, -w1_2, w1_2, w1_2, w2_2, w2_2, w1_2, -w1_2, -w2_2, w2_2, -w2_2, -w2_2]

        y_arr = [h1_2, h2_2, h2_2, h1_2, h1_2, -h1_2, -h1_2, h1_2, h2_2, -h2_2, -h1_2, -h1_2, -h2_2, -h2_2, -h2_2, h2_2]
        
        z_arr = [0, l, l, 0, 0, 0, 0, 0, l, l, 0, 0, l, l, l, l]

        return [x_arr, y_arr, z_arr]

