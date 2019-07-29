from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class SGuide(SPrim):
    def __init__(self, g_pos=(0., 0., 0.), w1=0, h1=0, w2=0, h2=0, l=0, 
                 R0=0, Qc=0, alpha=0, m=1, W=0, idx=0, ctx=0, max_bounces=50):
        self.g_pos  = np.array((g_pos[0], g_pos[1], g_pos[2], 0. ), dtype=clarr.vec.float3)
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
                                self.g_pos,
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
                                self.max_bounces).wait()