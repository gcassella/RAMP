from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class SPowder1(SPrim):
    def __init__(self, d_spacing=0.0, pack=0.0, vc=0.0, sigma_abs=0.0,
                 multiplicity=0, DW=0.0, F2=0.0, d_phi=180.0, idx=0, ctx=None):
        # F2 in barns         

        # Scattering cross section will be divided by v^2 in Kernel
        # Absorption cross section will be divided by v in Kernel
        q = (2 * np.pi) / d_spacing

        sigma_scat_numer = 4 * np.pi ** 3.0 * pack * multiplicity * F2 * 100 * DW
        sigma_scat_denom = vc ** 2.0 * q

        self.sigma_scat_v2 = np.float32(sigma_scat_numer / sigma_scat_denom)
        self.sigma_abs_v = np.float32(pack * sigma_abs / vc * 100)
        self.q = np.float32(q)
        self.idx = np.uint32(idx)
        self.d_phi = np.float32(d_phi)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'powder1.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.powder1(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.sigma_scat_v2,
                          self.sigma_abs_v,
                          self.q,
                          self.d_phi)