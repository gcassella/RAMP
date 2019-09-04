import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

from random import randint

import os
import re

class MCont():
    def __init__(self, ctx=None, spec_file="Let_Base.mcstas",
                 mod_dim=(0.0, 0.0), target_dim=(0.0, 0.0), target_dist=0.0,
                 E_min=0.0, E_max=0.0, T1=0.0, I1=0.0, T2=0.0, I2=0.0,
                 T3 = 0.0, I3 = 0.0):

        self.mod_dim = np.array((mod_dim[0], mod_dim[1]), dtype=clarr.vec.float2)
        self.target_dim = np.array((target_dim[0], target_dim[1]), dtype=clarr.vec.float2)
        self.target_dist = np.float32(target_dist)
        self.E_min = np.float32(E_min)
        self.E_max = np.float32(E_max)

        self.T1 = np.float32(T1)
        self.I1 = np.float32(I1)
        self.T2 = np.float32(T2)
        self.I2 = np.float32(I2)
        self.T3 = np.float32(T3)
        self.I3 = np.float32(I3)

        # self.calc_str_area(mod_dim, target_dim, target_dist)
        self.str_area = np.float32(1)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cont_mod.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
            
    def calc_str_area(self, mod_dim, target_dim, target_dist):
        n_steps = 50
        A = 0.0

        target_x = np.linspace(0, target_dim[0] / 2.0, num=n_steps)
        target_y = np.linspace(0, target_dim[1] / 2.0, num=n_steps)

        TX, TY = np.meshgrid(target_x, target_y)

        mod_x = np.linspace(-mod_dim[0] / 2.0, mod_dim[0] / 2.0, num=n_steps)
        mod_y = np.linspace(-mod_dim[1] / 2.0, mod_dim[1] / 2.0, num=n_steps)
        MX, MY = np.meshgrid(mod_x, mod_y)
        for mx in MX:
            for my in MY:
                A += np.sum(1 / (np.power(-1.0*TX + mx, 2.0) + np.power(-1.0*TY + my, 2.0) + target_dist**2.0))
        A *= (mod_dim[0] * mod_dim[1]) / (n_steps**4.0)
        A *= target_dim[0] * target_dim[1] * 10000

        self.str_area = np.float32(A)

    def gen_prg(self, queue, N, neutron_buf, intersection_buf):
        self.prg.generate_neutrons(queue, (N,), None, neutron_buf, intersection_buf,
                                   self.mod_dim,
                                   self.target_dim,
                                   self.target_dist,
                                   self.E_min,
                                   self.E_max,
                                   self.T1,
                                   self.I1,
                                   self.T2,
                                   self.I2,
                                   self.T3,
                                   self.I3,
                                   self.str_area,
                                   np.int32(N),
                                   np.int32(randint(0,4096)))
