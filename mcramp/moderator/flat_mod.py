import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

from random import randint

import os
import re

class MFlat():
    def __init__(self, ctx=None, spec_file="Let_Base.mcstas",
                 mod_dim=(0.0, 0.0), target_dim=(0.0, 0.0), target_dist=0.0,
                 E_min=0.0, E_max=0.0):

        self.mod_dim = np.array((mod_dim[0], mod_dim[1]), dtype=clarr.vec.float2)
        self.target_dim = np.array((target_dim[0], target_dim[1]), dtype=clarr.vec.float2)
        self.target_dist = np.float32(target_dist)
        self.E_min = np.float32(E_min)
        self.E_max = np.float32(E_max)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'flat_mod.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def gen_prg(self, queue, N, neutron_buf, intersection_buf):
        self.prg.generate_neutrons(queue, (N,), None, neutron_buf, intersection_buf,
                                   self.mod_dim,
                                   self.target_dim,
                                   self.target_dist,
                                   self.E_min,
                                   self.E_max,
                                   np.int32(N),
                                   np.int32(randint(0,4096)))
