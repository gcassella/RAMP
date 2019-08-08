from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os

class SCounter(SPrim):
    def __init__(self, idx=0, ctx=None,
                 filename=None):
        self.idx = idx
        
        self.counts = np.float32(0.0)
        self.filename = filename


        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'counter.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.counter(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx))
           
        neutrons = np.zeros((N, ), dtype=clarr.vec.float16)
        cl.enqueue_copy(queue, neutrons, neutron_buf).wait()

        counted = np.where((neutrons['s12'].astype(np.uint32) == self.idx))
        self.counts = np.sum(neutrons['s9'][counted]).astype(np.float32)

        print(self.counts)

        if self.filename:
            np.save(self.filename, self.counts)