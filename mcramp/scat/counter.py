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
        
        self.counts = np.zeros((2,), dtype=np.float32)
        self.filename = filename

        mf = cl.mem_flags
        self.counts_cl = cl.Buffer(ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.counts)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'counter.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.counter(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.counts_cl)
        
        cl.enqueue_copy(queue, self.counts, self.counts_cl)

        if self.filename:
            np.save(self.filename, self.counts)