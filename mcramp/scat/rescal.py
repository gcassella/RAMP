from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class SRescal(SPrim):
    def __init__(self, target=[0,0,0], E0=0, dE=0, focus_r=0, idx=0, ctx=None, **kwargs):
        self.idx = np.uint32(idx)
        self.ctx = ctx

        self.target = np.array((target[0], target[1], target[2], 0.),
                                 dtype=clarr.vec.float3)
        self.E0=np.float32(E0)
        self.dE=np.float32(dE)
        self.focus_r=np.float32(focus_r)

        self.qz_vals=np.array([])
        self.qy_vals=np.array([])
        self.qx_vals=np.array([])
        self.dE_vals=np.array([])
        self.p_vals =np.array([])

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rescal.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.events = np.zeros((N,), dtype=clarr.vec.float8)
        mf               = cl.mem_flags
        self.events_cl    = cl.Buffer(self.ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.events)

        self.neutron_buf = neutron_buf
        self.N = N

        self.prg.rescal(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.target,
                          self.E0,
                          self.dE,
                          self.focus_r,
                          self.events_cl)

    def data_reduce(self, queue):
        neutrons = np.zeros((self.N, ), dtype=clarr.vec.float16)
        
        cl.enqueue_copy(queue, self.events, self.events_cl)
        cl.enqueue_copy(queue, neutrons, self.neutron_buf).wait()

        eventlist_reduced = self.events[np.where(neutrons["s14"] > 0.0)]

        self.qz_vals=np.concatenate((self.qz_vals, eventlist_reduced["s2"]))
        self.qy_vals=np.concatenate((self.qy_vals, eventlist_reduced["s1"]))
        self.qx_vals=np.concatenate((self.qx_vals, eventlist_reduced["s0"]))
        self.dE_vals=np.concatenate((self.dE_vals, eventlist_reduced["s7"]))
        self.p_vals =np.concatenate((self.p_vals, eventlist_reduced["s6"]))

    def data(self, queue):
        return (self.qx_vals, self.qy_vals, self.qz_vals, self.dE_vals, self.p_vals)