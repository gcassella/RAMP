from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class MantidDetector(SPrim):
    def __init__(self, position=(0, 0, 0), y_binning=(0, 0, 0),
                 theta_binning=(0, 0, 0), idx=0, ctx=None, N=0):
        self.theta_binning  = theta_binning
        self.y_binning      = y_binning
        self.position       = position
        self.idx            = idx
        self.ctx            = ctx

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'detector.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.queue = queue
        self.N = N

        self.theta_histo = np.zeros((self.N, ), dtype=np.float32)
        self.y_histo = np.zeros((self.N, ), dtype=np.float32)
        
        mf               = cl.mem_flags
        self.theta_histo_cl    = cl.Buffer(self.ctx,
                                    mf.WRITE_ONLY,
                                    size=4*N)

        self.y_histo_cl    = cl.Buffer(self.ctx,
                                    mf.WRITE_ONLY,
                                    size=4*N)

        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.y_histo_cl,
                          self.position,
                          self.y_binning,
                          np.uint32(2),
                          np.uint32(1),
                          np.uint32(1))


        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.theta_histo_cl,
                          self.position,
                          self.theta_binning,
                          np.uint32(0),
                          np.uint32(0),
                          np.uint32(1))
    
    def get_histo(self):
        cl.enqueue_copy(self.queue, self.theta_histo, self.theta_histo_cl).wait()
        cl.enqueue_copy(self.queue, self.y_histo, self.y_histo_cl).wait()
        return (self.y_histo, self.theta_histo)

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, val):
        self._position = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

    @property
    def y_binning(self):
        return self._y_binning

    @property
    def theta_binning(self):
        return self._theta_binning

    @y_binning.setter
    def y_binning(self, val):
        self._y_binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

    
    @theta_binning.setter
    def theta_binning(self, val):
        self._theta_binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)