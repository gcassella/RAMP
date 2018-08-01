from .sprim import SPrim, Float3

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class QWDetector(SPrim):
    position   = Float3('position')
    ki_norm    = Float3('ki_norm')
    q_binning  = Float3('q_binning')
    de_binning = Float3('de_binning')

    def __init__(self, position=(0, 0, 0), ei=0, ki_norm=(0, 0, 1),
                  q_binning=(0, 0, 0), de_binning=(0, 0, 0), idx=0, ctx=None):
        self.position    = position
        self.idx         = idx

        self.ei          = np.float32(ei)
        self.ki_norm     = ki_norm
        self.q_binning   = q_binning
        self.de_binning  = de_binning

        self.num_q_bins    = np.ceil((q_binning[2] - q_binning[0])/q_binning[1]).astype(np.uint32)
        self.num_de_bins   = np.ceil((de_binning[2] - de_binning[0])/de_binning[1]).astype(np.uint32)
        self.histo       = np.zeros((self.num_q_bins*self.num_de_bins, ), dtype=np.float32)

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                    mf.WRITE_ONLY,
                                    self.histo.nbytes)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'qwdetector.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.histo_cl,
                          self.position,
                          self.ei,
                          self.ki_norm,
                          self.q_binning,
                          self.de_binning).wait()

        cl.enqueue_copy(queue, self.histo, self.histo_cl).wait()

    def reduce_histogram(self):
        q = np.linspace(self.q_binning['x'], self.q_binning['z'], self.num_q_bins)
        de = np.linspace(self.de_binning['x'], self.de_binning['z'], num=self.num_de_bins)

        Q, DE = np.meshgrid(q, de, indexing='ij')
        return (Q, DE, self.histo.reshape((self.num_q_bins, self.num_de_bins)))