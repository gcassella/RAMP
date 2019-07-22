from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class PSD2dMon(SPrim):
    def __init__(self, sample_pos=(0, 0, 0), shape="", axis1_binning=(0, 0, 0),
                 axis2_binning=(0, 0, 0), idx=0, ctx=None):
        
        shapes = {"plane" : 0, "banana": 1}

        self.axis1_binning = axis1_binning
        self.axis2_binning = axis2_binning
        self.sample_pos = sample_pos
        self.shape = np.uint32(shapes[shape])
        self.idx = np.uint32(idx)

        self.axis1_num_bins = np.uint32(np.ceil((axis1_binning[2] - axis1_binning[0]) / axis1_binning[1]))
        self.axis2_num_bins = np.uint32(np.ceil((axis2_binning[2] - axis2_binning[0]) / axis2_binning[1]))
        self.num_bins = np.uint32(self.axis1_num_bins * self.axis2_num_bins)
        self.histo = np.zeros((self.num_bins,), dtype=np.float32)
        self.histo2d = np.zeros((self.axis1_num_bins, self.axis2_num_bins))

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                    mf.WRITE_ONLY,
                                    self.histo.nbytes)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'psd2d_mon.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.histo_cl,
                          self.sample_pos,
                          self.axis1_binning,
                          self.axis2_binning,
                          self.axis1_num_bins,
                          self.axis2_num_bins,
                          self.shape)

        cl.enqueue_copy(queue, self.histo, self.histo_cl)
        self.histo2d = self.histo.reshape((self.axis1_num_bins, self.axis2_num_bins))

    @property
    def sample_pos(self):
        return self._sample_pos

    @sample_pos.setter
    def sample_pos(self, val):
        self._sample_pos = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

    @property
    def axis1_binning(self):
        return self._axis1_binning

    @axis1_binning.setter
    def axis1_binning(self, val):
        self._axis1_binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

    @property
    def axis2_binning(self):
        return self._axis2_binning

    @axis2_binning.setter
    def axis2_binning(self, val):
        self._axis2_binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)