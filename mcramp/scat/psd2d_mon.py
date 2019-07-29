from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os

class PSD2dMon(SPrim):
    def __init__(self, sample_pos=(0, 0, 0), shape="", axis1_binning=(0, 0, 0),
                 axis2_binning=(0, 0, 0), idx=0, ctx=None):
        
        shapes = {"plane" : 0, "banana": 1, "thetatof": 2}

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

        neutrons = np.zeros((N, ), dtype=clarr.vec.float16)
        cl.enqueue_copy(queue, neutrons, neutron_buf).wait()

        counted = np.where((neutrons['s14'] > 0) & (neutrons['s12'].astype(np.uint32) == self.idx))
        self.histo, _ = np.histogram(neutrons['s14'][counted], bins=range(self.num_bins + 1), weights=neutrons['s9'][counted])
        self.histo2d = self.histo.reshape((self.axis1_num_bins, self.axis2_num_bins))

        self.plot_histo()

    def plot_histo(self):
        plt.figure()
        x = np.linspace(self.axis1_binning['s0'], self.axis1_binning['s2'], num=self.axis1_num_bins)
        y = np.linspace(self.axis2_binning['s0'], self.axis2_binning['s2'], num=self.axis2_num_bins)

        X, Y = np.meshgrid(x, y)

        plt.pcolormesh(X, Y, self.histo2d.T, cmap='jet', shading='gouraud')
        plt.colorbar()

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