from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os

class PSD2d(SPrim):
    def __init__(self, shape="", axis1_binning=(0, 0, 0),
                 axis2_binning=(0, 0, 0), restore_neutron=False, idx=0, ctx=None,
                 filename=None, logscale = False):
        
        shapes = {"plane" : 0, "banana": 1, "thetatof": 2, "div" : 3, "divpos": 4}

        self.axis1_binning = axis1_binning
        self.axis2_binning = axis2_binning
        self.shape = np.uint32(shapes[shape])
        self.idx = np.uint32(idx)
        self.restore_neutron = np.uint32(1 if restore_neutron else 0)
        self.filename = filename
        self.logscale = logscale

        self.axis1_num_bins = np.uint32(np.ceil((axis1_binning[2] - axis1_binning[0]) / axis1_binning[1]))
        self.axis2_num_bins = np.uint32(np.ceil((axis2_binning[2] - axis2_binning[0]) / axis2_binning[1]))
        self.num_bins = np.uint32(self.axis1_num_bins * self.axis2_num_bins)
        self.histo = np.zeros((self.num_bins,), dtype=np.float32)
        self.histo2d = np.zeros((self.axis1_num_bins, self.axis2_num_bins))

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.histo)
        
        x = np.linspace(self.axis1_binning['s0'], self.axis1_binning['s2'], num=self.axis1_num_bins)
        y = np.linspace(self.axis2_binning['s0'], self.axis2_binning['s2'], num=self.axis2_num_bins)

        self.X, self.Y = np.meshgrid(x, y)
        self.Z = np.zeros(self.histo2d.T.shape)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'psd2d.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.histo_cl,
                          self.axis1_binning,
                          self.axis2_binning,
                          self.axis1_num_bins,
                          self.axis2_num_bins,
                          self.shape,
                          self.restore_neutron)

    def plot_histo(self, queue):
        cl.enqueue_copy(queue, self.histo, self.histo_cl).wait()
        self.Z = self.histo.reshape((self.axis1_num_bins, self.axis2_num_bins)).T

        plt.figure()
        Z = (np.log(self.Z + 1e-7) if self.logscale else self.Z)

        plt.pcolormesh(self.X, self.Y, Z, cmap='jet', shading='gouraud')
        plt.colorbar()

    def save_histo(self, queue):
        cl.enqueue_copy(queue, self.histo, self.histo_cl).wait()
        self.Z = self.histo.reshape((self.axis1_num_bins, self.axis2_num_bins)).T

        if self.filename:
            np.save(self.filename + 'X.dat', self.X)
            np.save(self.filename + 'Y.dat', self.Y)
            np.save(self.filename + 'Z.dat', self.Z)

    def get_histo(self):
        return (self.X, self.Y, self.Z)

    def sum_histo(self):
        self.Z += self.histo2d.T

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