from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os

class EMon(SPrim):
    def __init__(self, binning=(0, 0, 0), restore_neutron=False, idx=0, ctx=None,
                 filename=None):
        self.binning     = binning
        self.idx         = idx

        self.num_bins    = np.ceil((binning[2] - binning[0])/binning[1]).astype(np.uint32)
        self.histo = np.zeros((self.num_bins,), dtype=np.float64)
        self.restore = np.uint32(1 if restore_neutron else 0)
        self.filename = filename

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                    mf.WRITE_ONLY,
                                    self.histo.nbytes)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e_mon.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.histo_cl,
                          self.binning,
                          self.restore)

        neutrons = np.zeros((N, ), dtype=clarr.vec.double16)
        cl.enqueue_copy(queue, neutrons, neutron_buf).wait()

        counted = np.where((neutrons['s14'] > 0) & (neutrons['s12'].astype(np.uint32) == self.idx))
        self.histo, _ = np.histogram(neutrons['s14'][counted], bins=range(self.num_bins + 1), weights=neutrons['s9'][counted])

        self.plot_histo()

    def plot_histo(self):
        plt.figure()
        axis = np.linspace(self.binning['s0'], self.binning['s2'], num=self.num_bins)
        plt.plot(axis, self.histo)

        if self.filename:
            np.save(self.filename + 'X.dat', axis)
            np.save(self.filename + 'Y.dat', self.histo)

    @property
    def binning(self):
        return self._binning

    @binning.setter
    def binning(self, val):
        self._binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.double3)