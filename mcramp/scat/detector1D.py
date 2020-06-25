from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os
import datetime

class SDetector1D(SPrim):
    """
    Scattering kernel for a one axis monitor supporting a variety of axis variables.

    Parameters
    ----------
    binning : 3-tuple of floats
        Lower bin edge, bin size, and upper bin edge for axis variable
    restore_neutron : Boolean
        If False, neutron is terminated upon intersection with this component
    var : { "energy", "theta", "tof" }
        Quantity which defines the axis along which neutron weights are histogrammed
    filename : str or None
        Name of the file to which spectrum will be saved. No file saved if
        filename is None

    Methods
    -------
    Data
        Returns a 2-tuple of numpy arrays, the first containing the generated \
        binning axis and the second containing the histogrammed neutron weights in each\
        bin.
    Plot
        Displays a plot of histogrammed neutron weights as a function of neutron\
        energy
    Save
        Saves the histogram axis and histogrammed neutron weights in each energy bin to\
        numpy files "filename_X.dat" and "filename_Z.dat" if filename is not None.

    """

    def __init__(self, binning=(0, 0, 0), restore_neutron=False, idx=0, var="energy",
                 ctx=None, filename=None, **kwargs):

        var_dict = { "energy" : 0, "theta" : 1, "tof" : 2 }
        self.var_label_dict = { 0 : "Energy [meV]", 1 : "Theta [deg]", 2 : "Time-of-flight [us]" }
        self.var = np.uint32(var_dict[var])
        
        self.last_ran_datetime = datetime.datetime.now()
        self.last_copy_datetime = datetime.datetime.now()
        
        self.binning     = binning
        self.idx         = idx

        self.num_bins    = np.ceil((binning[2] - binning[0])/binning[1]).astype(np.uint32)
        self.histo = np.zeros((self.num_bins,), dtype=np.float32)
        self.restore = np.uint32(1 if restore_neutron else 0)
        self.filename = filename

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.histo)

        self.axis = np.linspace(self.binning['s0'], self.binning['s2'], num=self.num_bins)
        self.Z = np.zeros(self.histo.shape)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'detector1D.cl'), mode='r') as f:
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
                          self.restore,
                          self.var)

        self.last_ran_datetime = datetime.datetime.now()

    def plot(self, queue):
        self._cached_copy(queue)

        plt.figure()
        plt.plot(self.axis, self.histo)
        plt.ylabel("Intensity")
        plt.xlabel(self.var_label_dict[self.var])
        plt.tight_layout()

    def save(self, queue):
        self._cached_copy(queue)
        if self.filename:
            np.save(self.filename + 'X', self.axis)
            np.save(self.filename + 'Z', self.histo)

    def data(self, queue):
        self._cached_copy(queue)
        return (self.axis, self.histo)

    def get_histo(self):
        return (self.axis, self.histo)

    @property
    def binning(self):
        return self._binning

    @binning.setter
    def binning(self, val):
        self._binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)
    
    def _cached_copy(self, queue):        
        if self.last_ran_datetime > self.last_copy_datetime:
            cl.enqueue_copy(queue, self.histo, self.histo_cl).wait()
        
        self.last_copy_datetime = datetime.datetime.now()