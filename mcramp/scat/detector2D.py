from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os
import datetime

class SDetector2D(SPrim):
    """
    Scattering kernel for a two axis monitor supporting a variety of axis variables.

    Parameters
    ----------
    axis1_var, axis2_var : { "x", "y", "theta", "alpha", "tof", "divX", "divY" }
        Chooses the variables of each axis. These correspond to:
            - "x" : x coordinate of intersection point with the detector
            - "y" : y coordinate of intersection point with the detector
            - "theta" : Angle made by intersection point in xz plane with z axis, \
                typically the in-plane scattering angle
            - "alpha" : Angle made by intersection point in yz plane with z axis, \
                typically the out-of-plane scattering angle
            - "divX" : Horizontal divergence of neutron velocity
            - "divY" : Vertical divergence of neutron velocity
    axis1_binning : 3-tuple of floats
        Lower bin edge, bin size, and upper bin edge for axis 1
    axis2_binning : 3-tuple of floats
        Lower bin edge, bin size, and upper bin edge for axis 2
    restore_neutron : Boolean
        If False, neutron is terminated upon intersection with this component
    filename : str or None
        Name of the file to which histogram will be saved. No file saved if
        filename is None
    logscale : Boolean
        If True, histogram intensity is plotted on a logarithmic scale

    Methods
    -------
    Data
        Returns a 3-tuple of numpy arrays, the first two containing the binning\
        for axes 1 and 2, respectively, and the third containing the histogrammed\
        neutron weights for each bin
    Plot
        Displays a plot of the two dimensional histogram of neutron weights along\
        the chosen axes.
    Save
        Saves the binning along axes 1 and 2 and the histogrammed neutron weights\
        for each bin as filename_X, filename_Y, and filename_Z, respectively, if\
        filename is not None.

    """

    def __init__(self, shape="", axis1_binning=(0, 0, 0),
                 axis2_binning=(0, 0, 0), axis1_var="theta", axis2_var="tof",
                 restore_neutron=False, idx=0, ctx=None, filename=None, logscale = False, **kwargs):
        var_dict = {
            "x" : 0,
            "y" : 1,
            "theta" : 2,
            "alpha" : 3,
            "tof" : 4,
            "divX" : 5,
            "divY" : 6
        }

        self.var_label_dict = {
            0 : "x [m]",
            1 : "y [m]",
            2 : "Theta [deg]",
            3 : "Alpha [deg]",
            4 : "Time-of-flight [us]",
            5 : "Horizontal divergence [deg]",
            6 : "Vertical divergence [deg]"
        }
        
        self.last_ran_datetime = datetime.datetime.now()
        self.last_copy_datetime = datetime.datetime.now()

        self.axis1_var = np.uint32(var_dict[axis1_var])
        self.axis2_var = np.uint32(var_dict[axis2_var])
        self.axis1_binning = axis1_binning
        self.axis2_binning = axis2_binning
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

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'detector2D.cl'), mode='r') as f:
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
                          self.axis1_var,
                          self.axis2_var,
                          self.restore_neutron)

        self.last_ran_datetime = datetime.datetime.now()

    def plot(self, queue):
        self._cached_copy(queue)

        plt.figure()
        Z = (np.log(self.Z + 1e-7) if self.logscale else self.Z)

        plt.pcolormesh(self.X, self.Y, Z, cmap='jet', shading='gouraud')
        plt.xlabel(self.var_label_dict[self.axis1_var])
        plt.ylabel(self.var_label_dict[self.axis2_var])
        plt.colorbar()
        plt.tight_layout()

    def save(self, queue):
        self._cached_copy(queue)

        if self.filename:
            np.save(self.filename + 'X.dat', self.X)
            np.save(self.filename + 'Y.dat', self.Y)
            np.save(self.filename + 'Z.dat', self.Z)

    def data(self, queue):
        self._cached_copy(queue)
        return (self.X, self.Y, self.Z)

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

    def _cached_copy(self, queue):        
        if self.last_ran_datetime > self.last_copy_datetime:
            cl.enqueue_copy(queue, self.histo, self.histo_cl).wait()
            self.Z = self.histo.reshape((self.axis1_num_bins, self.axis2_num_bins)).T
        
        self.last_copy_datetime = datetime.datetime.now()
