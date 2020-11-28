from mcramp.scat import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os
import datetime

class SDetectorPix(SPrim):
    def __init__(self, shape="", y_binning=(0, 0, 0),
                 theta_binning=(0, 0, 0), tof_binning=(0, 0, 0),
                 restore_neutron=False, idx=0, ctx=None, filename=None, logscale = False, **kwargs):

        self.last_ran_datetime = datetime.datetime.now()
        self.last_copy_datetime = datetime.datetime.now()

        self.y_binning =  np.array((*(y for y in y_binning), 0.),
                                 dtype=clarr.vec.float3)
        self.theta_binning = np.array((*(th for th in theta_binning), 0.),
                                 dtype=clarr.vec.float3)
        self.tof_binning = np.array((*(tof for tof in tof_binning), 0.),
                                 dtype=clarr.vec.float3)
        self.idx = np.uint32(idx)
        self.restore_neutron = np.uint32(1 if restore_neutron else 0)
        self.filename = filename
        self.logscale = logscale

        self.y_num_bins = np.uint32(np.ceil((y_binning[2] - y_binning[0]) / y_binning[1]))
        self.theta_num_bins = np.uint32(np.ceil((theta_binning[2] - theta_binning[0]) / theta_binning[1]))
        self.tof_num_bins = np.uint32(np.ceil((tof_binning[2] - tof_binning[0]) / tof_binning[1]))

        self.det_num_bins = self.theta_num_bins*self.y_num_bins

        self.num_bins = np.uint32(self.det_num_bins * self.tof_num_bins)

        self.histo = np.zeros((self.num_bins,), dtype=np.float32)
        self.histo_err = np.zeros((self.num_bins,), dtype=np.float32)

        mf               = cl.mem_flags
        self.histo_cl    = cl.Buffer(ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.histo)
        self.histo_err_cl = cl.Buffer(ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.histo_err)                 
        
        x = np.arange(0, self.det_num_bins, 1)
        y = np.linspace(self.tof_binning['s0'], self.tof_binning['s2'], num=self.tof_num_bins)

        self.X, self.Y = np.meshgrid(x, y)
        self.Z = np.zeros((self.tof_num_bins, self.det_num_bins))
        self.E = np.zeros((self.tof_num_bins, self.det_num_bins))

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'detectorpix.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.histo_cl,
                          self.histo_err_cl,
                          self.y_binning,
                          self.theta_binning,
                          self.tof_binning,
                          self.y_num_bins,
                          self.theta_num_bins,
                          self.tof_num_bins,
                          self.restore_neutron)

        self.last_ran_datetime = datetime.datetime.now()

    def plot(self, queue):
        self._cached_copy(queue)

        plt.figure()
        Z = (np.log(self.Z + 1e-7) if self.logscale else self.Z)

        plt.pcolormesh(self.X, self.Y, Z, cmap='jet', shading='gouraud')
        plt.xlabel("Det num")
        plt.ylabel("TOF")
        plt.colorbar()
        plt.tight_layout()

    def save(self, queue):
        self._cached_copy(queue)

        if self.filename:
            np.save(self.filename + 'X', self.X)
            np.save(self.filename + 'Y', self.Y)
            np.save(self.filename + 'Z', self.Z)
            np.save(self.filename + 'E', self.E)

    def data(self, queue):
        self._cached_copy(queue)
        return (self.X, self.Y, self.Z, self.E)

    def get_histo(self):
        return (self.X, self.Y, self.Z, self.E)

    def sum_histo(self):
        self.Z += self.histo2d.T

    def _cached_copy(self, queue):        
        if self.last_ran_datetime > self.last_copy_datetime:
            cl.enqueue_copy(queue, self.histo, self.histo_cl).wait()          
            cl.enqueue_copy(queue, self.histo_err, self.histo_err_cl).wait()
            self.Z = self.histo.reshape((self.det_num_bins, self.tof_num_bins)).T
            self.E = self.histo_err.reshape((self.det_num_bins, self.tof_num_bins)).T
        
        self.last_copy_datetime = datetime.datetime.now()
