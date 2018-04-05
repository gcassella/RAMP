import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import simps

import os, csv, re

from time import time

import tkinter as tk
from PIL import Image, ImageTk

from enum import Enum

from abc import ABCMeta, abstractmethod

os.environ['PYOPENCL_CTX']='0:0'
os.environ['PYOPENCL_COMPILER_OUTPUT']='1'

class GSphere():
    def __init__(self, radius, position, idx, ctx):
        self.radius     = np.float32(radius)
        self.position   = position
        self.idx        = idx

        with open('geom/sphere.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_sphere(queue, (N,),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.position,
                                  self.radius)

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, val):
        self._position = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

class GBanana():
    def __init__(self, radius, position, height, mintheta, maxtheta, idx, ctx):
        self.radius     = np.float32(radius)
        self.height     = np.float32(height)
        self.mintheta   = np.float32(mintheta)
        self.maxtheta   = np.float32(maxtheta)
        self.position   = position
        self.idx        = idx

        with open('geom/banana.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_banana(queue, (N, ),
                                  None,
                                  neutron_buf,
                                  intersection_buf,
                                  iidx_buf,
                                  np.uint32(self.idx),
                                  self.position,
                                  self.radius,
                                  self.height,
                                  self.mintheta,
                                  self.maxtheta)

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, val):
        self._position = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

class SIsotropic():
    def __init__(self, fn, idx, ctx):
        q,w,rho,sigma_abs,sigma_scat,pw_cdf,pq_cdf,sqw = self._LoadSQW(fn)

        pq_cdf = pq_cdf.flatten()

        self.q = q
        self.w = w
        self.rho = rho
        self.sigma_abs = sigma_abs
        self.sigma_scat = sigma_scat
        self.pw_cdf = pw_cdf
        self.pq_cdf = pq_cdf.flatten()

        self.idx = idx

        self.q_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)
        self.w_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=w)
        self.pw_cdf_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pw_cdf)
        self.pq_cdf_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pq_cdf)

        with open('scat/isotropic.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf, histo_buf):
        self.prg.isotropic_scatter(queue, (N,),
                                   None,
                                   neutron_buf,
                                   intersection_buf,
                                   iidx_buf,
                                   np.uint32(self.idx),
                                   self.q_opencl,
                                   self.w_opencl,
                                   self.pw_cdf_opencl,
                                   self.pq_cdf_opencl,
                                   np.uint32(len(self.q)),
                                   np.uint32(len(self.w)),
                                   np.float32(self.rho),
                                   np.float32(self.sigma_abs),
                                   np.float32(self.sigma_scat))

    def _LoadSQW(self, fn):
        with open(fn, 'r') as fin:
            lines = fin.readlines()
            q = np.fromstring(lines[17], dtype=np.float32, sep=' ')
            w = np.fromstring(lines[20], dtype=np.float32, sep=' ')
            rho = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", lines[6])[0])
            sigma_abs = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", lines[9])[0])
            sigma_coh = np.float32(re.findall(r"[-+]?\d*\.\d+|\d+", lines[10])[0])

        sqw = np.loadtxt(fn, skiprows=23).astype(np.float32) + np.finfo(float).eps.astype(np.float32)
        sigma_scat = sigma_coh

        pw = simps(sqw.T, axis=1, x=q) / np.linalg.norm(sqw)
        pw_cdf = np.array([simps(pw[:i], x=w[:i]) for i in range(1,len(w)+1)], dtype=np.float32)
        pw_cdf = pw_cdf / np.max(pw_cdf)

        pq = np.array([sqw[:,i] / np.linalg.norm(sqw[:,i]) for i in range(0,len(w))])
        pq_cdf = np.array(
            [[simps(pq[j,:i], x=q[:i]) / simps(pq[j,:len(q)+1], x=q[:len(q)+1]) \
            for i in range(1,len(q)+1)] for j in range(0,len(w))], dtype=np.float32)

        return (q, w, rho, sigma_abs, sigma_scat, pw_cdf, pq_cdf, sqw)

class SPowder():
    def __init__(self, fn, idx, ctx):
        reflections = self._LoadLAZ(fn)

        self.nreflections = len(reflections)
        self.idx = idx
        mf = cl.mem_flags

        self.reflections_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=reflections)

        with open('scat/powder.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf, histo_buf):
        self.prg.powder_scatter(queue, (N,),
                                   None,
                                   neutron_buf,
                                   intersection_buf,
                                   iidx_buf,
                                   np.uint32(self.idx),
                                   self.reflections_opencl,
                                   np.uint32(self.nreflections))

    def _LoadLAZ(self, fn):
        reflections_temp = np.loadtxt(fn)

        reflections = np.empty((0,), dtype=clarr.vec.float3)

        sum_intensity = sum(reflections_temp[:,11])

        reflections = np.array([(ref[4], ref[5], ref[11]/sum_intensity, 0.) for ref in reflections_temp],
                               dtype=clarr.vec.float3 )

        return reflections

class SRandom():
    def __init__(self, idx, ctx):
        self.idx = idx

        with open('scat/random.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf, histo_buf):
        self.prg.random_scatter(queue, (N,),
                                None,
                                neutron_buf,
                                intersection_buf,
                                iidx_buf,
                                np.uint32(self.idx))

class MGaussian():
    def __init__(self, ctx, pos, norm, rad, t_dim, t_pos, lamb, dlamb):
        self.pos    = np.array((pos[0], pos[1], pos[2], 0.0), dtype=clarr.vec.float3)
        self.norm   = np.array((norm[0], norm[1], norm[2], 0.0), dtype=clarr.vec.float3)
        self.rad    = np.float32(rad)
        self.t_dim  = np.array((t_dim[0], t_dim[1]), dtype=clarr.vec.float2)
        self.t_pos  = np.array((t_pos[0], t_pos[1], t_pos[2], 0.0), dtype=clarr.vec.float3)
        self.lamb   = np.float32(lamb)
        self.dlamb  = np.float32(dlamb)

        with open('moderator/gaussian.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def gen_prg(self, queue, N, neutron_buf, intersection_buf):
        self.prg.generate_neutrons(queue, (N,), None, neutron_buf, intersection_buf,
                                   self.pos,
                                   self.norm,
                                   self.rad,
                                   self.t_dim,
                                   self.t_pos,
                                   self.lamb,
                                   self.dlamb)

class Component:
    def __init__(self, geom_kernel, scat_kernel):
        self.geom_kernel = geom_kernel
        self.scat_kernel = scat_kernel

class Detector:
    def __init__(self, position, binning, var, idx, ctx):
        self.binning     = binning
        self.var         = np.uint32(var)
        self.position    = position
        self.idx         = idx

        with open('scat/detector.cl', mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf, histo_buf):
        self.prg.detector(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          histo_buf,
                          self.position,
                          self.binning,
                          self.var)
    
    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, val):
        self._position = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

    @property
    def binning(self):
        return self._binning

    @binning.setter
    def binning(self, val):
        self._binning = np.array((val[0], val[1], val[2], 0.),
                                 dtype=clarr.vec.float3)

class Instrument:
    def __init__(self, source, components, ctx, queue):
        self.source = source
        self.components = components
        self.ctx = ctx
        self.queue = queue

        self.num_bins   = 1000

    def _initialize_buffers(self, N):
        self.neutrons           = np.zeros((N, ), dtype=clarr.vec.float16)
        self.histo        = np.zeros((self.num_bins, ), dtype=np.float32)
        self.intersections      = np.zeros((N, ), dtype=clarr.vec.float8)
        self.iidx               = np.zeros((N, ), dtype=np.uint32)

        self.ctx                = cl.create_some_context()
        self.queue              = cl.CommandQueue(ctx)
        mf                      = cl.mem_flags

        self.neutrons_cl        = cl.Buffer(ctx, 
                                        mf.READ_WRITE | mf.COPY_HOST_PTR, 
                                        hostbuf=self.neutrons)
        self.intersections_cl   = cl.Buffer(ctx,
                                        mf.READ_WRITE, 
                                        self.neutrons.nbytes)
        self.iidx_cl            = cl.Buffer(ctx,
                                        mf.WRITE_ONLY,
                                        self.iidx.nbytes)
        self.histo_cl     = cl.Buffer(ctx,
                                        mf.WRITE_ONLY,
                                        self.histo.nbytes)

        with open('scat/terminator.cl', 'r') as f:
            self.term_prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.getcwd()))

    def linear_sim(self, N):
        rtime = time()
        self._initialize_buffers(N)
        print('Buffer initialization completed in {} seconds'.format(time() - rtime))

        self.source.gen_prg(self.queue,
                            N,
                            self.neutrons_cl,
                            self.intersections_cl)

        rtime = time()

        for (idx, comp) in self.components.items():
            comp.geom_kernel.intersect_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl)
            comp.scat_kernel.scatter_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl,
                                           self.histo_cl)

        self.queue.finish()

        print('Ray tracing completed in {} seconds'.format(time() - rtime))

        cl.enqueue_copy(self.queue, self.histo, self.histo_cl)
        cl.enqueue_copy(self.queue, self.neutrons, self.neutrons_cl)

    def non_linear_sim(self, N, max_events):
        self._initialize_buffers(N)

        self.source.gen_prg(self.queue,
                            N,
                            self.neutrons_cl,
                            self.intersections_cl)
        events = 0

        while events < max_events:
            for (idx, comp) in self.components.items():
                comp.geom_kernel.intersect_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl)

            self.term_prg.terminate(queue, (N,), None, self.neutrons_cl, 
                                    self.intersections_cl)

            for (idx, comp) in self.components.items():
                comp.scat_kernel.scatter_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl,
                                           self.histo_idx_cl)

            events += 1

        self.queue.finish()

if __name__ == '__main__':
    N = 1000000
    N_bins=1000

    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Create instrument
    source = MGaussian(ctx, (0, 0, 0), (0, 0, 1), 0.01, (0.001, 0.001), (0, 0, 5), 2, 0.01)
    sample2 = Component(GSphere(0.1, (0., 0, 5), 2, ctx), SPowder('data/Y3Fe5O12_YIG.laz', 2, ctx))
    detector = Component(GBanana(1.1, (0., 0., 5.), 0.1, 0, np.pi, 3, ctx),
                         Detector((0., 0., 5.), (0, np.pi/N_bins/2, np.pi/2), 0, 3, ctx))
    comps = {1: sample2, 2: detector}

    rtime = time()

    inst = Instrument(source, comps, ctx, queue)
    inst.linear_sim(N)

    print(time() - rtime)

    plt.xlabel('Scattering angle / deg')
    plt.ylabel('Counts')
    plt.plot(np.arange(len(inst.histo))*180/2000,inst.histo)
    
    plt.show()

    


