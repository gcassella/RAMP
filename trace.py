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

def LoadSQW(fn):
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

def LoadLAZ(fn):
    reflections_temp = np.loadtxt(fn)

    reflections = np.empty((0,), dtype=clarr.vec.float3)

    sum_intensity = sum(reflections_temp[:,11])

    reflections = np.array([(ref[4], ref[5], ref[11]/sum_intensity, 0.) for ref in reflections_temp],
                           dtype=clarr.vec.float3 )

    return reflections

class Component:
    __metaclass__ = ABCMeta

    @abstractmethod
    def intersect(self):
        raise NotImplementedError

    @abstractmethod
    def scatter(self):
        raise NotImplementedError

class IsotropicSphere(Component):
    def __init__(self, pos, radius, sqw):
        q,w,rho,sigma_abs,sigma_scat,pw_cdf,pq_cdf,sqw=LoadSQW('He4_liq_coh.sqw')

    def intersect(self, neutrons, intersections):




if __name__ == '__main__':
    N = 10000000
    N_bins=250

    # Initialize neutron array with some random seed
    # neutron = (x,y,z,vx,vy,vz,sx,sy,sz,p,t,mc,_,_,_,abs)
    neutrons        = np.zeros((N, ), dtype=clarr.vec.float16)
    histo           = np.zeros((N,), dtype=np.uint32)
    histo_data      = np.zeros((N_bins+1,), dtype=np.float32)
    intersections   = np.zeros((N,), dtype=clarr.vec.float8)

    # OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    mf = cl.mem_flags

    # Initialize buffers
    neutrons_opencl      = cl.Buffer(ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=neutrons)
    neutrons_genned      = cl.Buffer(ctx, mf.READ_WRITE, neutrons.nbytes)
    intersections_opencl = cl.Buffer(ctx, mf.WRITE_ONLY, intersections.nbytes)

    histo_opencl         = cl.Buffer(ctx, mf.WRITE_ONLY, histo.nbytes)

    # Reflections buffer for powder scattering
    reflections = LoadLAZ("Y3Fe5O12_YIG.laz")
    reflections_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=reflections)

    # Sqw buffer for isotropic scattering
    q,w,rho,sigma_abs,sigma_scat,pw_cdf,pq_cdf,sqw=LoadSQW('He4_liq_coh.sqw')

    pq_cdf = pq_cdf.flatten()

    q_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)
    w_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=w)
    pw_cdf_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pw_cdf)
    pq_cdf_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pq_cdf)

    with open('trace.cl', 'r') as f:
        prg = cl.Program(ctx, f.read()).build(options=r'-I "C:\Users\wyz83832\Documents\opencl\include"')

    rtime = time()

    prg.generate_neutrons(queue, neutrons.shape, None, neutrons_opencl,
                          np.array((0., 0., 0., 0.), dtype=clarr.vec.float3),
                          np.array((0., 0., 1., 0.), dtype=clarr.vec.float3),
                          np.float32(0.01),
                          np.array((0.001, 0.001), dtype=clarr.vec.float2),
                          np.array((0., 0., 5., 0.), dtype=clarr.vec.float3),
                          np.float32(2),
                          np.float32(0.001))

    prg.reset_intersections(queue, neutrons.shape, None, intersections_opencl)


    prg.intersect_sphere(queue, neutrons.shape, None, neutrons_opencl,
                         intersections_opencl,
                         np.array((0., 0., 5., 0.), dtype=clarr.vec.float3),
                         np.float32(0.02))

    prg.isotropic_scatter(queue, neutrons.shape, None, neutrons_opencl,
                        intersections_opencl,
                        q_opencl,
                        w_opencl,
                        pw_cdf_opencl,
                        pq_cdf_opencl,
                        np.uint32(len(q)),
                        np.uint32(len(w)),
                        np.float32(rho),
                        sigma_abs,
                        np.float32(sigma_scat)
                        )
    
    prg.reset_intersections(queue, neutrons.shape, None, intersections_opencl)
    #prg.intersect_banana(queue, neutrons.shape, None, neutrons_opencl,
    #                     intersections_opencl,
    #                     np.array((0., 0., 5., 0.), dtype=clarr.vec.float3),
    #                     np.float32(1.1),
    #                     np.float32(.1),
    #                     np.float32(-np.pi),
    #                     np.float32(np.pi))
    prg.intersect_sphere(queue, neutrons.shape, None, neutrons_opencl,
                         intersections_opencl,
                         np.array((0., 0., 5., 0.), dtype=clarr.vec.float3),
                         np.float32(1.1))
    #prg.detector(queue, neutrons.shape, None, neutrons_opencl,
    #                      intersections_opencl,
    #                      histo_opencl,
    #                      np.array((0., 0., 5., 0.), dtype=clarr.vec.float3),
    #                      np.array((0,np.pi/N_bins,np.pi, 0.0), dtype=clarr.vec.float3),
    #                      np.uint32(0))
    prg.detector(queue, neutrons.shape, None, neutrons_opencl,
                          intersections_opencl,
                          histo_opencl,
                          np.array((0., 0., 5., 0.), dtype=clarr.vec.float3),
                          np.array((0,5/N_bins,5, 0.0), dtype=clarr.vec.float3),
                          np.uint32(3))

    queue.finish()

    cl.enqueue_copy(queue, neutrons, neutrons_opencl)
    cl.enqueue_copy(queue, histo, histo_opencl)

    print(time() - rtime)
    
    for (h,n) in zip(histo,neutrons):
        histo_data[h] += 1

    #plt.ylabel('Intensity')
    #plt.errorbar(np.arange(len(histo_data)),
    #             histo_data,
    #             yerr=(1/np.sqrt(N))*histo_data)
    #plt.show()

    H,xe,ye= np.histogram2d([n[14] for n in neutrons],
                            [n[15] for n in neutrons],
                            bins=[len(q),len(w)],
                            range=[[0.01,np.max(q)], [0.01,np.max(w)]])
    H=H.T
    fig, ax = plt.subplots()
    X,Y = np.meshgrid(xe, ye)
    ax.pcolormesh(X, Y, H)

    #H,_ = np.histogram([n[14] for n in neutrons], bins=len(q),
    #                    range=[0.01,5])

    #print(sum(H)/N)
    #plt.plot(q,H)
    #plt.plot(np.arange(len(histo_data)),histo_data)
    
    plt.show()

    


