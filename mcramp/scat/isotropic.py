import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

from scipy.integrate import simps

import os
import re

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

        mf = cl.mem_flags

        self.q_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)
        self.w_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=w)
        self.pw_cdf_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pw_cdf)
        self.pq_cdf_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pq_cdf)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'isotropic.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}\include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
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
                                   np.float32(self.sigma_scat)).wait()

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
