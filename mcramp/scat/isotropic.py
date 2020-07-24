from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

from scipy.integrate import simps

import os
import re

class SIsotropic(SPrim):
    """
    Scattering kernel for an isotropic scatterer which samples from a S(Q, w) distribution
    specified in a file format as set out by the
    `corresponding McStas component <http://mcstas.org/download/components/samples/Isotropic_Sqw.pure.html>`_.

    Parameters
    ----------
    fn : str
        Filename of S(Q, w) file
    
    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None
    """

    def __init__(self, fn_coh='', fn_inc='', fn_mag='', temperature=0, idx=0, ctx=None, **kwargs):
        self.temperature = np.float32(temperature)

        if fn_coh == '' and fn_inc == '':
            raise ValueError("Invalid filename combination")

        self.idx = np.uint32(idx)

        ############
        # COHERENT #
        ############

        mf = cl.mem_flags

        if not fn_coh == '':
            q,w,rho,sigma_abs,sigma_scat,pw_cdf,pq_cdf,sqw = self._LoadSQW(fn_coh)

            pq_cdf = pq_cdf.flatten()

            self.coh_args = (
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=w),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pw_cdf),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pq_cdf),
                np.uint32(len(q)),
                np.uint32(len(w)),
                np.float32(rho),
                np.float32(sigma_abs),
                np.float32(sigma_scat)
            )
        else:
            # Dummy values
            self.coh_args = (
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                np.uint32(0),
                np.uint32(0),
                np.float32(0),
                np.float32(0),
                np.float32(0)
            )

        ##############
        # INCOHERENT #
        ##############

        if not fn_inc == '':
            q,w,rho,sigma_abs,sigma_scat,pw_cdf,pq_cdf,sqw = self._LoadSQW(fn_inc)

            pq_cdf = pq_cdf.flatten()

            self.inc_args = (
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=w),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pw_cdf),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pq_cdf),
                np.uint32(len(q)),
                np.uint32(len(w)),
                np.float32(rho),
                np.float32(sigma_abs),
                np.float32(sigma_scat)
            )
        else:
            # Dummy values
            self.inc_args = (
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                np.uint32(0),
                np.uint32(0),
                np.float32(0),
                np.float32(0),
                np.float32(0)
            )

        ############
        # MAGNETIC #
        ############

        if not fn_mag == '':
            q,w,rho,sigma_abs,sigma_scat,pw_cdf,pq_cdf,sqw = self._LoadSQW(fn_mag)

            pq_cdf = pq_cdf.flatten()

            self.mag_args = (
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=w),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pw_cdf),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=pq_cdf),
                np.uint32(len(q)),
                np.uint32(len(w)),
                np.float32(rho),
                np.float32(sigma_abs),
                np.float32(sigma_scat)
            )
        else:
            # Dummy values
            self.mag_args = (
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array([0], dtype=np.float32)),
                np.uint32(0),
                np.uint32(0),
                np.float32(0),
                np.float32(0),
                np.float32(0)
            )

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'isotropic.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.isotropic_scatter(queue, (N,),
                                   None,
                                   neutron_buf,
                                   intersection_buf,
                                   iidx_buf,
                                   np.uint32(self.idx),
                                   *self.coh_args,
                                   *self.inc_args,
                                   *self.mag_args,
                                   self.temperature)

    def _LoadSQW(self, fn):
        header = {}
        flag = 0
        with open(fn, 'r') as fin:
            for i, line in enumerate(fin):
                data = re.findall(r"#\s+([A-za-z_]+)\s+([-+]?\d*\.\d+|\d+)", line)
                if not data == []:
                    header[data[0][0]] = data[0][1]

                if not line[0] == '#':
                    if flag == 0:
                        q = np.fromstring(line, dtype=np.float32, sep=' ')[:-1]
                        flag += 1
                        continue
                    elif flag == 1:
                        w = np.fromstring(line, dtype=np.float32, sep=' ')
                        flag += 1
                        continue
                    else:
                        break
            
            sqw = np.genfromtxt(
                fin.readlines()
            ).astype(np.float32) + np.finfo(float).eps.astype(np.float32)
                       
        sigma_abs = np.float32(header["sigma_abs"])
        sigma_coh = np.float32(header["sigma_coh"])
        rho = np.float32(header["V_rho"])

        sigma_scat = sigma_coh
        pw = simps(sqw.T, axis=1, x=q) / np.linalg.norm(sqw)
        pw_cdf = np.array([simps(pw[:i], x=w[:i]) for i in range(1,len(w)+1)], dtype=np.float32)
        pw_cdf = pw_cdf / np.max(pw_cdf)

        pq = np.array([sqw[:,i] / np.linalg.norm(sqw[:,i]) for i in range(0,len(w))])
        pq_cdf = np.array(
            [[simps(pq[j,:i], x=q[:i]) / simps(pq[j,:len(q)+1], x=q[:len(q)+1]) \
            for i in range(1,len(q)+1)] for j in range(0,len(w))], dtype=np.float32)

        return (q, w, rho, sigma_abs, sigma_scat, pw_cdf, pq_cdf, sqw)
