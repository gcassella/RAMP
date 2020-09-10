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
    fn_coh : str
        Filename of S(Q, w) file specifying coherent scattering (no polarisation effect)
    fn_inc : str
        Filename of S(Q, w) file specifying incoherent scattering (2/3rd spin flip)
    fn_mag : str
        Filename of S(Q, w) file specifying magnetic scattering (Halpern-Johnson pol rule)
    temperature : float
        Temperature of scattering sample for detailed balance
    transmit : int
        Flag which specifies if transmitted beam should be retained (0 = off, 1 = on)

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None
    """

    def __init__(self, fn_coh='', fn_inc='', fn_mag='', temperature=0, transmit=1, idx=0, ctx=None, **kwargs):
        self.temperature = np.float32(temperature)
        self.transmit = np.uint32(transmit)

        if fn_coh == '' and fn_inc == '' and fn_mag == '':
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
                                   self.temperature,
                                   self.transmit)

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
        pw = simps(q*sqw.T, axis=1, x=q)
        pw /= simps(pw, x=w)
        pq = np.array([q*sqw[:,i] / simps(sqw.T, axis=0, x=w) for i in range(0,len(w))])


        pw_cdf = cdf_from_pd(w, pw)
        pq_cdf = np.array([cdf_from_pd(q, pq[j,:]) for j in range(0,len(w))], dtype=np.float32)
        return (q, w, rho, sigma_abs, sigma_scat, pw_cdf, pq_cdf, sqw)

def cdf_from_pd(x, pd):
    unnormed = np.array([simps(pd[:i], x=x[:i]) for i in range(1,len(x)+1)], dtype=np.float32)
    return unnormed / unnormed[-1]