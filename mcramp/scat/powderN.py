from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class SPowderN(SPrim):
    """
    Scattering kernel for Bragg scattering powder sample. Recreates the
    functionality of the PowderN component in McStas. If a neutron satisfies
    the Bragg condition, it is scattered at an angle twotheta into a random angle
    on the Debye-Scherrer cone. Phi focusing is implemented to improve simulation
    performance.

    WARNING: Neutron weights are NOT verified to be physically accurate for the
    scattering from this component, however the lineshape is correct.

    Parameters
    ----------
    d_spacing : float array
        Lattice spacings corresponding to the Bragg powder lines in AA
    pack : float in range [0, 1]
        Packing fraction of the sample
    vc : float
        Volume of the sample unit cell in AA^3
    sigma_abs : float
        Absorption cross section of the sample at 2200 m/s in barns
    multiplicity : int array
        Multiplicity of the powder lines
    DW : float in range [0, 1]
        Debye-Waller factor
    F2 : float array
        Structure factors of the powder lines
    d_phi : float in range [0.0, 180.0]
        Max angle around Debye-Scherrer cone into which neutrons are scattered

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, d_spacing=[0.0], pack=0.0, vc=0.0, sigma_abs=0.0,
                 multiplicity=[0], DW=0.0, F2=[0.0], d_phi=180.0, transmit=1,
                 idx=0, ctx=None, **kwargs):
        # F2 in barns         

        mf = cl.mem_flags

        # Scattering cross section will be divided by v^2 in Kernel
        # Absorption cross section will be divided by v in Kernel
        q = np.array((2 * np.pi) / np.array(d_spacing), dtype=np.float32)

        sigma_scat_numer = np.array(4 * np.pi ** 3.0 * pack * np.array(multiplicity) * np.array(F2) * 100 * DW, dtype=np.float32)
        sigma_scat_denom = vc ** 2.0 * q

        sigma_scat_v2 = (sigma_scat_numer / sigma_scat_denom)
        self.sigma_abs_v = np.float32(pack * sigma_abs / vc * 100)
        self.idx = np.uint32(idx)
        self.d_phi = np.float32(d_phi)

        self.sigma_scat_v2_cl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=sigma_scat_v2)
        self.q_cl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=q)

        self.num_lines = np.uint32(len(d_spacing))

        self.transmit = np.uint32(transmit)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'powderN.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.powderN(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.sigma_scat_v2_cl,
                          self.sigma_abs_v,
                          self.q_cl,
                          self.d_phi,
                          self.num_lines,
                          self.transmit)