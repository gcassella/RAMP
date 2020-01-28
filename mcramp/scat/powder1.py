from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class SPowder1(SPrim):
    """
    Scattering kernel for single Bragg scattering powder sample. Recreates the
    functionality of the Powder1 component in McStas. If a neutron satisfies
    the Bragg condition, it is scattered at an angle twotheta into a random angle
    on the Debye-Scherrer cone. Phi focusing is implemented to improve simulation
    performance.

    WARNING: Neutron weights are NOT verified to be physically accurate for the
    scattering from this component, however the lineshape is correct.

    Parameters
    ----------
    d_spacing : float
        Lattice spacing corresponding to the Bragg powder line in AA
    pack : float in range [0, 1]
        Packing fraction of the sample
    vc : float
        Volume of the sample unit cell in AA^3
    sigma_abs : float
        Absorption cross section of the sample at 2200 m/s in barns
    multiplicity : int
        Multiplicity of the powder line
    DW : float in range [0, 1]
        Debye-Waller factor
    F2 : float
        Structure factor of the powder line
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

    def __init__(self, d_spacing=0.0, pack=0.0, vc=0.0, sigma_abs=0.0,
                 multiplicity=0, DW=0.0, F2=0.0, d_phi=180.0, idx=0, ctx=None,
                 **kwargs):
        # F2 in barns         

        # Scattering cross section will be divided by v^2 in Kernel
        # Absorption cross section will be divided by v in Kernel
        q = (2 * np.pi) / d_spacing

        sigma_scat_numer = 4 * np.pi ** 3.0 * pack * multiplicity * F2 * 100 * DW
        sigma_scat_denom = vc ** 2.0 * q

        self.sigma_scat_v2 = np.float32(sigma_scat_numer / sigma_scat_denom)
        self.sigma_abs_v = np.float32(pack * sigma_abs / vc * 100)
        self.q = np.float32(q)
        self.idx = np.uint32(idx)
        self.d_phi = np.float32(d_phi)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'powder1.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.powder1(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          self.idx,
                          self.sigma_scat_v2,
                          self.sigma_abs_v,
                          self.q,
                          self.d_phi)