from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import mcpl

def MCPL_to_RAMP(particle, i):
    neutron = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., np.float32(np.random.randint(0, 2**30)+i), 0., 0., 0., 0.]
    SE2V = 437.393377

    KE = particle.ekin*1e9
    direction = particle.direction
    polarisation = particle.polarisation
    position = particle.position

    neutron[0] = position[0]*1e-2
    neutron[1] = position[1]*1e-2
    neutron[2] = position[2]*1e-2

    neutron[3] = direction[0]*np.sqrt(KE)*SE2V
    neutron[4] = direction[1]*np.sqrt(KE)*SE2V
    neutron[5] = direction[2]*np.sqrt(KE)*SE2V

    neutron[6] = polarisation[0]
    neutron[7] = polarisation[1]
    neutron[8] = polarisation[2]

    neutron[9] = particle.weight

    neutron[10] = particle.time*1e-3

    return tuple(neutron)

class SMCPLIn(SPrim):
    """
    Scattering kernel for MCPLIn component - Loads neutron buffer with neutrons from
    MCPL file.

    Parameters
    ----------
    None

    Methods
    -------
    Data
        None
    Plot
        None
    Save
        None

    """

    def __init__(self, filename="", idx=0, ctx=0, **kwargs):
        self.filename = filename
        return

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        neutrons = np.zeros((N, ), dtype=clarr.vec.float16)
        myfile = mcpl.MCPLFile(self.filename)
        particles = myfile.particles
        M = myfile.nparticles
        
        i = 0
        for p in particles:
            if i > N and i > M:
                break
            neutrons[i] = MCPL_to_RAMP(p, i)
            i+=1

        while i < N:
            neutrons[i] = (0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.)
            i+=1

        cl.enqueue_copy(queue, neutron_buf, neutrons)