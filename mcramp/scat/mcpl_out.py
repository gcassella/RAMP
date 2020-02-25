from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

def write_hdr(file, count):
    # Write MCPL header
    file.write(b'MCPL')
    # Version - 3
    file.write(b'003')
    # Endianness, litle
    file.write(b'L')
    # Number of particles
    file.write(np.uint64(count))
    # Number of comments
    file.write(np.uint32(0))
    # Number of custom data blobs
    file.write(np.uint32(0))
    # User flags flag
    file.write(np.uint32(0))
    # Polarisation vectors
    file.write(np.uint32(1))
    # Single precision
    file.write(np.uint32(1))
    # PDG code field
    file.write(np.uint32(2112))
    # Length per particle
    file.write(np.uint32(44))
    # Universal weight - no
    file.write(np.uint32(0))
    # Source name
    file.write(np.uint32(4))
    file.write(str("RAMP").encode(encoding='ascii'))

def pack_vector(x, y, z):
    ux = 0.0
    uy = 0.0
    si = 1.0
    if x > y and x > z:
        ux = 1.0 / z
        uy = uy
        si = np.sign(x)
    elif y > x and y > z:
        ux = x
        uy = 1.0 / z
        si = np.sign(y)
    elif z > x and z > y:
        ux = ux
        uy = uy
        si = np.sign(z)

    return [ux, uy, si]

def KE_from_v(x, y, z):
    vmagn2 = x**2.0 + y**2.0 + z**2.0
    return 5.22703725e-15 * vmagn2

def write_particle(file, particle):
    # Polarisation vector
    file.write(np.float32(particle[6]))
    file.write(np.float32(particle[7]))
    file.write(np.float32(particle[8]))
    # Position vector
    file.write(np.float32(particle[0] * 1e2))
    file.write(np.float32(particle[1] * 1e2))
    file.write(np.float32(particle[2] * 1e2))
    # Packed direction vector and KE
    packed = pack_vector(np.float32(particle[3]), np.float32(particle[4]), np.float32(particle[5]))
    KE = KE_from_v(np.float32(particle[3]), np.float32(particle[4]), np.float32(particle[5]))
    file.write(np.float32(packed[0]))
    file.write(np.float32(packed[1]))
    file.write(np.float32(packed[2] * KE))
    # Time of flight
    file.write(np.float32(particle[10] * 1e3))
    # Particle weight
    file.write(np.float32(particle[9]))

class SMCPLOut(SPrim):
    """
    Scattering kernel for MCPLIn component - Dumps neutron buffer into MCPL file.

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
        cl.enqueue_copy(queue, neutrons, neutron_buf)
        queue.finish()
        count = 0
        open(self.filename, 'wb').close()
        with open(self.filename, 'rb+') as fh:
            for n in neutrons:
                if not n[15] == 1.0:
                    write_particle(fh, n)
                    count+=1

            fh.seek(0, 0)
            particle_list = fh.read()
            fh.seek(0, 0)
            write_hdr(fh, count)
            fh.write(particle_list)