from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

from stl import mesh

def get_bounds(your_mesh):
    all_vertices = np.concatenate((your_mesh.v0, your_mesh.v1, your_mesh.v2)).T
    x_max = np.max(all_vertices[0])
    y_max = np.max(all_vertices[1])
    z_max = np.max(all_vertices[2])

    x_min = np.min(all_vertices[0])
    y_min = np.min(all_vertices[1])
    z_min = np.min(all_vertices[2])

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    return (
        np.float32(x_min),
        np.float32(x_max),
        np.float32(y_min),
        np.float32(y_max),
        np.float32(z_min),
        np.float32(z_max)
    )

class GMesh(GPrim):
    """
    Geometry kernel for 'mesh' geometry.

    Parameters
    ----------
    filename: str
        Name of STL file specifying mesh vertices
    interior: bool
        Does this kernel correspond to the 'interior' of a component? Used for
        multiple scattering in multi mode. EXPERIMENTAL

    Notes
    -----
    Intersection 1 :
        First point of intersection with the mesh, 'entering' mesh.
    Intersection 2 :
        Second point of intersection with the mesh, 'exiting' mesh.

    Works even for non-convex meshes!

    Methods
    -------
    None
    """

    def __init__(self, filename='', interior=False, idx=0, ctx=None):
        self.idx        = np.uint32(idx)

        self.mesh = mesh.Mesh.from_file(filename)

        mf = cl.mem_flags

        self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max = get_bounds(self.mesh)
        self.points_cl = cl.Buffer(
            ctx, 
            mf.READ_ONLY | mf.COPY_HOST_PTR, 
            hostbuf=self.mesh.points.flatten().astype(np.float32)
        )
        self.num_tri = num_tri = np.uint32(len(self.mesh.points.flatten())/9)

        if interior:
            self.interior = np.uint32(1)
        else:
            self.interior = np.uint32(0)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mesh.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect(queue, (N, ),
                           None,
                           neutron_buf,
                           intersection_buf,
                           iidx_buf,
                           self.idx,
                           self.points_cl,
                           self.num_tri,
                           self.x_min,
                           self.x_max,
                           self.y_min,
                           self.y_max,
                           self.z_min,
                           self.z_max,
                           self.interior)