from .gprim import GPrim #pylint: disable=E0401

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os

class GPlane(GPrim):
    """
    Geometry kernel for 'plane' geometry.

    Parameters
    ----------
    width : float
        The width of the plane
    height : float
        The height of the plane
    orientation : {"xy", "yz"}
        The orientation of the plane. "xy" gives a plane normal to the z axis,
        "yz" gives a plane normal to the x axis.

    Notes
    -----
    Intersection 1 :
        Point of intersection with the plane
    Intersection 2 :
        Same as Intersection 1.

    Methods
    -------
    None
    """

    def __init__(self, width=0, height=0, idx=0, orientation="xy", ctx=None):
        orientations = {"xy": 0, "yz": 1}

        self.orientation_str = orientation
        self.orientation = np.uint32(orientations[orientation])
        self.width      = np.float32(width)
        self.height     = np.float32(height)
        self.idx        = np.uint32(idx)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plane.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


    def intersect_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.intersect_plane(queue, (N, ),
                                 None,
                                 neutron_buf,
                                 intersection_buf,
                                 iidx_buf,
                                 self.idx,
                                 self.width,
                                 self.height,
                                 self.orientation)

    def lines(self):
        x_arr = [-self.width / 2.0, 
                 self.width / 2.0, 
                 self.width / 2.0,
                 -self.width / 2.0,
                 -self.width / 2.0]

        y_arr = [self.height / 2.0, 
                 self.height / 2.0, 
                 -self.height / 2.0,
                 -self.height / 2.0,
                 self.height / 2.0]

        z_arr = [0.0, 0.0, 0.0, 0.0, 0.0]

        if self.orientation_str == 'xy':
            return [x_arr, y_arr, z_arr]
        elif self.orientation_str == 'yz':
            return [z_arr, y_arr, x_arr]
