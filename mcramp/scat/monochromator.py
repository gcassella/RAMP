from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class SMonochromator(SPrim):
    def __init__(self, slab_width=0.0, slab_height=0.0, gap=0.0, n_horizontal=1,
                 n_vertical=1, mosaic_horizontal=0.0, mosaic_vertical=0.0,
                 r0=1.0, d_spacing=0.0, radius_vertical=0.0, radius_horizontal=0.0,
                 idx=0, ctx=0):
                
        mf = cl.mem_flags

        gaussX = np.array([-0.987992518020485, -0.937273392400706, -0.848206583410427,
-0.724417731360170, -0.570972172608539, -0.394151347077563,
-0.201194093997435, 0, 0.201194093997435,
0.394151347077563, 0.570972172608539, 0.724417731360170,
0.848206583410427, 0.937273392400706, 0.987992518020485]).astype(np.float64)
        gaussW = np.array([0.030753241996117, 0.070366047488108, 0.107159220467172,
0.139570677926154, 0.166269205816994, 0.186161000115562,
0.198431485327111, 0.202578241925561, 0.198431485327111,
0.186161000115562, 0.166269205816994, 0.139570677926154,
0.107159220467172, 0.070366047488108, 0.030753241996117]).astype(np.float64)
        self.gausslen = np.uint32(len(gaussX))

        self.gaussX_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gaussX)
        self.gaussW_opencl = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=gaussW)

        self.idx = np.uint32(idx)
        self.slab_width = np.float64(slab_width)
        self.slab_height = np.float64(slab_height)
        self.gap = np.float64(gap)
        self.n_horizontal = np.uint32(n_horizontal)
        self.n_vertical = np.uint32(n_vertical)
        self.mosaic_horizontal = np.float64(mosaic_horizontal)
        self.mosaic_vertical = np.float64(mosaic_vertical)
        self.r0 = np.float64(r0)
        self.d_spacing = np.float64(d_spacing)
        self.radius_vertical = np.float64(radius_vertical)
        self.radius_horizontal = np.float64(radius_horizontal)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'monochromator.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.monochromator(queue, (N, ),
                                None,
                                neutron_buf,
                                intersection_buf,
                                iidx_buf,
                                self.idx,
                                self.slab_width,
                                self.slab_height,
                                self.gap,
                                self.n_horizontal,
                                self.n_vertical,
                                self.mosaic_horizontal,
                                self.mosaic_vertical,
                                self.r0,
                                self.d_spacing,
                                self.radius_vertical,
                                self.radius_horizontal,
                                self.gaussX_opencl,
                                self.gaussW_opencl,
                                self.gausslen)
