import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os, json, importlib

from time import time

class Component:
    def __init__(self, geom_kernel, scat_kernel):
        self.geom_kernel = geom_kernel
        self.scat_kernel = scat_kernel

class ExecutionBlock:
    def __init__(self, components, parent, linear, max_events):
        self.max_events = max_events
        self.linear = linear
        self.components = components
        self.parent = parent

    @staticmethod
    def fromJSON(block, parent, linear, events):
        # TODO add a check here that moderators are in a linear execution block

        comps = {}
        i = 1

        for comp in block.values():

            gk = getattr(importlib.import_module("mcramp"), comp['geom_kernel']['name'])
            gargs = {k: v for (k, v) in comp['geom_kernel'].items() if not k == 'name'}
            gargs['idx'] = i
            gargs['ctx'] = ctx

            sk = getattr(importlib.import_module("mcramp"), comp['scat_kernel']['name'])
            sargs = {k: v for (k, v) in comp['scat_kernel'].items() if not k == 'name'}
            sargs['idx'] = i
            sargs['ctx'] = ctx

            comps[str(i)] = Component(gk(**gargs), sk(**sargs))

        ex_block = ExecutionBlock(comps, parent, linear, events)

        return ex_block

    def execute(self, N):
        if self.linear:

            for (idx, comp) in self.components.items():
                comp.geom_kernel.intersect_prg(parent.queue,
                                               torun,
                                               parent.neutrons_cl,
                                               parent.intersections_cl,
                                               parent.iidx_cl)

                comp.scat_kernel.scatter_prg(parent.queue,
                                             torun,
                                             parent.neutrons_cl,
                                             parent.intersections_cl,
                                             parent.iidx_cl)

            parent.queue.finish()

        else:
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'scat/terminator.cl'), mode='r') as f:
                self.term_prg = cl.Program(parent.ctx,
                                           f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))

            events = 0

            while events < self.max_events:
                for (idx, comp) in self.components.items():
                    comp.geom_kernel.intersect_prg(parent.queue,
                                                   torun,
                                                   parent.neutrons_cl,
                                                   parent.intersections_cl,
                                                   parent.iidx_cl)

                self.term_prg.terminate(parent.queue, (torun, ), None, parent.neutrons_cl, parent.intersections_cl)

                for (idx, comp) in self.components.items():
                    comp.scat_kernel.scatter_prg(parent.queue,
                                                 torun,
                                                 parent.neutrons_cl,
                                                 parent.intersections_cl,
                                                 parent.iidx_cl)

                events += 1

class Instrument:
    def __init__(self, blocks, ctx, queue):
        self.ctx = ctx
        self.queue = queue
        self.blocks = blocks

        self.dev = self.ctx.devices[0]
        self.max_buf = int(1E7)

    @staticmethod
    def fromJSON(fn, ctx, queue):
        inst = json.load(open(fn, 'r'))

        blocks = []
        i = 1

        for block in inst.values():

            if "multi" in block:
                linear = False
                events = block["multi"]
            else:
                linear = True
                events = 1

            blocks.append(ExecutionBlock.fromJSON(block, self, linear, events))

        return Instrument(blocks, ctx, queue)

    def _initialize_buffers(self, N):
        self.neutrons = np.zeros((N, ), dtype=clarr.vec.float16)
        self.intersections = np.zeros((N, ), dtype=clarr.vec.float8)
        self.iidx = np.zeros((N, ), dtype=np.uint32)

        mf = cl.mem_flags

        self.neutrons_cl = cl.Buffer(self.ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.neutrons)
        self.intersections_cl = cl.Buffer(self.ctx,
                                          mf.READ_WRITE,
                                          self.neutrons.nbytes)
        self.iidx_cl = cl.Buffer(self.ctx,
                                 mf.WRITE_ONLY,
                                 self.iidx.nbytes)

    def execute(self, N):
        # TODO Add buffer chunking to allow neutrons in excess of max buf size

        self._initialize_buffers(N)

        for block in blocks:
            block.execute(N)

        return 0

