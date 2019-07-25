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
    def __init__(self, source, components, parent, linear, max_events):
        self.source = source
        self.max_events = max_events
        self.linear = linear
        self.components = components
        self.parent = parent

    @staticmethod
    def fromJSON(block, parent, linear, events):
        comps = {}
        i = 1

        source = None

        for (name, comp) in block.items():
            if name == "linear" or name == "multi":
                continue
            if "source" in comp:
                mk = getattr(importlib.import_module("mcramp"), comp['moderator_kernel']['name'])
                args = {k : v for (k,v,) in comp['moderator_kernel'].items() if not k == 'name'}
                args['ctx'] = parent.ctx

                source = mk(**args)
            else:
                gk = getattr(importlib.import_module("mcramp"), comp['geom_kernel']['name'])
                gargs = {k: v for (k, v) in comp['geom_kernel'].items() if not k == 'name'}
                gargs['idx'] = i
                gargs['ctx'] = parent.ctx

                sk = getattr(importlib.import_module("mcramp"), comp['scat_kernel']['name'])
                sargs = {k: v for (k, v) in comp['scat_kernel'].items() if not k == 'name'}
                sargs['idx'] = i
                sargs['ctx'] = parent.ctx

                comps[name] = Component(gk(**gargs), sk(**sargs))
            
            i += 1

        ex_block = ExecutionBlock(source, comps, parent, linear, events)

        return ex_block

    def execute(self, N):
        if self.linear:

            for (idx, comp) in self.components.items():
                comp.geom_kernel.intersect_prg(self.parent.queue,
                                               N,
                                               self.parent.neutrons_cl,
                                               self.parent.intersections_cl,
                                               self.parent.iidx_cl)

                comp.scat_kernel.scatter_prg(self.parent.queue,
                                             N,
                                             self.parent.neutrons_cl,
                                             self.parent.intersections_cl,
                                             self.parent.iidx_cl)

            self.parent.queue.finish()

        else:
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'scat/terminator.cl'), mode='r') as f:
                self.term_prg = cl.Program(self.parent.ctx,
                                           f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))

            events = 0

            while events < self.max_events:
                for (idx, comp) in self.components.items():
                    comp.geom_kernel.intersect_prg(self.parent.queue,
                                                   N,
                                                   self.parent.neutrons_cl,
                                                   self.parent.intersections_cl,
                                                   self.parent.iidx_cl)

                self.term_prg.terminate(self.parent.queue, (N,), None,
                                        self.parent.neutrons_cl,
                                        self.parent.intersections_cl)

                for (idx, comp) in self.components.items():
                    comp.scat_kernel.scatter_prg(self.parent.queue,
                                                 N,
                                                 self.parent.neutrons_cl,
                                                 self.parent.intersections_cl,
                                                 self.parent.iidx_cl)

                events += 1

class Instrument:
    def __init__(self, fn, ctx, queue):
        self.ctx = ctx
        self.queue = queue
        self.blocks = self.fromJSON(fn, ctx, queue)

        self.dev = self.ctx.devices[0]
        self.max_buf = int(1E7)

    def fromJSON(self, fn, ctx, queue):
        inst = json.load(open(fn, 'r'))

        blocks = []

        for block in inst.values():

            if "linear" in block:
                linear = block["linear"]

                if not linear:
                    events = block["multi"]
                else:
                    events = 1
            else:
                linear = True
                events = 1

            blocks.append(ExecutionBlock.fromJSON(block, self, linear, events))

        return blocks

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

        # Find source, should be in 1st block but save the headache and search for it
        # then generate initial neutron buffer
        i = 0
        for block in self.blocks:
            if block.source is not None:
                source_idx = i
                break
            
            i += 1
        
        self.blocks[source_idx].source.gen_prg(self.queue,
                                N,
                                self.neutrons_cl,
                                self.intersections_cl)

        for block in self.blocks:
            block.execute(N)

        return 0

