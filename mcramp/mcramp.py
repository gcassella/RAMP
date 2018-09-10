import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os, json, importlib

from time import time

class Component:
    def __init__(self, geom_kernel, scat_kernel):
        self.geom_kernel = geom_kernel
        self.scat_kernel = scat_kernel

class Instrument:
    def __init__(self, source, components, ctx, queue):
        self.source = source
        self.components = components
        self.ctx = ctx
        self.queue = queue

        self.dev = self.ctx.devices[0]
        self.max_buf = int(1E7)

    @staticmethod
    def fromJSON(fn, ctx, queue):
        inst = json.load(open(fn, 'r'))
        comps = {}
        sourceidx, mantididx, sampleidx = (0, 0, 0)
        i = 1

        for comp in inst.values():
            if "source" in comp:
                mk = getattr(importlib.import_module("mcramp"), comp['moderator_kernel']['name'])
                args = {k : v for (k,v,) in comp['moderator_kernel'].items() if not k == 'name'}
                args['ctx'] = ctx

                source = mk(**args)

                sourceidx = i
            else:
                gk = getattr(importlib.import_module("mcramp"), comp['geom_kernel']['name'])
                gargs = {k : v for (k,v) in comp['geom_kernel'].items() if not k == 'name'}
                gargs['idx'] = i
                gargs['ctx'] = ctx

                sk = getattr(importlib.import_module("mcramp"), comp['scat_kernel']['name'])
                sargs = {k : v for (k,v) in comp['scat_kernel'].items() if not k ==  'name'}
                sargs['idx'] = i
                sargs['ctx'] = ctx

                comps[str(i)] = Component(gk(**gargs), sk(**sargs))

                if comp['scat_kernel']['name'] == "MantidDetector":
                    mantididx = i

                if comp['scat_kernel']['name'] == "SPowder" or comp['scat_kernel']['name'] == "SIsotropic":
                    sampleidx = i

            i += 1

        inst = Instrument(source, comps, ctx, queue)
        inst.sourceidx = sourceidx
        inst.sampleidx = sampleidx
        inst.mantididx = mantididx

        return inst

    def _initialize_buffers(self, N):
        self.neutrons           = np.zeros((N, ), dtype=clarr.vec.float16)
        self.intersections      = np.zeros((N, ), dtype=clarr.vec.float8)
        self.iidx               = np.zeros((N, ), dtype=np.uint32)

        mf                      = cl.mem_flags

        self.neutrons_cl        = cl.Buffer(self.ctx, 
                                        mf.READ_WRITE | mf.COPY_HOST_PTR, 
                                        hostbuf=self.neutrons)
        self.intersections_cl   = cl.Buffer(self.ctx,
                                        mf.READ_WRITE, 
                                        self.neutrons.nbytes)
        self.iidx_cl            = cl.Buffer(self.ctx,
                                        mf.WRITE_ONLY,
                                        self.iidx.nbytes)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scat/terminator.cl'), mode='r') as f:
            self.term_prg = cl.Program(self.ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))


    def linear_sim(self, N):
        left = N

        rtime = time()

        if N > self.max_buf:
            self._initialize_buffers(self.max_buf)
        else:
            self._initialize_buffers(N)

        torun = 0

        while True:
            if left <= 0:
                break
            elif left >= self.max_buf:
                torun = self.max_buf
            else:
                torun = left

            self.source.gen_prg(self.queue,
                                torun,
                                self.neutrons_cl,
                                self.intersections_cl)

            for (idx, comp) in self.components.items():
                comp.geom_kernel.intersect_prg(self.queue, 
                                               torun, 
                                               self.neutrons_cl, 
                                               self.intersections_cl, 
                                               self.iidx_cl)

                comp.scat_kernel.scatter_prg(self.queue, 
                                               torun, 
                                               self.neutrons_cl, 
                                               self.intersections_cl, 
                                               self.iidx_cl)

            self.queue.finish()

            left -= self.max_buf

        print('Ray tracing completed in {} seconds'.format(time() - rtime))
        return time() - rtime

    def non_linear_sim(self, N, max_events):
        left = N

        if N > self.max_buf:
            self._initialize_buffers(self.max_buf)
        else:
            self._initialize_buffers(N)

        torun = 0

        while True:
            if left <= 0:
                break
            elif left >= self.max_buf:
                torun = self.max_buf
            else:
                torun = left

            self._initialize_buffers(torun)

            rtime = time()

            self.source.gen_prg(self.queue,
                                torun,
                                self.neutrons_cl,
                                self.intersections_cl)

            events = 0

            while events < max_events:
                for (idx, comp) in self.components.items():
                    comp.geom_kernel.intersect_prg(self.queue, 
                                               torun, 
                                               self.neutrons_cl, 
                                               self.intersections_cl, 
                                               self.iidx_cl)
                
                self.term_prg.terminate(self.queue, (torun, ), None, self.neutrons_cl, self.intersections_cl)

                for (idx, comp) in self.components.items():
                    comp.scat_kernel.scatter_prg(self.queue, 
                                               torun, 
                                               self.neutrons_cl, 
                                               self.intersections_cl, 
                                               self.iidx_cl)

                events += 1

            self.queue.finish()

            left -= self.max_buf

            rt_time = time() - rtime

        print("Raytracing took {} seconds".format(time() - rtime))

        return rt_time

    def visualize(self, fig=None, ax=None, **kwargs):
        if fig is None and ax is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
        elif fig is None and ax is not None:
            return ax
        elif fig is not None and ax is None:
            ax = fig.gca(projection='3d')

        ax.set_aspect('equal', 'box')

        for comp in self.components.values():
            lines = []
            lines += comp.scat_kernel.lines()
            lines += comp.geom_kernel.lines()

            for line in lines:
                ax.plot(*line, **kwargs)


        ax.set_axis_off()
        return ax

