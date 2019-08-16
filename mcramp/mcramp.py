import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os, json, importlib, re

from time import time

class Component:
    def __init__(self, geom_kernel, scat_kernel, restore_neutron, pos, rot=[0, 0, 0]):
        self.geom_kernel = geom_kernel
        self.scat_kernel = scat_kernel

        # Should neutron be terminated if it doesn't hit this component?
        self.restore_neutron = np.uint32(1) if restore_neutron else np.uint32(0)

        self.pos = np.array((pos[0], pos[1], pos[2], 0.), dtype=clarr.vec.double3)
        self.rot = np.array((rot[0], rot[1], rot[2], 0.), dtype=clarr.vec.double3)

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
                pos = comp['position'] if 'position' in comp else [0, 0, 0]
                rot = comp['rotation'] if 'rotation' in comp else[0, 0, 0]
                restore_neutron = comp['restore_neutron'] if 'restore_neutron' in comp else False
                
                if 'relative' in comp:
                    rcomp_name = comp['relative']
                    rcomp_pos = comps[rcomp_name].pos
                    rcomp_rot = comps[rcomp_name].rot

                    pos = frame_rotate(pos, [rcomp_rot['s0'], rcomp_rot['s1'], rcomp_rot['s2']])

                    pos = np.add(pos, [rcomp_pos['s0'], rcomp_pos['s1'], rcomp_pos['s2']])
                    rot = np.add(rot, [rcomp_rot['s0'], rcomp_rot['s1'], rcomp_rot['s2']])

                gk = getattr(importlib.import_module("mcramp"), comp['geom_kernel']['name'])
                gargs = {k: v for (k, v) in comp['geom_kernel'].items() if not k == 'name'}
                gargs['idx'] = i
                gargs['ctx'] = parent.ctx

                sk = getattr(importlib.import_module("mcramp"), comp['scat_kernel']['name'])
                sargs = {k: v for (k, v) in comp['scat_kernel'].items() if not k == 'name'}
                sargs['idx'] = i
                sargs['ctx'] = parent.ctx

                comps[name] = Component(gk(**gargs), sk(**sargs), restore_neutron, pos, rot)
            
            i += 1

        ex_block = ExecutionBlock(source, comps, parent, linear, events)

        return ex_block

    def execute(self, N):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'scat/terminator.cl'), mode='r') as f:
                self.term_prg = cl.Program(self.parent.ctx,
                                           f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'scat/frametransform.cl'), mode='r') as f:
                self.trans_prg = cl.Program(self.parent.ctx,
                                           f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))

        if self.linear:

            for (idx, comp) in self.components.items():
                self.trans_prg.transform(self.parent.queue, (N,), None,
                                         self.parent.neutrons_cl,
                                         comp.pos,
                                         comp.rot)

                comp.geom_kernel.intersect_prg(self.parent.queue,
                                               N,
                                               self.parent.neutrons_cl,
                                               self.parent.intersections_cl,
                                               self.parent.iidx_cl)

                self.term_prg.terminate(self.parent.queue, (N,), None,
                                        self.parent.neutrons_cl,
                                        self.parent.intersections_cl,
                                        getattr(comp, 'restore_neutron', np.uint32(0)))           

                comp.scat_kernel.scatter_prg(self.parent.queue,
                                             N,
                                             self.parent.neutrons_cl,
                                             self.parent.intersections_cl,
                                             self.parent.iidx_cl)

                self.trans_prg.untransform(self.parent.queue, (N,), None,
                                           self.parent.neutrons_cl,
                                           comp.pos,
                                           comp.rot)

                #cl.enqueue_copy(self.parent.queue, self.parent.neutrons, self.parent.neutrons_cl)
                #self.parent.queue.finish()

                #print(self.parent.neutrons)

            self.parent.queue.finish()

        else:
            events = 0

            while events < self.max_events:
                for (idx, comp) in self.components.items():
                    self.trans_prg.transform(self.parent.queue, (N,), None,
                                         self.parent.neutrons_cl,
                                         comp.pos,
                                         comp.rot)
                    comp.geom_kernel.intersect_prg(self.parent.queue,
                                                   N,
                                                   self.parent.neutrons_cl,
                                                   self.parent.intersections_cl,
                                                   self.parent.iidx_cl)
                    self.trans_prg.untransform(self.parent.queue, (N,), None,
                                           self.parent.neutrons_cl,
                                           comp.pos,
                                           comp.rot)

                self.term_prg.terminate(self.parent.queue, (N,), None,
                                        self.parent.neutrons_cl,
                                        self.parent.intersections_cl)

                for (idx, comp) in self.components.items():
                    self.trans_prg.transform(self.parent.queue, (N,), None,
                                         self.parent.neutrons_cl,
                                         comp.pos,
                                         comp.rot)
                    comp.scat_kernel.scatter_prg(self.parent.queue,
                                                 N,
                                                 self.parent.neutrons_cl,
                                                 self.parent.intersections_cl,
                                                 self.parent.iidx_cl)
                    self.trans_prg.untransform(self.parent.queue, (N,), None,
                                           self.parent.neutrons_cl,
                                           comp.pos,
                                           comp.rot)

                events += 1

class Instrument:
    def __init__(self, fn, ctx, queue, **kwargs):
        self.ctx = ctx
        self.queue = queue
        self.blocks = self.fromJSON(fn, ctx, queue, **kwargs)

        self.dev = self.ctx.devices[0]
        self.max_buf = int(1E7)

    def substitute_params(self, json_str, **kwargs):
        # FIXME: if a token contains another as a substring things get effed up
        pattern = r"(?:\$)(.*?)(?:\$)"
        tokens = re.findall(pattern, json_str)

        for t in tokens:
            subbed_t = t
            for key, value in kwargs.items():
                subbed_t = subbed_t.replace("{}".format(key), str(value))

            json_str = json_str.replace("${}$".format(t), str(eval(subbed_t)))

        return json_str

    def fromJSON(self, fn, ctx, queue, **kwargs):
        with open(fn, 'r') as f:
            json_str = f.read()
            json_str = self.substitute_params(json_str, **kwargs)

        inst = json.loads(json_str)

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
        self.neutrons = np.zeros((N, ), dtype=clarr.vec.double16)
        self.intersections = np.zeros((N, ), dtype=clarr.vec.double8)
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
        device = self.queue.get_info(cl.command_queue_info.DEVICE)
        max_mem_alloc_size = device.get_info(cl.device_info.MAX_MEM_ALLOC_SIZE)
        buf_max = int(max_mem_alloc_size / 8 / 26)

        i = 0
        for block in self.blocks:
            if block.source is not None:
                source_idx = i
                break

            i += 1
        
        if (N > buf_max):
            remaining = N
            while remaining > 0.0:
                torun = buf_max if remaining > buf_max else remaining
                remaining = remaining - torun
                self._initialize_buffers(torun)

                self.blocks[source_idx].source.gen_prg(self.queue,
                                        torun,
                                        self.neutrons_cl,
                                        self.intersections_cl)

                for block in self.blocks:
                    block.execute(torun)

def frame_rotate(vec, rot):
    x_a = rot[0]
    y_a = rot[1]
    z_a = rot[2]

    r00 = np.cos(y_a)*np.cos(z_a)
    r01 = -np.cos(x_a)*np.sin(z_a) + np.sin(x_a)*np.sin(y_a)*np.cos(z_a)
    r02 = np.sin(x_a)*np.sin(z_a) + np.cos(x_a)*np.sin(y_a)*np.cos(z_a)
    r10 = np.cos(y_a)*np.sin(z_a)
    r11 = np.cos(x_a)*np.cos(z_a) + np.sin(x_a)*np.sin(y_a)*np.sin(z_a)
    r12 = -np.sin(x_a)*np.cos(z_a) + np.cos(x_a)*np.sin(y_a)*np.sin(z_a)
    r20 = -np.sin(y_a)
    r21 = np.sin(x_a)*np.cos(y_a)
    r22 = np.cos(x_a) * np.cos(y_a)

    row1 = [r00, r01, r02]
    row2 = [r10, r11, r12]
    row3 = [r20, r21, r22]
    
    return [np.dot(row1, vec), np.dot(row2, vec), np.dot(row3, vec)]