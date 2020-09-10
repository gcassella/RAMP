import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os, json, importlib, re, sys

from time import time

DETECTOR_KERNELS = ["PSD2d", "EMon"]

def build_kernel(filename, ctx):
    with open(filename, mode='r') as f:
            prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(__file__)))

    return prg

# FIXME: KernelRef and Component should really just be dictionaries (basically what they are already)

class KernelRef:
    # Keeps track of the execution block and component that a kernel
    # belongs to in order to retrieve it's histogram at the end of execution
    #
    # Also tracks position and rotation for visualisation purposes
    def __init__(self, block, comp, pos, rot, vis):
        self.block = block
        self.comp = comp
        self.pos = pos
        self.rot = rot
        self.vis = vis
        
        self.comp_name = self.comp["name"]

class ExecutionBlock:
    def __init__(self, source, components, parent, linear, max_events):
        self.source = source
        self.max_events = max_events
        self.linear = linear
        self.components = components
        self.parent = parent

        self.trace_lines = []

    @staticmethod
    def fromJSON(block, parent, linear, events, idx):
        comps = {}
        i = 1

        source = None

        for (name, comp) in block.items():
            if name == "linear" or name == "multi":
                continue

            if "disable" in comp and comp["disable"]:
                continue
            
            if "source" in comp:
                mk = getattr(importlib.import_module("mcramp"), comp['moderator_kernel']['name'])
                args = {k : v for (k,v,) in comp['moderator_kernel'].items() if not k == 'name'}
                args['ctx'] = parent.ctx

                source = mk(**args)
            else:
                pos = comp['position'] if 'position' in comp else [0, 0, 0]
                rot = comp['rotation'] if 'rotation' in comp else [0, 0, 0]
                vis = comp['visualise'] if 'visualise' in comp else True
                restore_neutron = comp['restore_neutron'] if 'restore_neutron' in comp else False
                
                if 'relative' in comp:
                    rcomp_name = comp['relative']
                    rcomp_pos = comps[rcomp_name]["pos"]
                    rcomp_rot = comps[rcomp_name]["rot"]

                    pos = frame_rotate(pos, [rcomp_rot['s0'], rcomp_rot['s1'], rcomp_rot['s2']])

                    pos = np.add(pos, [rcomp_pos['s0'], rcomp_pos['s1'], rcomp_pos['s2']])
                    rot = np.add(rot, [rcomp_rot['s0'], rcomp_rot['s1'], rcomp_rot['s2']])

                try:
                    gk = getattr(importlib.import_module("mcramp"), comp['geom_kernel']['name'])
                except AttributeError:
                    gk = getattr(importlib.import_module(comp['geom_kernel']['name']), comp['geom_kernel']['name'])
                
                gargs = {k: v for (k, v) in comp['geom_kernel'].items() if not k == 'name'}
                gargs['idx'] = i
                gargs['ctx'] = parent.ctx

                try:
                    sk = getattr(importlib.import_module("mcramp"), comp['scat_kernel']['name'])
                except AttributeError:
                    sk = getattr(importlib.import_module(comp['scat_kernel']['name']), comp['scat_kernel']['name'])
                
                sargs = {k: v for (k, v) in comp['scat_kernel'].items() if not k == 'name'}
                sargs['idx'] = i
                sargs['ctx'] = parent.ctx
                sargs['inst'] = parent
                sargs['block'] = idx

                comps[name] = {
                    "name" : name,
                    "geom_kernel": gk(**gargs),
                    "scat_kernel": sk(**sargs),
                    "restore_neutron": np.uint32(1) if restore_neutron else np.uint32(0),
                    "pos": np.array((pos[0], pos[1], pos[2], 0.), dtype=clarr.vec.float3),
                    "rot": np.array((rot[0], rot[1], rot[2], 0.), dtype=clarr.vec.float3),
                    "vis": vis
                }

                parent.kernel_refs.append(KernelRef(idx, comps[name], pos, rot, vis))
            
            i += 1

        ex_block = ExecutionBlock(source, comps, parent, linear, events)

        return ex_block

    def execute(self, N, debug=False, trace=False):
        if self.linear:

            if trace:
                self.prev_locs = np.zeros(self.parent.neutrons.shape, dtype=np.dtype((np.float32, (3,))))

            for (_, comp) in self.components.items():
                self.parent.trans_prg.transform(self.parent.queue, (N,), None,
                                         self.parent.neutrons_cl,
                                         comp["pos"],
                                         comp["rot"])

                comp["geom_kernel"].intersect_prg(self.parent.queue,
                                               N,
                                               self.parent.neutrons_cl,
                                               self.parent.intersections_cl,
                                               self.parent.iidx_cl)

                self.parent.term_prg.terminate(self.parent.queue, (N,), None,
                                        self.parent.neutrons_cl,
                                        self.parent.intersections_cl,
                                        comp["restore_neutron"])

                if debug or (trace and comp["vis"]):
                    cl.enqueue_copy(self.parent.queue, self.parent.intersections, self.parent.intersections_cl)
                if debug:
                    print(self.parent.intersections)

                comp["scat_kernel"].scatter_prg(self.parent.queue,
                                             N,
                                             self.parent.neutrons_cl,
                                             self.parent.intersections_cl,
                                             self.parent.iidx_cl)

                self.parent.trans_prg.untransform(self.parent.queue, (N,), None,
                                           self.parent.neutrons_cl,
                                           comp["pos"],
                                           comp["rot"])
                
                if debug or (trace and comp["vis"]):
                    cl.enqueue_copy(self.parent.queue, self.parent.neutrons, self.parent.neutrons_cl)
                    self._get_trace_lines(offset=comp["pos"], rot=comp["rot"])
                if debug:
                    print(self.parent.neutrons)

        else:
            events = 0

            while events < self.max_events:
                for (idx, comp) in self.components.items():
                    self.parent.trans_prg.transform(self.parent.queue, (N,), None,
                                         self.parent.neutrons_cl,
                                         comp["pos"],
                                         comp["rot"])
                    comp["geom_kernel"].intersect_prg(self.parent.queue,
                                                   N,
                                                   self.parent.neutrons_cl,
                                                   self.parent.intersections_cl,
                                                   self.parent.iidx_cl)
                    self.parent.trans_prg.untransform(self.parent.queue, (N,), None,
                                           self.parent.neutrons_cl,
                                           comp["pos"],
                                           comp["rot"])

                self.parent.term_prg.terminate(self.parent.queue, (N,), None,
                                        self.parent.neutrons_cl,
                                        self.parent.intersections_cl)

                for (idx, comp) in self.components.items():
                    self.parent.trans_prg.transform(self.parent.queue, (N,), None,
                                         self.parent.neutrons_cl,
                                         comp["pos"],
                                         comp["rot"])
                    comp["scat_kernel"].scatter_prg(self.parent.queue,
                                                 N,
                                                 self.parent.neutrons_cl,
                                                 self.parent.intersections_cl,
                                                 self.parent.iidx_cl)
                    self.parent.trans_prg.untransform(self.parent.queue, (N,), None,
                                           self.parent.neutrons_cl,
                                           comp["pos"],
                                           comp["rot"])

                events += 1

    def _get_trace_lines(self, offset, rot):
        self.parent.queue.finish()
        non_terminated = np.transpose(np.where(self.parent.neutrons["s15"] == 0)).flatten()
        for idx in non_terminated:
            labframe_intersection = frame_rotate(
                [
                    self.parent.intersections["s0"][idx],
                    self.parent.intersections["s1"][idx],
                    self.parent.intersections["s2"][idx]
                ],
                [
                    rot["s0"],
                    rot["s1"],
                    rot["s2"]
                ]
            )

            coords = [
                [self.prev_locs[idx][0], labframe_intersection[0] + offset["s0"]],
                [self.prev_locs[idx][1], labframe_intersection[1] + offset["s1"]],
                [self.prev_locs[idx][2], labframe_intersection[2] + offset["s2"]]
            ]

            self.trace_lines.append(coords)
        
            self.prev_locs[idx] = [
                labframe_intersection[0] + offset["s0"],
                labframe_intersection[1] + offset["s1"],
                labframe_intersection[2] + offset["s2"]
            ]


class Instrument:
    """
    Core class which instantiates instrument definiton files, and provides public
    methods for executing simulations, and plotting, saving, and analyzing results

    Parameters
    ----------
    fn : str
        Instrument definition file name
    ctx : pyopencl.Context
        OpenCL context within which the simulation is to be executed
    queue : pyopencl.CommandQueue
        OpenCL command queue for enqueueing simulation kernels
    **kwargs
        Values for variables contained in the instrument definition file
    """

    def __init__(self, fn, ctx, queue, **kwargs):
        self.ctx = ctx
        self.queue = queue
        self.kernel_refs = []

        self.blocks = self._fromJSON(fn, ctx, queue, **kwargs)

        self.dev = self.ctx.devices[0]

    def data(self):
        """
        Returns a dictionary whose keys are the names of components in the instrument
        definition file, and whose values are the data arrays returned by the respective
        component
        """

        data = {}

        for d in self.kernel_refs:
            data[d.comp_name] = self.blocks[d.block].components[d.comp_name]["scat_kernel"].data(self.queue)

        return data

    def plot(self):
        """
        Calls the plotting function of each component and displays the output
        """
        
        from matplotlib.pyplot import show

        for d in self.kernel_refs:
            self.blocks[d.block].components[d.comp_name]["scat_kernel"].plot(self.queue)

        show()

    def save(self):
        """
        Saves the data array returned by each component's `data()` function to a numpy
        file "filename.npy", where filename is the value of the filename attribute in
        the instrument definition file
        """

        for d in self.kernel_refs:
            self.blocks[d.block].components[d.comp_name]["scat_kernel"].save(self.queue)

    def visualise(self, controls=True, xlim=None, ylim=None, zlim=None, focus=None, **kwargs):
        """
        Opens a plotting window containing orthogonal projections of the instrument
        geometry.
        """

        from .visualisation import Visualisation

        vis = Visualisation(self, controls=controls, xlim=xlim, ylim=ylim, zlim=zlim, focus=focus, **kwargs)
        vis.show()

    def _substitute_params(self, json_str, **kwargs):
        # FIXME: if a token contains another as a substring things get effed up
        pattern = r"(?:\$)(.*?)(?:\$)"
        tokens = re.findall(pattern, json_str)

        for t in tokens:
            subbed_t = t
            for key, value in kwargs.items():
                subbed_t = subbed_t.replace("{}".format(key), str(value))

            json_str = json_str.replace("${}$".format(t), str(eval(subbed_t)))

        # Split on all newlines

        split_str = json_str.split("\n")

        # Remove comments

        for line_num, line in enumerate(split_str):
            if '//' in line:
                comment_split = line.split("//", 1)
                split_str[line_num] = comment_split[0]

        # Merge and return

        json_str = '\n'.join(split_str)

        return json_str

    def _fromJSON(self, fn, ctx, queue, **kwargs):
        with open(fn, 'r') as f:
            json_str = f.read()
            json_str = self._substitute_params(json_str, **kwargs)


        inst = json.loads(json_str)

        blocks = []
        i = 0

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

            blocks.append(ExecutionBlock.fromJSON(block, self, linear, events, i))
            i += 1

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
                                          self.intersections.nbytes)
        self.iidx_cl = cl.Buffer(self.ctx,
                                 mf.WRITE_ONLY,
                                 self.iidx.nbytes)

    def execute(self, N, debug=False, trace=False):
        """
        Executes the instrument simulation

        Parameters
        ----------
        N : int
            Number of neutron trajectories to simulate
        """
        device = self.queue.get_info(cl.command_queue_info.DEVICE)
        max_mem_alloc_size = device.get_info(cl.device_info.MAX_MEM_ALLOC_SIZE)
        buf_max = int(max_mem_alloc_size / 8 / 32)

        i = 0
        for block in self.blocks:
            if block.source is not None:
                source_idx = i
                break

            i += 1
        
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'scat/terminator.cl'), mode='r') as f:
                self.term_prg = cl.Program(self.ctx,
                                           f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      'scat/frametransform.cl'), mode='r') as f:
                self.trans_prg = cl.Program(self.ctx,
                                           f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.abspath(__file__))))

        self._initialize_buffers(buf_max if N > buf_max else N)

        remaining = N
        while remaining > 0.0:
            torun = buf_max if remaining > buf_max else remaining
            remaining = remaining - torun

            self.blocks[source_idx].source.gen_prg(self.queue,
                                    torun,
                                    self.neutrons_cl,
                                    self.intersections_cl)

            for block in self.blocks:
                block.execute(torun, debug, trace)

            for d in self.kernel_refs:
                self.blocks[d.block].components[d.comp_name]["scat_kernel"].data_reduce(self.queue)

            self.queue.finish()

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