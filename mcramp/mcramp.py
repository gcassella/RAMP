import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

    def gen_idf(self):
        if self.mantididx == 0:
            return "NO MANTID DETECTOR IN INSTRUMENT - CHECK JSON"
        elif self.sourceidx == 0:
            return "NO SOURCE IN INSTRUMENT - CHECK JSON"

        import xml.etree.ElementTree as ET

        mdet = self.components[str(self.mantididx)]
        source = self.source
        sample = self.components[str(self.sampleidx)].geom_kernel

        theta_binning = mdet.scat_kernel.theta_binning
        n_thetabins = int(np.floor((theta_binning['z'] - theta_binning['x'])/(theta_binning['y'])))
        
        y_binning = mdet.scat_kernel.y_binning
        n_ybins = int(np.floor((y_binning['z'] - y_binning['x'])/(y_binning['y'])))

        radius = mdet.geom_kernel.radius

        idf_instrument = ET.Element("instrument",
                            attrib={"xmlns": "http://www.mantidproject.org/IDF/1.0",
                                    "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
                                    "xsi:schemaLocation": "http://www.mantidproject.org/IDF/1.0 http://schema.mantidproject.org/IDF/1.0/IDFSchema.xsd",
                                    "valid-from": "1900-01-31 23:59:59",
                                    "valid-to": "2100-01-31 23:59:59"})

        idf_defaults = ET.Element("defaults")

        idf_defaults_length = ET.SubElement(idf_defaults,
                                            "length",
                                            attrib = {"unit": "meter"})

        idf_defaults_angle = ET.SubElement(idf_defaults,
                                            "angle",
                                            attrib = {"unit": "radians"})

        idf_defaults_location = ET.SubElement(idf_defaults,
                                              "location",
                                              attrib = {"r": "0.0",
                                                        "t": "0.0",
                                                        "p": "0.0",
                                                        "ang": "0.0",
                                                        "axis-x": "0.0",
                                                        "axis-y": "0.0",
                                                        "axis-z": "1.0"})

        idf_defaults_along_beam = ET.SubElement(idf_defaults,
                                                "along-beam",
                                                attrib={"axis": "z"})
        
        idf_defaults_pointing_up = ET.SubElement(idf_defaults,
                                                 "pointing-up",
                                                 attrib={"axis": "y"})
        
        idf_defaults_handedness = ET.SubElement(idf_defaults,
                                                "handedness",
                                                attrib={"val": "right"})

        idf_defaults_origin = ET.SubElement(idf_defaults,
                                            "origin",
                                            attrib={"val": "beam"})

        idf_sample = ET.Element("component",
                                attrib={"type": "sample",
                                        "name": "sample"})

        idf_sample_location = ET.SubElement(idf_sample,
                                            "location",
                                            attrib={"x": str(sample.position['x']),
                                                    "y": str(sample.position['y']),
                                                    "z": str(sample.position['z'])})

        idf_sample_type = ET.Element("type",
                                     attrib={"name": "sample", "is": "SamplePos"})

        idf_source = ET.Element("component",
                                attrib={"type": "source",
                                        "name": "source"})

        idf_source_location = ET.SubElement(idf_source,
                                            "location",
                                            attrib={"x": str(source.pos['x']),
                                                    "y": str(source.pos['y']),
                                                    "z": str(source.pos['z'])})

        idf_source_type = ET.Element("type",
                                     attrib={"name": "source", "is": "Source"})

        idf_pixel = ET.Element("type", attrib={"is": "detector",
                                               "name": "pixel"})
        
        idf_mantid_detector = ET.Element("component",
                                         attrib={"type": "mantid_detector",
                                                 "name": "mantid_detector",
                                                 "idlist": "mantid_detector_list"})

        idf_mantid_detector_pos = ET.SubElement(idf_mantid_detector,
                                                "locations",
                                                attrib={"x": str(mdet.scat_kernel.position['x']),
                                                        "y": str(mdet.scat_kernel.position['y']+y_binning['x']),
                                                        "y-end": str(mdet.scat_kernel.position['y']+y_binning['z']),
                                                        "n-elements": "1",
                                                        "z": str(mdet.scat_kernel.position['z']),
                                                        "axis-x": "0.0",
                                                        "axis-y": "1.0",
                                                        "axis-z": "0.0"})

        idf_type_mantid_detector = ET.Element("type",
                                              attrib={"name": "mantid_detector"})
        
        idf_type_mantid_detector_pixels = ET.SubElement(idf_type_mantid_detector,
                                                        "component",
                                                        attrib={"type": "pixel"})

        idf_type_mantid_detector_pixels_location = ET.SubElement(idf_type_mantid_detector_pixels,
                                                                 "locations",
                                                                 attrib={"r": str(radius),
                                                                         "t": str(theta_binning['x']),
                                                                         "t-end": str(theta_binning['z']),
                                                                         "n-elements": str(n_thetabins*n_ybins),
                                                                         "axis-x" : "0.0",
                                                                         "axis-y": "1.0",
                                                                         "axis-z": "0.0"})

        idf_idlist = ET.Element("idlist",
                                attrib={"idname": "mantid_detector_list"})
        
        idf_idlist_id = ET.SubElement(idf_idlist,
                                      "id",
                                      attrib={"start": "0",
                                              "end": str(n_thetabins*n_ybins-1)})

        idf_instrument.append(idf_source)
        idf_instrument.append(idf_source_type)
        idf_instrument.append(idf_sample)
        idf_instrument.append(idf_sample_type)
        idf_instrument.append(idf_defaults)
        idf_instrument.append(idf_pixel)
        idf_instrument.append(idf_mantid_detector)
        idf_instrument.append(idf_type_mantid_detector)
        idf_instrument.append(idf_idlist)

        tree = ET.ElementTree(element=idf_instrument)
        tree.write('mcramp.xml')

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
            self.term_prg = cl.Program(self.ctx, f.read()).build(options=r'-I "{}\include"'.format(os.path.dirname(os.path.abspath(__file__))))


    def linear_sim(self, N):
        rtime = time()
        self._initialize_buffers(N)
        print('Buffer initialization completed in {} seconds'.format(time() - rtime))

        self.source.gen_prg(self.queue,
                            N,
                            self.neutrons_cl,
                            self.intersections_cl)

        rtime = time()

        for (idx, comp) in self.components.items():
            comp.geom_kernel.intersect_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl)
                                           
            comp.scat_kernel.scatter_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl)

        self.queue.finish()

        print('Ray tracing completed in {} seconds'.format(time() - rtime))

    def non_linear_sim(self, N, max_events):
        self._initialize_buffers(N)

        rtime = time()

        self.source.gen_prg(self.queue,
                            N,
                            self.neutrons_cl,
                            self.intersections_cl)

        events = 0

        while events < max_events:
            for (idx, comp) in self.components.items():
                comp.geom_kernel.intersect_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl)
            
            cl.enqueue_copy(self.queue, self.neutrons, self.neutrons_cl)
            print(self.neutrons)

            self.term_prg.terminate(self.queue, (N, ), None, self.neutrons_cl, self.intersections_cl)

            for (idx, comp) in self.components.items():
                comp.scat_kernel.scatter_prg(self.queue, 
                                           N, 
                                           self.neutrons_cl, 
                                           self.intersections_cl, 
                                           self.iidx_cl)

            events += 1

        self.queue.finish()

        print("Raytracing took {} seconds".format(time() - rtime))

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

