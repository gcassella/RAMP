from .sprim import SPrim

import numpy as np
import pyopencl as cl
import pyopencl.array as clarr
import matplotlib.pyplot as plt

import os
import datetime

class SCounter(SPrim):
    """
    Scattering kernel for Counter detector. Marks a neutron as detected for usage
    alongside the 'SRescal' component, and outputs a floating point number equal
    to the sum of the weights of detected neutrons.

    Parameters
    ----------
    filename : str or None
        Name of the file to which neutron counts will be saved. No file saved if
        filename is None

    Methods
    -------
    Data
        Returns the sum of detected neutron weights as a floating point number.
    Plot
        None
    Save
        Saves the sum of detected neutron weights as a npy file if the parameter\
        filename is not None.

    """

    def __init__(self, idx=0, ctx=None,
                 filename=None, **kwargs):
                 
        self.last_ran_datetime = datetime.datetime.now()
        self.last_copy_datetime = datetime.datetime.now()

        self.idx = idx
        
        self.counts = np.zeros((2,), dtype=np.float32)
        self.filename = filename

        mf = cl.mem_flags
        self.counts_cl = cl.Buffer(ctx,
                                     mf.READ_WRITE | mf.COPY_HOST_PTR,
                                     hostbuf=self.counts)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'counter.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def scatter_prg(self, queue, N, neutron_buf, intersection_buf, iidx_buf):
        self.prg.counter(queue, (N, ),
                          None,
                          neutron_buf,
                          intersection_buf,
                          iidx_buf,
                          np.uint32(self.idx),
                          self.counts_cl)

        self.last_ran_datetime = datetime.datetime.now()

    def data(self, queue):
        self._cached_copy(queue)
        return self.counts[0]

    def save(self, queue):
        if self.filename:
            self._cached_copy(queue)
            np.save(self.filename, self.counts[0])

    def _cached_copy(self, queue):        
        if self.last_ran_datetime > self.last_copy_datetime:
            cl.enqueue_copy(queue, self.counts, self.counts_cl).wait()
        
        self.last_copy_datetime = datetime.datetime.now()