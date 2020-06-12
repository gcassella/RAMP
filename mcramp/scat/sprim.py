import numpy as np
import pyopencl.array as clarr

from abc import ABCMeta, abstractmethod

class SPrim():
    """
    Primitive class from which all scattering kernels inherit.

    Attributes
    ----------
    None

    Methods
    -------
    lines()
        Returns a list of lines used for rendering the geometry of the kernel
        in matplotlib.
    data_reduce(queue)
        Called every time the scattering kernel executes - used to reduce
        component data between execution of buffer chunks.
    data(queue)
        Returns the component data.
    plot(queue)
        Plots the component data.
    save(queue)
        Saves the component data.

    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def lines(self):
        return [[0],
                [0],
                [0]]

    @abstractmethod
    def data_reduce(self, queue):
        return 0

    @abstractmethod
    def data(self, queue):
        return 0

    @abstractmethod
    def plot(self, queue):
        return 0

    @abstractmethod
    def save(self, queue):
        return 0