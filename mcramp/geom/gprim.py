import pyopencl.array as clarr
import numpy as np

from abc import ABCMeta, abstractmethod

class GPrim():
    """
    Primitive class from which all geometry kernels inherit.

    Attributes
    ----------
    None

    Methods
    -------
    lines()
        Returns a list of lines used for rendering the geometry of the kernel
        in matplotlib.

    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def lines(self):
        return [[0],
                [0],
                [0]]
