import pyopencl.array as clarr
import numpy as np

from abc import ABCMeta, abstractmethod

class GPrim():
    __metaclass__ = ABCMeta

    @abstractmethod
    def lines(self):
        return []
