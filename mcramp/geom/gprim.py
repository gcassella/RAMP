import pyopencl.array as clarr
import numpy as np

from abc import ABCMeta, abstractmethod

class GPrim():
    __metaclass__ = ABCMeta

    @property
    @abstractmethod
    def position(self):
        return self._position
    
    @position.setter
    def position(self, val):
        self._position = np.array((val[0], val[1], val[2], 0.),
                                    dtype=clarr.vec.float3)

    @abstractmethod
    def lines(self):
        return []
