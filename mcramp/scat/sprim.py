import numpy as np

from abc import ABCMeta, abstractmethod

class SPrim():
    __metaclass__ = ABCMeta

    @abstractmethod
    def lines(self):
        return []
