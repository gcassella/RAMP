import pyopencl as cl
from mcramp import Instrument

import numpy as np

import matplotlib.pyplot as plt
import os

from time import time

os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"

if __name__ == '__main__':
    N = int(1e6)
    
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument('MCPL_testback.json', ctx, queue, v_foc=1000.0, pha_offset = 222e-6)

    inst.execute(N)
    queue.finish()
    inst.plot()