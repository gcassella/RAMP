#import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl
from mcramp import Instrument

import matplotlib.pyplot as plt
import os

from time import sleep

os.environ["PYOPENCL_CTX"] = "0:1"
os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "0"

if __name__ == '__main__':
    N = int(1e7)
    
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument('LET_withchop.json', ctx, queue, v_foc=1000.0, pha_offset = 222e-6)
    inst.execute(N)

    #cl.enqueue_barrier(queue)
    cl.enqueue_copy(queue, inst.neutrons, inst.neutrons_cl)
    queue.finish()

    plt.show()