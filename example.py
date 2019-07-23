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
    N = int(1e6)
    
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument('inst_modtest.json', ctx, queue)
    inst.execute(N)

    #cl.enqueue_barrier(queue)
    #cl.enqueue_copy(queue, inst.neutrons, inst.neutrons_cl)
    queue.finish()

    detector=inst.blocks[0].components['3'].scat_kernel
    x=np.linspace(detector.axis1_binning['s0'],
                  detector.axis1_binning['s2'],
                  num=detector.axis1_num_bins)
    y=np.linspace(detector.axis2_binning['s0'],
                    detector.axis2_binning['s2'],
                    num=detector.axis2_num_bins)

    X, Y=np.meshgrid(x, y)
    plt.imshow(detector.histo2d)
    plt.show()