#import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl
from mcramp import Instrument

import matplotlib.pyplot as plt
import os

from time import sleep

os.environ["PYOPENCL_CTX"] = "0:0"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"

if __name__ == '__main__':
    N = int(1e4)
    
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument('inst_modtest.json', ctx, queue)
    inst.execute(N)

    #cl.enqueue_barrier(queue)
    #cl.enqueue_copy(queue, inst.neutrons, inst.neutrons_cl)
    queue.finish()

    print(inst.blocks[0].components['2'].scat_kernel.histo)
    plt.plot(np.linspace(0.0, 40.0, num=len(inst.blocks[0].components['2'].scat_kernel.histo)),
             inst.blocks[0].components['2'].scat_kernel.histo)

    data = np.loadtxt('data/energy_data.E', skiprows=27, unpack=True)
    total_I = np.sum(data[1])
    plt.errorbar(data[0], data[1], yerr=data[2])

    plt.show()