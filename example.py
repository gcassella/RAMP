#import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl
from mcramp import Instrument

import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    N = int(1e7)
    
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument.fromJSON('inst_TOF.json', ctx, queue)
    inst.non_linear_sim(N, 8)

    fig = plt.figure()
    # Plot resulting histogram
    ax = fig.add_subplot(1, 1, 1)

    histo = inst.components['4']
    ax.plot(np.arange(0.008,0.012,step=1e-5),
        histo.scat_kernel.histo, 'b-')
    plt.show()


    plt.show()
