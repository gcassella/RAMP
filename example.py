#import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl
from mcramp import Instrument

#from mpl_toolkits.mplot3d import Axes3D

import os

if __name__ == '__main__':
    N = 100000

    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument.fromJSON('inst.json', ctx, queue)
    inst.non_linear_sim(N, 4)

    #fig = plt.figure()

    #ax = fig.add_subplot(2, 1, 1, projection='3d')

    #inst.visualize(fig=fig, ax=ax, linewidth=1)
    #ax.set_xlim((-0.5, 0.5))
    #ax.set_ylim((-0.5, 0.5))
    #ax.set_zlim((4, 5))

    ### Plot resulting histogram
    #ax = fig.add_subplot(2, 1, 2)
    #ax.set_xlabel('Scattering angle / deg')
    #ax.set_ylabel('Counts')
    #histo = inst.components['4']
    #ax.plot(np.linspace(np.rad2deg(0.1), np.rad2deg(2.44), num=histo.scat_kernel.histo.shape[0]),
    #    histo.scat_kernel.histo, 'b-')
    #plt.show()