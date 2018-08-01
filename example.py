import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl
from mcramp import Instrument

import matplotlib

from mpl_toolkits.mplot3d import Axes3D

import os

os.environ['PYOPENCL_CTX']=':'
matplotlib.rcParams.update({'font.size': 12})

if __name__ == '__main__':
    N = int(1e7)
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument.fromJSON('inst_noguide.json', ctx, queue)
    inst_van = Instrument.fromJSON('inst_vanadium.json', ctx, queue)
    inst.non_linear_sim(N, 8)
    inst_van.linear_sim(N)

    fig = plt.figure()

    #ax = fig.add_subplot(1, 1, 1, projection='3d')
 
    #inst.visualize(fig=fig, ax=ax, linewidth=1)
    #ax.set_xlim((-0.5, 0.5))
    #ax.set_ylim((-0.5, 0.5))
    #ax.set_zlim((4, 5))
 
    # Plot resulting histogram
    ax = fig.add_subplot(1, 1, 1)

    histo = inst.components['3']
    histo_van = inst_van.components['3']
    #ax.plot(np.arange(0.008,0.012,step=1e-5),
    #    histo.scat_kernel.histo, 'b-')
    #plt.show()

    Q, DE, I = histo.scat_kernel.reduce_histogram()
    _, _, I_van = histo_van.scat_kernel.reduce_histogram()

    I_van[np.where(I_van<=0)] = 1

    plt.pcolormesh(Q, DE, (I/np.max(I)), vmin=0, vmax=2)
    plt.colorbar()
    plt.show()