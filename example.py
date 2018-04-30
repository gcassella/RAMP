import matplotlib.pyplot as plt
import numpy as np
import pyopencl as cl
from mcramp import Instrument

from mpl_toolkits.mplot3d import Axes3D

import os

if __name__ == '__main__':
    N = 10000000

    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument
    inst = Instrument.fromJSON('inst2.json', ctx, queue)
    inst.non_linear_sim(N, 4)

    #fig = plt.figure()
#
    #ax = fig.add_subplot(2, 1, 1, projection='3d')
#
    #inst.visualize(fig=fig, ax=ax, linewidth=1)
    #ax.set_xlim((-0.5, 0.5))
    #ax.set_ylim((-0.5, 0.5))
    #ax.set_zlim((4, 5))
#
    ### Plot resulting histogram
    #ax = fig.add_subplot(2, 1, 2)
    #ax.set_xlabel('Scattering angle / deg')
    #ax.set_ylabel('Counts')
    histo = inst.components['4']
    plt.plot(np.linspace(-0.7, 2.44, num=histo.scat_kernel.histo.shape[0]), histo.scat_kernel.histo, 'b-')
    plt.show()
    #YH, TH = histo.scat_kernel.get_histo()
    #cl.enqueue_copy(queue, inst.neutrons, inst.neutrons_cl)
    #plt.hist(YH[np.where(abs(YH) > 0)], bins=200)
    #plt.show()

    #plt.hist(TH[np.where(abs(TH) > 0)], bins=2000)
    #plt.show()
    #weights = np.array([n[9] for n in inst.neutrons[:]])
    #times = np.array([n[10] for n in inst.neutrons[:]])
#
    #indxs = np.where(abs(TH) > 0)
#
    #events = np.dstack((YH[indxs],
    #           TH[indxs],
    #           weights[indxs],
    #           times[indxs]))
##
    #print(events)