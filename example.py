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
    d_spacing=3.3539
    mono_q=2 * np.pi / d_spacing

    Ei=4.0
    Li=(6.626e-34) / np.sqrt(2 * 1.674929e-27 * Ei * 1.602e-22) * 1e10
    ki=2 * np.pi / Li  # focus on lambda 4.522AA
    

    A2=2 * np.arcsin(mono_q / 2.0 / ki)
    A1=A2 / 2

    moz=40
    
    twotheta=30.0
    phi=0.0
    deltaE=0.5

    Eo=Ei + deltaE
    Lo=(6.626e-34) / np.sqrt(2 * 1.674929e-27 * Eo * 1.602e-22) * 1e10
    ko=2 * np.pi / Lo
    
    A2_ana=2 * np.arcsin(mono_q / 2.0 / ko)
    A1_ana=A2_ana / 2
    
    inst=Instrument('monotest.json', ctx, queue, A1=A1, A2=A2, moz=moz, d_spacing=d_spacing,
                    twotheta = twotheta, phi=phi, dE=deltaE, A_ana1=A1_ana, A_ana2=A2_ana)
    inst.execute(N)

    #cl.enqueue_barrier(queue)
    cl.enqueue_copy(queue, inst.neutrons, inst.neutrons_cl)
    queue.finish()

    print(inst.neutrons)

    plt.show()