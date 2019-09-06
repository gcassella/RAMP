import numpy as np
import pyopencl as cl
from pyopencl.algorithm import copy_if
from mcramp import Instrument, frame_rotate

import matplotlib.pyplot as plt
import os

os.environ["PYOPENCL_NO_CACHE"] = "0"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"

if __name__ == '__main__':
    N = int(1e7)
    
    ########### OpenCL setup and internals ###########
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    ########### Instrument parameters ###########

    # Constants
    h=6.62607015e-34
    mn=1.674929e-27
    meV=1.602176634e-22
    overAA = 1e10

    # Sample parameters
    d_spacing = 6.0

    # Monochromator parameters, set for choosing Ki
    mono_d_spacing=3.3539
    mono_q=2 * np.pi / mono_d_spacing

    # Determine the Ei that will give us the desired point in reciprocal space at
    # a chosen value of twotheta
    Ei=5.0

    Li=h / np.sqrt(2 * mn * Ei * meV) * overAA
    ki=2 * np.pi / Li 
    Mono_angle=np.arcsin(mono_q / 2.0 / ki)

    # Analyzer parameters, set for choosing |Kf| to look at
    Ef=5.0
    Lf=h / np.sqrt(2 * mn * Ef * meV) * overAA
    kf=2 * np.pi / Lf
    
    Ana_angle=np.arcsin(mono_q / 2.0 / kf)
    
    sample_target_rot=1.3  # radians
    sample_target = frame_rotate([0.0, 0.0, 1.0], [0.0, sample_target_rot, 0.0])

    ########### Run instrument ###########

    inst=Instrument('TAS.json', ctx, queue, Mono_angle=Mono_angle, d_spacing=d_spacing, 
                    Ana_angle=Ana_angle, monodspacing=mono_d_spacing, sample_target=sample_target, sample_trot=sample_target_rot)
    inst.execute(N)

    ########### Get and plot results ###########

    inst.plot()

    data=inst.data()
    qx,qy,qz,dE,p=data['sample']

    fig=plt.figure()
    ax=fig.subplots(2, 3)
    ax[0, 0].hist2d(qx, qy, weights=p, bins=50)
    ax[0, 0].set_xlabel("Qx")
    ax[0, 0].set_ylabel("Qy")

    ax[0, 1].hist2d(qx, qz, weights=p, bins=50)
    ax[0, 1].set_xlabel("Qx")
    ax[0, 1].set_ylabel("Qz")

    ax[0, 2].hist2d(qy, qz, weights=p, bins=50)
    ax[0, 2].set_xlabel("Qy")
    ax[0, 2].set_ylabel("Qz")

    ax[1, 0].hist2d(qx, dE, weights=p, bins=50)
    ax[1, 0].set_xlabel("Qx")
    ax[1, 0].set_ylabel("dE")

    ax[1, 1].hist2d(qy, dE, weights=p, bins=50)
    ax[1, 1].set_xlabel("Qy")
    ax[1, 1].set_ylabel("dE")

    ax[1, 2].hist2d(qz, dE, weights=p, bins=50)
    ax[1, 2].set_xlabel("Qz")
    ax[1, 2].set_ylabel("dE")

    plt.tight_layout()
    plt.show()