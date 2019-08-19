import numpy as np
import pyopencl as cl
from mcramp import Instrument

import matplotlib.pyplot as plt
import os

os.environ["PYOPENCL_NO_CACHE"] = "0"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
os.environ["PYOPENCL_CTX"] = ":1"

if __name__ == '__main__':
    N = int(1e8)
    
    ## OpenCL setup and internals
    plat=cl.get_platforms()[0]
    devices=plat.get_devices()
    dev=devices[1]
    ctx = cl.Context(devices=[dev], 
        properties=[(cl.context_properties.PLATFORM, plat)])
    queue = cl.CommandQueue(ctx)

    ## Load and simulate instrument

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
    Ef=Ei
    Lf=h / np.sqrt(2 * mn * Ef * meV) * overAA
    kf=2 * np.pi / Lf
    
    Ana_angle=np.arcsin(mono_q / 2.0 / kf)

    # Determine nominal kf to look exactly at sample scattering point
    anangle=np.linspace(Ana_angle-0.02, Ana_angle+0.02, num=20)
    I=np.ones(anangle.shape)
    
    Lf=4 * np.pi * np.sin(anangle) / mono_q
    
    Ef=(h / np.sqrt(2 * mn * meV) / Lf * overAA)**2.0
    deltaE=Ei - Ef

    twotheta=2*np.arcsin(Li  / 2 / d_spacing)

   #for i, aa in enumerate(anangle):
   #    inst=Instrument('TAS.json', ctx, queue, Mono_angle=Mono_angle, d_spacing=d_spacing,
   #                    Ana_angle=aa, monodspacing=mono_d_spacing, twotheta=twotheta)
   #    inst.execute(N)
   #    counts = np.load("tascounts.npy")
   #    I[i] = counts
   #    plt.close('all')

   #plt.figure()
   #plt.xlabel("deltaE [meV]")
   #plt.ylabel("intensity [arb]")
   #plt.plot(deltaE, I)
   #plt.show()

    # Run once a final time to show all of the monitors off

    inst=Instrument('TAS.json', ctx, queue, Mono_angle=Mono_angle, d_spacing=d_spacing, 
                    Ana_angle=Ana_angle, monodspacing=mono_d_spacing, twotheta=twotheta)
    inst.execute(N)
    inst.plot()

    os.remove("tascounts.npy")

