import pyopencl as cl
from mcramp import Instrument

import numpy as np

import matplotlib.pyplot as plt
import os

os.environ["PYOPENCL_NO_CACHE"] = "1"
os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
os.environ["PYOPENCL_BUILD_OPTIONS"] = "-D CONFIG_USE_DOUBLE=1"

if __name__ == '__main__':
    N = int(1e9)
    
    ## OpenCL setup and internals
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    # Calculate TOF bins

    V2K = lambda v: v*(1.58825361e-3)
    E2K = lambda E: np.sqrt(E*(1.0/2.072))
    K2E = lambda K: 2.072*(K**2.0)
    E2V = lambda E: np.sqrt(E)*437.393377
    Ei = 5.15
    vi = E2V(Ei)
    tof_i = 1.0*1e6 / vi
    tof_f = 3.5*1e6 / vi

    vi_lo = E2V(1.0)
    vi_hi = E2V(5.65)

    R = 3.5

    tof_f_lo = R*1e6 / vi_lo
    tof_f_hi = R*1e6 / vi_hi

    inst = Instrument(
        'LET.json', 
        ctx, 
        queue, 
        Ei=Ei,
        sample=r"'delta.sqw'",
        tof_bin = [tof_i + tof_f_hi, 100, tof_i + tof_f_lo]
    )

    alu = {
        "d" : [2.338, 2.0248, 1.4317, 1.221],
        "F2" : [21.3, 16.105, 8.327, 6.22],
        "m" : [8, 6, 12, 24]
    }

    ## Execute Vanadium

    inst.execute(N)
    queue.finish()

    inst.plot()
    
    vanadium = inst.data()['Detector']
    
    vanadium_int = np.sum(vanadium[2], axis=0)
    plt.plot(vanadium_int)
    plt.show()
    exit()

    # Execute sample

    #inst = Instrument(
    #    'LET.json', 
    #    ctx, 
    #    queue, 
    #    Ei=Ei,
    #    sample=r"'He4_liq_coh.sqw'",
    #    tof_bin = [tof_i + tof_f_hi, 5, tof_i + tof_f_lo]
    #)

    ## Normalize sample data

    #inst.execute(N)
    #queue.finish()
    #sample = inst.data()['Detector']

    ##sample_Z = np.divide(sample[2], vanadium_int)
    #sample_Z = sample[2]

    ## Convert detector number to phi

    #detnum = sample[0]
    #tof = sample[1]

    #theta_binning = [-40, 1.0, 140]
    #num_theta_bins = round((theta_binning[-1] - theta_binning[0]) / theta_binning[1])

    #y_binning = [-2.0, 0.1, 2.0]
    #num_y_bins = round((y_binning[-1] - y_binning[0]) / y_binning[1])


    #R = 3.5

    #theta = (detnum % num_theta_bins) * theta_binning[1] + theta_binning[0]
    #y = (detnum // num_theta_bins) * y_binning[1] + y_binning[0]
    #x = R*np.sin(theta)
    #z = R*np.cos(theta)

    #phi = np.sign(x) * np.degrees(np.arctan2(np.sqrt(x**2.0 + y**2.0), z))

    ## Convert (phi, tof) to (Q, dE)

    #L1 = 1.0
    #L2 = np.sqrt(3.5**2.0 + y**2.0)

    #vi = E2V(Ei)
    #ki = E2K(Ei)
    #kf = V2K(L2 / ((tof - tof_i)*1e-6))
    #Ef = K2E(kf)

    #dE = Ei - Ef
    #Q = np.sqrt(ki**2.0 + kf**2.0 - 2*ki*kf*np.cos(phi))

    #plt.hist2d(
    #    dE.flatten(), 
    #    Q.flatten(),
    #    weights=sample_Z.flatten(),
    #    bins=(100, 100)
    #)

    #plt.show()


    