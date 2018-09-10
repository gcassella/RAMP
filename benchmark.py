import numpy as np
import pyopencl as cl
from mcramp import Instrument

import os

os.environ['PYOPENCL_CTX']="0:0"

if __name__ == '__main__':
    N = int(1e7)
    runs = 30

    ## Load and simulate instrument

    with open("pow.dat", "w") as pow_output:
        pow_output.write("Steps Run Time[s]\n")
        pow_times = []
        for i in range(runs):
            for j in range(2,8):
                ctx = cl.create_some_context()
                queue = cl.CommandQueue(ctx)

                pow_inst = Instrument.fromJSON('inst_powder.json', ctx, queue)

                pow_time = pow_inst.non_linear_sim(N, 2)
                pow_times.append(pow_time)

                pow_output.write("{} {} {}\n".format(j, i, pow_time))
        pow_output.write("Device used: {}".format(ctx.get_info(cl.context_info.DEVICES)))

    with open("TOF.dat", "w") as TOF_output:
        TOF_output.write("Steps Run Time[s]\n")
        TOF_times = []
        for i in range(runs):
            for j in range(2,8):
                ctx = cl.create_some_context()
                queue = cl.CommandQueue(ctx)

                TOF_inst = Instrument.fromJSON('inst_TOF.json', ctx, queue)

                TOF_time = TOF_inst.non_linear_sim(N, 5)
                TOF_times.append(TOF_time)

                TOF_output.write("{} {} {}\n".format(j, i, TOF_time))

        TOF_output.write("Device used: {}".format(ctx.get_info(cl.context_info.DEVICES)))
