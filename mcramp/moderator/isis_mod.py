import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

from random import randint

import os
import re

class MISIS():
    def __init__(self, ctx=None, spec_file="Let_Base.mcstas",
                 mod_dim=(0.0, 0.0), target_dim=(0.0, 0.0), target_dist=0.0, E_min=0.0, E_max=0.0):

        self.ctx = ctx

        self.mod_dim = np.array((mod_dim[0], mod_dim[1]), dtype=clarr.vec.float2)
        self.target_dim = np.array((target_dim[0], target_dim[1]), dtype=clarr.vec.float2)
        self.target_dist = np.float32(target_dist)
        self.E_min = np.float32(E_min)
        self.E_max = np.float32(E_max)
        self.load_spectrum(spec_file)
        self.calc_str_area(mod_dim, target_dim, target_dist)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'isis_mod.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    def calc_fraction(self, lower, upper, minimum, maximum):
        # Calculates the fraction of the bin [lower, upper] contained in the range
        # [min, max]

        bin_range = upper - lower
        frac = (upper - minimum) / bin_range if (minimum > lower) else 1.0
        frac -= (upper - maximum) / bin_range if (maximum < upper) else 0.0

        return frac

    @staticmethod
    def omega(a, b, d):
        al = a / (2*d)
        be = b / (2*d)
        return 4*np.arccos(np.sqrt((1+al**2 + be**2) / ((1+al**2)*(1+be**2))))

    def calc_str_area(self, mod_dim, target_dim, target_dist):
        n_steps = 50
        A = 0.0

        target_x = np.linspace(0, target_dim[0] / 2.0, num=n_steps)
        target_y = np.linspace(0, target_dim[1] / 2.0, num=n_steps)

        TX, TY = np.meshgrid(target_x, target_y)

        mod_x = np.linspace(-mod_dim[0] / 2.0, mod_dim[0] / 2.0, num=n_steps)
        mod_y = np.linspace(-mod_dim[1] / 2.0, mod_dim[1] / 2.0, num=n_steps)
        MX, MY = np.meshgrid(mod_x, mod_y)
        for mx in MX:
            for my in MY:
                A += np.sum(1 / (np.power(-1.0*TX + mx, 2.0) + np.power(-1.0*TY + my, 2.0) + target_dist**2.0))
        A *= (mod_dim[0] * mod_dim[1]) / (n_steps**4.0)
        A *= target_dim[0] * target_dim[1] * 10000


        self.str_area = np.float32(A)

    def load_spectrum(self, spec_file):
        with open(spec_file, 'r') as mod_file:
            lines = mod_file.read()

            # Regex parses data values in line format FLOAT FLOAT (FLOAT)
            pattern = r"(?:\d*\.\d*[e][+\-]?\d+ \d*\.\d*[e][+\-]?\d+ \(\d*\.\d*[e][+\-]?\d+\) )"
            result = re.findall(pattern, lines)

            # Regex parses data values in line format FLOAT FLOAT FLOAT
            pattern = r"(?:[+\-]?\d*\.\d* [+\-]?\d*\.\d* [+\-]?\d*\.\d*)"
            rdum = re.findall(pattern, lines)
            rdum = np.array([r.split() for r in rdum[1:]]).astype(np.float32) # Skip "points at"

            mod_file.seek(0)

            line = ''
            offset_linenum = 0
            while 'TimeOffset' not in line:
                line = mod_file.readline()
                offset_linenum += 1
            self.time_offset = np.float32(line.split()[-1])

            time_linenum = offset_linenum
            while 'time' not in line:
                line = mod_file.readline()
                time_linenum += 1

            total_linenum = time_linenum
            while 'total' not in line:
                line = mod_file.readline()
                total_linenum += 1

            mod_file.seek(0)

            e_uppers = []
            e_lowers = []

            for line in mod_file:
                if 'energy bin' in line:
                    pattern = r"(?:\d*\.\d*[e][+\-]?\d+)"
                    bounds = re.findall(pattern, line)
                    lower_bound = float(bounds[0])
                    upper_bound = float(bounds[1])

                    e_uppers.append(upper_bound * 1e9)
                    e_lowers.append(lower_bound * 1e9)

        num_ener_bins = len(e_lowers) + 1
        num_time_bins = total_linenum - time_linenum - 1

        data=np.array([r.strip().replace('(', '').replace(')', '') for r in result])
        data=np.reshape(data, (int(np.size(data) / num_time_bins), num_time_bins))

        e_min_idx=np.searchsorted(e_lowers, self.E_min) # Index of the first bin larger than self.E_min
        e_max_idx=np.searchsorted(e_lowers, self.E_max) # Index of the first bin larger than self.E_max

        e_bins=[] # Edges of the energy bins
        if self.E_min == e_lowers[0]:
            e_min_idx += 1

        e_bins.append(self.E_min)
        for e in e_lowers[e_min_idx:e_max_idx]:
            e_bins.append(e)

        e_bins.append(self.E_max)
        e_bins = np.array(e_bins)
        num_ener_bins = len(e_bins)

        fractions = [self.calc_fraction(e_lowers[i], e_uppers[i], self.E_min, self.E_max) for i in range(e_min_idx - 1, e_max_idx)]
        
        t_bins = np.array([l.split(' ')[0] for l in data[0]]).astype(float) / 1e8 # Centres of time bins

        i_vals = np.zeros(data.shape)

        for i in range(data.shape[0]):
            data_chunk = np.array([l.split(' ')[1] for l in data[i]]).astype(float)

            i_vals[i] = data_chunk

        i_vals_cropped = i_vals[e_min_idx - 1:e_max_idx,:]
        fractions = np.repeat(np.array([[f] for f in fractions]), i_vals_cropped.shape[1], axis=1)
        i_vals_cropped = np.multiply(i_vals_cropped, fractions)

        flux = np.cumsum(i_vals_cropped.flatten()).astype(np.float32)
        e_int = np.insert((flux.reshape(i_vals_cropped.shape)[:,-1]), 0, 0.0).astype(np.float32)

        total = np.sum(i_vals_cropped)

        print("Energy range: {} to {} meV\nTime range {} to {} s".format(e_bins[0],
                                                                 e_bins[-1],
                                                                 t_bins[0],
                                                                 t_bins[-1]))

        mf = cl.mem_flags

        self.num_time_bins = np.int32(num_time_bins)
        self.num_ener_bins = np.int32(num_ener_bins)

        self.flux = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=flux)
        self.time_bins = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=t_bins.astype(np.float32))
        self.ener_bins = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=e_bins.astype(np.float32))
        self.e_int = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=e_int)
        self.total = np.float32(total)

    def gen_prg(self, queue, N, neutron_buf, intersection_buf):
        self.prg.generate_neutrons(queue, (N,), None, neutron_buf, intersection_buf,
                                   self.mod_dim,
                                   self.target_dim,
                                   self.target_dist,
                                   self.E_min,
                                   self.E_max,
                                   self.num_time_bins,
                                   self.num_ener_bins,
                                   self.flux,
                                   self.time_bins,
                                   self.ener_bins,
                                   self.e_int,
                                   self.total,
                                   self.str_area,
                                   self.time_offset,
                                   np.int32(N),
                                   np.int32(randint(0,4096)))
