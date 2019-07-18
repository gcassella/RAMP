import numpy as np
import pyopencl as cl
import pyopencl.array as clarr

import os
import re

class MISIS():
    def __init__(self, ctx=None, pos=(0.0, 0.0, 0.0), spec_file="Let_Base.mcstas",
                 mod_dim=(0.0, 0.0), target_dim=(0.0, 0.0), target_dist=0.0, E_min=0.0, E_max=0.0):

        self.ctx = ctx
                 
        self.pos    = np.array((pos[0], pos[1], pos[2], 0.0), dtype=clarr.vec.float3)
        self.mod_dim = np.array((mod_dim[0], mod_dim[1]), dtype=clarr.vec.float2)
        self.target_dim = np.array((target_dim[0], target_dim[1]), dtype=clarr.vec.float2)
        self.target_dist = np.float32(target_dist)
        self.E_min = np.float32(E_min)
        self.E_max = np.float32(E_max)
        self.load_spectrum(spec_file)
        self.calc_str_area(mod_dim, target_dim, target_dist)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'isis_mod.cl'), mode='r') as f:
            self.prg = cl.Program(ctx, f.read()).build(options=r'-I "{}/include"'.format(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

    @staticmethod
    def omega(a, b, d):
        al = a / (2*d)
        be = b / (2*d)
        return 4*np.arccos(np.sqrt((1+al**2 + be**2) / ((1+al**2)*(1+be**2))))

    def calc_str_area(self, mod_dim, target_dim, target_dist):
        n_steps = 50
        A = 0.0

        target_x = np.linspace(0, target_dim[0] / 2.0, num = 50)
        target_y = np.linspace(0, target_dim[1] / 2.0, num = 50)
        
        TX, TY = np.meshgrid(target_x, target_y)
        
        mod_x = np.linspace(-mod_dim[0] / 2.0, mod_dim[0] / 2.0, num = 50)
        mod_y = np.linspace(-mod_dim[1] / 2.0, mod_dim[1] / 2.0, num = 50)
        MX, MY = np.meshgrid(mod_x, mod_y)
        for mx in MX:
            for my in MY:
                A += np.sum(1 / (np.power(-1.0*TX + mx, 2.0) + np.power(-1.0*TY + my, 2.0) + target_dist**2.0))
        A *= (mod_dim[0] * mod_dim[1]) / (50**4.0)
        A *= target_dim[0] * target_dim[1] * 10000

        
        self.str_area = np.float32(A)
        
        print(self.str_area)

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
    
            num_time_bins = total_linenum - time_linenum - 1
    
            mod_file.seek(0)
    
            e_bins = []
            e_axis = []
    
            for line in mod_file:
                if 'energy bin' in line:
                    pattern = r"(?:\d*\.\d*[e][+\-]?\d+)"
                    bounds = re.findall(pattern, line)
                    lower_bound = float(bounds[0])
                    upper_bound = float(bounds[1])
    
                    bin_centre = (upper_bound + lower_bound) / 2.0
    
                    e_axis.append(bin_centre)
                    e_bins.append(lower_bound * 1e9)
    
            e_bins.append(upper_bound)

        data = np.array([r.strip().replace('(', '').replace(')', '') for r in result])
        data = np.reshape(data, (int(np.size(data) / num_time_bins), num_time_bins))

        t_axis = np.array([l.split(' ')[0] for l in data[0]]).astype(float) / 1e8
        e_axis = np.array(e_axis) * 1e9

        t_vals, e_vals = np.meshgrid(t_axis, e_axis)

        i_vals = np.zeros(e_vals.shape)

        for i in range(e_vals.shape[0]):
            data_chunk = np.array([l.split(' ')[1] for l in data[i]]).astype(float)

            i_vals[i] = data_chunk

        time_bins = np.array([l.split(' ')[0] for l in data[1]]).astype(float) / 1e8

        e_min = self.E_min  # [meV]
        e_max = self.E_max  # [meV]

        t_min = 0.0  # [s]
        t_max = 0.5  # [s]

        e_min_idx = (np.abs(e_axis - e_min)).argmin()
        e_max_idx = (np.abs(e_axis - e_max)).argmin()

        e_min = e_axis[e_min_idx]
        e_max = e_axis[e_max_idx]

        t_min_idx = (np.abs(t_axis - t_min)).argmin()
        t_max_idx = (np.abs(t_axis - t_max)).argmin()

        t_min = t_axis[t_min_idx]
        t_max = t_axis[t_max_idx]

        print("Energy range: {} to {} meV\nTime range {} to {} s".format(e_vals[e_min_idx, 0],
                                                                         e_vals[e_max_idx, 0],
                                                                         t_vals[0, t_min_idx],
                                                                         t_vals[0, t_max_idx]))

        i_vals_cropped = np.array(i_vals[e_min_idx:e_max_idx, t_min_idx:t_max_idx], dtype=np.float32)
        t_bins_cropped = np.array(time_bins[t_min_idx:t_max_idx], dtype=np.float32)
        e_bins_cropped = np.array(e_bins[e_min_idx:e_max_idx], dtype=np.float32)

        num_time_bins = len(t_bins_cropped)
        num_ener_bins = len(e_bins_cropped)

        EInt = np.zeros((e_bins_cropped.shape))
        EInt[0] = 0.0
        Flux = np.zeros(i_vals_cropped.shape)
        Tsum = 0.0

        for i in range(0, e_max_idx - e_min_idx):
            EInt[i] = Tsum
            for j in range(0, t_max_idx - t_min_idx):
                Tsum += i_vals_cropped[i, j]
                Flux[i, j] = Tsum

        Total = Tsum

        mf = cl.mem_flags

        self.num_time_bins = np.int32(num_time_bins)
        self.num_ener_bins = np.int32(num_ener_bins)

        self.flux = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Flux.flatten().astype(np.float32))
        self.time_bins = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=t_bins_cropped)
        self.ener_bins = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=e_bins_cropped)
        self.e_int = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=EInt.astype(np.float32))
        self.total = np.float32(Total)

        print(self.total)

    def gen_prg(self, queue, N, neutron_buf, intersection_buf):
        self.prg.generate_neutrons(queue, (N,), None, neutron_buf, intersection_buf,
                                   self.pos,
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
                                   np.int32(N))
