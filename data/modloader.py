#%%

import os, re

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import simps

os.chdir(r"C:\Users\Gino\Documents\RAMP\RAMP\data")

with open('Let_Base.mcstas', 'r') as mod_file:
    lines = mod_file.read()

    # Regex: (?:\d*\.\d*[e][+\-]?\d+ \d*\.\d*[e][+\-]?\d+ \(\d*\.\d*[e][+\-]?\d+\) [\n]) will capture lines
    pattern = r"(?:\d*\.\d*[e][+\-]?\d+ \d*\.\d*[e][+\-]?\d+ \(\d*\.\d*[e][+\-]?\d+\) )"
    result = re.findall(pattern, lines)

    mod_file.seek(0)

    line = ''
    time_linenum = 0
    while not 'time' in line:
        line = mod_file.readline()
        time_linenum += 1

    total_linenum = time_linenum
    while not 'total' in line:
        line = mod_file.readline()
        total_linenum += 1

    time_bins = total_linenum - time_linenum - 1

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
            e_bins.append(lower_bound)

    e_bins.append(upper_bound)        

data = np.array([r.strip().replace('(','').replace(')','') for r in result])
data = np.reshape(data, (int(np.size(data) / time_bins), time_bins))

t_axis = np.array([l.split(' ')[0] for l in data[0]]).astype(float) / 1e8
e_axis = np.array(e_axis) * 1e9

t_vals, e_vals = np.meshgrid(t_axis, e_axis)

i_vals = np.zeros(e_vals.shape)

for i in range(e_vals.shape[0]):
    data_chunk = np.array([l.split(' ')[1] for l in data[i]]).astype(float)
    
    i_vals[i] = data_chunk

#%%

#
# - Determine which energy bins lie in crop range
# - Determine which time bins lie in crop range
# - Take corresponding slice from the numpy array
#

e_min = 0.0  # [meV]
e_max = 40.0 # [meV]

t_min = 0.0  # [s]
t_max = 0.001  # [s]

e_min_idx = (np.abs(e_axis - e_min)).argmin()
e_max_idx = (np.abs(e_axis - e_max)).argmin()

t_min_idx = (np.abs(t_axis - t_min)).argmin()
t_max_idx = (np.abs(t_axis - t_max)).argmin()

print("Energy range: {} to {} meV\nTime range {} to {} meV".format(e_vals[e_min_idx,0], e_vals[e_max_idx,0], t_vals[0,t_min_idx], t_vals[0,t_max_idx]))

i_vals_cropped = i_vals[e_min_idx:e_max_idx, t_min_idx:t_max_idx]
t_vals_cropped = t_vals[e_min_idx:e_max_idx, t_min_idx:t_max_idx]
e_vals_cropped = e_vals[e_min_idx:e_max_idx, t_min_idx:t_max_idx]
t_bins_cropped = time_bins[t_min_idx:t_max_idx]
e_bins_cropped = e_bins[e_min_idx:e_max_idx]

# from this point onward everything needs to be ported to be carried out inside of a kernel - just prototyping in python first to understand the process

(M, N) = i_vals_cropped.shape
i_vals_reshaped = np.reshape(i_vals_cropped, (M * N, ))

total_flux = sum(i_vals_reshaped)
i_vals_normed = i_vals_reshaped / total_flux

i_vals_cdf = np.cumsum(i_vals_normed)

num_e_bins = e_max_idx - e_min_idx
num_t_bins = t_max_idx - t_min_idx

vals = []

deviate = np.random.random_sample()
idx = (np.abs(i_vals_cdf - deviate)).argmin()
e_idx = idx // num_t_bins
t_idx = idx % num_t_bins  

#sample = e_vals_cropped[n, m]

# TODO check sampling is accurate - compare with LET mod spectrum?

# TODO implement polynomial interpolation around bins
#   - some kind of linear interpolation in time axis? double check
#   - polynomial interpolation around energy axis, find intensity for n points around chosen point, linearly interpolate energy 'coordinate', polyfit the intensity at n points, obtain intensity at energy coordinate

deviate = np.random.random_sample()

time_spread = t_bins_cropped[t_idx + 1] - t_bins_cropped[t_idx]
time_fluxdiff = i_vals_cropped[e_idx, t_idx + 1] - i_vals_cropped[e_idx, t_idx]

something = total_flux * deviate - i_vals_cropped[e_idx, t_idx]
something /= time_fluxdiff
time_val = t_bins_cropped[t_idx]*something # ask about this?

deviate = np.random.random_sample()

# TODO from here down is wrong - E_int is the integrated intensity in the energy bins - modify code to match what is in McStas and try to understand it

energy_spread = e_bins[e_idx + 1] - e_bins[e_idx]

e_val = e_bins[e_idx] + deviate * energy_spread

# energy point chosen, now interpolate to desired intensity

interpol_start = e_idx - 3 if (e_idx > 3) else 0
interpol_end = e_idx + 3 if (len(e_bins_cropped) - 1 - e_idx) > 3 else len(e_bins_cropped) - 1

interpol_i_vals = i_vals_cropped[interpol_start:interpol_end, t_idx]
interpol_e_vals = e_bins[interpol_start:interpol_end]

p = np.polyfit(interpol_e_vals, interpol_i_vals, deg=6)
res = p[0]*e_val