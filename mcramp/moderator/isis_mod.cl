#include "rand.h"

int find_closest_in_cdf(__global float* arr, int arrlen, float val) {
  float nearest = fabs(arr[0] - val);
  int nearest_idx = 0;
  for (int i = 1; i < arrlen; i++) {
    if (fabs(arr[i] - val) <= nearest) {
      nearest = fabs(arr[i] - val);
      nearest_idx = i;
    } else {
      break;
    }
  }

  return nearest_idx;
}

__kernel void generate_neutrons(__global float16* neutrons,
    __global float8* intersections, float3 const pos, float2 const mod_dim,
    float2 const target_dim, float const target_dist, float const E_min,
    float const E_max, int const num_time_bins, int const num_ener_bins,
    __global float* flux, __global float* time_bins, __global float* ener_bins,
    global float* E_int, float const total, float const str_area,
    float const time_offset, int const num_sim) {

  int global_addr, idx, Epnt, Tpnt, interpol_start, interpol_end;
  float deviate, time_val, time_range, time_spread, R, ener_val, ener_spread, 
  accumulator, Pj, vel, Dx, Dy, displacement;
  float16 neutron;
  

  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];

  deviate = total * rand(&neutron, global_addr);

  idx = find_closest_in_cdf(flux, num_time_bins * num_ener_bins, deviate);

  Epnt = floor((float)idx / (float)num_time_bins);
  Tpnt = fmod((float)idx, (float)num_time_bins);

  // FIXME: find a better method for sampling the full range of bins so the entire
  // specified range is actually sampled
  if(Tpnt >= num_time_bins-1 || Epnt >= num_ener_bins - 1) {
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;
  }

  time_val = time_bins[Tpnt];

  time_range = time_bins[Tpnt + 1] - time_bins[Tpnt];
  time_spread = flux[Epnt * num_time_bins + Tpnt + 1] - flux[Epnt * num_time_bins + Tpnt];
  R = deviate - flux[Epnt * num_time_bins + Tpnt];
  R /= time_spread;
  time_val += time_range * R;

  ener_val = ener_bins[Epnt];
  ener_spread = E_int[Epnt+1] - E_int[Epnt];

  interpol_start = (Epnt > 3) ? Epnt - 3 : 0;
  interpol_end = ((num_ener_bins - 1 - Epnt) > 3) ? Epnt + 3 : (num_ener_bins) - 1;

  deviate = E_int[Epnt] + ener_spread * rand(&neutron, global_addr);
  accumulator = 0.0;
  for (int i = interpol_start; i <= interpol_end; i++) {
    Pj = 1.0;
    for (int j = interpol_start; j <= interpol_end; j++) {
      if (j != i) {
        Pj *= (deviate - E_int[j]) / (E_int[i] - E_int[j]);
      } else {
        Pj *= 1;
      }
    }
    Pj *= ener_bins[i];

    accumulator += Pj;
  }
  ener_val = accumulator;

  neutron.s0 = pos.x + mod_dim.x*(0.5 - rand(&neutron, global_addr));
  neutron.s1 = pos.y + mod_dim.y*(0.5 - rand(&neutron, global_addr));
  neutron.s2 = pos.z;

  vel = 438.01*sqrt(ener_val);
  Dx = target_dim.x*(0.5 - rand(&neutron, global_addr)) - neutron.s0;
  Dy = target_dim.y*(0.5 - rand(&neutron, global_addr)) - neutron.s1;
  displacement = sqrt(Dx*Dx + Dy*Dy + target_dist*target_dist);

  neutron.s3 = vel*Dx/displacement;
  neutron.s4 = vel*Dy/displacement;
  neutron.s5 = vel*target_dist/displacement;

  // Initialize weight
  // TODO: with buffer chunking this will exaggerate the intensity by the number
  // of buffer chunks (assuming equally sized chunks). Must find a way to fix!
  neutron.s9 = str_area*total*6.2415093e+12 / num_sim;

  // Initialize time
  neutron.sa = time_val;

  // Revive terminated neutrons
  neutron.sf = 0.;

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}