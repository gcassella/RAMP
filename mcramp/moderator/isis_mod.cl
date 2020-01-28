#include "rand.h"
#include "consts.h"

int find_closest_in_cdf(__global float* arr, int arrlen, float val) {
  // Find index of last element in arr which is less than val
  for (int i = 0; i < arrlen; i++) {
    if ((arr[i] - val) > 0.0f) {
      return i - 1;
    }
  }

  return 0;
}

__kernel void generate_neutrons(__global float16* neutrons,
    __global float8* intersections, float2 const mod_dim,
    float2 const target_dim, float const target_dist, float const E_min,
    float const E_max, int const num_time_bins, int const num_ener_bins,
    __global float* flux, __global float* time_bins, __global float* ener_bins,
    global float* E_int, float const total, float const str_area,
    float const time_offset, int const num_sim, int const seed) {

  int global_addr, idx, Epnt, Tpnt, interpol_start, interpol_end;
  float deviate, time_val, time_range, time_spread, R, ener_val, ener_spread, 
  accumulator, Pj, vel, Dx, Dy;
  float16 neutron;
  
  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];
  
  // Clear old neutron data and seed

  neutron = (float16){0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,(float)seed+neutron.sb,0.0f,0.0f,0.0f,0.0f};

  deviate = total * rand(&neutron, global_addr);

  idx = find_closest_in_cdf(flux, num_time_bins * num_ener_bins, deviate);

  Epnt = find_closest_in_cdf(E_int, num_ener_bins, deviate);
  Tpnt = fmod((float)idx, (float)num_time_bins);

  if(Tpnt >= num_time_bins-1) {
    Tpnt -= 1;
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
  interpol_end = ((num_ener_bins - Epnt - 1) > 3) ? Epnt + 3 : (num_ener_bins - 1);
  
  deviate = E_int[Epnt] + ener_spread * rand(&neutron, global_addr);
  accumulator = 0.0f;
  for (int i = interpol_start; i <= interpol_end; i++) {
    Pj = 1.0f;
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

  neutron.s0 = mod_dim.x*(0.5f - rand(&neutron, global_addr));
  neutron.s1 = mod_dim.y*(0.5f - rand(&neutron, global_addr));
  neutron.s2 = 0.0f;

  vel = SE2V*sqrt(ener_val);
  Dx = target_dim.x*(0.5f - rand(&neutron, global_addr)) - neutron.s0;
  Dy = target_dim.y*(0.5f - rand(&neutron, global_addr)) - neutron.s1;
  
  neutron.s345 = vel*normalize((float3)( Dx, Dy, target_dist ));

  // Initialize weight
  // TODO: with buffer chunking this will exaggerate the intensity by the number
  // of buffer chunks (assuming equally sized chunks). Must find a way to fix!
  neutron.s9 = str_area*total*6.2415093e+12f / num_sim;

  // Initialize time
  neutron.sa = time_val;

  // Revive terminated neutrons
  neutron.sf = 0.f;

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}