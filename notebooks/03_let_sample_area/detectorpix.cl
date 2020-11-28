#include "consts.h"

void atomicAdd_g_f(volatile __global float *addr, float val)
{
    union {
        unsigned int u32;
        float        f32;
    } next, expected, current;
	current.f32    = *addr;
    do {
	   expected.f32 = current.f32;
        next.f32     = expected.f32 + val;
		current.u32  = atomic_cmpxchg( (volatile __global unsigned int *)addr, 
                            expected.u32, next.u32);
    } while( current.u32 != expected.u32 );
}

int find_idx(float val, float3 binning) {
    int idx;
    if(binning.s0<=val && val<=binning.s2) {    
        idx = floor((val -  binning.s0) / binning.s1);
    } else {
        idx = -1;
    }

    return idx;
}



__kernel void detector(__global float16 *neutrons,
                       __global float8 *intersections,
                        __global uint *iidx,
                       uint const comp_idx, 
                       volatile __global float *histogram,
                       volatile __global float *histogram_err,
                       float3 const y_binning,
                       float3 const theta_binning,
                       float3 const tof_binning,
                       uint const y_num_bins,
                       uint const theta_num_bins,
                       uint const tof_num_bins,
                       uint const restore_neutron)
{
  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];

  int this_iidx, theta_idx, y_idx, det_num, flattened_idx, tof_idx;
  float y_val, theta_val, tof_val;
  this_iidx = iidx[global_addr];

  /* Check we are scattering from the intersected component -------------- */
  if (!(this_iidx == comp_idx))
      return;

  /* Check termination flag ---------------------------------------------- */
  if (neutron.sf > 0.f) 
      return;

  /* Perform scattering here --------------------------------------------- */

  theta_val = degrees(atan2(
        intersection.s4,
        intersection.s6
  ));

  y_val = intersection.s5;

  theta_idx = find_idx(theta_val, theta_binning);
  y_idx = find_idx(y_val, y_binning);

  det_num = y_idx * theta_num_bins + theta_idx;

  tof_val = (1.0e6f)*(neutron.sa + intersection.s7);
  tof_idx = find_idx(tof_val, tof_binning);

  if (!((theta_idx == -1) || (y_idx == -1)  || (tof_idx == -1))) {
    flattened_idx = det_num * tof_num_bins + tof_idx;
    atomicAdd_g_f(&histogram[flattened_idx], (float)neutron.s9);
    atomicAdd_g_f(&histogram_err[flattened_idx], (float)neutron.s9*neutron.s9);
  }

  if (restore_neutron == 0) {
    neutron.s012 = intersection.s456;
    neutron.sa += intersection.s7;
    neutron.sf = 1.f;
  }

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
