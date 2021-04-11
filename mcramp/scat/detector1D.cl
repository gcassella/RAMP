#include "consts.h"

inline void atomicAdd_g_f(volatile __global float *addr, float val)
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



__kernel void detector(__global float16 *neutrons,
                       __global float8 *intersections, __global uint *iidx,
                       uint const comp_idx, volatile __global float *histogram, 
                       volatile __global float *histogram_err,
                       float3 const binning, uint const restore_neutron, uint const var)
{

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];
  float var_val, min_var, step_var, max_var;

  uint this_iidx, idx;
  this_iidx = iidx[global_addr];

  if (!(this_iidx == comp_idx))
  {
      return;
  }

  if (NEUTRON_DIE  > 0.f)
  {
      return;
  }

  min_var = binning.s0;
  step_var = binning.s1;
  max_var = binning.s2;

  if (var == 0)
    var_val = VS2E*pow(length(NEUTRON_VEL), 2.0f);
  else if (var == 1)
    var_val = degrees(atan2(INTERSECTION_X2, INTERSECTION_Z2));
  else if (var == 2)
    var_val = (1.0e6f)*(NEUTRON_TOF + INTERSECTION_T2);
  else if (var == 3)
    var_val = 2*M_PI / (V2K*length(NEUTRON_VEL));

  if(min_var<=var_val && var_val<=max_var) {    
    idx = floor((var_val -  min_var) / step_var);
    atomicAdd_g_f(&histogram[idx], (float)NEUTRON_P);
    atomicAdd_g_f(&histogram_err[idx], (float)NEUTRON_P*NEUTRON_P);
  }
  
  
  if (restore_neutron == 0) {
    NEUTRON_POS = INTERSECTION_POS2;
    NEUTRON_TOF += INTERSECTION_T2;
    NEUTRON_DIE  = 1.f;
  }

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
