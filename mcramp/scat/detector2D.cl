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
        idx = round((val -  binning.s0) / binning.s1);
    } else {
        idx = -1;
    }

    return idx;
}



__kernel void detector(__global float16 *neutrons,
                       __global float8 *intersections, __global uint *iidx,
                       uint const comp_idx, volatile __global float *histogram,
                       volatile __global float *histogram_err,
                       float3 const axis1_binning, float3 const axis2_binning, 
                       uint const axis1_numbins, uint const axis2_numbins,
                       uint const axis1_var, uint const axis2_var,
                       uint const restore_neutron)
{

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];

  int this_iidx, axis1_idx, axis2_idx, flattened_idx;
  float axis1_val, axis2_val;
  this_iidx = iidx[global_addr];

  /* Check we are scattering from the intersected component -------------- */
  if (!(this_iidx == comp_idx))
      return;

  /* Check termination flag ---------------------------------------------- */
  if (NEUTRON_DIE  > 0.f) 
      return;

  /* Perform scattering here --------------------------------------------- */

  // This code snippet obviously violates DRY but any kind of solution that may
  // be more elegant and non-repetitive will ultimately only act to obscure
  // meaning, so I have chosen to stick with the naieve approach.

  // Find axis 1 bin
  switch(axis1_var) {
    case 0 :
      axis1_val = INTERSECTION_X2;
      break;
    case 1 :
      axis1_val = INTERSECTION_Y2;
      break;
    case 2 :
      axis1_val = degrees(atan2(
        INTERSECTION_X2,
        INTERSECTION_Z2
      ));
      break;
    case 3 :
      axis1_val = degrees(atan2(
        INTERSECTION_Y2, 
        sqrt(INTERSECTION_X2*INTERSECTION_X2 + INTERSECTION_Z2*INTERSECTION_Z2)
      ));
      break;
    case 4 :
      axis1_val = sign(INTERSECTION_Y2)*degrees(acos(
        INTERSECTION_X2/
        sqrt(INTERSECTION_X2*INTERSECTION_X2 + INTERSECTION_Y2*INTERSECTION_Y2)
      ));
      break;
    case 5 :
      axis1_val = (1.0e6f)*(NEUTRON_TOF + INTERSECTION_T2);
      break;
    case 6 :
      axis1_val = degrees(atan2(NEUTRON_VX, INTERSECTION_Y2));
      break;
    case 7 :
      axis1_val = degrees(atan2(NEUTRON_VY, INTERSECTION_Y2));
      break;
    case 8 :
      axis1_val = 2*M_PI / (V2K*length(NEUTRON_VEL));
      break;
    case 9 :
      axis1_val = VS2E*pow(length(NEUTRON_VEL), 2.0f);
      break;
    case 10 :
      axis1_val = sign(INTERSECTION_X2)*degrees(atan2(
        sqrt(INTERSECTION_X2*INTERSECTION_X2 + INTERSECTION_Y2*INTERSECTION_Y2),
        INTERSECTION_Z2
      ));
    default:
      break;
  }

  axis1_idx = find_idx(axis1_val, axis1_binning);

  // Find axis 2 bin
  switch(axis2_var) {
    case 0 :
      axis2_val = INTERSECTION_X2;
      break;
    case 1 :
      axis2_val = INTERSECTION_Y2;
      break;
    case 2 :
      axis2_val = degrees(atan2(
        INTERSECTION_X2,
        INTERSECTION_Z2
      ));
      break;
    case 3 :
      axis2_val = degrees(atan2(
        INTERSECTION_Y2, 
        sqrt(INTERSECTION_X2*INTERSECTION_X2 + INTERSECTION_Z2*INTERSECTION_Z2)
      ));
      break;
    case 4 :
      axis2_val = sign(INTERSECTION_Y2)*degrees(acos(
        INTERSECTION_X2/
        sqrt(INTERSECTION_X2*INTERSECTION_X2 + INTERSECTION_Y2*INTERSECTION_Y2)
      ));
      break;
    case 5 :
      axis2_val = (1.0e6f)*(NEUTRON_TOF + INTERSECTION_T2);
      break;
    case 6 :
      axis2_val = degrees(atan2(NEUTRON_VX, NEUTRON_VZ));
      break;
    case 7 :
      axis2_val = degrees(atan2(NEUTRON_VY, NEUTRON_VZ));
      break;
    case 8 :
      axis2_val = 2*M_PI / (V2K*length(NEUTRON_VEL));
      break;
    case 9 :
      axis2_val = VS2E*pow(length(NEUTRON_VEL), 2.0f);
      break;
    case 10 :
      axis2_val = sign(INTERSECTION_X2)*degrees(atan2(
        sqrt(INTERSECTION_X2*INTERSECTION_X2 + INTERSECTION_Y2*INTERSECTION_Y2),
        INTERSECTION_Z2
      ));
    default:
      break;
  }

  axis2_idx = find_idx(axis2_val, axis2_binning);


  if (!((axis1_idx == -1) || (axis2_idx == -1))) {
    flattened_idx = axis1_idx * axis2_numbins + axis2_idx;
    atomicAdd_g_f(&histogram[flattened_idx], (float)neutron.s9);
    atomicAdd_g_f(&histogram_err[flattened_idx], (float)neutron.s9*neutron.s9);
  }

  if (restore_neutron == 0) {
    NEUTRON_POS= intersection.s456;
    NEUTRON_TOF += intersection.s7;
    NEUTRON_DIE  = 1.f;
  }

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
