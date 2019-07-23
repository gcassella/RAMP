#ifndef M_PI
#define M_PI 3.14159265
#endif

void AtomicAdd(volatile global float *source, const float operand) {
    union {
        unsigned int intVal;
        float floatVal;
    } newVal;
    union {
        unsigned int intVal;
        float floatVal;
    } prevVal;
 
    do {
        prevVal.floatVal = *source;
        newVal.floatVal = prevVal.floatVal + operand;
    } while (atomic_cmpxchg((volatile global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}


int find_idx(float val, float3 binning) {
    uint idx;
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
                       float3 const pos, float3 const axis1_binning, 
                       float3 const axis2_binning, uint const axis1_numbins,
                       uint const axis2_numbins, uint const shape)
{

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];

  int this_iidx, axis1_idx, axis2_idx, flattened_idx;
  this_iidx = iidx[global_addr];

  if (!(this_iidx == comp_idx))
  {
      return;
  }

  if (neutron.sf > 0.)
  {
      return;
  }

  // Find axis 1 bin

  if (shape == 0) { // Plane detector, axis1 is x
    float x_val = intersection.s4 - pos.s0;
    axis1_idx = find_idx(x_val, axis1_binning);
  } else if (shape == 1 || shape == 2) { // Banana detector, axis1 is 2theta
    float3 sample_to_det = intersection.s456 - pos;
    float theta_val = degrees(atan2(sample_to_det.s0, sample_to_det.s2));
    axis1_idx = find_idx(theta_val, axis1_binning);
  }

  // Find axis 2 bin

  if (shape == 0) { // Plane detector, axis2 is y
    float y_val = intersection.s5 - pos.s1;
    axis2_idx = find_idx(y_val, axis2_binning);
  } else if (shape == 1) { // Banana detector, axis2 is alpha
    float3 sample_to_det = intersection.s456 - pos;
    float alpha_val = degrees(atan2(sample_to_det.s1, sample_to_det.s2));
    axis2_idx = find_idx(alpha_val, axis2_binning);
  } else if (shape == 2) {
    axis2_idx = find_idx(neutron.sa+intersection.s7, axis2_binning);
  }


  if (!((axis1_idx == -1) || (axis2_idx == -1))) {
    flattened_idx = axis1_idx * axis2_numbins + axis2_idx;
    neutron.se = flattened_idx;
  } else {
    neutron.se = -1;
  }
  // AtomicAdd(&histogram[flattened_idx], neutron.s9);

  iidx[global_addr] = 0;
  neutron.s012 = intersection.s456;
  neutron.sa += intersection.s7;
  neutron.sc = comp_idx;
  neutron.sf = 1.;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
