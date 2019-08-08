#include "consts.h"

void AtomicAdd(volatile global double *source, const double operand) {
    union {
        unsigned int intVal;
        double doubleVal;
    } newVal;
    union {
        unsigned int intVal;
        double doubleVal;
    } prevVal;
 
    do {
        prevVal.doubleVal = *source;
        newVal.doubleVal = prevVal.doubleVal + operand;
    } while (atomic_cmpxchg((volatile global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}


int find_idx(double val, double3 binning) {
    uint idx;
    if(binning.s0<=val && val<=binning.s2) {    
        idx = round((val -  binning.s0) / binning.s1);
    } else {
        idx = -1;
    }

    return idx;
}

__kernel void detector(__global double16 *neutrons,
                       __global double8 *intersections, __global uint *iidx,
                       uint const comp_idx, volatile __global double *histogram,
                       double3 const axis1_binning, double3 const axis2_binning, 
                       uint const axis1_numbins, uint const axis2_numbins,
                       uint const shape, uint const restore_neutron)
{

  uint global_addr = get_global_id(0);
  double16 neutron = neutrons[global_addr];
  double8 intersection = intersections[global_addr];

  int this_iidx, axis1_idx, axis2_idx, flattened_idx;
  this_iidx = iidx[global_addr];

  /* Check we are scattering from the intersected component -------------- */
  if (!(this_iidx == comp_idx))
      return;

  /* Check termination flag ---------------------------------------------- */
  if (neutron.sf > 0.) 
      return;

  /* Perform scattering here --------------------------------------------- */

  // Find axis 1 bin

  if (shape == 0 || shape == 4) { // Plane detector, axis1 is x
    double x_val = intersection.s4;
    axis1_idx = find_idx(x_val, axis1_binning);
  } else if (shape == 1 || shape == 2) { // Banana detector, axis1 is 2theta
    double3 sample_to_det = intersection.s456;
    double theta_val = degrees(atan2(sample_to_det.s0, sample_to_det.s2));
    axis1_idx = find_idx(theta_val, axis1_binning);
  } else if (shape == 3) {
    axis1_idx = find_idx(degrees(atan2(neutron.s3, neutron.s5)), axis2_binning);
  }

  // Find axis 2 bin

  if (shape == 0) { // Plane detector, axis2 is y
    double y_val = intersection.s5;
    axis2_idx = find_idx(y_val, axis2_binning);
  } else if (shape == 1) { // Banana detector, axis2 is alpha
    double3 sample_to_det = intersection.s456;
    double alpha_val = degrees(atan2(sample_to_det.s1, sample_to_det.s2));
    axis2_idx = find_idx(alpha_val, axis2_binning);
  } else if (shape == 2) {
    axis2_idx = find_idx(neutron.sa+intersection.s7, axis2_binning);
  } else if (shape == 3) {
    axis2_idx = find_idx(degrees(atan2(neutron.s4, neutron.s5)), axis2_binning);
  } else if (shape == 4) {
    axis2_idx = find_idx(degrees(atan2(neutron.s3, neutron.s5)), axis2_binning);
  }


  if (!((axis1_idx == -1) || (axis2_idx == -1))) {
    flattened_idx = axis1_idx * axis2_numbins + axis2_idx;
    neutron.se = flattened_idx;
  } else {
    neutron.se = -1;
  }

  if (restore_neutron == 0) {
    neutron.s012 = intersection.s456;
    neutron.sa += intersection.s7;
    neutron.sf = 1.;
  }

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
