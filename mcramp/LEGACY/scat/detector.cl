inline void AtomicAdd(volatile __global double *source, double const operand) {
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
    } while (atomic_cmpxchg((volatile __global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

__kernel void detector(__global double16* neutrons,
  __global double8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global double* histogram, double3 const det_pos,
  double3 const binning, uint const var) {

  uint global_addr = get_global_id(0);
  
  double16 neutron;
  double8 intersection;
  double3 planediff_padded, planediff_cross;
  double2 planediff;
  double varval;
  double minvar, maxvar, stepvar;
  uint idx;

  uint this_iidx = iidx[global_addr];
  
  neutron = neutrons[global_addr];
  intersection = intersections[global_addr];

  if(!(this_iidx == comp_idx)) {
    return;
  }

  if(neutron.sf > 0.) {
    return;
  }

  minvar = binning.s0;
  stepvar = binning.s1;
  maxvar = binning.s2;

  switch(var) {
    case 0:
      planediff = intersection.s46 - det_pos.s02;
      planediff_padded = (double3)( planediff.s0, 0.0f, planediff.s1 );
      planediff_cross = cross(normalize(planediff_padded), (double3)( 0.0f, 1.0f, 0.0f ));
      varval = acos(dot(normalize(planediff), (double2)( 0.0f, 1.0f )))
                * sign(planediff_cross.s2);
      break;
    case 1:
      varval = neutron.sa+intersection.s7;
      break;
    case 2:
      varval = intersection.s5 - det_pos.s1;
      break;
    default:
      break;
  }

  if(minvar<=varval && varval<=maxvar) {    
    idx = round((varval -  minvar) / stepvar);
    AtomicAdd(&histogram[idx], neutron.s9);
  }

  iidx[global_addr] = 0;
  neutron.s012 = intersection.s456;
  neutron.sa += intersection.s7;
  neutron.sc = comp_idx;
  neutron.sf = 1.;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}