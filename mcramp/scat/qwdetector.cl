inline void AtomicAdd(volatile __global float *source, float const operand) {
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
    } while (atomic_cmpxchg((volatile __global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

__kernel void detector(__global float16* neutrons,
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global float* histogram, float3 const det_pos,
  float const ei, float3 const ki_norm, 
  float3 const q_binning, float3 const de_binning) {

  uint global_addr = get_global_id(0);
  uint this_iidx = iidx[global_addr];
  
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];

  if(!(this_iidx == comp_idx)) {
    return;
  }

  if(neutron.sf > 0.) {
    return;
  }

  float3 planediff_padded, planediff_cross, ki, kf, kf_norm, Q;
  float2 planediff;
  float q_minvar, q_maxvar, q_stepvar, q_nbins, q_idx,
        de_minvar, de_maxvar, de_stepvar, de_nbins, de_idx;
  float ki_mag, kf_mag, de, ef, Q_mag;
  uint idx;
  uint missed = 0;

  ki_mag = 0.69297*sqrt(ei);
  ki = ki_norm * ki_mag;

  kf_norm = normalize(neutron.s345);
  kf_mag = 1.58818*pow(10., -3.)*length(neutron.s345);
  kf = kf_norm * kf_mag;

  ef = 2.082*kf_mag*kf_mag;
  de = ei - ef;

  Q = ki - kf;
  Q_mag = length(Q);

  q_minvar  = q_binning.s0;
  q_stepvar = q_binning.s1;
  q_maxvar  = q_binning.s2;
  q_nbins   = round((q_maxvar - q_minvar) / q_stepvar);
  
  de_minvar  = de_binning.s0;
  de_stepvar = de_binning.s1;
  de_maxvar  = de_binning.s2;
  de_nbins   = round((de_maxvar - de_minvar) / de_stepvar);

  if(q_minvar<=Q_mag && Q_mag<=q_maxvar) {    
    q_idx = round((Q_mag - q_minvar) / q_stepvar);
  } else {
    missed = 1;
  }

  if(de_minvar<=de && de<=de_maxvar) {    
    de_idx = round((de - de_minvar) / de_stepvar);
  } else {
    missed = 1;
  }

  if (missed == 0) {
    idx = q_idx*de_nbins + de_idx;
    AtomicAdd(&histogram[idx], neutron.s9);
  }

  iidx[global_addr] = 0;

  neutron.s012 = intersection.s456;
  neutron.sa += intersection.s7;
  neutron.sc = comp_idx;
  neutron.sf = 1.;

  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}