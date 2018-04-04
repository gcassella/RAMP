__kernel void detector(__global float16* neutrons,
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global uint* histogram, float3 const det_pos,
  float3 const binning, uint const var) {

  uint global_addr = get_global_id(0);
  
  float16 neutron;
  float8 intersection;
  float varval;
  float minvar, maxvar, stepvar;
  uint idx;

  uint this_iidx = iidx[global_addr];

  minvar = binning.s0;
  stepvar = binning.s1;
  maxvar = binning.s2;

  neutron = neutrons[global_addr];
  intersection = intersections[global_addr];

  if(!(this_iidx == comp_idx)) {
    return;
  }

  if(neutron.sf > 0.) {
    return;
  }

  if (length(intersection.s4567) != 0.) {
    switch(var) {
      case 0:
        varval = acos(dot(normalize(intersection.s46 - det_pos.s02), (float2)( 0.0f, 1.0f )));
        break;
      case 1:
        varval = intersection.s7;
        break;
      default:
        break;
    }

    if(minvar<=varval && varval<=maxvar) {
      idx = round((varval - minvar) / stepvar);
      neutron.s012a = intersection.s4567;

      histogram[global_addr] = idx;
    }

    neutron.sf = 1.;
  }

  neutrons[global_addr] = neutron;
}