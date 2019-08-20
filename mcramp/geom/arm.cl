__kernel void arm(__global float16* neutrons, 
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];
    float8 intersection = intersections[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.) {
        return;
    }

    // Have to add some intersection to avoid neutron termination
    intersection.s456 = (float3){ 1.0, 1.0, 1.0 };
    iidx[global_addr] = comp_idx;

    intersections[global_addr] = intersection;
    neutrons[global_addr] = neutron;
}