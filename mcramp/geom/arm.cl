__kernel void arm(__global double16* neutrons, 
  __global double8* intersections, __global uint* iidx,
  uint const comp_idx) {

    uint global_addr = get_global_id(0);

    double16 neutron = neutrons[global_addr];
    double8 intersection = intersections[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.) {
        return;
    }

    // Have to add some intersection to avoid neutron termination
    intersection.s456 = (double3){ 1.0, 1.0, 1.0 };
    iidx[global_addr] = comp_idx;

    intersections[global_addr] = intersection;
    neutrons[global_addr] = neutron;
}