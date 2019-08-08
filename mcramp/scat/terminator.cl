__kernel void terminate(__global double16* neutrons,
    __global double8* intersections, uint const restore_neutron) {

    uint global_addr = get_global_id(0);

    double16 neutron = neutrons[global_addr];
    double8 intersection = intersections[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.) {
        return;
    }

    /* Didn't intersect anything, terminate */
    if (length(intersection.s456) == 0. && restore_neutron == 0) {
        neutron.sf = 1.;
    }

    neutrons[global_addr] = neutron;
}