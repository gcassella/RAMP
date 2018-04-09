__kernel void terminate(__global float16* neutrons,
    __global float8* intersections) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];
    float8 intersection = intersections[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.) {
        return;
    }

    /* Didn't intersect anything, terminate */
    if (length(intersection.s456) == 0.) {
        neutron.sf = 1.;
    }

    neutrons[global_addr] = neutron;
}