__kernel void terminate(__global float16* neutrons,
    __global float8* intersections, uint const restore_neutron) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];
    float8 intersection = intersections[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.0f) {
        return;
    }

    /* Didn't intersect anything, terminate */
    if (length(intersection.s456) == 0.0f && restore_neutron == 0) {
        neutron.sf = 1.f;
    }

    neutrons[global_addr] = neutron;
}