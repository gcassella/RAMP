#include "consts.h"

__kernel void terminate(__global float16* neutrons,
    __global float8* intersections, uint const restore_neutron) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];
    float8 intersection = intersections[global_addr];

    /* Already terminated? */
    if (NEUTRON_DIE  > 0.0f) {
        return;
    }

    /* Didn't intersect anything, terminate */
    if (length(intersection.s456) == 0.0f && restore_neutron == 0) {
        NEUTRON_DIE  = 1.f;
    }

    neutrons[global_addr] = neutron;
}