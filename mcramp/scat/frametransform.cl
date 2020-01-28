#include "geom.h"

__kernel void transform(__global float16* neutrons,
    float3 const pos, float3 const rot) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.f) {
        return;
    }

    neutron.s012 -= pos;

    neutron.s012 = frame_derotate(neutron.s012, rot);
    neutron.s345 = frame_derotate(neutron.s345, rot);

    neutrons[global_addr] = neutron;
}

__kernel void untransform(__global float16* neutrons,
    float3 const pos, float3 const rot) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.f) {
        return;
    }

    neutron.s012 = frame_rotate(neutron.s012, rot);
    neutron.s012 += pos;
    neutron.s345 = frame_rotate(neutron.s345, rot);

    neutrons[global_addr] = neutron;
}