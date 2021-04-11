#include "geom.h"
#include "consts.h"

__kernel void transform(__global float16* neutrons,
    float3 const pos, float3 const rot) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];

    /* Already terminated? */
    if (NEUTRON_DIE  > 0.f) {
        return;
    }

    NEUTRON_POS -= pos;

    NEUTRON_POS = frame_derotate(NEUTRON_POS, rot);
    NEUTRON_VEL = frame_derotate(NEUTRON_VEL, rot);
    NEUTRON_POL = frame_derotate(NEUTRON_POL, rot);

    neutrons[global_addr] = neutron;
}

__kernel void untransform(__global float16* neutrons,
    float3 const pos, float3 const rot) {

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];

    /* Already terminated? */
    if (NEUTRON_DIE  > 0.f) {
        return;
    }

    NEUTRON_POS= frame_rotate(NEUTRON_POS, rot);
    NEUTRON_POS+= pos;
    NEUTRON_VEL = frame_rotate(NEUTRON_VEL, rot);
    NEUTRON_POL = frame_rotate(NEUTRON_POL, rot);

    neutrons[global_addr] = neutron;
}