#include "consts.h"

__kernel void intersect_plane(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, float const width, 
    float const height, uint const orientation) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection     = intersections[global_addr];

    /* Check termination flag */
    if (NEUTRON_DIE  > 0.f) 
        return;


    /* Perform raytracing here */

    float3 vel = NEUTRON_VEL;
    float3 pos = NEUTRON_POS;

    float t, x, y;
    
    if (orientation == 0) {
        t = (-pos.s2) / vel.s2;
        x = pos.s0 + t*vel.s0;
        y = pos.s1 + t*vel.s1;
    } else if (orientation == 1) {
        t = (-pos.s0) / vel.s0;
        x = pos.s2 + t*vel.s2;
        y = pos.s1 + t*vel.s1;
    }

    if ((fabs(x) < width / 2.0f) && (fabs(y) < height / 2.0f)
        && t < INTERSECTION_T1
        && t < INTERSECTION_T2
        && t > 0.0f
        && dot(vel, (float3)( 0.0f, 0.0f, 1.0f )) > 0.f) {

        INTERSECTION_POS1 = pos + t*vel;
        INTERSECTION_POS2 = pos + t*vel;
        INTERSECTION_T1   = t;
        INTERSECTION_T2   = t;

        iidx[global_addr] = comp_idx;
    }

    /* ----------------------- */

    /* Update global memory */
    intersections[global_addr] = intersection;
    neutrons[global_addr]      = neutron;
}