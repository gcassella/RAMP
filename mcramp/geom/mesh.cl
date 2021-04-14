#include "consts.h"

__kernel void intersect(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, __global float* points,
    uint const num_tri, float const x1, float const x2,
    float const y1, float const y2, float const z1,
    float const z2, uint const interior) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection     = intersections[global_addr];

    /* Check termination flag */
    if (NEUTRON_DIE  > 0.f) 
        return;

    /* Perform raytracing here */

    // Raytrace against bounding box

    float3 box_ll, box_ur, tA, tB;
    float t1, t2, t3, t4, t5, t6, tmin, tmax;
    uint num_intersections = 0;

    box_ll = (float3){ x1, y1, z1 };
    box_ur = (float3){ x2, y2, z2 };
    tA = (box_ll - NEUTRON_POS) / NEUTRON_VEL;
    tB = (box_ur - NEUTRON_POS) / NEUTRON_VEL;

    t1 = tA.s0;
    t2 = tB.s0;
    t3 = tA.s1;
    t4 = tB.s1;
    t5 = tA.s2;
    t6 = tB.s2;
    tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    if (tmax < 0.0f || tmin > tmax) {
        return; // Missed bounding box
    }

    // Raytrace against triangles

    float3 v0, v1, v2, e1, e2, h, s, q, norm;
    const float eps = 0.000000001f;
    float a, f, u, v, t;

    float8 out = (float8){ 0.0f, 0.0f, 0.0f, 1e8f, 0.0f, 0.0f, 0.0f, 1e8f };

    for(uint i=0;i<num_tri;i++) {
        // Loop over triangles

        v0 = (float3){ 
            points[i*9+0],
            points[i*9+1],
            points[i*9+2]
        };
        v1 = (float3){ 
            points[i*9+3],
            points[i*9+4],
            points[i*9+5]
        };
        v2 = (float3){ 
            points[i*9+6],
            points[i*9+7],
            points[i*9+8]
        };

        e1 = v1 - v0;
        e2 = v2 - v0;

        h = cross(NEUTRON_VEL, e2);
        a = dot(e1, h);

        if (fabs(a) < eps)
            continue;
        
        f = 1.0f / a;
        s = NEUTRON_POS - v0;
        u = f*dot(s, h);

        if (u < 0.0 || u > 1.0) {
            continue;
        }

        q = cross(s, e1);
        v = f*dot(NEUTRON_VEL, q);

        if (v < 0.0 || u + v > 1.0) {
            continue;
        }

        t = f * dot(e2, q);

        if (t > eps) {
            num_intersections++;
        }

        norm = cross(e1, e2);

        if ((t > eps) && (t < out.s3) && (dot(norm, NEUTRON_VEL) < 0)) {
            out.s012 = NEUTRON_POS + NEUTRON_VEL*t;
            out.s3 = t;
        } 
        
        if ((t > eps) && (t < out.s7) && (dot(norm, NEUTRON_VEL) > 0)) {
            out.s456 = NEUTRON_POS + NEUTRON_VEL*t;
            out.s7 = t;
        }
    }

    if ((interior == 1) && (num_intersections%2 == 1)) {
        out.s4567 = out.s0123;
    } else if ((interior == 1) && (num_intersections%2 == 0)) {
        return;
    }

    if(out.s3 < INTERSECTION_T1) {
        intersection = out;  
        iidx[global_addr] = comp_idx;
    }

    /* ----------------------- */

    /* Update global memory */
    intersections[global_addr] = intersection;
    neutrons[global_addr]      = neutron;
}