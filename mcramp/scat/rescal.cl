#include "rand.h"
#include "geom.h"
#include "consts.h"

__kernel void rescal(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, float3 const target, float const E0,
    float const dE, float const focus_r, __global float8* eventlist) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection     = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component -------------- */
    if (!(this_iidx == comp_idx))
        return;

    /* Check termination flag ---------------------------------------------- */
    if (neutron.sf > 0.0f) 
        return;

    /* Perform scattering here --------------------------------------------- */

    // Choose random scattering point inside sample

    float3 path = intersection.s456 - intersection.s012;
    neutron.s012 = intersection.s012 + rand(&neutron, global_addr)*path;

    float8 ev;
    float3 target_perp, target_perp_perp, target_final, ki, kf;
    float l2, costhetamax, theta, phi, vf, Ef;

    ki = V2K*neutron.s345;

    l2 = length(target)*length(target);
    costhetamax = sqrt(l2 / (focus_r*focus_r + l2));

    theta = acos(1.0f - rand(&neutron, global_addr) * (1.0f - costhetamax));
    phi = rand(&neutron, global_addr)*2.0f*M_PI;

    if (target.s0 == 0.0f && target.s2 == 0.0f) {
        target_perp = (float3){ 1.0f, 0.0f, 0.0f };
    } else {
        target_perp = (float3){ -target.s2, 0.0f, target.s0 };
    }

    target_perp_perp = cross(target, target_perp);

    target_final = target;
    rotate_about_axis(theta, target_perp_perp, &target_final);
    rotate_about_axis(phi, target, &target_final);
    
    target_final = normalize(target_final);

    Ef = E0 + dE*(1.0f - 2.0f*rand(&neutron, global_addr));
    vf = SE2V*sqrt(Ef);

    neutron.s345 = target_final*vf;
    kf = V2K*neutron.s345;

    ev.s012 = ki - kf;
    ev.s6 = neutron.s9;
    ev.s7 = VS2E*pow(length(K2V*ki), 2.0f) - Ef;

    eventlist[global_addr] = ev;
    /* ----------------------- */

    /* Update global memory and reset intersection */
    
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;
    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}