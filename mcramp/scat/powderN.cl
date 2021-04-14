#include "rand.h"
#include "geom.h"
#include "consts.h"

__kernel void powderN(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, __global float* sigma_scat_v2, float const sigma_abs_v,
    __global float* q, float const d_phi, uint const num_lines,
    uint const transmit) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection     = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component -------------- */
    if (!(this_iidx == comp_idx))
        return;

    /* Check termination flag ---------------------------------------------- */
    if (NEUTRON_DIE  > 0.f) 
        return;

    /* Perform scattering here --------------------------------------------- */

    float vel, sigma_scat, sigma_abs, p_scat, p_inter, full_path_length, phi,
            arg, twotheta, pen_depth, sum_sigma_scat, u, v, q_chosen, deltaT;
    float3 beam_para, beam_perp, path;
    uint i, attempts;
    // bool scattered = false; 

    // NEUTRON_TOF += intersection.s3;

    if (all(INTERSECTION_POS1 == INTERSECTION_POS2)) {
        path = INTERSECTION_POS2 - NEUTRON_POS;
    } else {
        path = INTERSECTION_POS2 - INTERSECTION_POS1;
    }

    full_path_length = length(path);
    vel = length(neutron.s345);
    sigma_abs = sigma_abs_v * 2200.0 / vel;

    // Calculate scattering probabilites

    sum_sigma_scat = 0.0f;
    for(i=0;i<num_lines;i++) {
        sigma_scat = sigma_scat_v2[i] / (vel * vel * V2K * V2K);
        sum_sigma_scat += sigma_scat;
    }

    p_scat = sum_sigma_scat / (sigma_abs + sum_sigma_scat);
    p_inter = 1 - exp(-(sum_sigma_scat + sigma_abs)*full_path_length);

    if (rand(&neutron, global_addr) > p_inter) {
        // Transmitted, return without modifying
        // neutron state, but multiply by weight factor
        NEUTRON_P *= 1 - p_inter;
        NEUTRON_POS= (intersection.s456+0.0001f*normalize(path));
        NEUTRON_TOF += intersection.s7;

        if (transmit == 0)
          NEUTRON_DIE  = 1.0f;

        neutron.sc = comp_idx;
        iidx[global_addr] = 0;
        neutrons[global_addr] = neutron;
        intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                             0.0f, 0.0f, 0.0f, 100000.0f );
        return;
    }

    // Select scattering line

    attempts = 0;

    while(true) {
        u = rand(&neutron, global_addr);
        v = 0.0f;
        i = 0;
        while(v < u) {
            sigma_scat = sigma_scat_v2[i] / (vel * vel * V2K * V2K);
            q_chosen = q[i];
            v += sigma_scat / sum_sigma_scat;

            i++;
        }

        arg = q_chosen * K2V / (2.0f * vel);

        if(arg < 1.0f) {
            break;
        } else if (arg > 1.0f && attempts > 2*num_lines) { // No bragg reflection possible, transmit
            NEUTRON_DIE  = 1;
            iidx[global_addr] = 0;
            neutrons[global_addr] = neutron;
            intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
            return;
        }

        attempts++;
    }

    twotheta = 2.0f*asin(arg);

    pen_depth = -log(1.0f - rand(&neutron, global_addr)*(1.0f - p_inter)) / (sigma_abs + sum_sigma_scat);

    beam_para = normalize(neutron.s345);
    beam_perp = cross(beam_para, (float3)(0.0f, 1.0f, 0.0f));

    if(length(beam_perp) < 1e-3f) { // Beam just happens to be along 0 1 0
        beam_perp = normalize(cross(beam_para, (float3)(1.0f, 0.0f, 0.0f)));
    } else {
        beam_perp = normalize(beam_perp);
    }

    if (d_phi)
      { 
        phi = (2.0f*rand(&neutron, global_addr) - 1.0f)*(d_phi * M_PI / 180.0f);
        if (rand(&neutron, global_addr) > 0.5f) {
            // go to other side of beam
            phi += M_PI;
        }

        phi += M_PI / 2.0f;
      }
      else
        phi = 2.0f * M_PI * rand(&neutron, global_addr);

    rotate_about_axis(phi, beam_para, &beam_perp);

    rotate_about_axis(twotheta, beam_perp, &beam_para);

    if (all(INTERSECTION_POS1 == intersection.s456)) {
        NEUTRON_POS += pen_depth * path;
        NEUTRON_TOF += length(pen_depth * path) / length(neutron.s345);
    } else {
        NEUTRON_POS= INTERSECTION_POS1 + pen_depth * path;   
        NEUTRON_TOF += INTERSECTION_T2 - intersection.s3;
    }
    
    NEUTRON_VEL = vel*normalize(beam_para);

    NEUTRON_P *= (1 - exp(-(sum_sigma_scat + sigma_abs)*pen_depth)) * sum_sigma_scat / (sum_sigma_scat + sigma_abs);

    /* ----------------------- */

    /* Update global memory and reset intersection */
    
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;
    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}