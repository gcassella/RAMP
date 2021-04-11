#include "rand.h"
#include "geom.h"
#include "consts.h"

__kernel void powder1(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, float const sigma_scat_v2, float const sigma_abs_v,
    float const q, float const d_phi) {

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
            arg, twotheta, pen_depth;
    float3 beam_para, beam_perp, path;
    // bool scattered = false; 

    NEUTRON_TOF += intersection.s3;

    path = INTERSECTION_POS2 - intersection.s012;

    full_path_length = length(path);
    vel = length(neutron.s345);
    sigma_scat = sigma_scat_v2 / (vel * vel * V2K * V2K);
    sigma_abs = sigma_abs_v * 2200.0 / vel;

    p_scat = sigma_scat / (sigma_abs + sigma_scat);
    p_inter = 1 - exp(-(sigma_scat + sigma_abs)*full_path_length);

    arg = q * K2V / (2.0f * vel);
    if (arg > 1.0f) { // No bragg reflection possible
        
        iidx[global_addr] = 0;
        neutron.sc = comp_idx;
        neutrons[global_addr]      = neutron;
        intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
        return;
    }

    twotheta = 2.0f*asin(arg);

    NEUTRON_P *= p_scat; // Adjust weight to model absorption

    pen_depth = -log(1.0f - rand(&neutron, global_addr)*(1.0f - p_inter)) / (sigma_abs + sigma_scat);

    beam_para = normalize(NEUTRON_VEL);
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

    NEUTRON_POS= INTERSECTION_POS1 + pen_depth * path;
    NEUTRON_VEL = vel*normalize(beam_para);

    NEUTRON_TOF += INTERSECTION_T2 - intersection.s3;
    NEUTRON_P *= full_path_length*sigma_scat*exp(-(sigma_scat+sigma_abs)*full_path_length);

    /* ----------------------- */

    /* Update global memory and reset intersection */
    
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;
    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}