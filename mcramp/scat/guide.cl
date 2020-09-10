/*
* The functionality of this component has been entirely ripped
* (to the extent of copy paste) from the McStas guide component.
*/

#include "ref.h"
#include "consts.h"

__kernel void guide_scatter(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx,
    real_t const w1, real_t const h1,
    real_t const w2, real_t const h2, real_t const l,
    real_t const R0, real_t const Qc, real_t const alpha,
    real_t const m, real_t const W, uint const max_bounces) {

    uint global_addr        = get_global_id(0);
    real16_t neutron         = convert_real16_t(neutrons[global_addr]);
    real8_t intersection = convert_real8_t(intersections[global_addr]);
    
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component */
    if (!(this_iidx == comp_idx)) {
        return;
    }

    /* Check termination flag */
    if (neutron.sf > 0.)  {
        return;
    }

    /* Perform scattering here */

    // Propogate neutron to guide entrance first

    NEUTRON_POS = INTERSECTION_POS2;
    NEUTRON_TOF += INTERSECTION_T2;

    real3_t pos, vel;
    real_t t1, t2, av, ah, bv, bh, cv1, cv2, ch1, ch2, d, ww, hh, whalf, hhalf;
    real_t vdotn_v1, vdotn_v2, vdotn_h1, vdotn_h2, q, nlen2, refl;

    uint i = 0;

    ww = 0.5*(w2 - w1); hh = 0.5*(h2 - h1);
    whalf = 0.5*w1; hhalf = 0.5*h1;

    while(true && i < max_bounces) {
        pos = NEUTRON_POS;
        vel = NEUTRON_VEL;

        av = l*vel.s0; bv = ww*vel.s2;
        ah = l*vel.s1; bh = hh*vel.s2;
        vdotn_v1 = bv + av;
        vdotn_v2 = bv - av;
        vdotn_h1 = bh + ah;
        vdotn_h2 = bh - ah;

        cv1 = -whalf*l - pos.s2*ww;
        cv2 = pos.s0*l;
        ch1 = -hhalf*l - pos.s2*hh;
        ch2 = pos.s1*l;

        t1 = (l - pos.s2)/vel.s2;
        i = 0;
        if(vdotn_v1 < 0.0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
        {
          t1 = t2;
          i = 1;
        }
        if(vdotn_v2 < 0.0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
        {
          t1 = t2;
          i = 2;
        }
        if(vdotn_h1 < 0.0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
        {
          t1 = t2;
          i = 3;
        }
        if(vdotn_h2 < 0.0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
        {
          t1 = t2;
          i = 4;
        }

        if (i == 0) {
            break;
        }

        NEUTRON_POS += t1*vel;
        NEUTRON_TOF += t1;

        switch(i)
        {
          case 1:                   /* Left vertical mirror */
            nlen2 = l*l + ww*ww;
            q = V2K*(-2.0)*vdotn_v1/sqrt(nlen2);
            d = 2.0*vdotn_v1/nlen2;
            NEUTRON_VX = NEUTRON_VX - d*l;
            NEUTRON_VZ = NEUTRON_VZ - d*ww;
            break;
          case 2:                   /* Right vertical mirror */
            nlen2 = l*l + ww*ww;
            q = V2K*(-2.0)*vdotn_v2/sqrt(nlen2);
            d = 2.0*vdotn_v2/nlen2;
            NEUTRON_VX = NEUTRON_VX + d*l;
            NEUTRON_VZ = NEUTRON_VZ - d*ww;
            break;
          case 3:                   /* Lower horizontal mirror */
            nlen2 = l*l + hh*hh;
            q = V2K*(-2.0)*vdotn_h1/sqrt(nlen2);
            d = 2.0*vdotn_h1/nlen2;
            NEUTRON_VY = NEUTRON_VY - d*l;
            NEUTRON_VZ = NEUTRON_VZ - d*hh;
            break;
          case 4:                   /* Upper horizontal mirror */
            nlen2 = l*l + hh*hh;
            q = V2K*(-2.0)*vdotn_h2/sqrt(nlen2);
            d = 2.0*vdotn_h2/nlen2;
            NEUTRON_VY = NEUTRON_VY + d*l;
            NEUTRON_VZ = NEUTRON_VZ - d*hh;
            break;
        }

        refl = reflectivity_func(q, R0, Qc, alpha, m, W);
        NEUTRON_P *= refl;

        // Delete extremely small reflections to avoid numerical issues
        if (refl < 1e-12 || refl > 1.0f)
          NEUTRON_P = 0.0;

        i++;
    }

    if (i >= max_bounces)
      NEUTRON_DIE = 1.0;

    /* ----------------------- */

    /* Update global memory and reset intersection */
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;

    neutrons[global_addr]      = convert_float16(neutron);
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}