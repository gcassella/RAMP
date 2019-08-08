/*
* The functionality of this component has been entirely ripped
* (to the extent of copy paste) from the McStas guide component.
*/

#include "ref.h"
#include "consts.h"

__kernel void guide_scatter(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx,
    double const w1, double const h1,
    double const w2, double const h2, double const l,
    double const R0, double const Qc, double const alpha,
    double const m, double const W, uint const max_bounces) {

    uint global_addr        = get_global_id(0);
    double16 neutron         = neutrons[global_addr];
    double8 intersection = intersections[global_addr];
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

    neutron.s012 = intersection.s456;
    neutron.sa += intersection.s7;

    double3 pos, vel;
    double t1, t2, av, ah, bv, bh, cv1, cv2, ch1, ch2, d, ww, hh, whalf, hhalf;
    double vdotn_v1, vdotn_v2, vdotn_h1, vdotn_h2, q, nlen2;

    uint i = 0;

    ww = 0.5*(w2 - w1); hh = 0.5*(h2 - h1);
    whalf = 0.5*w1; hhalf = 0.5*h1;

    while(true) {
        pos = neutron.s012;
        vel = neutron.s345;

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
        if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
        {
          t1 = t2;
          i = 1;
        }
        if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
        {
          t1 = t2;
          i = 2;
        }
        if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
        {
          t1 = t2;
          i = 3;
        }
        if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
        {
          t1 = t2;
          i = 4;
        }

        if (i == 0) {
            break;
        }

        neutron.s012 += t1*vel;
        neutron.sa += t1;

        switch(i)
        {
          case 1:                   /* Left vertical mirror */
            nlen2 = l*l + ww*ww;
            q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
            d = 2*vdotn_v1/nlen2;
            neutron.s3 = neutron.s3 - d*l;
            neutron.s5 = neutron.s5 - d*ww;
            break;
          case 2:                   /* Right vertical mirror */
            nlen2 = l*l + ww*ww;
            q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
            d = 2*vdotn_v2/nlen2;
            neutron.s3 = neutron.s3 + d*l;
            neutron.s5 = neutron.s5 - d*ww;
            break;
          case 3:                   /* Lower horizontal mirror */
            nlen2 = l*l + hh*hh;
            q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
            d = 2*vdotn_h1/nlen2;
            neutron.s4 = neutron.s4 - d*l;
            neutron.s5 = neutron.s5 - d*hh;
            break;
          case 4:                   /* Upper horizontal mirror */
            nlen2 = l*l + hh*hh;
            q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
            d = 2*vdotn_h2/nlen2;
            neutron.s4 = neutron.s4 + d*l;
            neutron.s5 = neutron.s5 - d*hh;
            break;
        }

        neutron.s9  *= reflectivity_func(q, R0, Qc, alpha, m, W);
    }


    /* ----------------------- */

    /* Update global memory and reset intersection */
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;

    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}