#include "rand.h"
#include "geom.h"
#include "consts.h"

#ifndef GAUSS
#define GAUSS(x,mean,rms) \
  (exp(-((x)-(mean))*((x)-(mean))/(2*(rms)*(rms)))/(sqrt(2*M_PI)*(rms)))
#endif

double randnorm(double16* neutron, uint global_addr) {

  // Samples normal distn by box muller transform

  double u1, u2;

  u1 = rand(neutron, global_addr);
  u2 = rand(neutron, global_addr);

  return sqrt(-2*log(u1))*cos(2*M_PI*u2);
}

__kernel void monochromator(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx, double const slab_width, double const slab_height, 
    double const gap, uint const n_horizontal, uint const n_vertical, 
    double const mosaic_horizontal, double const mosaic_vertical,
    double const r0, double const d_spacing, double const radius_vertical, 
    double const radius_horizontal, __global double* gaussx,
    __global double* gaussw, uint const gausslen) {

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

    neutron.s012 = intersection.s456;
    neutron.sa += intersection.s7;

    double mos_rms_y, mos_rms_z, mos_rms_max, Q, tilt_horizontal, tilt_vertical,
          zmin, zmax, ymin, ymax, ratio, Q_order, row, col, q0, q0x, theta, delta,
          p_reflect;
    double3 ki, ku;

    mos_rms_y = MIN2RAD*mosaic_horizontal / sqrt(8.0*log(2.0));
    mos_rms_z = MIN2RAD*mosaic_vertical / sqrt(8.0*log(2.0));
    mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

    Q = 2 * M_PI / d_spacing;

    zmax = ((n_horizontal*(slab_width+gap))-gap)/2.0;
    zmin = -zmax;
    ymax = ((n_vertical*(slab_height+gap))-gap)/2.0;
    ymin = -ymax;

    if (neutron.s2 > zmin && neutron.s2 < zmax && neutron.s1 > ymin && neutron.s1 < ymax) {
      col = ceil((neutron.s2 - zmin)/(slab_width + gap));
      row = ceil((neutron.s1 - ymin)/(slab_height + gap));

      tilt_horizontal = radius_horizontal ? asin((col - (double)(n_horizontal+1)/2.0) * (slab_width + gap) / radius_horizontal) : 0.0;
      tilt_vertical = radius_vertical ? -asin((row - (double)(n_vertical+1)/2.0)*(slab_height + gap) / radius_vertical) : 0.0;

      // Transform to slab frame MAKE SURE TO LEAVE SLAB FRAME WHEN WE'RE DONE

      double center_z = zmin + (col - 0.5)*(slab_width + gap) - gap / 2.0;
      double center_y = ymin + (row - 0.5)*(slab_height + gap) - gap / 2.0;
      double3 slab_pos = (double3){ 0.0, center_y, center_z };
      double3 slab_rot = (double3){ 0.0, tilt_horizontal, tilt_vertical };
      neutron.s345 = frame_derotate(neutron.s345, slab_rot);
      neutron.s012 -= slab_pos;
      neutron.s012 = frame_derotate(neutron.s012, slab_rot);

      if (fabs(neutron.s2) <= slab_width / 2.0 && fabs(neutron.s1) <= slab_height / 2.0) {
        // Didn't hit the gap

        ki = V2K*neutron.s345;

        // Determine order of Bragg scattering
        ratio = -2*ki.s0 / Q;
        Q_order = floor(ratio + 0.5);

        if (Q_order == 0.0) Q_order = ratio < 0 ? -1 : 1;
        if (Q_order < 0) Q_order = -Q_order;

        // Ensure scattering is possible

        ku = normalize(ki);

        if (Q_order > 2*length(ki) / Q) Q_order-= 1;

        // TODO: add higher order restriction here

        q0 = Q_order*Q;
        q0x = ratio < 0 ? -q0 : q0;
        theta = asin(q0 / (2*length(ki)));

        delta = asin(fabs(ku.s0)) - theta;

        p_reflect = fabs(r0) * exp(-ki.s2*ki.s2/length(ki.s12) * delta*delta/(2*mos_rms_y*mos_rms_y)) * exp(-ki.s1*ki.s1 / length(ki.s12) * delta*delta/(2*mos_rms_z*mos_rms_z));

        if (rand(&neutron, global_addr) <= p_reflect) { // reflect
          double3 b, a, c1, q;
          double phi, cos_2theta, k_sin_2theta, cos_phi, sin_phi;
          double total, w, mos_sample;

          cos_2theta = cos(2*theta);
          k_sin_2theta = length(ki)*sin(2*theta);

          b = normalize(cross(ki, (double3){ q0x, 0.0, 0.0 }));
          b *= k_sin_2theta;

          a = cross(b, ku);
          total = 0;

          mos_sample = mos_rms_max / cos(theta);
          c1 = ki*(cos_2theta - 1.0);

          for(uint i = 0; i<100; i++) {
            w = 5*mos_sample;
            cos_phi = cos(w);
            sin_phi = sin(w);
            q = c1 + cos_phi*a + sin_phi*b;
            q.s1 = q.s1 / mos_rms_z;
            q.s2 = q.s2 / mos_rms_y;

            if (q.s2*q.s2 + q.s1*q.s1 < (25.0/(2.0/3.0))*(q.s0*q.s0))
              break;
            mos_sample *= (2.0 / 3.0);
          }

          for(uint i = 0; i < gausslen; i++) {
            phi = w*gaussx[i];
            cos_phi = cos(phi);
            sin_phi = sin(phi);

            q = c1 + cos_phi*a + sin_phi*b;

            p_reflect = GAUSS((q.s2/q.s0),0,mos_rms_y)*GAUSS((q.s1/q.s0),0,mos_rms_z);
            total += gaussw[i]*p_reflect;
          }

          total *= w;

          phi = mos_sample*randnorm(&neutron, global_addr);

          cos_phi = cos(phi);
          sin_phi = sin(phi);

          q = c1 + cos_phi*a + sin_phi*b;

          p_reflect = GAUSS((q.s2/q.s0),0,mos_rms_y)*GAUSS((q.s1/q.s0),0,mos_rms_z);

          neutron.s345 = K2V*(ki + q);

          p_reflect /= total*GAUSS(phi, 0, mos_sample);
          // if (p_reflect <= 0) // absorb
          if (p_reflect > 1) p_reflect = 1;

          neutron.s9 *= p_reflect;
        } else {
          // transmit neutron
        }
      } else {
        // hit gap, transmit neutron
      }

      // Time to leave slab frame

      neutron.s012 = frame_rotate(neutron.s012, slab_rot);
      neutron.s012 += slab_pos;
      neutron.s345 = frame_rotate(neutron.s345, slab_rot);
    } else {
      // missed mono, should never happen in ramp
    }

    /* ----------------------- */

    /* Update global memory and reset intersection */
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;

    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}