#include "rand.h"
#include "geom.h"
#include "consts.h"

#ifndef GAUSS
#define GAUSS(x,mean,rms) \
  (exp(-((x)-(mean))*((x)-(mean))/(2.0f*(rms)*(rms)))/(sqrt(2.0f*M_PI)*(rms)))
#endif

float randnorm(float16* neutron, uint global_addr) {

  // Samples normal distn by box muller transform

  float u1, u2;

  u1 = rand(neutron, global_addr);
  u2 = rand(neutron, global_addr);

  return sqrt(-2.0f*log(u1))*cos(2.0f*M_PI*u2);
}

__kernel void monochromator(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, float const slab_width, float const slab_height, 
    float const gap, uint const n_horizontal, uint const n_vertical, 
    float const mosaic_horizontal, float const mosaic_vertical,
    float const r0, float const d_spacing, float const radius_vertical, 
    float const radius_horizontal, __global float* gaussx,
    __global float* gaussw, uint const gausslen) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component */
    if (!(this_iidx == comp_idx)) {
        return;
    }

    /* Check termination flag */
    if (neutron.sf > 0.f)  {
        return;
    }

    /* Perform scattering here */

    neutron.s012 = intersection.s456;
    neutron.sa += intersection.s7;

    float mos_rms_y, mos_rms_z, mos_rms_max, Q, tilt_horizontal, tilt_vertical,
          zmin, zmax, ymin, ymax, ratio, Q_order, row, col, q0, q0x, theta, delta,
          p_reflect;
    float3 ki, ku;

    mos_rms_y = MIN2RAD*mosaic_horizontal / sqrt(8.0f*log(2.0f));
    mos_rms_z = MIN2RAD*mosaic_vertical / sqrt(8.0f*log(2.0f));
    mos_rms_max = mos_rms_y > mos_rms_z ? mos_rms_y : mos_rms_z;

    Q = 2.0f * M_PI / d_spacing;

    zmax = ((n_horizontal*(slab_width+gap))-gap)/2.0f;
    zmin = -zmax;
    ymax = ((n_vertical*(slab_height+gap))-gap)/2.0f;
    ymin = -ymax;

    if (neutron.s2 > zmin && neutron.s2 < zmax && neutron.s1 > ymin && neutron.s1 < ymax) {
      col = ceil((neutron.s2 - zmin)/(slab_width + gap));
      row = ceil((neutron.s1 - ymin)/(slab_height + gap));

      tilt_horizontal = radius_horizontal > 0.0f ? asin(((float)col - (float)(n_horizontal+1)/2.0f) * (slab_width + gap) / radius_horizontal) : 0.0f;
      tilt_vertical = radius_vertical > 0.0f ? -asin(((float)row - (float)(n_vertical+1)/2.0f)*(slab_height + gap) / radius_vertical) : 0.0f;

      // Transform to slab frame MAKE SURE TO LEAVE SLAB FRAME WHEN WE'RE DONE

      float center_z = zmin + (col - 0.5f)*(slab_width + gap) - gap / 2.0f;
      float center_y = ymin + (row - 0.5f)*(slab_height + gap) - gap / 2.0f;
      float3 slab_pos = (float3){ 0.0f, center_y, center_z };
      float3 slab_rot = (float3){ 0.0f, tilt_horizontal, tilt_vertical };
      neutron.s345 = frame_derotate(neutron.s345, slab_rot);
      neutron.s012 -= slab_pos;
      neutron.s012 = frame_derotate(neutron.s012, slab_rot);

      if (fabs(neutron.s2) <= slab_width / 2.0f && fabs(neutron.s1) <= slab_height / 2.0f) {
        // Didn't hit the gap

        ki = V2K*neutron.s345;

        // Determine order of Bragg scattering
        ratio = -2.0f*ki.s0 / Q;
        Q_order = floor(ratio + 0.5f);

        if (Q_order == 0.0f) Q_order = ratio < 0.0f ? -1.0f : 1.0f;
        if (Q_order < 0.0f) Q_order = -Q_order;

        // Ensure scattering is possible

        ku = normalize(ki);

        if (Q_order > 2.0f*length(ki) / Q) Q_order-= 1.0f;

        // TODO: add higher order restriction here

        q0 = Q_order*Q;
        q0x = ratio < 0.0f ? -q0 : q0;
        theta = asin(q0 / (2.0f*length(ki)));

        delta = asin(fabs(ku.s0)) - theta;

        p_reflect = fabs(r0) * exp(-ki.s2*ki.s2/(length(ki.s12)*length(ki.s12)) * delta*delta/(2.0f*mos_rms_y*mos_rms_y)) * exp(-ki.s1*ki.s1 / (length(ki.s12)*length(ki.s12)) * delta*delta/(2.0f*mos_rms_z*mos_rms_z));

        if (rand(&neutron, global_addr) <= p_reflect) { // reflect
          float3 b, a, c1, q;
          float phi, cos_2theta, k_sin_2theta, cos_phi, sin_phi;
          float total, w, mos_sample;

          cos_2theta = cos(2.0f*theta);
          k_sin_2theta = length(ki)*sin(2.0f*theta);

          b = normalize(cross(ki, (float3){ q0x, 0.0f, 0.0f }));
          b *= k_sin_2theta;

          a = cross(b, ku);
          total = 0.0f;

          mos_sample = mos_rms_max / cos(theta);
          c1 = ki*(cos_2theta - 1.0f);

          for(uint i = 0; i<100; i++) {
            w = 5.0f*mos_sample;
            cos_phi = cos(w);
            sin_phi = sin(w);
            q = c1 + cos_phi*a + sin_phi*b;
            q.s1 = q.s1 / mos_rms_z;
            q.s2 = q.s2 / mos_rms_y;

            if (q.s2*q.s2 + q.s1*q.s1 < (25.0f/(2.0f/3.0f))*(q.s0*q.s0))
              break;
            mos_sample *= (2.0f / 3.0f);
          }

          for(uint i = 0; i < gausslen; i++) {
            phi = w*gaussx[i];
            cos_phi = cos(phi);
            sin_phi = sin(phi);

            q = c1 + cos_phi*a + sin_phi*b;

            p_reflect = GAUSS((q.s2/q.s0),0.0f,mos_rms_y)*GAUSS((q.s1/q.s0),0.0f,mos_rms_z);
            total += gaussw[i]*p_reflect;
          }

          total *= w;

          phi = mos_sample*randnorm(&neutron, global_addr);

          cos_phi = cos(phi);
          sin_phi = sin(phi);

          q = c1 + cos_phi*a + sin_phi*b;

          p_reflect = GAUSS((q.s2/q.s0),0.0f,mos_rms_y)*GAUSS((q.s1/q.s0),0.0f,mos_rms_z);

          neutron.s345 = K2V*(ki + q);

          p_reflect /= total*GAUSS(phi, 0.0f, mos_sample);
          // if (p_reflect <= 0) // absorb
          if (p_reflect > 1.0f) p_reflect = 1.0f;

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
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}