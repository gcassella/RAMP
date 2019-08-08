#include "rand.h"
#include "geom.h"

__kernel void powder_scatter(__global double16* neutrons,
  __global double8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global double3* reflections, uint const num_reflections,
  double const sigma_abs, double const Vc) {

  // Choose a reflection
  uint global_addr = get_global_id(0);
  double16 neutron = neutrons[global_addr];
  double8 intersection = intersections[global_addr];
  uint this_iidx = iidx[global_addr];

  if(!(this_iidx == comp_idx)) {
    return;
  }

  if (neutron.sf == 1.) {
    return;
  }

  uint NReflection, attempts, scattered;
  double3 reflection, perp, normvel, path;
  double x, y, z, Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, arg, q_p;
  double vel, q_v, alpha, sigma_s, sigma_a, sigma_tot, mu, ki;



  normvel = normalize(neutron.s345);
  vel = length(neutron.s345);

  attempts = 0;
  scattered = 0;

  path = intersection.s456 - intersection.s012;
  ki = 1.583*pow(10.,-3.)*length(neutron.s345);

  normvel = normalize(neutron.s345);

  // Monte carlo choice to see if our neutron scatters

  neutron.s012 = intersection.s012 + rand(&neutron, global_addr)*(intersection.s456 - intersection.s012);
  neutron.sa += intersection.s3;

  do {
    NReflection = floor(rand(&neutron, global_addr)*num_reflections);
    
    reflection = reflections[NReflection];
    // Check if reflection is okay?
    q_v = 3.97*pow(10., 3.) / reflection.s0;
    q_p = 1.58*pow(10., -3.) * q_v;

    arg = q_v / (2*vel);

    sigma_a = sigma_abs * 2200 * 100 / Vc / length(neutron.s345);
    sigma_s = 100*4*M_PI*M_PI*M_PI*reflection.s2*10*reflection.s1*reflection.s1 / (2*Vc*Vc*ki*ki*q_p);
    sigma_tot = sigma_a + sigma_s;

    mu = sigma_tot;
    
    if (arg < 1.0) {
      // Okay, scatter
      // Generate an arbitrary perpendicular vector to the velocity
      // x*vx + y*vy + z*vz = (vx+vy)+z*vz = 0
      // => z = -(vx+vy)/vz

      if (rand(&neutron, global_addr) < exp(-mu*length(path.s012))) {
        // Transmitted and attenuated
        neutron.s9*=exp(-mu*length(path.s012));


        neutron.sc = comp_idx;
        iidx[global_addr] = 0;
        neutrons[global_addr] = neutron;
        intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                             0.0f, 0.0f, 0.0f, 100000.0f );
        neutrons[global_addr] = neutron;
        return;
      }

      x = 1.;
      y = 1.;
      z = -(neutron.s3+neutron.s4)/(neutron.s5);
      perp = normalize((double3)( x, y, z ));
      // construct rotation matrix to randomly rotate the scattering
      // vector about the velocity

      alpha = 2*M_PI*rand(&neutron, global_addr);
      
      rotate_about_axis(alpha, normvel, (&perp));

      alpha = 2*asin(arg);

      rotate_about_axis(alpha, perp, (&normvel));

      neutron.s345 = vel*normvel;
      neutron.s9 *= length(path.s012)*sigma_s*exp(-sigma_a*length(path.s012));
      neutron.se = this_iidx;
      break;

    } else {
      // Not allowed, try again
      attempts += 1;
    }
  } while (attempts < 100);

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}