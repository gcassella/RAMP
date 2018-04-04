#include "rand.h"

__kernel void powder_scatter(__global float16* neutrons,
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global float3* reflections, uint const num_reflections) {

  // Choose a reflection
  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];
  uint this_iidx = iidx[global_addr];

  if(!(this_iidx == comp_idx)) {
    return;
  }

  if (neutron.sf == 1.) {
    return;
  }

  uint NReflection, attempts;
  float3 reflection, perp, normvel;
  float x, y, z, Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, arg;
  float vel, q_v, alpha;

  normvel = normalize(neutron.s345);
  vel = length(neutron.s345);

  attempts = 0;
  
  if (length(intersection.s0123) != 0.) {
    neutron.s012a = (intersection.s0123 + intersection.s4567) / 2;

    do {
      NReflection = floor(rand(&neutron, global_addr)*num_reflections);
      

      reflection = reflections[NReflection];

      // Check if reflection is okay?
      q_v = 3.97*pow(10., 3.) / reflection.s1;

      arg = q_v / (2*vel);

      if (arg < 1.0) {
        // Okay, scatter
        // Generate an arbitrary perpendicular vector to the velocity
        // x*vx + y*vy + z*vz = (vx+vy)+z*vz = 0
        // => z = -(vx+vy)/vz

        x = 1.;
        y = 1.;
        z = -(neutron.s3+neutron.s4)/(neutron.s5);
        perp = normalize((float3)( x, y, z ));
        // construct rotation matrix to randomly rotate the scattering
        // vector about the velocity

        alpha = 2*M_PI*rand(&neutron, global_addr);
        

        Rxx = cos(alpha)+normvel.s0*normvel.s0*(1-cos(alpha));
        Rxy = normvel.s0*normvel.s1*(1-cos(alpha))-normvel.s2*sin(alpha);
        Rxz = normvel.s0*normvel.s2*(1-cos(alpha))+normvel.s1*sin(alpha);
        Ryx = normvel.s1*normvel.s0*(1-cos(alpha))+normvel.s2*sin(alpha);
        Ryy = cos(alpha)+normvel.s1*normvel.s1*(1-cos(alpha));
        Ryz = normvel.s1*normvel.s2*(1-cos(alpha))-normvel.s0*sin(alpha);
        Rzx = normvel.s2*normvel.s0*(1-cos(alpha))-normvel.s1*sin(alpha);
        Rzy = normvel.s2*normvel.s1*(1-cos(alpha))+normvel.s0*sin(alpha);
        Rzz = cos(alpha)+normvel.s2*normvel.s2*(1-cos(alpha));

        perp = (float3)( Rxx*perp.s0+Rxy*perp.s1+Rxz*perp.s2, 
                         Ryx*perp.s0+Ryy*perp.s1+Ryz*perp.s2, 
                         Rzx*perp.s0+Rzy*perp.s1+Rzz*perp.s2 );

        // now construct rotation matrix to rotate the velocity 2theta
        // about the scattering vector
        alpha = 2*asin(arg);

        Rxx = cos(alpha)+perp.s0*perp.s0*(1-cos(alpha));
        Rxy = perp.s0*perp.s1*(1-cos(alpha))-perp.s2*sin(alpha);
        Rxz = perp.s0*perp.s2*(1-cos(alpha))+perp.s1*sin(alpha);
        Ryx = perp.s1*perp.s0*(1-cos(alpha))+perp.s2*sin(alpha);
        Ryy = cos(alpha)+perp.s1*perp.s1*(1-cos(alpha));
        Ryz = perp.s1*perp.s2*(1-cos(alpha))-perp.s0*sin(alpha);
        Rzx = perp.s2*perp.s0*(1-cos(alpha))-perp.s1*sin(alpha);
        Rzy = perp.s2*perp.s1*(1-cos(alpha))+perp.s0*sin(alpha);
        Rzz = cos(alpha)+perp.s2*perp.s2*(1-cos(alpha));

        normvel = (float3)( Rxx*normvel.s0+Rxy*normvel.s1+Rxz*normvel.s2, 
                            Ryx*normvel.s0+Ryy*normvel.s1+Ryz*normvel.s2, 
                            Rzx*normvel.s0+Rzy*normvel.s1+Rzz*normvel.s2 );

        neutron.s345 = vel*normvel;
        neutron.s9 *= reflection.s2;

        break;

      } else {
        // Not allowed, try again
        attempts += 1;
      }
    } while (attempts < 100);
  }

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 0.0f );
}