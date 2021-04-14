#include "consts.h"

__kernel void intersect(__global float16* neutrons, 
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx, float const radius, float const height) {

  uint global_addr = get_global_id(0);

  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float2 plane_pos, plane_vel;
  float theta2, a, b, c, t1, t2, quotient;

  if (NEUTRON_DIE  > 0.f) {
    return;
  }


  plane_pos = neutron.s02;
  plane_vel = neutron.s35;

  a = dot(plane_vel, plane_vel);
  b = 2.0f*dot(plane_pos, plane_vel);
  c = dot(plane_pos, plane_pos) - radius*radius;

  quotient = b*b-4.0f*a*c;

  t1 = (-b - sqrt(quotient)) / (2.f*a);
  t2 = (-b + sqrt(quotient)) / (2.f*a);

  intersection = intersections[global_addr];

  theta2 = acos(dot(normalize(plane_pos+t2*plane_vel), (float2)( 0.0f, 1.0f )));
  if ((quotient > 0.0f) &&
      ((-height/2.0f) < (NEUTRON_Y+ t2*NEUTRON_VY)) &&
      ((NEUTRON_Y+ t2*NEUTRON_VY) < (height/2))) {
  
    if (t1 < INTERSECTION_T1 && t1 > 0.0f && t2 > 0.0f) {
        
        INTERSECTION_POS1 = NEUTRON_POS+ t1*NEUTRON_VEL;
        INTERSECTION_T1   = t1;
      
        iidx[global_addr] = comp_idx;
        

        INTERSECTION_POS2 = NEUTRON_POS+ t2*NEUTRON_VEL;
        INTERSECTION_T2   = t2;
    }
  }

  neutrons[global_addr] = neutron;
  intersections[global_addr] = intersection;
}
