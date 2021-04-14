#include "consts.h"

__kernel void intersect_sphere(__global float16* neutrons, 
      __global float8* intersections, __global uint* iidx,
      uint const comp_idx, float const sphere_radius) {

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float3 pos, vel;
  float s_rad, a, b, c, d1, d2, quotient;
 
  if (NEUTRON_DIE  > 0.f) {
    return;
  }

  pos = NEUTRON_POS;
  vel = NEUTRON_VEL;

  s_rad = sphere_radius;

  intersection = intersections[global_addr];

  a = dot(vel, vel);
  b = 2.0f*(dot(vel, (pos)));
  c = dot((pos), (pos)) - s_rad*s_rad;
  quotient = b*b - 4.0f*a*c;

  if (quotient > 0.0f) {
    d1 = (-b - sqrt(quotient)) / (2.f*a);
    d2 = (-b + sqrt(quotient)) / (2.f*a);
    if (d1 < 0.0 && d2 < INTERSECTION_T1) {
      INTERSECTION_POS1 = pos + d2*vel;
      INTERSECTION_T1 = d2;
      iidx[global_addr] = comp_idx;
      INTERSECTION_POS2 = pos + d2*vel;
      INTERSECTION_T2 = d2;
    }
  }

  intersections[global_addr] = intersection;
  neutrons[global_addr] = neutron;
}