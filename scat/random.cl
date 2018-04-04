#include "rand.h"

__kernel void random_scatter(__global float16* neutrons, 
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx) {

  uint global_addr = get_global_id(0);
  float16 neutron;
  float8 intersection;
  float u1, u2, vel;

  uint this_iidx = iidx[global_addr];

  if(!(this_iidx == comp_idx)) {
    return;
  }
  

  neutron = neutrons[global_addr];
  intersection = intersections[global_addr]; 

  neutron.se = comp_idx;

  if (neutron.sf > 0.) {
    return;
  }

  if (length(intersection.s0123) != 0.) {

    neutron.s012a = (intersection.s0123 + intersection.s4567) / 2;

    vel = length(neutron.s345);

    u1 = rand(&neutron, global_addr) * 2 * M_PI;
    neutron.sb += 1;
    u2 = acos(2*rand(&neutron, global_addr)-1);
    neutron.sb += 1;

    neutron.s3 = vel*cos(u1)*sin(u2);
    neutron.s4 = vel*sin(u1)*sin(u2);
    neutron.s5 = vel*cos(u2);
  }

  neutrons[global_addr] = neutron;
}
