__kernel void intersect_annulus(__global float16* neutrons, 
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx, float3 const annulus_pos,
  float const annulus_inner_radius, float const annulus_outer_radius,
  float const annulus_height) {

  uint global_addr = get_global_id(0);

  float16 neutron = neutrons[global_addr];
  float8 intersection;

  if (neutron.sf > 0.) {
    return;
  }

  if (neutron.sc == comp_idx) {
    return;
  }

  // First ray trace against the outer radius

  float2 plane_pos, plane_vel;
  float a, b, c, quotient, t1, t2, t3, t4;

  plane_pos = neutron.s02 - annulus_pos.s02;
  plane_vel = neutron.s35;

  a = dot(plane_vel, plane_vel);
  b = 2*dot(plane_pos, plane_vel);
  c = dot(plane_pos, plane_pos) - annulus_outer_radius*annulus_outer_radius;

  quotient = b*b-4*a*c;

  t1 = (-b - sqrt(quotient)) / (2.*a);
  t2 = (-b + sqrt(quotient)) / (2.*a);

  // Now ray trace against the inner radius

  c = dot(plane_pos, plane_pos) - annulus_inner_radius*annulus_inner_radius;
  quotient = b*b-4*a*c;

  t3 = (-b - sqrt(quotient)) / (2.*a);
  t4 = (-b + sqrt(quotient)) / (2.*a);

  intersection = intersections[global_addr];

  if (t1 > 0. && t1 < intersection.s3) {                  // Outside annulus
    intersection.s012 = neutron.s012 + t1*neutron.s345;
    intersection.s3   = t1;
    intersection.s456 = neutron.s012 + t3*neutron.s345;
    intersection.s7   = t3;
    iidx[global_addr] = comp_idx;
  } else if (t1 < 0. && t4 > 0. && t4 < intersection.s3) { // Inside annulus
    intersection.s012 = neutron.s012 + t4*neutron.s345;
    intersection.s3   = t4;
    intersection.s456 = neutron.s012 + t2*neutron.s345;
    intersection.s7   = t2;
    iidx[global_addr] = comp_idx;
  }

  neutrons[global_addr] = neutron;
  intersections[global_addr] = intersection;
}
