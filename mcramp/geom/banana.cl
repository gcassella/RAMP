__kernel void intersect_banana(__global float16* neutrons, 
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx, float3 const banana_pos,
  float const banana_radius, float const banana_height,
  float const mintheta, float const maxtheta) {

  uint global_addr = get_global_id(0);

  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float2 plane_pos, plane_vel;
  float theta2, a, b, c, t1, t2, quotient;

  if (neutron.sf > 0.) {
    return;
  }

  if (neutron.sc == comp_idx) {
    return;
  }

  plane_pos = neutron.s02 - banana_pos.s02;
  plane_vel = neutron.s35;

  a = dot(plane_vel, plane_vel);
  b = 2*dot(plane_pos, plane_vel);
  c = dot(plane_pos, plane_pos) - banana_radius*banana_radius;

  quotient = b*b-4*a*c;

  t1 = (-b - sqrt(quotient)) / (2.*a);
  t2 = (-b + sqrt(quotient)) / (2.*a);

  intersection = intersections[global_addr];

  theta2 = acos(dot(normalize(plane_pos+t2*plane_vel), (float2)( 0.0f, 1.0f )));
  if ((quotient > 0) &&
      ((banana_pos.s1 - banana_height/2) < (neutron.s1 + t2*neutron.s4)) &&
      ((neutron.s1 + t2*neutron.s4) < (banana_pos.s1 + banana_height/2)) &&
      (mintheta < theta2 < maxtheta)) {
  
    if (t2 < intersection.s3 && t1 < 0 && t2 > 0) {
        
        intersection.s012 = neutron.s012 + t1*neutron.s345;
        intersection.s3   = t1;
      
        iidx[global_addr] = comp_idx;
        

        intersection.s456 = neutron.s012 + t2*neutron.s345;
        intersection.s7   = t2;
    }
  }

  neutrons[global_addr] = neutron;
  intersections[global_addr] = intersection;
}
