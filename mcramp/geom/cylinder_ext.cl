__kernel void intersect(__global float16* neutrons, 
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx, float const radius, float const height) {

  uint global_addr = get_global_id(0);

  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float2 plane_pos, plane_vel;
  float theta2, a, b, c, t1, t2, quotient;

  if (neutron.sf > 0.f) {
    return;
  }

  if (neutron.sc == comp_idx) {
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
      ((-height/2.0f) < (neutron.s1 + t2*neutron.s4)) &&
      ((neutron.s1 + t2*neutron.s4) < (height/2))) {
  
    if (t1 < intersection.s3 && t1 > 0.0f && t2 > 0.0f) {
        
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
