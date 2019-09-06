__kernel void intersect_box(__global float16* neutrons, 
      __global float8* intersections, __global uint* iidx,
      uint const comp_idx, float const width, float const height,
      float const depth) {

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float3 pos, vel;
 
  if (neutron.sf > 0.f) {
    return;
  }

  if (neutron.sc == comp_idx) {
    return;
  }

  pos = neutron.s012;
  vel = neutron.s345;

  float3 A, B, tA, tB;
  float t1, t2, t3, t4, t5, t6, tmin, tmax;

  A = (float3){ -width/2.0f, -height/2.0f, -depth/2.0f };
  B = (float3){ width/2.0f, height/2.0f, depth/2.0f };

  tA = (A - pos) / vel;
  tB = (B - pos) / vel;

  t1 = tA.s0;
  t2 = tB.s0;
  t3 = tA.s1;
  t4 = tB.s1;
  t5 = tA.s2;
  t6 = tB.s2;
  tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
  tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

  if (tmax < 0.0f) {
    // Hit the box in the past
  } else if (tmin > tmax) {
    // No intersection
  } else {    
    intersection.s012 = pos + tmin*vel;
    intersection.s3 = tmin;
    intersection.s456 = pos + tmax*vel;
    intersection.s7 = tmax;

    iidx[global_addr] = comp_idx;
  }

  intersections[global_addr] = intersection;
  neutrons[global_addr] = neutron;
}