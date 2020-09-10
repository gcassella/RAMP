__kernel void beamstop(__global float16 *neutrons, __global float8 *intersections, 
                        __global uint *iidx, uint const comp_idx)
{

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  //float8 intersection = intersections[global_addr];
 
  uint this_iidx;
  this_iidx = iidx[global_addr];

  if (!(this_iidx == comp_idx))
  {
      return;
  }

  if (neutron.sf > 0.f)
  {
      return;
  }

  neutron.sf = 1.0f;

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
