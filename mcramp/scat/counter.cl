__kernel void counter(__global double16 *neutrons,
                       __global double8 *intersections, __global uint *iidx,
                       uint const comp_idx)
{

  uint global_addr = get_global_id(0);
  double16 neutron = neutrons[global_addr];
  double8 intersection = intersections[global_addr];
 
  uint this_iidx, idx;
  this_iidx = iidx[global_addr];

  if (!(this_iidx == comp_idx))
  {
      return;
  }

  if (neutron.sf > 0.)
  {
      return;
  }

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
