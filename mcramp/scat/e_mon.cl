void atomicAdd_g_f(volatile __global float *addr, float val)
{
    union {
        unsigned int u32;
        float        f32;
    } next, expected, current;
	current.f32    = *addr;
    do {
	   expected.f32 = current.f32;
        next.f32     = expected.f32 + val;
		current.u32  = atomic_cmpxchg( (volatile __global unsigned int *)addr, 
                            expected.u32, next.u32);
    } while( current.u32 != expected.u32 );
}



__kernel void detector(__global double16 *neutrons,
                       __global double8 *intersections, __global uint *iidx,
                       uint const comp_idx, volatile __global float *histogram,
                       double3 const binning, uint const restore_neutron)
{

  uint global_addr = get_global_id(0);
  double16 neutron = neutrons[global_addr];
  double8 intersection = intersections[global_addr];
  double ener_val, min_var, step_var, max_var;

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

  min_var = binning.s0;
  step_var = binning.s1;
  max_var = binning.s2;

  ener_val = pow(length(neutron.s345) / 438.01, 2.0);

  if(min_var<=ener_val && ener_val<=max_var) {    
    idx = round((ener_val -  min_var) / step_var);
    atomicAdd_g_f(&histogram[idx], (float)neutron.s9);
    neutron.se = idx;
  } else {
    neutron.se = -1;
  }
  
  
  if (restore_neutron == 0) {
    neutron.s012 = intersection.s456;
    neutron.sa += intersection.s7;
    neutron.sf = 1.;
  }

  iidx[global_addr] = 0;
  neutron.sc = comp_idx;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

  neutrons[global_addr] = neutron;
}
