#include "consts.h"
#include "rand.h"

__kernel void chopper(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx, double const slit_width,
    double const radius, double const freq, uint const n_slits,
    double const phase, double const jitter) {

  uint global_addr        = get_global_id(0);
  double16 neutron         = neutrons[global_addr];
  double8 intersection = intersections[global_addr];
  uint this_iidx          = iidx[global_addr];

  /* Check we are scattering from the intersected component */
  if (!(this_iidx == comp_idx)) {
      return;
  }

  /* Check termination flag */
  if (neutron.sf > 0.)  {
      return;
  }

  /* Perform scattering here */

  neutron.s012 = intersection.s456;
  neutron.sa += intersection.s7;

  double Tg, toff, thi_chop,thi_neut,thi_diff,ang_diff,width,pha_error;

  Tg = 2*M_PI / (double)n_slits;
  
	/*pha_error is the error in the phasing of the chopper in radians*/
  pha_error=jitter*2.0*(rand(&neutron, global_addr)-0.5); 
  thi_chop=freq*(neutron.sa - fabs(phase))+pha_error;
  thi_neut=atan2(neutron.s0, neutron.s1+radius);
  thi_diff=fabs(thi_chop-thi_neut); 
  ang_diff=fmod(thi_diff,Tg);
  if (fabs(ang_diff-Tg) < ang_diff)
	  ang_diff=fabs(ang_diff-Tg);

  width=atan2(slit_width, neutron.s1+radius);

  /* does neutron hit the slit? */
  if (ang_diff>width/2.0) 
    neutron.sf = 1.0;
  

  /* ----------------------- */

  /* Update global memory and reset intersection */
  iidx[global_addr] = 0;
  neutron.sc = comp_idx;

  neutrons[global_addr]      = neutron;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                       0.0f, 0.0f, 0.0f, 100000.0f );

}