#include "rand.h"
#include "geom.h"

__kernel void isotropic_scatter(__global double16* neutrons,
  __global double8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global double* q, __global double* w,
  __global double* pw_cdf, __global double* pq_cdf,
  uint const qsamples, uint const wsamples,
  double const rho, double const sigma_abs,
  double const sigma_scat) {

  uint global_addr    = get_global_id(0);
  uint w_index;
  double16 neutron     = neutrons[global_addr];
  double8 intersection = intersections[global_addr];
  uint this_iidx = iidx[global_addr];
 
  if(!(this_iidx == comp_idx)) {
    return;
  }

  if (neutron.sf > 0.) {
    return;
  }

  double3 path, perp, normvel;
  double sigma_tot, path_length, ki, sigma_s, sigma_a,
        mu, eta, u, v, Q, omega, kf, arg, theta, x, y, z, alpha,
        Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, mindiff, TOF;

  path = intersection.s456 - intersection.s012;
  ki = 1.583*pow(10.,-3.)*length(neutron.s345);

  normvel = normalize(neutron.s345);

  sigma_s = sigma_scat / (2.*ki*ki);
  sigma_a = sigma_abs * 2200. / length(neutron.s345);
  sigma_tot = sigma_a + sigma_s;

  mu = rho*sigma_tot*100.;


  // Monte carlo choice to see if our neutron scatters
  if (rand(&neutron, global_addr) < exp(-mu*length(path.s012))) {
    // Transmitted, return without modifying
    // neutron state, but multiply by weight factor
    neutron.s9 *= 1.0 - sigma_s / sigma_tot;
    neutron.s012 = (intersection.s456+0.01f*normalize(path));
    neutron.sa += intersection.s7;

    neutron.sc = comp_idx;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;
  } else {
    // Scattered, multiply by weight factor
    // to model absorption
    neutron.s9 *= sigma_s / sigma_tot;
  }

  // Monte carlo choice to find scattering point along
  // path through sample

  eta = rand(&neutron, global_addr);
  path_length = (-1./mu)*log(1 - eta*(1 - exp(-mu*length(path.s012))));

  // Monte carlo choice to find values of Q and w from
  // cumulative distribution functions using inverse
  // transform sampling

  u = rand(&neutron, global_addr);
  v = rand(&neutron, global_addr);

  mindiff=1.;

  for(uint i=0;i<wsamples;i++) {
    if (pw_cdf[i] == 1.) {
      break;
    }

    if (fabs(pw_cdf[i] - v) < mindiff) {
      mindiff = fabs(pw_cdf[i] - v);
      omega = w[i] + (rand(&neutron, global_addr) - 0.5f)*(w[i+1] - w[i]);
      w_index = i;
    } 
  }

  mindiff=1.;

  for(uint i=0;i<qsamples;i++) {
    if (pq_cdf[w_index*qsamples + i] == 1.) {
      break;
    }

    if (fabs(pq_cdf[w_index*qsamples + i] - u) < mindiff) {
      mindiff = fabs(pq_cdf[w_index*qsamples + i] - u);
      Q = q[i] + (rand(&neutron, global_addr) - 0.5f)*(q[i+1] - q[i]);
    } 
  }

  // Test conservation laws can be satisfied

  if (ki*ki < (0.693*omega)) {
    // No real valued kf possible
    
    neutron.sf = 1;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;
  }

  if (rand(&neutron, global_addr) > 0.5) {
    kf = sqrt(ki*ki - 0.693*omega);
  } else {
    kf = sqrt(ki*ki + 0.693*omega);
  }

  arg = (ki*ki + kf*kf - Q*Q)/(2*ki*kf);

  if (fabs(arg) > 1) {
    // Unphysical scattering direction
    neutron.sf = 1;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;

  } else {
    theta = acos(arg);
  }

  // Rotate neutron wavevector theta around
  // random vector perpendicular to it, first we
  // construct the perpendicular vector and rotate
  // it randomly on a unit circle

  x = 1.;
  y = 1.;
  z = -(neutron.s3+neutron.s4)/(neutron.s5);
  perp = normalize((double3)( x, y, z ));
  // construct rotation matrix to randomly rotate the scattering
  // vector about the velocity
  alpha = 2*M_PI*rand(&neutron, global_addr);

  rotate_about_axis(alpha, normvel, (&perp));
  
  alpha = theta;

  rotate_about_axis(alpha, perp, (&normvel));

  TOF = path_length / length(neutron.s345);

  neutron.s012 = intersection.s012 + (path_length)*path;
  neutron.sa   += intersection.s3 + TOF;
  
  neutron.s345 = normvel*kf/(double)(1.583*pow(10.,-3.));

  neutron.sc = comp_idx;
  iidx[global_addr] = 0;
  neutrons[global_addr] = neutron;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}