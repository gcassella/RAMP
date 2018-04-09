#include "rand.h"

__kernel void isotropic_scatter(__global float16* neutrons,
  __global float8* intersections, __global uint* iidx,
  uint const comp_idx,
  __global float* q, __global float* w,
  __global float* pw_cdf, __global float* pq_cdf,
  uint const qsamples, uint const wsamples,
  float const rho, float const sigma_abs,
  float const sigma_scat) {

  uint global_addr    = get_global_id(0);
  uint w_index;
  float16 neutron     = neutrons[global_addr];
  float8 intersection = intersections[global_addr];
  uint this_iidx = iidx[global_addr];
 
  if(!(this_iidx == comp_idx)) {
    return;
  }

  if (neutron.sf > 0.) {
    return;
  }

  float3 path, perp, normvel;
  float sigma_tot, path_length, ki, sigma_s, sigma_a,
        mu, eta, u, v, Q, omega, kf, arg, theta, x, y, z, alpha,
        Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, mindiff, TOF;

  path = intersection.s456 - intersection.s012;
  ki = 1.583*pow(10.,-3.)*length(neutron.s345);

  normvel = normalize(neutron.s345);

  sigma_s = sigma_scat / (2*ki*ki);
  sigma_a = sigma_abs * 2200 / length(neutron.s345);
  sigma_tot = sigma_a + sigma_s;

  mu = rho*sigma_tot*100.;


  // Monte carlo choice to see if our neutron scatters

  if (rand(&neutron, global_addr) < exp(-mu*length(path.s012))) {
    // Transmitted, return without modifying
    // neutron state
    neutron.sc = comp_idx;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    neutrons[global_addr] = neutron;
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
      omega = w[i];
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
      Q = q[i];
    } 
  }

  // Test conservation laws can be satisfied

  if (ki*ki < (0.007706*omega)) {
    // No real valued kf possible
    return;
  }

  if (rand(&neutron, global_addr) > 0.5) {
    kf = sqrt(ki*ki - 0.007706*omega);
  } else {
    kf = -sqrt(ki*ki - 0.007706*omega);
  }

  arg = (ki*ki + kf*kf - Q*Q)/(2*ki*kf);

  if (fabs(arg) > 1) {
    // Unphysical scattering direction
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
  perp = normalize((float3)( x, y, z ));
  // construct rotation matrix to randomly rotate the scattering
  // vector about the velocity
  alpha = 2*M_PI*rand(&neutron, global_addr);
  
  Rxx = cos(alpha)+normvel.s0*normvel.s0*(1-cos(alpha));
  Rxy = normvel.s0*normvel.s1*(1-cos(alpha))-normvel.s2*sin(alpha);
  Rxz = normvel.s0*normvel.s2*(1-cos(alpha))+normvel.s1*sin(alpha);
  Ryx = normvel.s1*normvel.s0*(1-cos(alpha))+normvel.s2*sin(alpha);
  Ryy = cos(alpha)+normvel.s1*normvel.s1*(1-cos(alpha));
  Ryz = normvel.s1*normvel.s2*(1-cos(alpha))-normvel.s0*sin(alpha);
  Rzx = normvel.s2*normvel.s0*(1-cos(alpha))-normvel.s1*sin(alpha);
  Rzy = normvel.s2*normvel.s1*(1-cos(alpha))+normvel.s0*sin(alpha);
  Rzz = cos(alpha)+normvel.s2*normvel.s2*(1-cos(alpha));

  perp = (float3)( Rxx*perp.s0+Rxy*perp.s1+Rxz*perp.s2, 
                   Ryx*perp.s0+Ryy*perp.s1+Ryz*perp.s2, 
                   Rzx*perp.s0+Rzy*perp.s1+Rzz*perp.s2 );
  
  alpha = theta;

  Rxx = cos(alpha)+perp.s0*perp.s0*(1-cos(alpha));
  Rxy = perp.s0*perp.s1*(1-cos(alpha))-perp.s2*sin(alpha);
  Rxz = perp.s0*perp.s2*(1-cos(alpha))+perp.s1*sin(alpha);
  Ryx = perp.s1*perp.s0*(1-cos(alpha))+perp.s2*sin(alpha);
  Ryy = cos(alpha)+perp.s1*perp.s1*(1-cos(alpha));
  Ryz = perp.s1*perp.s2*(1-cos(alpha))-perp.s0*sin(alpha);
  Rzx = perp.s2*perp.s0*(1-cos(alpha))-perp.s1*sin(alpha);
  Rzy = perp.s2*perp.s1*(1-cos(alpha))+perp.s0*sin(alpha);
  Rzz = cos(alpha)+perp.s2*perp.s2*(1-cos(alpha));
  
  normvel = (float3)( Rxx*normvel.s0+Rxy*normvel.s1+Rxz*normvel.s2, 
                      Ryx*normvel.s0+Ryy*normvel.s1+Ryz*normvel.s2, 
                      Rzx*normvel.s0+Rzy*normvel.s1+Rzz*normvel.s2 );

  TOF = path_length / length(neutron.s345);

  neutron.s012 = intersection.s012 + path_length*path;
  neutron.sa   += intersection.s3 + TOF;
  
  neutron.s345 = normvel*kf/(1.583*pow(10.,-3.));

  neutron.sc = comp_idx;
  iidx[global_addr] = 0;
  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}