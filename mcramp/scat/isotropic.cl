#include "rand.h"
#include "geom.h"
#include "consts.h"

float2 sample_distn(float16* neutron, uint global_addr, __global float* q, __global float* w, __global float* pw_cdf, 
  __global float* pq_cdf, uint const qsamples, uint const wsamples) {
    
  float mindiff, u, v, omega, Q;
  uint w_index;
  
  u = rand(neutron, global_addr);
  v = rand(neutron, global_addr);
  mindiff=1.;

  for(uint i=0;i<wsamples;i++) {
    if (pw_cdf[i] == 1.) {
      break;
    }

    if (fabs(pw_cdf[i] - v) < mindiff) {
      mindiff = fabs(pw_cdf[i] - v);
      omega = (i > 0 ? w[i-1] : 0.0f) + (rand(neutron, global_addr))*(w[i] - (i > 0 ? w[i-1] : 0.0f));
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
      Q = (i > 0 ? q[i-1] : 0.0f) + (rand(neutron, global_addr))*(q[i] - (i > 0 ? q[i-1] : 0.0f));
    } 
  }

  return (float2){ omega, Q };
}

__kernel void isotropic_scatter(
  __global float16* neutrons,
  __global float8* intersections,
  __global uint* iidx,
  uint const comp_idx,
  __global float* coh_q, 
  __global float* coh_w,
  __global float* coh_pw_cdf, 
  __global float* coh_pq_cdf,
  uint const coh_qsamples, 
  uint const coh_wsamples,
  float const coh_rho, 
  float const coh_sigma_abs,
  float const coh_sigma_scat, 
  __global float* inc_q, 
  __global float* inc_w,
  __global float* inc_pw_cdf, 
  __global float* inc_pq_cdf,
  uint const inc_qsamples, 
  uint const inc_wsamples,
  float const inc_rho, 
  float const inc_sigma_abs,
  float const inc_sigma_scat, 
  __global float* mag_q, 
  __global float* mag_w,
  __global float* mag_pw_cdf, 
  __global float* mag_pq_cdf,
  uint const mag_qsamples, 
  uint const mag_wsamples,
  float const mag_rho, 
  float const mag_sigma_abs,
  float const mag_sigma_scat, 
  float const temperature,
  uint const transmit) {

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
  float sigma_tot, path_length, ki, sigma_s, sigma_a, p_trans,
        mu, eta, u, v, Q, omega, kf, arg, theta, x, y, z, alpha,
        Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, mindiff, TOF;

  uint flag = 0;

  if (all(intersection.s012 == intersection.s456))
    path = INTERSECTION_POS2 - NEUTRON_POS;
  else
    path = INTERSECTION_POS2 - INTERSECTION_POS1;

  ki = V2K*length(neutron.s345);

  normvel = normalize(neutron.s345);

  sigma_s = (coh_sigma_scat + inc_sigma_scat + mag_sigma_scat) / (2.*ki*ki);
  sigma_a = (coh_sigma_abs + inc_sigma_abs + mag_sigma_abs) * 2200. / length(neutron.s345);
  sigma_tot = sigma_a + sigma_s;

  mu = (coh_rho + inc_rho + mag_rho)*sigma_tot*100.;

  // Monte carlo choice to see if our neutron scatters
  p_trans = exp(-mu*length(path.s012));
  if (rand(&neutron, global_addr) < p_trans) {
    // Transmitted, return without modifying
    // neutron state, but multiply by weight factor
    neutron.s9 *= p_trans;
    neutron.s012 = (intersection.s456+0.01f*normalize(path));
    neutron.sa += intersection.s7;

    if (transmit == 0)
      neutron.sf = 1.0f;

    neutron.sc = comp_idx;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;
  } else {
    // Scattered, multiply by weight factor
    // to model absorption
    neutron.s9 *= (1 - p_trans) * sigma_s / sigma_tot;
  }

  // Monte carlo choice to find scattering point along
  // path through sample

  eta = rand(&neutron, global_addr);
  path_length = (-1./mu)*log(1 - eta*(1 - exp(-mu*length(path.s012))));

  // decide between incoherent and coherent

  u = rand(&neutron, global_addr);
  u *= sigma_s;

  if (u < (coh_sigma_scat / (2.*ki*ki))) {
    flag = 1;
  } else if (u > ((coh_sigma_scat) / (2.*ki*ki)) &&
             u < ((coh_sigma_scat + inc_sigma_scat) / (2.*ki*ki)))  {
    flag = 2;
  } else if (u > ((coh_sigma_scat + inc_sigma_scat) / (2.*ki*ki))) {
    flag = 3;
  }

  // Monte carlo choice to find values of Q and w from
  // cumulative distribution functions using inverse
  // transform sampling

  float2 sample = (float2){ 0.0f, 0.0f };

  if (flag == 1) {
    sample = sample_distn(
      &neutron, 
      global_addr, 
      coh_q, 
      coh_w, 
      coh_pw_cdf, 
      coh_pq_cdf, 
      coh_qsamples, 
      coh_wsamples
    );
  } else if (flag == 2) {
    sample = sample_distn(
      &neutron, 
      global_addr, 
      inc_q, 
      inc_w, 
      inc_pw_cdf, 
      inc_pq_cdf, 
      inc_qsamples, 
      inc_wsamples
    );
  } else if (flag == 3) {
    sample = sample_distn(
      &neutron, 
      global_addr, 
      mag_q, 
      mag_w, 
      mag_pw_cdf, 
      mag_pq_cdf, 
      mag_qsamples, 
      mag_wsamples
    );
  }
  
  omega = sample.s0;
  Q = sample.s1;

  // Test conservation laws can be satisfied

  if (ki*ki < (E2KS*omega)) {
    // No real valued kf possible
    
    neutron.sf = 1;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;
  }

  if (rand(&neutron, global_addr) > 0.5*exp(-omega / (kB*temperature))) {
    kf = sqrt(ki*ki - E2KS*omega);
    neutron.se = omega;
  } else {
    kf = sqrt(ki*ki + E2KS*omega); 
    neutron.se = -omega;
  }

  arg = (ki*ki + kf*kf - Q*Q)/(2*ki*kf);

  if (fabs(arg) > 1) {
    // Unphysical scattering direction
    neutron.sf = 1;
    iidx[global_addr] = 0;
    neutrons[global_addr] = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
    return;

  } else {
    theta = acos(arg);
  }

  // Rotate neutron wavevector theta around
  // random vector perpendicular to it, first we
  // construct the perpendicular vector and rotate
  // it randomly on a unit circle

  float3 oldvel = neutron.s345;

  x = 1.;
  y = 1.;
  z = -(neutron.s3+neutron.s4)/(neutron.s5);
  perp = normalize((float3)( x, y, z ));
  // construct rotation matrix to randomly rotate the scattering
  // vector about the velocity
  alpha = 2*M_PI*rand(&neutron, global_addr);

  rotate_about_axis(alpha, normvel, (&perp));
  
  alpha = theta;

  rotate_about_axis(alpha, perp, (&normvel));

  TOF = path_length / length(neutron.s345);

  neutron.sd = Q;

  neutron.s012 = intersection.s012 + (path_length)*path;
  neutron.sa   += intersection.s3 + TOF;
  
  neutron.s345 = normvel*kf*K2V;

  // Modify beam polarisation based on scattering
  float3 q_norm = normalize((oldvel - neutron.s345)*V2K);

  if (flag == 2) {
    // Incoherent, model spin flip 2/3
    neutron.s678 *= -1.0f/3.0f;
  } else if (flag == 3) {
    // Magnetic, model directional polarisation
    // according to halpern-johnson eqn
    
    neutron.s678 = -q_norm*dot(neutron.s678, q_norm);
  }

  neutron.sc = comp_idx;
  iidx[global_addr] = 0;
  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}