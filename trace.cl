#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

#include <Random123/threefry.h>
#include <Random123/u01fixedpt.h>

__kernel void intersect_sphere(__global float16* neutrons, 
      __global float8* intersections, float3 const sphere_pos, 
      float const sphere_radius) {

  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float3 pos, vel, s_pos;
  float s_rad, a, b, c, d1, d2, quotient;

  pos = neutron.s012;
  vel = neutron.s345;

  s_pos = sphere_pos;
  s_rad = sphere_radius;

  intersection = intersections[global_addr];

  a = dot(vel, vel);
  b = 2.0*(dot(vel, (pos - s_pos)));
  c = dot((pos - s_pos), (pos - s_pos)) - s_rad*s_rad;
  quotient = b*b - 4.0*a*c;

  if (quotient > 0) {
    d1 = (-b - sqrt(quotient)) / (2.*a);
    d2 = (-b + sqrt(quotient)) / (2.*a);
    if (d1 < intersection.s3) {
      if (d1>0) {
        intersection.s012 = pos + d1*vel;
        intersection.s3 = d1 + neutron.sa;
      }

      if (d2>0) {
        intersection.s456 = pos + d2*vel;
        intersection.s7 = d2 + neutron.sa;
      }
    }
  }

  intersections[global_addr] = intersection;
  neutrons[global_addr] = neutron;
}

/*
__kernel void intersect_plane(__global float16* neutrons,
  __global float8* intersections, float3 const plane_a,
  float3 const plane_b, float3 const plane_c) {
  
  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float3 p01, p02;
  float det, t, u, v;

  intersection = intersections[global_addr];

  p01 = plane_b - plane_a;
  p02 = plane_c - plane_a;

  det = -dot(neutron.s345, cross(p01, p02));
  
  if(det == 0.0) {
    // Does not intersect

    intersection = (float8)( 0.0f, 0.0f, 0.0f, 0.0f,
                             0.0f, 0.0f, 0.0f, 0.0f );
  } else {
    // Intersects with infinite plane
    t = (1/det)*dot(cross(p01, p02), neutron.s012 - plane_a);

    u = (1/det)*dot(cross(p02, -neutron.s345), neutron.s012 - plane_a);
    v = (1/det)*dot(cross(-neutron.s345, p01), neutrons.s012 - plane_a);

    if (0<=u<=1 && 0<=v<=1) {
      // Intersects with finite plane
      intersection.s0123 = (float4)( 0.0f, 0.0f, 0.0f, 0.0f );

      intersection.s456 = neutron.s012 + t*neutron.s345;
      intersection.s7 = t;
    } else {
      // Does not intersect

      intersection = (float8)( 0.0f, 0.0f, 0.0f, 0.0f,
                               0.0f, 0.0f, 0.0f, 0.0f );
    }
  }

  intersections[global_addr] = intersection;
}
*/

__kernel void intersect_banana(__global float16* neutrons, 
  __global float8* intersections, float3 const banana_pos,
  float const banana_radius, float const banana_height,
  float const mintheta, float const maxtheta) {

  uint global_addr = get_global_id(0);

  float16 neutron = neutrons[global_addr];
  float8 intersection;
  float2 plane_pos, plane_vel;
  float theta, a, b, c, t1, t2, quotient;

  plane_pos = neutron.s02 - banana_pos.s02;
  plane_vel = neutron.s35;

  a = dot(plane_vel, plane_vel);
  b = 2*dot(plane_pos, plane_vel);
  c = dot(plane_pos, plane_pos) - banana_radius*banana_radius;

  quotient = b*b-4*a*c;

  t1 = (-b - sqrt(quotient)) / (2.*a);
  t2 = (-b + sqrt(quotient)) / (2.*a);

  intersection = intersections[global_addr];

  theta = acos(dot(normalize(plane_pos+t2*plane_vel), (float2)( 0.0f, 1.0f )));

  if ((quotient > 0) &&
      ((banana_pos.s1 - banana_height/2) < (neutron.s1 + t2*neutron.s4)) &&
      ((neutron.s1 + t2*neutron.s4) < (banana_pos.s1 + banana_height/2)) &&
      (mintheta < theta < maxtheta)) {
    
    if (t1 < intersection.s3) {
      if (t1 > 0) {   
        intersection.s012 = neutron.s012 + t1*neutron.s345;
        intersection.s3   = neutron.sa + t1;
      }

      if (t2 > 0) {
        intersection.s456 = neutron.s012 + t2*neutron.s345;
        intersection.s7   = neutron.sa + t2;
      }
    }
  }

  intersections[global_addr] = intersection;
}

float rand(float16* neutron, uint tid) {
  float counter = (*neutron).sb;

  threefry4x32_key_t k = {{tid, 0xdecafbad, 0xfacebead, 0x12345678}};
  threefry4x32_ctr_t c = {{0, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};

  union {
    threefry4x32_ctr_t c;
    int4 i;
  } u;

  c.v[0]+=counter;
  u.c = threefry4x32(c, k);

  (*neutron).sb+=1;

  return u01fixedpt_closed_closed_32_24(u.i.x);
}

__kernel void generate_neutrons(__global float16* neutrons,
    float3 const source_pos, float3 const source_normal,
    float const source_radius, float2 const target_dimensions,
    float3 const target_pos, float const lambda,
    float const dlambda) {

  uint global_addr;

  float chi, x, y, r2, tx, ty, u1, u2, wl, vel;
  float3 dir, target_norm;
  float16 neutron;

  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];

  // Choose a point on the moderator to emit neutrons
  chi = 2.*M_PI_F*(rand(&neutron, global_addr));
  
  r2 = pow(source_radius * rand(&neutron, global_addr), (float)2.);
  

  x = sqrt(r2) * cos(chi);
  y = sqrt(r2) * sin(chi);

  neutron.s0 = x;
  neutron.s1 = y;
  neutron.s2 = source_pos.s2;

  // Choose a point on the target to aim towrad
  tx = (2*rand(&neutron, global_addr)-1)/2*target_dimensions.s0;
  
  ty = (2*rand(&neutron, global_addr)-1)/2*target_dimensions.s1;
  

  target_norm = normalize(source_pos - target_pos);

  // Calculate emission point -> target point vector
  dir = (tx*cross(target_norm, (float3)( 0.0f, 1.0f, 0.0f )) + 
    ty*(float3)( 0.0f, 1.0f, 0.0f )) + target_pos - neutron.s012;

  // Generate a normally distributed wavelength

  u1 = rand(&neutron, global_addr);
  
  u2 = rand(&neutron, global_addr);
  
  wl = sqrt(-2.0 * log(u1)) * cos(2 * M_1_PI * u2) * dlambda + lambda;
  vel = 3.97*pow(10.,-7.) / (wl*pow(10.,-10.));

  dir = normalize(dir)*vel;

  neutron.s3 = dir.s0;
  neutron.s4 = dir.s1;
  neutron.s5 = dir.s2;

  // Initialize weight
  neutron.s9 = 1;

  neutrons[global_addr] = neutron;
}

__kernel void reset_intersections(__global float8* intersections) {
  uint global_addr = get_global_id(0);

  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 0.0f );
}

__kernel void detector(__global float16* neutrons,
  __global float8* intersections,
  __global uint* histogram, float3 const det_pos,
  float3 const binning, uint const var) {

  uint global_addr = get_global_id(0);
  
  float16 neutron;
  float8 intersection;
  float varval;
  float minvar, maxvar, stepvar;
  uint idx;

  minvar = binning.s0;
  stepvar = binning.s1;
  maxvar = binning.s2;

  neutron = neutrons[global_addr];
  intersection = intersections[global_addr];

  if (length(intersection.s4567) != 0.) {
    switch(var) {
      case 0:
        varval = acos(dot(normalize(intersection.s46 - det_pos.s02), (float2)( 0.0f, 1.0f )));
        varval = clamp(varval, minvar, maxvar);
        break;
      case 1:
        varval = clamp(intersection.s7, minvar, maxvar);
        break;
      case 2:
        varval = clamp(neutron.se, minvar, maxvar);
        break;
      case 3:
        varval = clamp(neutron.sf, minvar, maxvar);
        break;
      default:
        break;
    }

    idx = round((varval - minvar) / stepvar);
    neutron.s012a = intersection.s4567;

    histogram[global_addr] = idx;
  }
}

__kernel void random_scatter(__global float16* neutrons, 
  __global float8* intersections, float3 const pos,
  float const radius) {

  uint global_addr = get_global_id(0);
  float16 neutron;
  float8 intersection;
  float u1, u2, vel;

  neutron = neutrons[global_addr];
  intersection = intersections[global_addr];

  if (length(intersection.s0123) != 0.) {

    neutron.s012a = (intersection.s0123 + intersection.s4567) / 2;

    vel = length(neutron.s345);

    u1 = rand(&neutron, global_addr) * 2 * M_PI;
    neutron.sb += 1;
    u2 = acos(2*rand(&neutron, global_addr)-1);
    neutron.sb += 1;

    neutron.s3 = vel*cos(u1)*sin(u2);
    neutron.s4 = vel*sin(u1)*sin(u2);
    neutron.s5 = vel*cos(u2);
  }

  neutrons[global_addr] = neutron;
}

__kernel void powder_scatter(__global float16* neutrons,
  __global float8* intersections, __global float3* reflections, 
  uint const num_reflections) {

  // Choose a reflection
  uint global_addr = get_global_id(0);
  float16 neutron = neutrons[global_addr];
  float8 intersection = intersections[global_addr];

  if (neutron.sf == 1.) {
    return;
  }

  uint NReflection, attempts;
  float3 reflection, perp, normvel;
  float x, y, z, Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, arg;
  float vel, q_v, alpha;

  normvel = normalize(neutron.s345);
  vel = length(neutron.s345);

  attempts = 0;
  
  if (length(intersection.s0123) != 0.) {
    neutron.s012a = (intersection.s0123 + intersection.s4567) / 2;

    do {
      NReflection = floor(rand(&neutron, global_addr)*num_reflections);
      

      reflection = reflections[NReflection];

      // Check if reflection is okay?
      q_v = 3.97*pow(10., 3.) / reflection.s1;

      arg = q_v / (2*vel);

      if (arg < 1.0) {
        // Okay, scatter
        // Generate an arbitrary perpendicular vector to the velocity
        // x*vx + y*vy + z*vz = (vx+vy)+z*vz = 0
        // => z = -(vx+vy)/vz

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

        // now construct rotation matrix to rotate the velocity 2theta
        // about the scattering vector
        alpha = 2*asin(arg);

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

        neutron.s345 = vel*normvel;
        neutron.s9 *= reflection.s2;

        break;

      } else {
        // Not allowed, try again
        attempts += 1;
      }
    } while (attempts < 100);
  }

  neutrons[global_addr] = neutron;
}

__kernel void isotropic_scatter(__global float16* neutrons,
  __global float8* intersections,
  __global float* q, __global float* w,
  __global float* pw_cdf, __global float* pq_cdf,
  uint const qsamples, uint const wsamples,
  float const rho, float const sigma_abs,
  float const sigma_scat) {

  uint global_addr    = get_global_id(0);
  uint w_index;
  float16 neutron     = neutrons[global_addr];
  float8 intersection = intersections[global_addr];

  float4 path;
  float3 perp, normvel;
  float sigma_tot, path_length, ki, sigma_s, sigma_a,
        mu, eta, u, v, Q, omega, kf, arg, theta, x, y, z, alpha,
        Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz, mindiff;

  path = intersection.s4567 - intersection.s0123;
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
    neutron.s9 = 0.;
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

  neutron.s012a = intersection.s0123 + path_length*path;
  neutron.s345 = normvel*kf/(1.583*pow(10.,-3.));

  neutron.se = Q;
  neutron.sf = omega;

  neutrons[global_addr] = neutron;
}