#include "rand.h"

__kernel void generate_neutrons(__global float16* neutrons,
    __global float8* intersections,
    float3 const source_pos, float3 const source_normal,
    float const source_radius, float2 const target_dimensions,
    float3 const target_pos, float const T) {

  uint global_addr;

  float chi, x, y, r2, tx, ty, u11, u12, u21, u22, u31, u32, wl, u1, u2, u3, prefac, vx, vy, vz;
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

  u11 = rand(&neutron, global_addr);
  u12 = rand(&neutron, global_addr);
  u1 = sqrt(-2*log(u11))*cos(2*M_PI*u12);
  u21 = rand(&neutron, global_addr);
  u22 = rand(&neutron, global_addr);
  u2 = sqrt(-2*log(u21))*cos(2*M_PI*u22);
  u31 = rand(&neutron, global_addr);
  u32 = rand(&neutron, global_addr);
  u3 = sqrt(-2*log(u31))*cos(2*M_PI*u32);

  prefac = 64.3;
  
  vx = prefac*sqrt(T)*u1;
  vy = prefac*sqrt(T)*u2;
  vz = prefac*sqrt(T)*u3;

  dir = normalize(dir)*length((float3)( vx, vy, vz ));

  neutron.s3 = dir.s0;
  neutron.s4 = dir.s1;
  neutron.s5 = dir.s2;

  // Initialize weight
  neutron.s9 = 1.;

  // Initialize time
  neutron.sa = 0.;

  // Revive terminated neutrons
  neutron.sf = 0.;

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}