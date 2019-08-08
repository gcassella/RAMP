#include "rand.h"

__kernel void generate_neutrons(__global double16* neutrons,
    __global double8* intersections,
    double3 const source_pos, double3 const source_normal,
    double const source_radius, double2 const target_dimensions,
    double3 const target_pos, double const E,
    double const dE) {

  uint global_addr;

  double chi, x, y, r2, tx, ty, u1, u2, calcE, vel;
  double3 dir, target_norm;
  double16 neutron;

  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];

  // Choose a point on the moderator to emit neutrons
  chi = 2.*M_PI_F*(rand(&neutron, global_addr));
  
  r2 = pow(source_radius * rand(&neutron, global_addr), (double)2.);
  

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
  dir = (tx*cross(target_norm, (double3)( 0.0f, 1.0f, 0.0f )) + 
    ty*(double3)( 0.0f, 1.0f, 0.0f )) + target_pos - neutron.s012;

  // Generate a normally distributed wavelength

  u1 = rand(&neutron, global_addr);
  
  u2 = rand(&neutron, global_addr);
  
  calcE = (E + (sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2)) * dE);

  vel = 438.01*sqrt(calcE);

  dir = normalize(dir)*vel;

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
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}