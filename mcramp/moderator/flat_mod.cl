#include "rand.h"
#include "consts.h"

__kernel void generate_neutrons(__global float16* neutrons,
    __global float8* intersections, float2 const mod_dim,
    float2 const target_dim, float const target_dist, float const E_min,
    float const E_max, int const num_sim, int const seed) {

  // FIXME: make sure the emission intensities are correct for this

  int global_addr;
  float16 neutron;
  
  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];

  // Clear old neutron data and seed

  neutron = (float16){0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,(float)seed+neutron.sb,0.0f,0.0f,0.0f,0.0f};

  float E_range, E_val, vel, Dx, Dy;

  E_range = (E_max - E_min);
  E_val = E_min + E_range*rand(&neutron, global_addr);

  neutron.s9 = 1.0f / num_sim;

  neutron.s0 = mod_dim.x*(0.5f - rand(&neutron, global_addr));
  neutron.s1 = mod_dim.y*(0.5f - rand(&neutron, global_addr));
  neutron.s2 = 0.0f;

  vel = SE2V*sqrt(E_val);
  Dx = target_dim.x*(0.5f - rand(&neutron, global_addr)) - neutron.s0;
  Dy = target_dim.y*(0.5f - rand(&neutron, global_addr)) - neutron.s1;
  
  neutron.s345 = vel*normalize((float3)( Dx, Dy, target_dist ));

  // Revive terminated neutrons
  neutron.sf = 0.f;

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}