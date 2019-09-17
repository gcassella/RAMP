#include "rand.h"
#include "consts.h"

float maxwell_energy_distn(float E, float T) {
  float rootE = sqrt(E);
  float overKT = 1.0f / (k_B * T);

  return 2.0f*rootE/sqrt(M_PI)*pow(overKT, (float)(3.0f/2.0f))*exp(-E*overKT);
}

__kernel void generate_neutrons(__global float16* neutrons,
    __global float8* intersections, float2 const mod_dim,
    float2 const target_dim, float const target_dist, float const E_min,
    float const E_max, float const T1, float const I1, float const T2,
    float const I2, float const T3, float const I3, float const str_area,
    int const num_sim, int const seed) {

  // FIXME: make sure the emission intensities are correct for this

  int global_addr;
  float16 neutron;
  
  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];

  // Clear old neutron data and seed

  neutron = (float16){0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,(float)seed+neutron.sb,0.0f,0.0f,0.0f,0.0f};

  // Generate an energy val by linearly sampling the range, then weight
  // from a joint Maxwellian distribution according to T1/T2/T3

  float E_range, E_val, E_val_Joules, M_intensity, vel, Dx, Dy, distmax_T1, distmax_T2, distmax_T3;

  E_range = (E_max - E_min);
  E_val = E_min + E_range*rand(&neutron, global_addr);
  E_val_Joules = 1.6021765e-22f * E_val;

  distmax_T1 = maxwell_energy_distn(k_B*T1, T1);
  distmax_T2 = maxwell_energy_distn(k_B*T2, T2);
  distmax_T3 = maxwell_energy_distn(k_B*T3, T3);

  M_intensity = I1 * maxwell_energy_distn(E_val_Joules, T1) / distmax_T1;
  M_intensity += I2 * maxwell_energy_distn(E_val_Joules, T2) / distmax_T2;
  M_intensity += I3 * maxwell_energy_distn(E_val_Joules, T3) / distmax_T3;

  neutron.s9 = M_intensity / num_sim;

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