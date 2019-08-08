#include "rand.h"
#include "consts.h"

double maxwell_energy_distn(double E, double T) {
  double rootE = sqrt(E);
  double overKT = 1.0 / (k_B * T);

  return 2*rootE/sqrt(M_PI)*pow(overKT, (double)(3.0/2.0))*exp(-E*overKT);
}

__kernel void generate_neutrons(__global double16* neutrons,
    __global double8* intersections, double2 const mod_dim,
    double2 const target_dim, double const target_dist, double const E_min,
    double const E_max, double const T1, double const I1, double const T2,
    double const I2, double const T3, double const I3, double const str_area,
    int const num_sim) {

  // FIXME: make sure the emission intensities are correct for this

  int global_addr;
  double16 neutron;
  
  global_addr = get_global_id(0);
  neutron = neutrons[global_addr];

  // Generate an energy val by linearly sampling the range, then weight
  // from a joint Maxwellian distribution according to T1/T2/T3

  double E_range, E_val, E_val_Joules, M_intensity, vel, Dx, Dy, distmax_T1, distmax_T2, distmax_T3;

  E_range = (E_max - E_min);
  E_val = E_min + E_range*rand(&neutron, global_addr);
  E_val_Joules = 1.6021765e-22 * E_val;

  distmax_T1 = maxwell_energy_distn(k_B*T1, T1);
  distmax_T2 = maxwell_energy_distn(k_B*T2, T2);
  distmax_T3 = maxwell_energy_distn(k_B*T3, T3);

  M_intensity = I1 * maxwell_energy_distn(E_val_Joules, T1) / distmax_T1;
  M_intensity += I2 * maxwell_energy_distn(E_val_Joules, T2) / distmax_T2;
  M_intensity += I3 * maxwell_energy_distn(E_val_Joules, T3) / distmax_T3;

  neutron.s9 = M_intensity / num_sim;

  neutron.s0 = mod_dim.x*(0.5 - rand(&neutron, global_addr));
  neutron.s1 = mod_dim.y*(0.5 - rand(&neutron, global_addr));
  neutron.s2 = 0.0;

  vel = 437.393377*sqrt(E_val);
  Dx = target_dim.x*(0.5 - rand(&neutron, global_addr)) - neutron.s0;
  Dy = target_dim.y*(0.5 - rand(&neutron, global_addr)) - neutron.s1;
  
  neutron.s345 = vel*normalize((double3)( Dx, Dy, target_dist ));

  // Revive terminated neutrons
  neutron.sf = 0.;

  neutrons[global_addr] = neutron;
  intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );
}