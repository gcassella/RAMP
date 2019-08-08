__kernel void collimator(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx, double const length, double const slope_H,
    double const slope_V, double const transmission) {

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

    double phi;

    if (slope_H > 0.0) {
        phi = fabs(neutron.s3 / neutron.s5);
        if (phi > slope_H) {
            neutron.sf = 1.0
        } else {
            neutron.s9 *= transmission*(1.0 - phi / slope_H);
        }
    }

    if (slope_V > 0.0) {
        phi = fabs(neutron.s4 / neutron.s5);
        if (phi > slope_V) {
            neutron.sf = 1.0;
        } else {
            neutron.s9 *= transmission*(1.0 - phi / slope_V);
        }
    }

    /* ----------------------- */

    /* Update global memory and reset intersection */
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;

    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}