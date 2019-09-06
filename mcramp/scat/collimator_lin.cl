__kernel void collimator(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx, float const length, float const slope_H,
    float const slope_V, float const transmission) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component */
    if (!(this_iidx == comp_idx)) {
        return;
    }

    /* Check termination flag */
    if (neutron.sf > 0.f)  {
        return;
    }

    /* Perform scattering here */

    float phi;

    if (slope_H > 0.0f) {
        phi = fabs(neutron.s3 / neutron.s5);
        if (phi > slope_H) {
            neutron.sf = 1.0f
        } else {
            neutron.s9 *= transmission*(1.0f - phi / slope_H);
        }
    }

    if (slope_V > 0.0f) {
        phi = fabs(neutron.s4 / neutron.s5);
        if (phi > slope_V) {
            neutron.sf = 1.0f;
        } else {
            neutron.s9 *= transmission*(1.0f - phi / slope_V);
        }
    }

    /* ----------------------- */

    /* Update global memory and reset intersection */
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;

    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}