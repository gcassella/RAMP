#include "rand.h"
#include "geom.h"

__kernel void vanadium_scatter(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection     = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component */
    if (!(this_iidx == comp_idx))
        return;

    /* Check termination flag */
    if (neutron.sf > 0.) 
        return;

    /* Perform scattering here */

    float3 path, perp, normvel;
    float vel, phi, theta, x, y, z, q2ki, alpha;

    q2ki       = rand(&neutron, global_addr);

    vel     = length(neutron.s345);
    normvel = normalize(neutron.s345);

    x = 1.;
    y = 1.;
    z = -(neutron.s3+neutron.s4)/(neutron.s5);
    perp = normalize((float3)( x, y, z ));

    // construct rotation matrix to randomly rotate the scattering
    // vector about the velocity

    alpha = 2*M_PI*rand(&neutron, global_addr);
    
    rotate_about_axis(alpha, normvel, (&perp));

    alpha = 2*asin(q2ki);

    rotate_about_axis(alpha, perp, (&normvel));

    neutron.s345 = (vel)*normvel;

    path         = intersection.s456 - intersection.s012;
    neutron.s012 = intersection.s012 + rand(&neutron, global_addr)*(path);
    neutron.sa  += intersection.s7;

    /* ----------------------- */

    /* Update global memory and reset intersection */
    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}