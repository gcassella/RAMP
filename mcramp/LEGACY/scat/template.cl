#include "rand.h"

__kernel void scatter_template(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx) {

    uint global_addr        = get_global_id(0);
    double16 neutron         = neutrons[global_addr];
    double8 intersection     = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component */
    if (!(this_iidx == comp_idx))
        return;

    /* Check termination flag */
    if (neutron.sf > 0.) 
        return;

    /* Perform scattering here */

    /* ----------------------- */

    /* Update global memory and reset intersection */
    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (double8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}