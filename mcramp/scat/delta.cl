#include "rand.h"
#include "geom.h"

#ifndef K2V
#define K2V 629.622368
#endif

#ifndef V2K
#define V2K 1.58825361e-3
#endif

#ifndef V2SE
#define V2SE 5.22703725e-6
#endif

#ifndef SE2V
#define SE2V 437.393377 
#endif

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

__kernel void delta(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx, uint const comp_idx,
    float const twotheta, float const phi, float const deltaE) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersection     = intersections[global_addr];
    uint this_iidx          = iidx[global_addr];

    /* Check we are scattering from the intersected component -------------- */
    if (!(this_iidx == comp_idx))
        return;

    /* Check termination flag ---------------------------------------------- */
    if (neutron.sf > 0.) 
        return;

    /* Perform scattering here --------------------------------------------- */

    float ki, kf, Ei, Ef, vel, path_length, scat_point, vel_f;
    float3 Q_v, path;

    scat_point = rand(&neutron, global_addr);
    neutron.s012 = intersection.s012 + scat_point*(intersection.s456 - intersection.s012);
    neutron.sa += intersection.s3 + scat_point*(intersection.s7 - intersection.s3);

    vel = length(neutron.s345);

    Ei = pow(vel / 438.01f, 2.0f);
    Ef = Ei + deltaE;
    vel_f = 438.01f*sqrt(Ef);

    Q_v = vel_f*(float3){ -sin(twotheta)*cos(phi), -sin(twotheta)*sin(phi), (vel / vel_f) - cos(twotheta) };

    neutron.s345 -= Q_v;

    /* ----------------------- */

    /* Update global memory and reset intersection */
    
    iidx[global_addr] = 0;
    neutron.sc = comp_idx;
    neutrons[global_addr]      = neutron;
    intersections[global_addr] = (float8)( 0.0f, 0.0f, 0.0f, 100000.0f,
                                         0.0f, 0.0f, 0.0f, 100000.0f );

}