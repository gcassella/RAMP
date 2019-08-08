__kernel void intersect_plane(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx, double const width, 
    double const height, uint const orientation) {

    uint global_addr        = get_global_id(0);
    double16 neutron         = neutrons[global_addr];
    double8 intersection     = intersections[global_addr];

    /* Check termination flag */
    if (neutron.sf > 0.) 
        return;

    if (neutron.sc == comp_idx) {
        return;
    }

    /* Perform raytracing here */

    double3 vel = neutron.s345;
    double3 pos = neutron.s012;

    double t, x, y;
    
    if (orientation == 0) {
        t = (-pos.s2) / vel.s2;
        x = pos.s0 + t*vel.s0;
        y = pos.s1 + t*vel.s1;
    } else if (orientation == 1) {
        t = (-pos.s0) / vel.s0;
        x = pos.s2 + t*vel.s2;
        y = pos.s1 + t*vel.s1;
    }

    if ((fabs(x) < width / 2) && (fabs(y) < height / 2)
        && t < intersection.s3
        && t < intersection.s7
        && t > 0
        && dot(vel, (double3)( 0.0f, 0.0f, 1.0f )) > 0.) {

        intersection.s012 = pos + t*vel;
        intersection.s456 = pos + t*vel;
        intersection.s3   = t;
        intersection.s7   = t;

        iidx[global_addr] = comp_idx;
    }

    /* ----------------------- */

    /* Update global memory */
    intersections[global_addr] = intersection;
    neutrons[global_addr]      = neutron;
}