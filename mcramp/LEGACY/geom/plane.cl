__kernel void intersect_plane(__global double16* neutrons,
    __global double8* intersections, __global uint* iidx,
    uint const comp_idx, double3 const g_pos, double const width, 
    double const height) {

    /* Geometry kernel for linear components such as guides which have an
     * entrance rectangle and an exit rectangle. It is assumed that
     * we only care about neutrons who pass through the entrance rectangle.
     * It is possible to relax this requirement but I don't care to.
     *
     * The plane normal points along the z axis
     */

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

    double t = (g_pos.s2 - pos.s2) / vel.s2;
    double x = pos.s0 + t*vel.s0;
    double y = pos.s1 + t*vel.s1;

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