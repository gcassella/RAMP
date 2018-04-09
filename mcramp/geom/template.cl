__kernel void intersect_template(__global float16* neutrons,
    __global float8* intersections, __global uint* iidx,
    uint const comp_idx) {

    uint global_addr        = get_global_id(0);
    float16 neutron         = neutrons[global_addr];
    float8 intersections    = intersections[global_addr];

    /* Check termination flag */
    if (neutron.sf > 0.) 
        return;

    /* Perform raytracing here */

    /* ----------------------- */

    /* Update global memory */
    intersections[global_addr] = intersection;
    neutrons[global_addr]      = neutron;
}