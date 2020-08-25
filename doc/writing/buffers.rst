Conventions: Buffer structure
=============================

The neutron buffer
------------------

The neutron buffer holds a ``float16`` entry for each simulated trajectory. The elements of these entries correspond to

 - ``(x, y, z)``, the neutron position. In kernels, this is given in the component frame.
 - ``(vx, vy, vz)``, the neutron velocity. Again, in the component frame.
 - ``(px, py, pz)``, the neutron spin. As of the current update, unused.
 - ``p``, the neutron weighting i.e. how many counts this trajectory will contribute to a detector.
 - ``t``, time of flight.
 - ``mc``, a counter that is incremented upon each Monte Carlo decision, used in random number generation.
 - ``prevcomp``, the index of the previous component hit. Prevents spurious multiple collisions due to rounding errors.
 - ``(_, _)``, two unused entries. May be used as custom flags / storage by kernels but this should be made extremely explicit in documentation to avoid conflicts.
 - ``T``, termination flag. Trajectory is no longer simulated by kernels once this flag is set not equal to zero.

Altering the data in the neutron buffer is the role of scattering kernels. These entires have assosciated convenient
macros defined in the header file ``consts.h`` as follows

.. code-block:: C

    neutron.s0 == NEUTRON_X     // neutron x coordinate
    neutron.s1 == NEUTRON_Y     // neutron y coordinate
    neutron.s2 == NEUTRON_Z     // neutron z coordinate
    neutron.s012 == NEUTRON_POS // float3 of neutron coordinates (x, y, z)
    neutron.s3 == NEUTRON_VX    // neutron x velocity
    neutron.s4 == NEUTRON_VY    // neutron y velocity
    neutron.s5 == NEUTRON_VZ    // neutron z velocity
    neutron.s345 == NEUTRON_VEL // float3 of neutron velocity (vx, vy, vz)
    neutron.s6 == NEUTRON_PX    // x component of neutron polarisation
    neutron.s7 == NEUTRON_PY    // y component of neutron polarisation
    neutron.s8 == NEUTRON_PZ    // z component of neutron polarisation
    neutron.s678 == NEUTRON_POL // float3 of neutron polarisation (px, py, pz)
    neutron.s9 == NEUTRON_P     // monte carlo weighting
    neutron.sa == NEUTRON_TOF   // neutron time-of-flight
    neutron.sf == NEUTRON_DIE   // neutron termination flag

For these macros to function as intended, the neutron at ``global_addr`` in the neutron buffer must be assigned to a
``float16`` named ``neutron``, i.e. the following code should appear at the start of a kernel that uses them

.. code-block:: C

    uint global_addr = get_global_id(0);
    float16 neutron = neutrons[global_addr];

The intersection buffer
-----------------------

The intersection buffer holds a ``float8`` entry for each simulated trajectory. The elements of these entries correspond to

 - ``(x1, y1, z1)``, the coordinates (in component frame) of the earliest intersection with the instrument geometry.
 - ``t1``, the time of flight from the neutron position at the previous scattering step to ``(x1, y1, z1)``.
 - ``(x2, y2, z2)``, the coordinates (in component frame) of the second intersection point with the component intersected at ``(x1, y1, z1)``, e.g. where the neutron leaves a component.
 - ``t2``, the time of flight from the neutron position at the previous scattering step to ``(x2, y2, z2)``.

The units of distance are meters and the units of time are seconds. As for the neutron buffer, the following macros are defined in ``consts.h``

.. code-block:: C

    intersection.s0 == INTERSECTION_X1 // x coordinate of first intersection
    intersection.s1 == INTESRECTION_Y1 // y coordinate of first intersection
    intersection.s2 == INTERSECTION_Z1 // z coordinate of first intersection
    intersection.s3 == INTERSECTION_T1 // tof to first intersection
    intersection.s4 == INTERSECTION_X2 // x coordinate of second intersection
    intersection.s5 == INTERSECTION_Y2 // y coordinate of second intersection
    intersection.s6 == INTERSECTION_Z2 // z coordinate of second intersection
    intersection.s7 == INTERSECTION_T2 // tof to second intersection
    intersection.s012 == INTERSECTION_POS1 // float3 of coordinates of first intersection
    intersection.s456 == INTERSECTION_POS2 // float3 of coordinates of second intersection

These macros require that the intersection buffer entry is retrieved into a ``float8`` named ``intersection`` at the start of a kernel that uses them

.. code-block:: C

    uint global_addr = get_global_id(0);
    float8 intersection = intersections[global_addr];

Altering the data in the intersection buffer is the role of geometry kernels.

The index buffer
----------------

The index buffer holds a ``uint`` entry for each simulated trajectory. Each entry corresponds to an index for the earliest intersected component at each scattering step for each trajectory. This allows the execution of the appropriate kernel for each trajectory.

Buffer manipulation conventions
-------------------------------

At the end of each scattering kernel, the intersection buffer entries are reset to ```(0.0, 0.0, 0.0, 1e5, 0.0, 0.0, 0.0, 1e5)```. At each scattering step, if a trajectory retains this intersection entry after all the geometry kernels have been executed (i.e., it does not intersect with the instrument geometry) it is terminated.

Geometry kernels corresponding to a 'flat' geometry such as a plane, where there is only a single point of intersection, should store identical entries in the first and second halves of the intersection buffer entry. See the ``GPlane`` kernel for reference.

Because negative time intersections are typically discarded, 'interior' geometry of shapes, such as the interior face of a sphere, should be treated as flat geometries are. That is, if one wishes to intersect with the interior of a sphere, the kernel should carry out a ray-sphere intersection and store the second (positive time) intersection with the sphere in both halves of the intersection buffer entry.