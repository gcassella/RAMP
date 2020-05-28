Conventions: Buffer structure
=============================

The neutron buffer
------------------

The neutron buffer holds a ``float16`` entry for each simulated trajectory. The elements of these entries correspond to

 - ``(x, y, z)``, the neutron position. In kernels, this is given in the component frame.
 - ``(vx, vy, vz)``, the neutron velocity. Again, in the component frame.
 - ``(sx, sy, sz)``, the neutron spin. As of the current update, unused.
 - ``p``, the neutron weighting i.e. how many counts this trajectory will contribute to a detector.
 - ``t``, time of flight.
 - ``mc``, a counter that is incremented upon each Monte Carlo decision, used in random number generation.
 - ``prevcomp``, the index of the previous component hit. Prevents spurious multiple collisions due to rounding errors.
 - ``(_, _)``, two unused entries. May be used as custom flags / storage by kernels but this should be made extremely explicit in documentation to avoid conflicts.
 - ``T``, termination flag. Trajectory is no longer simulated by kernels once this flag is set not equal to zero.

Altering the data in the neutron buffer is the role of scattering kernels.

The intersection buffer
-----------------------

The intersection buffer holds a ``float8`` entry for each simulated trajectory. The elements of these entries correspond to

 - ``(x1, y1, z1)``, the coordinates (in component frame) of the earliest intersection with the instrument geometry.
 - ``t1``, the time of flight from the neutron position at the previous scattering step to ``(x1, y1, z1)``.
 - ``(x2, y2, z2)``, the coordinates (in component frame) of the second intersection point with the component intersected at ``(x1, y1, z1)``, e.g. where the neutron leaves a component.
 - ``t2``, the time of flight from the neutron position at the previous scattering step to ``(x2, y2, z2)``.

The units of distance are meters and the units of time are seconds.

Altering the data in the intersection buffer is the role of geometry kernels.

The index buffer
----------------

The index buffer holds a ``uint`` entry for each simulated trajectory. Each entry corresponds to an index for the earliest intersected component at each scattering step for each trajectory. This allows the execution of the appropriate kernel for each trajectory.

Buffer manipulation conventions
-------------------------------

At the end of each scattering kernel, the intersection buffer entries are reset to ```(0.0, 0.0, 0.0, 1e5, 0.0, 0.0, 0.0, 1e5)```. At each scattering step, if a trajectory retains this intersection entry after all the geometry kernels have been executed (i.e., it does not intersect with the instrument geometry) it is terminated.

Geometry kernels corresponding to a 'flat' geometry such as a plane, where there is only a single point of intersection, should store identical entries in the first and second halves of the intersection buffer entry. See the ``GPlane`` kernel for reference.

Because negative time intersections are typically discarded, 'interior' geometry of shapes, such as the interior face of a sphere, should be treated as flat geometries are. That is, if one wishes to intersect with the interior of a sphere, the kernel should carry out a ray-sphere intersection and store the second (positive time) intersection with the sphere in both halves of the intersection buffer entry.