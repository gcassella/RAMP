Instrument definition file format
=================================

Execution blocks
----------------

Top level objects in RAMP instrument definition files are execution blocks. These \
blocks are used to separate sections of the instrument into linear and non-linear \
execution. The name of execution block objects is arbitrary.

In linear execution blocks, each successive component is executed in turn a single \
time. This is useful when simulating, for example, guide sections, where each neutron \
is expected to proceed in an orderly fashion from one guide component to the next.

In non-linear execution blocks, the geometry kernel of every component is executed \
and the earliest intersected component for each neutron is stored. Then, only the \
scattering kernel for the earliest intersected component acts upon the neutron. This \
is useful when simulating, for example, sample environments, where there is no guarantee \
that neutrons will intersect components in the order they are specified in the instrument \
definition file.

The following is an example of a truncated instrument definition file containing \
two blocks: a linearly executed guide section and a non-linearly executed sample \
environment::

 {
     "guides" : {
         "linear" : true,

         ...

     },

     "sample_env" : {
         "linear" : false

         ...

     }
 }

Components
----------

Below the level of execution blocks come component definitions. An execution block \
contains several components. 

The position and rotation of a component are specified using the `position` and \
`rotation` attributes of the component object. Both attributes are specified as \
an array of three floating point values. For `position` these are the Cartesian \
coordinates of the component in meters and for `rotation` these are the rotation \
angles about the Cartesian axes with the component at the origin, in the order \
z, y, x.

The position and rotation of a component can also be specified relative to a component \
declared earlier in the instrument definition file using the `relative` attribute. \
This is best demonstrated by example - the following is a truncated instrument definition \
file demonstrating a basic monochromator setup using relative rotations::

 {
    "all" : {
        "linear" : true,
        "neutron_source" : {
            "source": true,
            "position" : [0.0, 0.0, 0.0],
            
            ...

        },

        "Mono_arm" : {
            "position" : [0.0, 0.0, 1.0],
            "rotation" : [0.0, 1.0, 0.0],
            "relative" : "neutron_source"
            
            ...

        },

        "Mono_out" : {
            "position" : [0.0, 0.0, 0.0],
            "rotation" : [0.0, 1.0, 0.0],
            "relative" : "Mono_arm",
            
            ...

        },

        "mono" : {
            "position" : [0.0, 0.0, 0.0],
            "relative" : "Mono_arm",
            
            ...

        },

        "sample" : {
            "position" : [0.0, 0.0, 1.0],
            "relative" : "Mono_out",

            ...

        }
    }
 }

Kernels
-------

All of the calculations in RAMP are handled by OpenCL kernels - programs which run \
on OpenCL capable devices. There are three classes of kernel in a RAMP simulation:

 - Moderator kernels
 - Geometry kernels
 - Scattering kernels

Moderator kernels
~~~~~~~~~~~~~~~~~

Typically an instrument will contain a single component which executes a moderator \
kernel, to generate the neutrons at the beginning of the simulation. There is a \
special component level attribute which must be specified for neutron sources: the \
`source` attribute should be set to `true`. For example, the following component \
defines an ISIS style moderator using the `MISIS` moderator kernel::

 "mod" : {
            "source": true,
            "position" : [0.0, 0.0, 0.0],
            "moderator_kernel": {
                "name": "MISIS",
                
                ...

            }
        }

Geometry and scattering kernels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The kernels which dictate how a component interacts with neutrons are specified \
below the level of component objects in the kernel objects `geom_kernel` and `scat_kernel`.

Each component which is not a source should contain a `geom_kernel` and `scat_kernel` \
object. The `name` attribute of these objects specify which kernel the component \
should use. The other attributes of the kernel objects are used to specify the parameters \
of the scattering kernel. For example, one would specify the radius of a spherical \
geometry kernel, or the lattice spacing of a monochromator scattering kernel.

The following is an example of the definition of a flat monochromator::

 "mono" : {
     "position" : [0.0, 0.0, 0.0],
     "geom_kernel" : {
         "name": "GPlane",
         "width": 0.10,
         "height": 0.10,
         "orientation": "yz"
     },
     "scat_kernel": {
         "name": "SMonochromator",
         "slab_width" : 0.1,
         "slab_height" : 0.1,
         "mosaic_horizontal" : 40,
         "mosaic_vertical" : 40,
         "r0" : 0.9,
         "d_spacing" : 3.53,
         "radius_vertical" : 0.0
     }
 }

Variables
---------

It is often inconvenient to edit the instrument definition file every time a component \
parameter such as the angle of a monochromator must be adjusted. To remedy this, \
RAMP supports a non-standard notation for its JSON files to allow variable to be \
set directly from the Python script. Variable names surrounded by \$ signs inside \
instrument definition files will be substituted for by keyword arguments provided \
when the instrument is instantiated in Python.

For example, if the following component was specified inside an instrument definition \
file `inst.json`::

 "mod" : {
     "source": true,
     "position" : [0.0, 0.0, 0.0],
     "moderator_kernel": {
         "name": "MISIS",
         "spec_file": "Let_Base.mcstas",
         "mod_dim": [0.04, 0.09],
         "target_dim": [0.04, 0.09],
         "target_dist": 1.7,
         "E_min": $emin$,
         "E_max": $emax$
     }
 },

The moderator attributes `E_min` and `E_max` could be set when the instrument is \
instantiated to 1.0 and 9.0, respectively, in Python via::

 inst = Instrument('inst.json', ctx, queue, emin=1.0, emax=9.0)

The variable syntax also supports basic arithmetic. After the variable names have \
been substituted for the values specified, the resulting expression within the \$ \
signs is evaluated as a Python expression. For example, if the instrument definition \
file `inst.json` were to contain two choppers with the same constant phase offset \
but different initial phases, this could be specified as follows::

 "Chopper1" : {
     "position" : [0.0, 0.0, 5.0],
     "geom_kernel" : {
         "name" : "GPlane",
         "width" : 0.5,
         "height" : 0.5
     },
     "scat_kernel" : {
         "name" : "SChopper",
         "radius": 0.5,
         "freq" : 314.1,
         "n_slits" : 6,
         "jitter" : 7e-7,
         "slit_width" : 0.04,
         "phase" : $initial_pha_chop1 + pha_offset$
     }
 },

 "Chopper2" : {
     "position" : [0.0, 0.0, 10.0],
     "geom_kernel" : {
         "name" : "GPlane",
         "width" : 0.5,
         "height" : 0.5
     },
     "scat_kernel" : {
         "name" : "SChopper",
         "radius": 0.5,
         "freq" : -314.1,
         "n_slits" : 6,
         "jitter" : 7e-7,
         "slit_width" : 0.04,
         "phase" : $initial_pha_chop2 + pha_offset$
     }
 }

and in the Python script::

 inst = Instrument(
     'inst.json', 
     ctx, 
     queue, 
     initial_pha_chop1 = 0.1,
     initial_pha_chop2 = 0.7,
     pha_offset = 55.0e-3
 )

`NOTE: once variables have been added to an instrument definition file it is no \
longer a strictly valid JSON file, and many programs that interpret JSON files will \
no longer properly load the instrument definition file.`