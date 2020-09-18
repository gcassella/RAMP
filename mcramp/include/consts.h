/* Scientific constants */

#ifndef k_B
#define k_B (1.38066e-23f)
#endif

#ifndef M_PI
#define M_PI (3.14159265f)
#endif

#ifndef V2K
#define V2K (1.58825361e-3f)
#endif

#ifndef K2V
#define K2V (629.622368f)
#endif

#ifndef MIN2RAD
#define MIN2RAD (M_PI/(180.0f*60.0f))
#endif

#ifndef VS2E
#define VS2E     5.22703725e-6f     /* Convert (v[m/s])**2 to E[meV] */
#endif

#ifndef SE2V
#define SE2V     437.393377f        /* Convert sqrt(E)[meV] to v[m/s] */
#endif

#ifndef E2KS
#define E2KS     (1.0f/2.072f)
#endif

#ifndef kB
#define kB       8.6173*10e-2f
#endif


// The following is taken without modification from
// https://streamhpc.com/blog/2013-10-17/writing-opencl-code-single-double-precision/

#if CONFIG_USE_DOUBLE

#if defined(cl_khr_fp64)  // Khronos extension available?
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#elif defined(cl_amd_fp64)  // AMD extension available?
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#endif

#endif // CONFIG_USE_DOUBLE

#if defined(DOUBLE_SUPPORT_AVAILABLE)

// double
typedef double real_t;
typedef double2 real2_t;
typedef double3 real3_t;
typedef double4 real4_t;
typedef double8 real8_t;
typedef double16 real16_t;

#define convert_real_t(X)   convert_double(X)
#define convert_real2_t(X)   convert_double2(X)
#define convert_real3_t(X)   convert_double3(X)
#define convert_real4_t(X)   convert_double4(X)
#define convert_real8_t(X)   convert_double8(X)
#define convert_real16_t(X)   convert_double16(X)

#else

// float
typedef float real_t;
typedef float2 real2_t;
typedef float3 real3_t;
typedef float4 real4_t;
typedef float8 real8_t;
typedef float16 real16_t;

#define convert_real_t(X)   convert_float(X)
#define convert_real2_t(X)   convert_float2(X)
#define convert_real3_t(X)   convert_float3(X)
#define convert_real4_t(X)   convert_float4(X)
#define convert_real8_t(X)   convert_float8(X)
#define convert_real16_t(X)   convert_float16(X)

#endif

#ifndef NEUTRON_MACROS

#define NEUTRON_MACROS
#define NEUTRON_X   neutron.s0
#define NEUTRON_Y   neutron.s1
#define NEUTRON_Z   neutron.s2
#define NEUTRON_POS neutron.s012
#define NEUTRON_VX  neutron.s3
#define NEUTRON_VY  neutron.s4
#define NEUTRON_VZ  neutron.s5
#define NEUTRON_VEL neutron.s345
#define NEUTRON_PX  neutron.s6
#define NEUTRON_PY  neutron.s7
#define NEUTRON_PZ  neutron.s8
#define NEUTRON_POL neutron.s678
#define NEUTRON_P   neutron.s9
#define NEUTRON_TOF neutron.sa
#define NEUTRON_DIE neutron.sf

#define INTERSECTION_X1     intersection.s0
#define INTERSECTION_Y1     intersection.s1
#define INTERSECTION_Z1     intersection.s2
#define INTERSECTION_T1     intersection.s3
#define INTERSECTION_X2     intersection.s4
#define INTERSECTION_Y2     intersection.s5
#define INTERSECTION_Z2     intersection.s6
#define INTERSECTION_T2     intersection.s7
#define INTERSECTION_POS1   intersection.s012
#define INTERSECTION_POS2   intersection.s456

#endif