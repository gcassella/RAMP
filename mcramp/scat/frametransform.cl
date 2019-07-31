__kernel void transform(__global float16* neutrons,
    float3 const pos, float3 const rot) {
        
    // Rotates neutrons in the xz plane

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.) {
        return;
    }

    float x_a = rot.s0;
    float y_a = rot.s1;
    float z_a = rot.s2;

    float r00, r01, r02, r10, r11, r12, r20, r21, r22;

    r00 = cos(y_a)*cos(z_a);
    r01 = -cos(x_a)*sin(z_a) + sin(x_a)*sin(y_a)*cos(z_a);
    r02 = sin(x_a)*sin(z_a) + cos(x_a)*sin(y_a)*cos(z_a);
    r10 = cos(y_a)*sin(z_a);
    r11 = cos(x_a)*cos(z_a) + sin(x_a)*sin(y_a)*sin(z_a);
    r12 = -sin(x_a)*cos(z_a) + cos(x_a)*sin(y_a)*sin(z_a);
    r20 = -sin(y_a);
    r21 = sin(x_a)*cos(y_a);
    r22 = cos(x_a)*cos(y_a);

    float3 row1, row2, row3;

    row1 = (float3){ r00, r01, r02 };
    row2 = (float3){ r10, r11, r12 };
    row3 = (float3){ r20, r21, r22 };

    neutron.s012 -= pos;

    float x, y, z, vx, vy, vz;

    x = dot(row1, neutron.s012);
    y = dot(row2, neutron.s012);
    z = dot(row3, neutron.s012);
    vx = dot(row1, neutron.s345);
    vy = dot(row2, neutron.s345);
    vz = dot(row3, neutron.s345);

    neutron.s012 = (float3){ x, y, z };
    neutron.s345 = (float3){ vx, vy, vz};

    neutrons[global_addr] = neutron;
}

__kernel void untransform(__global float16* neutrons,
    float3 const pos, float3 const rot) {
        
    // Rotates neutrons in the xz plane

    uint global_addr = get_global_id(0);

    float16 neutron = neutrons[global_addr];

    /* Already terminated? */
    if (neutron.sf > 0.) {
        return;
    }

    float x_a = rot.s0;
    float y_a = rot.s1;
    float z_a = rot.s2;

    float r00, r01, r02, r10, r11, r12, r20, r21, r22;

    r00 = cos(y_a)*cos(z_a);
    r01 = -cos(x_a)*sin(z_a) + sin(x_a)*sin(y_a)*cos(z_a);
    r02 = sin(x_a)*sin(z_a) + cos(x_a)*sin(y_a)*cos(z_a);
    r10 = cos(y_a)*sin(z_a);
    r11 = cos(x_a)*cos(z_a) + sin(x_a)*sin(y_a)*sin(z_a);
    r12 = -sin(x_a)*cos(z_a) + cos(x_a)*sin(y_a)*sin(z_a);
    r20 = -sin(y_a);
    r21 = sin(x_a)*cos(y_a);
    r22 = cos(x_a)*cos(y_a);

    float3 row1, row2, row3;

    row1 = (float3){ r00, r10, r20 };
    row2 = (float3){ r01, r11, r21 };
    row3 = (float3){ r02, r12, r22 };

    float x, y, z, vx, vy, vz;

    x = dot(row1, neutron.s012);
    y = dot(row2, neutron.s012);
    z = dot(row3, neutron.s012);
    vx = dot(row1, neutron.s345);
    vy = dot(row2, neutron.s345);
    vz = dot(row3, neutron.s345);

    neutron.s012 = (float3){ x, y, z };
    neutron.s012 += pos;
    neutron.s345 = (float3){ vx, vy, vz};

    neutrons[global_addr] = neutron;
}