void rotate_about_axis(float angle,
    float3 axis, float3* vect) {

    float Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz;

    Rxx = cos(angle)+axis.s0*axis.s0*(1-cos(angle));
    Rxy = axis.s0*axis.s1*(1-cos(angle))-axis.s2*sin(angle);
    Rxz = axis.s0*axis.s2*(1-cos(angle))+axis.s1*sin(angle);
    Ryx = axis.s1*axis.s0*(1-cos(angle))+axis.s2*sin(angle);
    Ryy = cos(angle)+axis.s1*axis.s1*(1-cos(angle));
    Ryz = axis.s1*axis.s2*(1-cos(angle))-axis.s0*sin(angle);
    Rzx = axis.s2*axis.s0*(1-cos(angle))-axis.s1*sin(angle);
    Rzy = axis.s2*axis.s1*(1-cos(angle))+axis.s0*sin(angle);
    Rzz = cos(angle)+axis.s2*axis.s2*(1-cos(angle));

    (*vect) = (float3)( Rxx*(*vect).s0+Rxy*(*vect).s1+Rxz*(*vect).s2, 
                     Ryx*(*vect).s0+Ryy*(*vect).s1+Ryz*(*vect).s2, 
                     Rzx*(*vect).s0+Rzy*(*vect).s1+Rzz*(*vect).s2 );
}

float3 frame_rotate(float3 vec, float3 rot) {
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

    float x, y, z;

    x = dot(row1, vec);
    y = dot(row2, vec);
    z = dot(row3, vec);

    vec = (float3){ x, y, z };

    return vec;
}

float3 frame_derotate(float3 vec, float3 rot) {
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

    float x, y, z;

    x = dot(row1, vec);
    y = dot(row2, vec);
    z = dot(row3, vec);

    vec = (float3){ x, y, z };

    return vec;
}