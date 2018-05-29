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