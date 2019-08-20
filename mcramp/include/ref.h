float reflectivity_func(float q, float R0,
    float Qc, float alpha, float m,
    float W) {

    /* Simpler parametrization from Henrik Jacobsen uses these values that depend on m only.
       float m_value=m*0.9853+0.1978;
       float W=-0.0002*m_value+0.0022;
       float alpha=0.2304*m_value+5.0944;
       float beta=-7.6251*m_value+68.1137; 
       If W and alpha are set to 0, use Henrik's approach for estimating these parameters
       and apply the formulation:
       arg = R0*0.5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc));

       Above lifted from the mcstas reflectivity implementation
    */  

    float beta, arg, R = 0.;

    if (W==0. && alpha==0.) {
        m     = m*0.9853+0.1978;
        W     = -0.0002*m+0.0022;
        alpha = 0.2304*m+5.0944;
        beta  = -7.6251*m+68.1137;

        if (m<=3) {
            alpha = m;
            beta = 0;
        }
    }

    arg = W > 0 ? (q - m*Qc)/W : 11;

    if ( arg > 10 || m <= 0 || Qc <= 0 || R0 <= 0) {
        R = 0.;
        return R;
    }

    if (m < 1) {
        Qc *= m;
        m   = 1;
    }

    if (q <= Qc) {
        R = R0;
        return R;
    }

    R = R0*0.5*(1 - tanh(arg))*(1 - alpha*(q - Qc) + 
                 beta*(q - Qc)*(q - Qc));
    
    return R;
}