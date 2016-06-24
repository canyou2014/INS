#include "INS_transformation.h"

INS_transformation::INS_transformation()
{

}
void INS_transformation::state_update(){
    Vector3d w_tb = wn*ts;
    Matrix4d ou = Matrix4d::Zero();
    ou(0,1) = -w_tb(0)*0.5; ou(0,2) = -w_tb(1)*0.5; ou(0,3) = -w_tb(2)*0.5;
    ou(1,0) = w_tb(0)*0.5; ou(1,2) = w_tb(2)*0.5; ou(1,3) = -w_tb(1)*0.5;
    ou(2,0) = w_tb(1)*0.5; ou(2,1) = -w_tb(2)*0.5; ou(2,3) = w_tb(0)*0.5;
    ou(3,0) = w_tb(2)*0.5; ou(3,1) = w_tb(1)*0.5; ou(3,2) = -w_tb(0)*0.5;
    double v = w_tb.norm();
    if (v != 0)
    {
        quat = ( cos(v/2)* I4 + 2/v*sin(v/2)*ou ) * quat;
        quat.normalized();
    }
    quat2cbn(quat, cbn);
    an = cbn * ar + g_n;
    Matrix<double, 6, 6> A = Matrix<double, 6, 6>::Zero();
    A(0,3) = ts; A(1,4) = ts; A(2,5) = ts;
    Matrix<double, 6, 3> B = Matrix<double, 6, 3>::Zero();
    B << ts*ts/2 * I3, ts * I3;
    x_state = A * x_state + B * an;
}
inline void INS_transformation::quat2cbn(Matrix<double, 4, 1>& quaternion, Matrix3d& dcm){

    dcm << (quaternion(0)*quaternion(0)+quaternion(1)*quaternion(1)-quaternion(2)*quaternion(2)-quaternion(3)*quaternion(3)),
            (2*(quaternion(1)*quaternion(2)-quaternion(0)*quaternion(3))),
            (2*(quaternion(1)*quaternion(3)+quaternion(0)*quaternion(2))),
            (2*(quaternion(1)*quaternion(2)+quaternion(0)*quaternion(3))),
            (quaternion(0)*quaternion(0)-quaternion(1)*quaternion(1)+quaternion(2)*quaternion(2)-quaternion(3)*quaternion(3)),
            (2*(quaternion(2)*quaternion(3)-quaternion(0)*quaternion(1))),
            (2*(quaternion(1)*quaternion(3)-quaternion(0)*quaternion(2))),
            (2*(quaternion(2)*quaternion(3)+quaternion(0)*quaternion(1))),
            (quaternion(0)*quaternion(0)-quaternion(1)*quaternion(1)-quaternion(2)*quaternion(2)+quaternion(3)*quaternion(3));

    double normq=quaternion(0)*quaternion(0)+quaternion(1)*quaternion(1)+quaternion(2)*quaternion(2)+quaternion(3)*quaternion(3);
    if (normq != 0) dcm /= sqrt(normq);
    else dcm = I3;

}
inline void INS_transformation::cbn2quat(Matrix3d& dcm, Matrix<double, 4, 1>& quaternion){
    double T = 1 + dcm(0,0) + dcm(1,1) + dcm(2,2);
    double S = 0.0;
    if (T > 1e-8)
    {
        S = 0.5 / sqrt(T);
        quaternion(0) = 0.25 / S;
        quaternion(1) = ( dcm(2,1) - dcm(1,2) ) * S;
        quaternion(2) = ( dcm(0,2) - dcm(2,0) ) * S;
        quaternion(3) = ( dcm(1,0) - dcm(0,1) ) * S;
    }
    else if ( (dcm(0,0) > dcm(1,1)) && (dcm(0,0) > dcm(2,2)))
    {
        S = sqrt( 0 + dcm(0,0) - dcm(1,1) - dcm(2,2)) * 2;
        quaternion(0) = (dcm(2,1) - dcm(1,2)) / S;
        quaternion(1) = 0.25 * S;
        quaternion(2) = (dcm(0,1) + dcm(1,0)) / S;
        quaternion(3) = (dcm(0,2) + dcm(2,0)) / S;
    }
    else if (dcm(1,1) > dcm(2,2))
    {
        S = sqrt( 0 + dcm(1,1) - dcm(0,0) - dcm(2,2) ) * 2;
        quaternion(0) = (dcm(0,2) - dcm(2,0)) / S;
        quaternion(1) = (dcm(0,1) + dcm(1,0)) / S;
        quaternion(2) = 0.25 * S;
        quaternion(3) = (dcm(1,2) + dcm(2,1)) / S;
    }
    else
    {
        S = sqrt( 0 + dcm(2,2) - dcm(0,0) - dcm(1,1) ) * 2;
        quaternion(0) = (dcm(1,0) - dcm(0,1)) / S;
        quaternion(1) = (dcm(0,2) + dcm(2,0)) / S;
        quaternion(2) = (dcm(1,2) + dcm(2,1)) / S;
        quaternion(3) = 0.25 * S;
    }
}


INS_transformation::~INS_transformation()
{
    //dtor
}
