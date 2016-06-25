#include "INS_transformation.h"


INS_transformation::INS_transformation()
{
    a_bias.setZero();
    w_bias.setZero();
    g=9.8116;
    ts=0.00401;
    sigma_a = 0.03;
    sigma_g = math_radians * (0.43333);
    window_size = 3;
    gamma = 30000;
    P_pos = 0.00001;
    P_vel = 0.00001;
    P_att = math_radians * (0.1);
    P_acc_bias = 0.005;
    P_gyro_bias = math_radians * (0.003);
    Q_acc_bias_noise=0.0000001;
    Q_gyro_bias_noise=math_radians * (0.0000001);
    Q_acc = 0.16;
    Q_gyro =math_radians * (1.6);

    P.block(9, 9, 3, 3) << I3*(P_acc_bias*P_acc_bias);
    P.block(12, 12, 3, 3) << I3*(P_gyro_bias*P_gyro_bias);
    P.block(0, 0, 3, 3) << I3*(P_pos*P_pos);
    P.block(3, 3, 3, 3) << I3*(P_vel*P_vel);
    P.block(6, 6, 3, 3) << I3*(P_att*P_att);
    Q.block(9, 9, 3, 3) << I3*(Q_gyro_bias_noise*Q_gyro_bias_noise);
    Q.block(0, 0, 3, 3) << I3*(Q_acc*Q_acc);
    Q.block(3, 3, 3, 3) << I3*(Q_gyro*Q_gyro);
    Q.block(6, 6, 3, 3) << I3*(Q_acc_bias_noise*Q_acc_bias_noise);
    R *= R_vel;
    F.setZero();
    G.setZero();
    x_state.setZero();

}
void INS_transformation::state_update(Matrix<double, 3, 1> & ar, Matrix<double, 3, 1> & wr){
    Vector3d w_tb = wr*ts;
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
inline void INS_transformation::correct_state(){
    x_state += dx.block<6, 1>(0, 0);
    Matrix3d fii; fii << 0, -dx(8), dx(7), dx(8), 0, -dx(6), -dx(7), dx(6), 0;
    cbn = (I3 - fii) * cbn;
    cbn2quat(cbn, quat);
}
inline void INS_transformation::state_transmision_matrix(){
    Matrix3d St; St << 0, -an(2), an(1), an(2), 0, -an(0), -an(1), an(0);
    //G << O3, O3, O3, O3, cbn, O3, O3, O3, O3, -cbn, O3, O3, O3, O3, I3, O3, O3, O3, O3, I3;
    G.block(3, 0, 3, 3) << cbn; G.block(6, 3, 3, 3) << -cbn;
    G.block(9, 6, 3, 3) << I3; G.block(12, 3, 9, 3) << -cbn;
    F.block(0, 3, 3, 3) << I3; G.block(3, 9, 3, 3) << cbn;
    F.block(3, 6, 3, 3) << St; G.block(6, 12, 3, 3) << -cbn;
    F = I15 + F * ts;
    G *= ts;
}
void INS_transformation::bias_calibration(Matrix<double, 3, 1> & ar, Matrix<double, 3, 1> & wr){}
void INS_transformation::start_transformation(std::string file_name){
    ifstream imu_file(file_name);
    ofstream in;
    in.open("result.txt");
    string line;

    if( !imu_file.good() ){
        cerr << "no imu file found at" << file_name;
        return;
    }
    getline(imu_file, line);
    int count = 0;

    while(!imu_file.eof()){
        count ++;
        if( !getline(imu_file, line)){
            cout << "Finished.";
            return;
        }
        stringstream stream(line);
        string s;
        int calibr_count = int(calibr_time / ts);
        getline(stream, s, ',');
        in << s << ',';
        for (int j = 0; j < 3; ++j){
            getline(stream, s, ',');
            wr(j) = stof(s);
        }
        for (int j = 0; j < 3; ++j){
            getline(stream, s, ',');
            ar(j) = stof(s);
        }
        if( count < calibr_count){
            bias_calibration(ar, wr);
        }
        else
        {
            ar -= a_bias; wr -= a_bias;
            Kalman_transformation();
            in << x_state(0) << ',' << x_state(1) << ',' << x_state(2)
                << ',' << cbn(0) << ',' << cbn(1) << ',' << cbn(2) << ',' << cbn(3) << ',' << cbn(4)
                << ',' << cbn(5) << ',' << cbn(6) << ',' << cbn(7) << ',' << cbn(8);
        }
    }


}
void INS_transformation::start_transformation(){}
void INS_transformation::Kalman_transformation(){
    ar += dx.block(9, 0, 3, 1); wr += dx.block(12, 0, 3, 1);
    state_update(ar, wr);
    state_transmision_matrix();
    P = F * P * F.transpose() + G * Q * G.transpose();
    P = (P + P.transpose()) * 0.5;
    zupt_detector();
    if(zupt_state){
        K = P.block(0, 3, 15, 3) * (P.block(3, 3 ,3, 3) + R).inverse();
        dx = K * (-x_state.block(3, 0, 3, 1));
        correct_state();
        Matrix15d tempK = Matrix15d::Zero();
        tempK.block(0, 3, 15, 3)  = K;
        P -= tempK*P;
        P = (P + P.transpose()) * 0.5;

    }
}
void INS_transformation::zupt_detector(){}
INS_transformation::~INS_transformation()
{

}
