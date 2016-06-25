#ifndef INS_TRANSFORMATION_H
#define INS_TRANSFORMATION_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
using namespace Eigen;
using Matrix15d = Matrix<double, 15, 15>;
using Matrix12d = Matrix<double ,12, 12>;
using Vector6d = Matrix<double, 6, 1>;
class INS_transformation
{
    public:
        INS_transformation();
        void set_motion(int motion);
        void start_transformation();
        void start_transformation(std::string file);
        void state_update(Matrix<double, 3, 1> & ar, Matrix<double, 3, 1> & wr);
        virtual ~INS_transformation();
        void Kalman_transformation();
    protected:
    private:
        inline void att2cbn(Matrix<double, 3, 1>& rot, Matrix3d& dcm);
        inline void cbn2quat(Matrix3d& dcm, Matrix<double, 4, 1>& quaternion);
        inline void quat2cbn(Matrix<double, 4, 1>& quaternion, Matrix3d& dcm);
        inline void correct_state();
        inline void state_transmision_matrix();
        void zupt_detector();
        void bias_calibration(Matrix<double, 3, 1> & ar, Matrix<double, 3, 1> & wr);
        Matrix3d cbn;
        Matrix3d R = Matrix3d::Identity();
        const Matrix3d I3 = Matrix3d::Identity();
        Matrix15d P, F;
        const Matrix15d I15 = Matrix15d::Identity();
        const Matrix3d O3 = Matrix3d::Zero();
        Matrix12d Q;
        const double math_radians = 0.01745329252;
        Matrix<double, 15, 12> G;
        Matrix<double, 15, 1> dx;
        Matrix<double, 3, 15> H;
        Matrix<double, 3, 1> ar, wr, an, wn, att;
        Matrix<double, 4, 1> quat;
        Vector6d x_state;
        const Matrix4d I4 = Matrix4d::Identity();
        double g, ts, sigma_a, sigma_g, window_size, gamma, R_vel,
                P_pos, P_vel, P_att, P_acc_bias, P_gyro_bias, Q_acc_bias_noise,
                Q_gyro_bias_noise, Q_acc, Q_gyro;
        const Matrix<double, 3, 1> g_n{0, 0, g};
        const int calibr_time = 3;
        Matrix<double ,3 ,1> a_bias, w_bias;
        int zupt_state;
        Matrix<double ,15, 3> K;
};

#endif // INS_TRANSFORMATION_H
