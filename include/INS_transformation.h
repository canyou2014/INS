#ifndef INS_TRANSFORMATION_H
#define INS_TRANSFORMATION_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>
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
        void state_update();
        virtual ~INS_transformation();

    protected:
    private:
        inline void att2cbn(Matrix<double, 3, 1>& rot, Matrix3d& dcm);
        inline void cbn2quat(Matrix3d& dcm, Matrix<double, 4, 1>& quaternion);
        inline void quat2cbn(Matrix<double, 4, 1>& quaternion, Matrix3d& dcm);
        Matrix3d cbn, R;
        const Matrix3d I3 = Matrix3d::Identity();
        Matrix15d P, F;
        const Matrix15d I15 = Matrix15d::Identity();
        Matrix12d Q;
        Matrix<double, 3, 15> H;
        Matrix<double, 3, 1> ar, wr, an, wn, att;
        Matrix<double, 4, 1> quat;
        Vector6d x_state;
        const Matrix4d I4 = Matrix4d::Identity();
        double g, ts, sigma_a, sigma_g, window_size, gama, R_vel,
                P_pos, P_vel, P_att, P_acc_bias, P_gyro_bias, Q_acc_bias_noise,
                Q_gyro_bias_noise, Q_acc, Q_gyro;
        const Matrix<double, 3, 1> g_n{0, 0, g};


};

#endif // INS_TRANSFORMATION_H
