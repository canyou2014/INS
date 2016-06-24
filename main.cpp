#include <iostream>
#include <Eigen/Dense>
using namespace std;

int main()
{
    Eigen::Matrix3d m3;
    m3.setIdentity();
    cout << m3;
    cout << "Hello world!" << endl;
    return 0;
}
