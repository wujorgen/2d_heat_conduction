#include <Eigen/Dense>

using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;

void EigenPtr(ArrayXXd* testarr)
{
    (*testarr)(0, all) = 1234;
}