#include <Eigen/Dense>

#include "ProblemInfo.hpp"

using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;

void ApplyUBoundary(ArrayXXd& u, const BoundaryConditions BC)
{
    // left
    u(seq(1, last-1), 0) = BC.U_L;
    // right
    u(seq(1, last-1), last) = BC.U_R;
    // top
    u(0, all) = BC.U_T * 2 - u(1, all);
    // bottom
    u(last, all) = BC.U_B * 2 - u(last-1, all);
}

void ApplyVBoundary(ArrayXXd& v, const BoundaryConditions BC)
{
    // top 
    v(0, seq(1, last-1)) = BC.V_T;
    // bottom
    v(last, seq(1, last-1)) = BC.V_B;
    // left
    v(all, 0) = BC.V_L * 2 - v(all, 1);
    // right
    v(all, last) = BC.V_R * 2 - v(all, last-1);
}

void NoPressureGradientAtBoundary(ArrayXXd& p)
{
    // top
    p(0, all) = p(1, all);
    // bottom
    p(last, all) = p(last-1, all);
    // left
    p(all, 0) = p(all, 1);
    // right
    p(all, last) = p(all, last-1);
}
