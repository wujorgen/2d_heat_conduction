#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "momentum_eqn.hpp"
#include "corrections.hpp"
#include "boundary.hpp"

using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;

int main()
{
    int NX = 6;
    double LX = 1.0;

    int NY = 6;
    double LY = 1.0;

    double dx = LX / (NX - 1);
    double dy = LY / (NY - 1);

    // x and y locations
    VectorXd x = VectorXd::LinSpaced(NX, 0, LX);
    VectorXd y = VectorXd::LinSpaced(NY, 0, LY);

    // note that the matrices are indexed as:
    // 0-------j
    // |
    // | 
    // i

    // while the grid is indexed as
    // Y
    // | 
    // | 
    // 0-------X

    // The grid defines vertices.

    // PRESSURE ghost grids: (NY + 1, NX + 1)
    // The pressure ghost grid is fully staggered from the vertices.
    ArrayXXd p = ArrayXXd::Ones(NY + 1, NX + 1);
    ArrayXXd p_star = ArrayXXd::Ones(NY + 1, NX + 1);
    ArrayXXd p_corr = ArrayXXd::Ones(NY + 1, NX + 1);
    ArrayXXd p_b = ArrayXXd::Ones(NY + 1, NX + 1);

    // u grid: (NY + 1, NX)
    // staggered only in y-direction
    // note that row 0 (top) and row NY + 1 (bottom) must be zero, as those are ghost velocities.
    ArrayXXd u = ArrayXXd::Zero(NY + 1, NX);
    ArrayXXd u_star = ArrayXXd::Zero(NY + 1, NX);
    ArrayXXd u_corr = ArrayXXd::Zero(NY + 1, NX);
    ArrayXXd d_e = ArrayXXd::Zero(NY + 1, NX);

    // v grid: (NY, NX + 1)
    // staggered only in x-direction
    // note that col 0 (left) and col NX + 1 (right) must be zero, as those are ghost velocities.
    ArrayXXd v = ArrayXXd::Zero(NY, NX + 1);
    ArrayXXd v_star = ArrayXXd::Zero(NY, NX + 1);
    ArrayXXd v_corr = ArrayXXd::Zero(NY, NX + 1);
    ArrayXXd d_n = ArrayXXd::Zero(NY, NX + 1);

    // Boundary Conditions
    double U_LID = 1;
    u(0, all) = U_LID;

    // Fluid Properties
    double mu = 1;  // 0.0010518  # dynamic viscosity, Pa*s
    double rho = 1000;  // 1000  # density, kg/m^3

    // Convergence
    double alpha = 0.8;
    double alpha_p = 0.5;

    cout << p << endl;

    EigenPtr(&p);

    cout << "asdfasdfasdf" << endl;
    cout << p << endl;

    return 0;
}