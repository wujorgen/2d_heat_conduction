#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "ProblemInfo.hpp"
#include "MomentumEquations.hpp"

using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;

void EigenRef(ArrayXXd& testarr)
{
    // order matters. even tho this function is declared in header,
    //      it must be defined before use.
    testarr(0, all) = 1234;
}

void SIMPLE(const BoundaryConditions BC, const GridInfo Mesh)
{
    // PRESSURE ghost grids: (NY + 1, NX + 1)
    // The pressure ghost grid is fully staggered from the vertices.
    ArrayXXd p = ArrayXXd::Ones(Mesh.NY + 1, Mesh.NX + 1);
    ArrayXXd p_star = ArrayXXd::Ones(Mesh.NY + 1, Mesh.NX + 1);
    ArrayXXd p_corr = ArrayXXd::Ones(Mesh.NY + 1, Mesh.NX + 1);
    ArrayXXd p_b = ArrayXXd::Ones(Mesh.NY + 1, Mesh.NX + 1);

    // u grid: (NY + 1, NX)
    // staggered only in y-direction
    // note that row 0 (top) and row NY + 1 (bottom) must be zero, as those are ghost velocities.
    ArrayXXd u = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    ArrayXXd u_star = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    ArrayXXd u_corr = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    ArrayXXd d_e = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);

    // v grid: (NY, NX + 1)
    // staggered only in x-direction
    // note that col 0 (left) and col NX + 1 (right) must be zero, as those are ghost velocities.
    ArrayXXd v = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    ArrayXXd v_star = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    ArrayXXd v_corr = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    ArrayXXd d_n = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);

    int itr = 0;
    int maxitr = 1000;
    double error = 1;
    double ethresh = 1e-4;

    vector<double> errors;

    // TODO: timestep loop around this, with conditional to disable transient term
    while (error > ethresh && itr < maxitr)
    {
        // calc u-momentum
        CalcUStar(u, u_star, d_e);

        // apply u-momentum boundary conditions

        // call v-momentum
        CalcVStar(v, v_star, d_n);

        // apply v-momentum boundary conditions

        // zero out pressure corrections

        // calculate pressure corrections

        // apply pressure corrections to pressure field

        // apply pressure boundary conditions

        // correct u-momentum

        // apply u-momentum boundary conditions

        // correct v-momentum

        // apply v-momentum boundary conditions

        // check for convergence
        itr++;
    }
    
    cout << p << endl;
    EigenRef(p);
    cout << p << endl;
}

