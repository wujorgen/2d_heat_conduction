#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "ProblemInfo.hpp"
#include "MomentumEquations.hpp"
#include "Boundary.hpp"
#include "Corrections.hpp"

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

void SIMPLE(const BoundaryConditions& BC, const GridInfo& Mesh, const ProblemInfo& Problem)
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
    int maxitr = 3;
    double error = 1;
    double ethresh = 1e-4;

    vector<double> errors;

    // TODO: timestep loop around this, with conditional to disable transient term
    while (error > ethresh && itr < maxitr)
    {
        // calc u-momentum
        CalcUStar(u_star, u, v, p, d_e, Mesh, Problem); cout << u_star << endl;

        // apply u-momentum boundary conditions
        ApplyUBoundary(u_star, BC);

        // call v-momentum
        CalcVStar(v_star, u, v, p, d_n, Mesh, Problem);

        // apply v-momentum boundary conditions
        ApplyVBoundary(v_star, BC);

        // zero out pressure corrections
        p_corr.setZero();
        p_b.setZero();

        // calculate pressure corrections
        CalcPressureCorrection(p_corr, p_b, u_star, v_star, d_e, d_n, Mesh);

        // apply pressure corrections to pressure field
        ApplyPressureCorrection(p, p_corr, Problem);

        // apply pressure boundary conditions
        NoPressureGradientAtBoundary(p);

        // correct u-momentum
        ApplyUCorrection(u, u_star, d_e, p_corr, Problem);

        // apply u-momentum boundary conditions
        ApplyUBoundary(u, BC);

        // correct v-momentum
        ApplyVCorrection(v, v_star, d_n, p_corr, Problem);

        // apply v-momentum boundary conditions
        ApplyVBoundary(v, BC);

        // check for convergence
        cout << p_b.matrix().norm() << endl;
        itr++;
    }
    
    //cout << p << endl;
    //EigenRef(p);
    //cout << p << endl;
    //cout << u_star << endl;
}

