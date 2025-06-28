#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include "Boundary.hpp"
#include "SIMPLE.hpp"

using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;

int main()
{
    int NX = 101;
    double LX = 0.5;

    int NY = 101;
    double LY = 0.3;

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

    // Boundary Conditions
    double U_LID = 2;

    // Fluid Properties
    double mu = 1;  // 0.0010518  # dynamic viscosity, Pa*s
    double rho = 50;  // 1000  # density, kg/m^3

    // Convergence
    double alpha = 0.5;
    double alpha_p = 0.3;

    ProblemInfo Problem;
    Problem.mu = mu;
    Problem.rho = rho;
    Problem.relax = alpha;
    Problem.relaxp = alpha_p;

    GridInfo Mesh;
    Mesh.NX = NX;
    Mesh.NY = NY;
    Mesh.LX = LX;
    Mesh.LY = LY;
    Mesh.dx = dx;
    Mesh.dy = dy;
    Mesh.x = x;
    Mesh.y = y;

    BoundaryConditions BC;
    BC.U_T = U_LID;

    SIMPLE(BC, Mesh, Problem);

    return 0;
}