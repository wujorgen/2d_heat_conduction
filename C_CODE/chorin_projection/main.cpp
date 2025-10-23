#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <iomanip>

using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;


ArrayXXd CentralDifferenceX(const ArrayXXd& f, const double& ElementLength)
{
    ArrayXXd diff(f.rows(), f.cols());
    diff.setZero();
    diff(seq(1, last-1), seq(1, last-1)) = (
        f(seq(1, last-1), seq(2, last)) - f(seq(1, last-1), seq(0, last-2))
    ) / (2.0 * ElementLength);

    return diff;
}

ArrayXXd CentralDifferenceY(const ArrayXXd& f, const double& ElementLength)
{
    ArrayXXd diff(f.rows(), f.cols());
    diff.setZero();
    diff(seq(1, last-1), seq(1, last-1)) = (
        f(seq(2, last), seq(1, last-1)) - f(seq(0, last-2), seq(1, last-1))
    ) / (2.0 * ElementLength);

    return diff;
}

ArrayXXd Laplace(const ArrayXXd& f, const double& ElementLength)
{
    ArrayXXd diff(f.rows(), f.cols());
    diff.setZero();
    diff(seq(1, last-1), seq(1, last-1)) = (
        f(seq(1, last-1), seq(0, last-2)) +
        f(seq(0, last-2), seq(1, last-1)) +
        f(seq(1, last-1), seq(2, last)) + 
        f(seq(2, last), seq(1, last-1)) -
        4 * f(seq(1, last-1), seq(1, last-1))
    ) / (pow(ElementLength, 2));

    return diff;
}

bool writeArrayToCSV(const ArrayXXd& array, 
                     const string& filename,
                     char delimiter = ',',
                     int precision = 6) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return false;
    }
    
    // Set precision for floating point output
    file << std::fixed << std::setprecision(precision);
    
    // Write array data row by row
    for (int row = 0; row < array.rows(); ++row) {
        for (int col = 0; col < array.cols(); ++col) {
            file << array(row, col);
            
            // Add delimiter except for last column
            if (col < array.cols() - 1) {
                file << delimiter;
            }
        }
        // Add newline after each row
        file << '\n';
    }
    
    file.close();
    return true;
}

int main()
{
    int N_POINTS = 41;
    double DOMAIN_SIZE = 1.0;
    int N_ITERATIONS = 500;  // for explicit timestep
    double TIME_STEP_LENGTH = 0.001;
    double KINEMATIC_VISCOSITY = 0.1;
    double DENSITY = 1.0;
    double HORIZONTAL_VELOCITY_TOP = 2.0;  // in the positive x direction
    int N_PRESSURE_POISSON_ITERATIONS = 100;

    double ELEMENT_LENGTH = DOMAIN_SIZE / (N_POINTS - 1);

    // x and y locations
    VectorXd x = VectorXd::LinSpaced(N_POINTS, 0, DOMAIN_SIZE);
    VectorXd y = VectorXd::LinSpaced(N_POINTS, 0, DOMAIN_SIZE);

    // meshgrid
    ArrayXXd X(N_POINTS, N_POINTS);
    ArrayXXd Y(N_POINTS, N_POINTS);
    for (int i = 0; i < N_POINTS; ++i) {
        for (int j = 0; j < N_POINTS; ++j) {
            X(i, j) = x(j);
            Y(i, j) = y(i);
        }
    }

    ArrayXXd u_prev(X.rows(), X.cols());
    ArrayXXd v_prev(X.rows(), X.cols());
    ArrayXXd p_prev(X.rows(), X.cols());
    u_prev.setZero();
    v_prev.setZero();
    p_prev.setZero();

    ArrayXXd u_next(X.rows(), X.cols());
    ArrayXXd v_next(X.rows(), X.cols());
    ArrayXXd p_next(X.rows(), X.cols());
    u_next.setZero();
    v_next.setZero();
    p_next.setZero();

    ArrayXXd u_tent(X.rows(), X.cols());
    ArrayXXd v_tent(X.rows(), X.cols());
    u_tent.setZero();
    v_tent.setZero();

    ArrayXXd d_u_prev__d_x(X.rows(), X.cols());
    ArrayXXd d_u_prev__d_y(X.rows(), X.cols());
    ArrayXXd d_v_prev__d_x(X.rows(), X.cols());
    ArrayXXd d_v_prev__d_y(X.rows(), X.cols());

    ArrayXXd laplace__u_prev(X.rows(), X.cols());
    ArrayXXd laplace__v_prev(X.rows(), X.cols());

    ArrayXXd d_u_tent__d_x(X.rows(), X.cols());
    ArrayXXd d_v_tent__d_y(X.rows(), X.cols());

    ArrayXXd rhs(X.rows(), X.cols());

    ArrayXXd d_p_next__d_x(X.rows(), X.cols());
    ArrayXXd d_p_next__d_y(X.rows(), X.cols());

    for (int itr = 0; itr < N_ITERATIONS; itr++)
    {
        // cout << "itr = " << itr << endl;
        d_u_prev__d_x = CentralDifferenceX(u_prev, ELEMENT_LENGTH);
        d_u_prev__d_y = CentralDifferenceY(u_prev, ELEMENT_LENGTH);
        d_v_prev__d_x = CentralDifferenceX(v_prev, ELEMENT_LENGTH);
        d_v_prev__d_y = CentralDifferenceY(v_prev, ELEMENT_LENGTH);
        laplace__u_prev = Laplace(u_prev, ELEMENT_LENGTH);
        laplace__v_prev = Laplace(v_prev, ELEMENT_LENGTH);

        // Solve tentative velocity field w/o pressure gradient
        u_tent = (
            u_prev
            + TIME_STEP_LENGTH * (
                -(u_prev * d_u_prev__d_x + v_prev * d_u_prev__d_y)
                + KINEMATIC_VISCOSITY * laplace__u_prev
            )
        );
        v_tent = (
            v_prev
            + TIME_STEP_LENGTH * (
                -(u_prev * d_v_prev__d_x + v_prev * d_v_prev__d_y)
                + KINEMATIC_VISCOSITY * laplace__v_prev
            )
        );

        // Enforce boundary conditions
        // Homogenous dirichlet BC everywhere, except for horizontal velocity at the top
        u_tent(0,   all ) = 0.0;  // bottom 
        u_tent(all, 0   ) = 0.0;  // left
        u_tent(all, last) = 0.0;  // right
        u_tent(last, all) = HORIZONTAL_VELOCITY_TOP;  // top
        v_tent(0,    all) = 0.0;  // bottom 
        v_tent(all, 0   ) = 0.0;  // left
        v_tent(all, last) = 0.0;  // right
        v_tent(last, all) = 0.0;  // top

        // terms for divergence of field
        d_u_tent__d_x = CentralDifferenceX(u_tent, ELEMENT_LENGTH);
        d_v_tent__d_y = CentralDifferenceY(v_tent, ELEMENT_LENGTH);

        // pressure correction from solving pressure poisson eqn
        rhs = (DENSITY / TIME_STEP_LENGTH) * (d_u_tent__d_x + d_v_tent__d_y);
        // cout << "Max RHS: " << rhs.abs().maxCoeff() << endl;
        //p_prev.setZero();
        p_next.setZero();
        for (int itr_p = 0; itr_p < N_PRESSURE_POISSON_ITERATIONS; itr_p++)
        {
            p_next.setZero();
            p_next(seq(1, last-1), seq(1, last-1)) = 0.25 * (
                p_prev(seq(1, last-1), seq(2, last)) +  // right
                p_prev(seq(1, last-1), seq(0, last-2)) +  // left
                p_prev(seq(2, last), seq(1, last-1)) +  // top
                p_prev(seq(0, last-2), seq(1, last-1)) -  // bottom
                ELEMENT_LENGTH * ELEMENT_LENGTH * rhs(seq(1, last-1), seq(1, last-1))
            );

            // Pressure BCs
            p_next(all, last) = p_next(all, last - 1);  // right Neumann
            p_next(0, all) = p_next(1, all);            // bottom Neumann
            p_next(all, 0) = p_next(all, 1);            // left Neumann
            p_next(last, all).setZero();                // top Dirichlet

            ArrayXXd p_prev = p_next;
        }


        // pressure gradient
        d_p_next__d_x = CentralDifferenceX(p_next, ELEMENT_LENGTH);
        d_p_next__d_y = CentralDifferenceY(p_next, ELEMENT_LENGTH);

        // Velocity corrections (for incompressible fluid!)
        u_next = (
            u_tent - TIME_STEP_LENGTH / DENSITY * d_p_next__d_x
        );
        v_next = (
            v_tent - TIME_STEP_LENGTH / DENSITY * d_p_next__d_y
        );

        // Enforce boundary conditions (again)
        // Homogenous dirichlet BC everywhere, except for horizontal velocity at the top
        u_next(0,   all ) = 0.0;  // bottom 
        u_next(all, 0   ) = 0.0;  // left
        u_next(all, last) = 0.0;  // right
        u_next(last, all) = HORIZONTAL_VELOCITY_TOP;  // top
        v_next(0,    all) = 0.0;  // bottom 
        v_next(all, 0   ) = 0.0;  // left
        v_next(all, last) = 0.0;  // right
        v_next(last, all) = 0.0;  // top

        // 
        u_prev = u_next;
        v_prev = v_next;
        p_prev = p_next;

    }

    writeArrayToCSV(p_next, "test_pressure_sol.csv");
    writeArrayToCSV(u_next, "test_u_sol.csv");
    writeArrayToCSV(v_next, "test_v_sol.csv");

    return 0;
}