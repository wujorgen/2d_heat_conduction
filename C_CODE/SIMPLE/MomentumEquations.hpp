#ifndef MOMENTUMEQUATIONS
#define MOMENTUMEQUATIONS

#include "ProblemInfo.hpp"
#include <Eigen/Dense>

void CalcUStar(
    Eigen::ArrayXXd&, // u_star
    Eigen::ArrayXXd&, // u
    Eigen::ArrayXXd&, // v
    Eigen::ArrayXXd&, // p
    Eigen::ArrayXXd&, // d_e
    const GridInfo,
    const ProblemInfo
);

void CalcVStar(
    Eigen::ArrayXXd&, // v_star
    Eigen::ArrayXXd&, // u
    Eigen::ArrayXXd&, // v
    Eigen::ArrayXXd&, // p
    Eigen::ArrayXXd&, // d_n
    const GridInfo,
    const ProblemInfo
);

#endif