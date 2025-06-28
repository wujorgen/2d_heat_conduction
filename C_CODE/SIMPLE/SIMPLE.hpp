#ifndef SIMPLE_HPP
#define SIMPLE_HPP

#ifndef PROBLEMINFO
#include "ProblemInfo.hpp"
#endif

using Eigen::ArrayXXd;
using Eigen::VectorXd;

void SIMPLE(const BoundaryConditions&, const GridInfo&, const ProblemInfo&);

void EigenRef(ArrayXXd&);

#endif