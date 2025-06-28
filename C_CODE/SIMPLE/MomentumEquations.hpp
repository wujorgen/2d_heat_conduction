#ifndef MOMENTUMEQUATIONS
#define MOMENTUMEQUATIONS

#include <Eigen/Dense>

void CalcUStar(Eigen::ArrayXXd&, Eigen::ArrayXXd&, Eigen::ArrayXXd&);

void CalcVStar(Eigen::ArrayXXd&, Eigen::ArrayXXd&, Eigen::ArrayXXd&);

#endif