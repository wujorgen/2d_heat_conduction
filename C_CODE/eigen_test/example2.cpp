#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

int main()
{
	MatrixXd m = MatrixXd::Random(3,3);
	m = (m + MatrixXd::Constant(3,3,1.2)) * 50;
	std::cout << "m =" << std::endl << m << std::endl;
	// to set at compile time:
	// Vector3d v(1,2,3);
	// to set at run time:
	VectorXd v(3);
	v << 1, 2, 3;
	std::cout << "m * v = " << std::endl << m * v << std::endl;

	return 0;
}
