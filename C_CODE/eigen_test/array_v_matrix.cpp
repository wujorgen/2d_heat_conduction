#include <Eigen/Dense>
#include <iostream>

using Eigen::ArrayXXf;
using Eigen::VectorXd;
using Eigen::seq;
using Eigen::all;
using Eigen::last;
using namespace std;

int main()
{
    VectorXd vec(4);
    ArrayXXf arr(4,4);
    ArrayXXf brr(4,4);
    arr << 1, 2, 3, 4,
           2, 4, 6, 8,
           3, 6, 9, 12,
           4, 8, 12,16;
    vec.setConstant(2);
    //brr.setOnes();
    brr.setConstant(3);
    cout << "arr = " << endl << arr << endl;

    cout << "arr but sliced = " << endl << arr(seq(1, last-1), seq(1, last-1)) << endl;

    cout << "brr = " << endl << brr << endl;

    cout << "arr + brr =" << endl << arr + brr << endl;
    cout << "arr * brr (as array type) =" << endl << arr * brr << endl;
    cout << "arr @ brr (as matrix type) =" << endl << arr.matrix() * brr.matrix() << endl;
    cout << "brr @ arr (as matrix type) =" << endl << brr.matrix() * arr.matrix() << endl;
    cout << "size of arr is " << arr.size() << endl;
    cout << "shape of arr is " << arr.rows() << ", " << arr.cols() << endl;
    cout << "vector vec is " << vec << endl;
    cout << "shape of vec is " << vec.size() << endl;
    cout << VectorXd::LinSpaced(41, 0, 1) << endl;
    return 0;
}