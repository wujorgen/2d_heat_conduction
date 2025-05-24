#include <Eigen/Dense>
#include <iostream>

using Eigen::ArrayXXf;
using Eigen::seq;
using Eigen::all;
using namespace std;

int main()
{
    ArrayXXf arr(4,4);
    ArrayXXf brr(4,4);
    arr << 1, 2, 3, 4,
           2, 4, 6, 8,
           3, 6, 9, 12,
           4, 8, 12,16;
    //brr.setOnes();
    brr.setConstant(3);
    cout << "arr = " << endl << arr << endl;

    cout << "arr but sliced = " << endl << arr(seq(0, 2), seq(0, 2)) << endl;

    cout << "brr = " << endl << brr << endl;

    cout << "arr + brr =" << endl << arr + brr << endl;
    cout << "arr * brr (as array type) =" << endl << arr * brr << endl;
    cout << "arr @ brr (as matrix type) =" << endl << arr.matrix() * brr.matrix() << endl;
    cout << "brr @ arr (as matrix type) =" << endl << brr.matrix() * arr.matrix() << endl;

    return 0;
}