#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>

using Eigen::ArrayXXd;
using namespace std;

void WriteArrayXXdToCSV(const ArrayXXd& arr, const string& fname)
{
    ofstream file(fname);

    for (int i = 0; i < arr.rows(); i++)
    {
        for (int j = 0; j < arr.cols(); j++)
        {
            file << arr(i,j);
            if (j < arr.cols() - 1)
                file << ", ";
        }
        file << endl;
    }

    file.close();
}