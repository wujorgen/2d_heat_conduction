#ifndef UTILITIES
#define UTILITIES

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <string>

template<typename T>
void WriteVectorToCSV(const std::vector<T>& vec, const std::string& fname)
{
    std::ofstream file(fname);

    // file << "itr, error" << std::endl;

    for (size_t i = 0; i < vec.size(); i++)
    {
        file << vec[i] << std::endl;
        //file << i << ", " << vec[i] << std::endl;
    }

    file.close();
};

void WriteArrayXXdToCSV(const Eigen::ArrayXXd&, const std::string&);

#endif