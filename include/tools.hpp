#ifndef TOOLS_H_
#define TOOLS_H_

#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <loos.hpp>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

constexpr unsigned int str2int(const char* str, int h);

double deduceMass(const std::string name);

Eigen::Vector3d gcoord_to_eigenv(const loos::GCoord in);

Eigen::Matrix3d gmatrix_to_eigenm(const loos::GMatrix in);

Eigen::Matrix3d inertia_tensor(const loos::AtomicGroup group);

void writexyz(const loos::AtomicGroup group, std::ostream &output);

// Eigen::Matrix3d align(loos::AtomicGroup ref, loos::AtomicGroup mol, Eigen::Matrix<double,Eigen::Dynamic,3> w);
#endif

