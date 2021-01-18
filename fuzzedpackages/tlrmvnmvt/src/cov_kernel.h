#ifndef COV_KERNEL_H
#define COV_KERNEL_H
#include <RcppEigen.h>
#include <functional>

std::function<double(double)> exp_kernel(double beta);
std::function<double(double)> matern_kernel(double rho, double nu);

#endif
