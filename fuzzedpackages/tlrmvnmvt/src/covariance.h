#ifndef COVARIANCE_H
#define COVARIANCE_H
#include <functional>
#include <vector>
#include <RcppEigen.h>
#include "tlr.h"
void tlr_aca_covM(const Eigen::MatrixXd &geom, std::vector<Eigen::MatrixXd> &B,
        std::vector<TLRNode> &UV, std::function<double(double)> kernel, int m, 
	const std::vector<int> &idxVec, double epsl, int allocSz);
void tlr_aca_covM(const Eigen::MatrixXd &covM, std::vector<Eigen::MatrixXd> &B,
        std::vector<TLRNode> &UV, int m, double epsl, int allocSz);
Eigen::MatrixXd dense_covM(const Eigen::MatrixXd &geom, 
        std::function<double(double)> kernel);
#endif
