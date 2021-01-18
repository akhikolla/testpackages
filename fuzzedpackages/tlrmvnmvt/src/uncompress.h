#ifndef UNCOMPRESS_H
#define UNCOMPRESS_H
#include <RcppEigen.h>
#include <vector>
#include "tlr.h"
Eigen::MatrixXd uncompress(const std::vector<Eigen::MatrixXd> &B, const 
	std::vector<TLRNode> &UV, bool symm);
Eigen::MatrixXd uncompress(const std::vector<Eigen::MatrixXd> &B, const
        std::vector<Eigen::MatrixXd> &U, const std::vector<Eigen::MatrixXd> &V, 
        bool symm);
#endif
