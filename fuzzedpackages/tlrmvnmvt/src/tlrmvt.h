#ifndef TLRMVT_H
#define TLRMVT_H
#include <RcppEigen.h>
#include <vector>
#include "tlr.h"
int tlrmvt(int N, double nu, const Eigen::VectorXd &mu, 
        const std::vector<Eigen::MatrixXd> &B,
        const std::vector<TLRNode> &UV, const Eigen::VectorXd &a1, 
        const Eigen::VectorXd &b1, double &v, double &e, int ns, int &scaler_in,
        double *workDbl, int lworkDbl, int *workInt, int lworkInt);
#endif

