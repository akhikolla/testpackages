#ifndef TLRMVN_H
#define TLRMVN_H
#include <RcppEigen.h>
int tlrmvn(int N, const std::vector<Eigen::MatrixXd> &B,
        const std::vector<TLRNode> &UV, const Eigen::VectorXd &a1,
        const Eigen::VectorXd &b1, double &v, double &e, int ns, int &scaler_in,
        double *workDbl, int lworkDbl, int *workInt, int lworkInt);
#endif
