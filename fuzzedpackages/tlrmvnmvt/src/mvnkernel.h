#ifndef MVNKERNEL_H
#define MVNKERNEL_H
#include <RcppEigen.h>
#include <vector>
int mvndns(int n,  int N, const Eigen::MatrixXd& L, const Eigen::MatrixXd &x,
       const Eigen::MatrixXd &a, const Eigen::MatrixXd &b,
       Eigen::ArrayXd &p, Eigen::MatrixXd &y, int yOffset, std::vector<int>& scaler,
       Eigen::MatrixXd &workM, Eigen::VectorXi &workVecI);
int mvndns(int n,  int N, const Eigen::MatrixXd& L, const Eigen::MatrixXd &x,
       const Eigen::MatrixXd &a, const Eigen::MatrixXd &b,
       double *p, double *y, int ldy, int *scaler,
       double *workDbl, int lworkDbl, int *workInt, int lworkInt);
#endif
