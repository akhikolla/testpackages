#ifndef SVD_H
#define SVD_H
#include <RcppEigen.h>
int svd(Eigen::MatrixXd &A , int ARowNum , int AColNum , Eigen::MatrixXd &U ,
        Eigen::MatrixXd &VT , Eigen::VectorXd &S , double *work , int lwork);
#endif
