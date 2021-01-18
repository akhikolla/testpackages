#ifndef QR_H
#define QR_H
#include <RcppEigen.h>
int qr(Eigen::MatrixXd &A , int ARowNum , int AColNum , Eigen::MatrixXd &R ,
        Eigen::VectorXd &tau , double *work , int lwork);
#endif
