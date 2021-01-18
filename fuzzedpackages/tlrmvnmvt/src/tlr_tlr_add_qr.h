#ifndef TLR_TLR_ADD_QR_H
#define TLR_TLR_ADD_QR_H
#include <RcppEigen.h>
#include "tlr.h"
//void tlr_tlr_add_qr(TLRNode &chgNode , const TLRNode &refNode , double epsl , double
//        *work , int lwork , Eigen::MatrixXd &tmpR1 , Eigen::MatrixXd &tmpR2 , Eigen::
//	VectorXd &tau , Eigen::MatrixXd &tmpU , Eigen::MatrixXd &tmpVT);

void tlr_tlr_add_qr(TLRNode &chgNode, const double *U, const double *V, int k2, double
        epsl, double *work, int lwork);
#endif
