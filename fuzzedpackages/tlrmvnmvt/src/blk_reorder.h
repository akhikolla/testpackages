#ifndef BLK_REORDER_H
#define BLK_REORDER_H
#include <RcppEigen.h>
#include <vector>
void blk_reorder(std::vector<Eigen::MatrixXd> &B, double *a, double *b, double *p,
        double *y, int *oldIdx, double *workDbl, int lworkDbl, int *workInt, int
        lworkInt);
#endif
