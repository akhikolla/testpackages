#ifndef RECUR_BLK_REORDER_H
#define RECUR_BLK_REORDER_H
#include <RcppEigen.h>
#include <vector>
#include "tlr.h"
int recur_blk_reorder(std::vector<Eigen::MatrixXd> &B, std::vector<TLRNode> &UV,
        double *a, double *b, double *p, double *y, int *oldIdx, double epsl,
        double *workDbl, int lworkDbl, int *workInt, int lworkInt);
#endif
