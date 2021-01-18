#ifndef MISC_H
#define MISC_H
#include <RcppEigen.h>
#include <vector>
#include <algorithm>

std::vector<int> primes(int n);
void primes(int n, int sz, int *primeVec);
void mean_std (const std::vector<double> &v, double &mean, double &stdev);
void mean_std (int n, const double *v, double &mean, double &stdev);
void reorder_row(Eigen::MatrixXd &A, const int *idx, int *workInt, int lworkInt, int 
	ncol = -1);

/*
        array[i] = array[idx[i]]
        workInt should be have dimensions no smaller than n
	added workInt checker
        2019/09/15
*/
template<class ForwardIt>
void reorder(ForwardIt array, const int *idx, int n, int *workInt, int lworkInt)
{
	if(lworkInt < n) Rcpp::stop("Dimension of workInt is insufficient\n");
        std::copy_n(idx, n, workInt);
        for(int i = 0; i < n; ++i)
        {
                int oldIdx = i;
                while(workInt[oldIdx] != i)
                {
                        std::swap(array[oldIdx], array[workInt[oldIdx]]);
                        std::swap(oldIdx, workInt[oldIdx]);
                }
                workInt[oldIdx] = oldIdx;
        }
}

#endif
