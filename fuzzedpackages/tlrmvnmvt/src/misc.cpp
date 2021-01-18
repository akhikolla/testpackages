#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <functional>
#include "misc.h"

using namespace std;
using namespace Eigen;

void mean_std (const std::vector<double> &v, double &mean, double &stdev)
{
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        mean = sum / v.size();

        std::vector<double> diff(v.size());
        std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        stdev = std::sqrt(sq_sum / v.size());
}

void mean_std (int n, const double *v, double &mean, double &stdev)
{
        double sum = std::accumulate(v, v+n, 0.0);
        mean = sum / (double) n;

        std::vector<double> diff(n);
        std::transform(v, v+n, diff.begin(), [mean](double x) { return x - mean; });
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        stdev = std::sqrt(sq_sum / (double) n);
}

std::vector<int> primes(int n)
{
        std::vector<int> v;
        if(n > 2)
        {
                v.push_back(2);
                for(int i = 3; i <= n; i++)
                {
                        int sqroot = sqrt(i);
                        bool prime = true;
                        for(auto j = v.begin(); j != v.end(); j++)
                        {
                                if(*(j) > sqroot)
                                        break;
                                if(i % *(j) == 0)
                                {
                                        prime = false;
                                        break;
                                }
                        }
                        if(prime)
                                v.push_back(i);
                }
        }
        return v;
}

/*
	exit on i == n || idx == sz
	2019/09/25
*/
void primes(int n, int sz, int *primeVec)
{
	int idx = 0;
        if(n > 2 && sz > 0)
        {
                primeVec[idx] = 2;
		idx++;
		if(idx == sz)
			return;
                for(int i = 3; i <= n; i++)
                {
                        int sqroot = sqrt(i);
                        bool prime = true;
                        for(int j = 0; j < idx; j++)
                        {
                                if(primeVec[j] > sqroot)
                                        break;
                                if(i % primeVec[j] == 0)
                                {
                                        prime = false;
                                        break;
                                }
                        }
                        if(prime)
			{
                                primeVec[idx] = i;
				idx++;
				if(idx == sz)
					return;
			}
                }
        }
}

/*
	workInt should be no smaller than A.rows()
	added lworkInt check
	2019/09/15
*/
void reorder_row(Eigen::MatrixXd &A, const int *idx, int *workInt, int lworkInt, int 
	ncol)
{
	int nrow = A.rows();
	if(lworkInt < nrow) Rcpp::stop("Dimension of workInt is insufficient\n");
	if(ncol < 0) ncol = A.cols();
	for(int j = 0; j < ncol; ++j)
	{
		double *ptrCol = A.col(j).data();
		reorder(ptrCol, idx, nrow, workInt, lworkInt);
	}
}
