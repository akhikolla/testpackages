#include <RcppEigen.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R.h>
#include <algorithm>
#include <cmath>	
#include <ctime>
#include "mvnkernel.h"
#include "misc.h"

using namespace std;
using namespace Eigen;

/*
	lworkDbl should be no smaller than 9*n*N + n + ns + 17*N
	lworkInt should be no smaller than 4*N + n + ns + 1
	2019/09/25
	Removed the step of extracting the exponent of p to avoid extracting the 
		exponent of 0
	2019/09/26
*/
int mvt(int N, double nu, const VectorXd &mu, const MatrixXd& L, 
	const VectorXd &a1, const VectorXd &b1, 
	double& v, double& e, int ns, int &scaler_in, double *workDbl, int lworkDbl,
	int *workInt, int lworkInt)
{
        int n = L.rows();
	// dbl work
	double *pN = workDbl; // 2*N
	double *y = pN + 2*N; // 2*n*N
	double *values = y + 2*n*N; // ns
	Map<MatrixXd> a(values + ns, n, N << 1); // 2*n*N
	Map<MatrixXd> b(a.data()+2*n*N, n, N << 1); // 2*n*N
	Map<MatrixXd> qN(b.data()+2*n*N, n+1, N); // (n+1)*N
	Map<MatrixXd> x(qN.data()+(n+1)*N, n, N << 1); // 2*n*N
	Map<VectorXd> xr(x.data()+2*n*N, n); // n
	Map<RowVectorXd> s(xr.data()+n, N << 1); // 2*N
	double *subworkDbl = s.data()+2*N; // 6*(2*N)
	int lsubworkDbl = 12*N;
	int memDbl = 2*N + 2*n*N + ns + 2*n*N + 2*n*N + (n+1)*N + 2*n*N + n + 
		2*N + 12*N; // 9*n*N + n + ns + 17*N
	if(lworkDbl < memDbl)
		Rcpp::stop("Insufficient memory for mvt \n");
	// int work
	int *prime = workInt; 
	int *scaler_N = prime + n+1;
	int *scaler_ns = scaler_N + 2*N;
	int *subworkInt = scaler_ns + ns;
	int lsubworkInt = 2*N;
	int memInt = n+1 + 2*N + ns + 2*N;
	if(lworkInt < memInt)
		Rcpp::stop("Insufficient memory for mvt \n");
	// generate n+1 prime numbers
	primes(5*(n+1)*log((double)(n+1)+1)/4, n+1, prime);
	// qN = prime * one2N
	transform(prime, prime+n+1, qN.data(), [](int x){return 
		sqrt((double) x);});
	for(int i = 1; i < qN.cols(); i++)
		transform(qN.data(), qN.data()+n+1, qN.col(i-1).data(), qN.col(i).
			data(), [](double x1,double x2){return x1 + x2;});
	// MC
        for(int i = 0; i < ns; i++)
        {
		fill(scaler_N, scaler_N+2*N, 0);
		GetRNGstate();
                double sr = unif_rand();
		PutRNGstate();
		for(int j = 0; j < N; j++) s(j) = qN(0,j) + sr;
		for_each(s.data(), s.data()+N, [](double &x){x = abs(2.0 * 
			(x - int(x)) - 1.0);});
		transform(s.data(), s.data()+N, s.data()+N, [](double &x){return 
			1 - x;});
		// chisq quantile
		for_each(s.data(), s.data()+(N<<1), [&nu](double &x){x = 
			sqrt(R::qchisq(x, nu, 0, 0)/nu);});
		// x with R RNG
		GetRNGstate();
		for_each(xr.data(), xr.data()+n, [](double &x){x = unif_rand();});
		PutRNGstate();
		for(int j = 0; j < N; j++)
			transform(qN.col(j).data()+1, qN.col(j).data()+n+1, xr.data(),
				x.col(j).data(), [](double &xFix, double &xr){
				return xFix + xr;});
		for_each(x.data(), x.data()+n*N, [](double &x){x = abs(2.0 * 
			(x - int(x)) - 1.0);});
		transform(x.data(), x.data()+n*N, x.data()+n*N, [](double &x){return 
			1 - x;});
		// scale a, b
                a.noalias() = a1 * s;
                b.noalias() = b1 * s;
		a.colwise() -= mu;
		b.colwise() -= mu;
		// kernel
                mvndns(n, N<<1, L, x, a, b, pN, y, n, scaler_N, subworkDbl, 
			lsubworkDbl, subworkInt, lsubworkInt);
		// scale
                int scaler_max = *(std::max_element(scaler_N, scaler_N+2*N));
                std::transform(scaler_N, scaler_N+2*N, pN, pN,
                        [scaler_max](int scaler, double base){return scalbn(base,
                        scaler - scaler_max);});
		// store
		values[i] = accumulate(pN, pN+(N<<1), 0.0) / (double) (N<<1);
                scaler_ns[i] = scaler_max;
        }
	int scaler_max = *(std::max_element(scaler_ns, scaler_ns+ns));
        std::transform(scaler_ns, scaler_ns+ns, values, values, 
		[scaler_max](int scaler, double base){return scalbn(base,
                scaler - scaler_max);});
        mean_std(ns, values, v, e);
        e = e / sqrt( (double) ns);
	scaler_in = scaler_max;

        return 0;
}

