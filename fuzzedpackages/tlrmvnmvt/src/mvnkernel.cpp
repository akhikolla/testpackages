#include <RcppEigen.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include <chrono>
#include "mvphi.h"

#define MIN_PROB_EXPO -1E3

using namespace std;
using namespace Eigen;

typedef std::chrono::time_point<std::chrono::steady_clock> TimeStamp;

/*
	workM should be of dimension (k, N), where k >= 6
	workVecI should be of dimension N
	p should be of dimension N
	2019/09/03
*/
int mvndns(int n,  int N, const Eigen::MatrixXd& L, const Eigen::MatrixXd &x,
       const Eigen::MatrixXd &a, const Eigen::MatrixXd &b,
       Eigen::ArrayXd &p, Eigen::MatrixXd &y, int yOffset, std::vector<int>& scaler,
       Eigen::MatrixXd &workM, Eigen::VectorXi &workVecI)
{
	p.setOnes();
	int stride = workM.rows();
	workM.row(5).setZero();
        // iter through int var
        for (int i = 0; i < n; i++) {
                if (i > 0) {
			workM.row(5) = x.row(i-1).cwiseProduct(workM.row(4));
			workM.row(2).noalias() += workM.row(5);
		        lc_vdCdfNormInv(N, workM.row(2).data(), workM.row(5).data(),
				stride);
		        y.row(yOffset+i-1) = workM.row(5);
			workM.row(5).noalias() = L.row(i).segment(0,i) * y.block(
				yOffset, 0, i, N);
		}
                double ct = L(i,i);
		workM.row(0) = a.row(i) - workM.row(5); // ai
		workM.row(1) = b.row(i) - workM.row(5); // bi
		workM.row(0).noalias() = workM.row(0) / ct;
		workM.row(1).noalias() = workM.row(1) / ct;
		// c and d
                lc_vdCdfNorm(N, workM.row(0).data(), workM.row(2).data(), stride); 
                lc_vdCdfNorm(N, workM.row(1).data(), workM.row(3).data(), stride); 
                workM.row(4) = workM.row(3) - workM.row(2); // dc
		// cumularte p 
                p = p * workM.row(4).transpose().array();
                // scale
		std::transform(p.data(), p.data()+N, workVecI.data(), []
                        (double val){return ilogb(val);});
                std::transform(workVecI.data(), workVecI.data()+N, p.data(), p.data(),
                        [](int expo, double base){return scalbn(base, -expo);});
                std::transform(workVecI.data(), workVecI.data()+N, scaler.begin(),
                        scaler.begin(), [](int scaler1, int scaler2){return
                        scaler1 + scaler2;});
        }
	workM.row(5) = x.row(n-1).cwiseProduct(workM.row(4));
	workM.row(2).noalias() += workM.row(5);
        lc_vdCdfNormInv(N, workM.row(2).data(), workM.row(5).data(), stride);
        y.row(yOffset+n-1) = workM.row(5);

        return 0;
}

/*
	workM should be of dimension (k, N), where k >= 6
	workVecI should be of dimension N
	p should be of dimension N
	2019/09/03
	workM is mapped from workDbl and workInt
	workVecI is substituted by a pointer directly
	p is substituted by a pointer
	y is substituted by a pointer, no need for offset but need for leading dim
	scaler is substituted by a pointer
	2019/09/25
	Added checking abnormal exponent after each iteration
	2019/09/26
*/
int mvndns(int n,  int N, const Eigen::MatrixXd& L, const Eigen::MatrixXd &x,
       const Eigen::MatrixXd &a, const Eigen::MatrixXd &b,
       double *p, double *y, int ldy, int *scaler,
       double *workDbl, int lworkDbl, int *workInt, int lworkInt)
{
	if(lworkDbl < 6*N)
		Rcpp::stop("Dimension of workDbl is insufficient\n");
	if(lworkInt < N)
		Rcpp::stop("Dimension of workInt is insufficient\n");
	if(ldy < n)
		Rcpp::stop("The leading dimension of y is smaller than n\n");

	fill(p, p+N, 1.0);
	fill(workDbl + 5*N, workDbl + 6*N, 0.0);
	double alpha, beta;
	int stride = 1;
	const double *aPtr = a.data();
	const double *bPtr = b.data();
	const double *xPtr = x.data();
	const double *LPtr = L.data();
	double *dblPtr1, *dblPtr2;
	Map<RowVectorXd> s(workDbl+5*N, N);
	Map<MatrixXd, 0, OuterStride<>> yMap(y, n, N, OuterStride<>(ldy));
        // iter through int var
        for (int i = 0; i < n; i++) {
                if (i > 0) {
			dblPtr1 = workDbl + 4*N;
			dblPtr2 = workDbl + 5*N;
			for(int j = 0; j < N; j++) *(dblPtr2+j) = *(dblPtr1+j) * 
				*(xPtr+i-1+j*n);
			transform(workDbl+2*N, workDbl+3*N, workDbl+5*N, workDbl+2*N,
				[](double x1, double x2){return x1 + x2;});
		        lc_vdCdfNormInv(N, workDbl+2*N, workDbl+5*N);
			F77_CALL(dcopy)(&N, workDbl+5*N, &stride, y+i-1, &ldy);
//			alpha = 1.0;
//			beta = 0.0;
//			F77_CALL(dgemv)("T", &i, &N, &alpha, y, &ldy, LPtr+i, 
//				&n, &beta, workDbl+5*N, &stride);
			s.noalias() = L.row(i).segment(0,i) * yMap.block(0,0,i,N);
		}
                double ct = *(LPtr+i+i*n);
		for(int j = 0; j < N; j++) *(workDbl+j) = *(aPtr+i+j*n) - 
			*(workDbl+5*N+j); // ai
		for(int j = 0; j < N; j++) *(workDbl+N+j) = *(bPtr+i+j*n) - 
			*(workDbl+5*N+j); // bi
		for_each(workDbl, workDbl+2*N, [ct](double &x){x = x / ct;}); // scale
                lc_vdCdfNorm(N, workDbl, workDbl+2*N); // c
                lc_vdCdfNorm(N, workDbl+N, workDbl+3*N); // d
		transform(workDbl+3*N, workDbl+4*N, workDbl+2*N, workDbl+4*N, []
			(double &x, double &y){return x-y;}); // dc
		// cumulate p 
		transform(p, p+N, workDbl+4*N, p, [](double x, double y)
			{return x*y;});
                // scale
		std::transform(p, p+N, workInt, [](double val){return ilogb(val);});
		std::transform(workInt, workInt+N, p, p, [](int e, double base)
			{return e < MIN_PROB_EXPO ? 0.0 : base;});
		std::for_each(workInt, workInt+N, [](int &e){e = e < 
			MIN_PROB_EXPO ? 0 : e;});
                std::transform(workInt, workInt+N, p, p, [](int expo, double base)
			{return scalbn(base, -expo);});
                std::transform(workInt, workInt+N, scaler,
                        scaler, [](int scaler1, int scaler2){return
                        scaler1 + scaler2;});
        }
	dblPtr1 = workDbl + 4*N;
	dblPtr2 = workDbl + 5*N;
	for(int j = 0; j < N; j++) *(dblPtr2+j) = *(dblPtr1+j) * *(xPtr+n-1+j*n);
	transform(workDbl+2*N, workDbl+3*N, workDbl+5*N, workDbl+2*N,
		[](double x1, double x2){return x1 + x2;});
        lc_vdCdfNormInv(N, workDbl+2*N, workDbl+5*N);
	F77_CALL(dcopy)(&N, workDbl+5*N, &stride, y+n-1, &ldy);

        return 0;
}
