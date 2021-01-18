// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/Print.h>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <functional>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mvphi.h"
#include "morton.h"
#include "mvnkernel.h"
#include "misc.h"
#include "aca.h"
#include "cholesky.h"
#include "covariance.h"
#include "mvn.h"
#include "uncompress.h"
#include "tlr.h"
#include "tlrmvn.h"
#include "tlr_tlr_add_qr.h"
#include "uni_reorder.h"
#include "blk_reorder.h"
#include "recur_blk_reorder.h"
#include "mvt.h"
#include "tlrmvt.h"
#include "cov_kernel.h"

using namespace std;
using namespace Eigen;

typedef std::chrono::time_point<std::chrono::steady_clock> TimeStamp;

// [[Rcpp::export]]
Rcpp::List mvn_internal(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd 
	covM, bool useLog2, int N)
{
	// tmp var
	int n = covM.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeChol, timeInt;
	// scale a, b, covM
	{
		VectorXd diagVec = covM.diagonal();
		diagVec.noalias() = diagVec.unaryExpr([](double x){return 
			1.0/sqrt(x);});
		auto diagM = diagVec.asDiagonal();
		a.noalias() = diagM * a;
		b.noalias() = diagM * b;
		covM.noalias() = diagM * covM;
		covM.noalias() = covM * diagM;
	}
	// mem alloc
	int lworkDbl = 9*n*N + n + ns + 14*N;
	int lworkInt = max(4*N + n + ns + 1, 2*n);
	double *workDbl;
	int *workInt;
	workDbl = new double[lworkDbl];
	workInt = new int[lworkInt];
	// uni reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *subworkDbl = bcp + n;
	int lsubworkDbl = 6*n;
	int *idx = workInt;
	int *subworkInt = idx + n;
	int lsubworkInt = n;
	copy(a.data(), a.data()+n, acp);
	copy(b.data(), b.data()+n, bcp);
	iota(idx, idx+n, 0);
	fail = uni_reorder(n, covM.data(), covM.rows(), acp, bcp, v, y, idx, 
		subworkDbl, lsubworkDbl);
	if(fail)
		Rcpp::stop("Cholesky failed. Please check the positive definiteness "
			"of the input covariance matrix\n");
	reorder(a.data(), idx, n, subworkInt, lsubworkInt);
	reorder(b.data(), idx, n, subworkInt, lsubworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the mvn function
	start = std::chrono::steady_clock::now();
	mvn(N, covM, a, b, v, e, ns, scaler, workDbl, lworkDbl, workInt, lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}
}

// [[Rcpp::export]]
Rcpp::List mvn_internal2(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd 
	geom, int kernelType, Eigen::VectorXd para, double nugget, bool useLog2, 
	int N)
{
	// tmp var
	int n = geom.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeCovM, timeChol, timeInt;
	// build matern corr mat (btm half)
	start = std::chrono::steady_clock::now();
	function<double(double)> kernel;
	switch(kernelType){
		case 1:
			kernel = matern_kernel(para[1], para[2]);
			break;
		default:
			Rcpp::stop("Invalid kernel type\n");
			break;
	}
	MatrixXd covM = dense_covM(geom, kernel);
	nugget = nugget / para[0] / para[0];
	if(nugget != 0.0) for(int i = 0; i < n; i++) covM(i,i) += nugget;
	if(para[0] != 1.0) {a.noalias() = a / para[0]; b.noalias() = b / para[0];}
	end = std::chrono::steady_clock::now();
	timeCovM = std::chrono::duration<double>(end - start).count();
	// mem alloc
	int lworkDbl = 9*n*N + n + ns + 14*N;
	int lworkInt = max(4*N + n + ns + 1, 2*n);
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// uni reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *subworkDbl = bcp + n;
	int lsubworkDbl = 6*n;
	int *idx = workInt;
	int *subworkInt = idx + n;
	int lsubworkInt = n;
	copy(a.data(), a.data()+n, acp);
	copy(b.data(), b.data()+n, bcp);
	iota(idx, idx+n, 0);
	fail = uni_reorder(n, covM.data(), covM.rows(), acp, bcp, v, y, idx, 
		subworkDbl, lsubworkDbl);
	if(fail)
		Rcpp::stop("Cholesky failed. Please check the positive definiteness "
			"of the input covariance matrix\n");
	reorder(a.data(), idx, n, subworkInt, lsubworkInt);
	reorder(b.data(), idx, n, subworkInt, lsubworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the mvn function
	start = std::chrono::steady_clock::now();
	mvn(N, covM, a, b, v, e, ns, scaler, workDbl, lworkDbl, workInt, lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Building covariance matrix") = timeCovM,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Building covariance matrix") = timeCovM,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}
}

// [[Rcpp::export]]
Rcpp::List tlrmvn_internal(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd 
	covM, bool useLog2, int m, double epsl, int N)
{
	// tmp var
	int n = covM.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeCovM, timeChol, timeInt;
	// scale a, b, covM
	{
		VectorXd diagVec = covM.diagonal();
		diagVec.noalias() = diagVec.unaryExpr([](double x){return 
			1.0/sqrt(x);});
		auto diagM = diagVec.asDiagonal();
		a.noalias() = diagM * a;
		b.noalias() = diagM * b;
		covM.noalias() = diagM * covM;
		covM.noalias() = covM * diagM;
	}
	// build tlr covM, extend a and b if necessary, resize covM after
	start = std::chrono::steady_clock::now();
	vector<MatrixXd> B;
	vector<TLRNode> UV;
	int allocSz = max(m >> 2, 20);
	double epslACA = epsl / (double) m;
	tlr_aca_covM(covM, B, UV, m, epslACA, allocSz);
	int lastBlkDim = n % m;
	if(lastBlkDim > 0)
	{
		VectorXd aCp = a;
		VectorXd bCp = b;
		a.resize(n + m - lastBlkDim);
		b.resize(n + m - lastBlkDim);
		a.segment(0, n) = aCp;
		b.segment(0, n) = bCp;
		a.segment(n, m - lastBlkDim).setConstant(-10.0);
		b.segment(n, m - lastBlkDim).setConstant(10.0);
	}
	covM.resize(0,0);
	end = std::chrono::steady_clock::now();
	timeCovM = std::chrono::duration<double>(end - start).count();
	// mem alloc
	if(lastBlkDim > 0) n += m - lastBlkDim;
	int lworkDbl = max(19*m*m+16*m+4*n, (5*n + 4*m + 19)*N + ns + m);
	int lworkInt = 4*N + 2*n + 2*m + ns + 1;
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// recur blk reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *p = bcp + n;
	double *subworkDbl = p + n;
	int lsubworkDbl = 19*m*m + 16*m;
	if(subworkDbl + lsubworkDbl > workDbl + lworkDbl)
		Rcpp::stop("Memory overflow\n");
	int *idx = workInt;
	int *subworkInt = idx + n;
	int lsubworkInt = max(2*m, n);
	if(subworkInt + lsubworkInt > workInt + lworkInt)
		Rcpp::stop("Memory overflow\n");
	copy(a.data(), a.data()+n, acp);
	copy(b.data(), b.data()+n, bcp);
	iota(idx, idx+n, 0);
	fail = recur_blk_reorder(B, UV, acp, bcp, p, y, idx, epsl, subworkDbl, 
		lsubworkDbl, subworkInt, lsubworkInt);
	if(fail)
		Rcpp::stop("TLR Cholesky failed. Either the original covariance "
			"matrix is not positive definite or the error in Cholesky "
			"factorization is significant\n");
	reorder(a.data(), idx, n, subworkInt, lsubworkInt);
	reorder(b.data(), idx, n, subworkInt, lsubworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the tlrmvn function
	start = std::chrono::steady_clock::now();
	tlrmvn(N, B, UV, a, b, v, e, ns, scaler, workDbl, lworkDbl, workInt, 
		lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// compute average rank
	int rkSum = 0;
	for(auto &tile : UV) rkSum += tile.crtColNum;
	int rkAvg = rkSum / UV.size();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}
}

// [[Rcpp::export]]
Rcpp::List tlrmvn_internal2(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd 
	geom, int kernelType, Eigen::VectorXd para, double nugget, bool useLog2, 
	int m, double epsl, int N)
{
	// tmp var
	int n = geom.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeCovM, timeChol, timeInt;
	// build tlr covM, extend a and b if necessary, scale a and b
	start = std::chrono::steady_clock::now();
	vector<MatrixXd> B;
	vector<TLRNode> UV;
	int allocSz = max(m >> 2, 20);
	epsl = epsl * sqrt((double) m);
	double epslACA = epsl / (double) m;
	function<double(double)> kernel;
        switch(kernelType){
                case 1:
                        kernel = matern_kernel(para[1], para[2]);
                        break;
                default:
                        Rcpp::stop("Invalid kernel type\n");
                        break;
        }
	if(geom.cols() != 2) Rcpp::stop("The geometry for `mvn.tlr` should be 2D, "
		"each row represents one location\n");
	if(geom.maxCoeff() > 1.0 || geom.minCoeff() < 0.0) Rcpp::stop("The geometry "
		"for `mvn.tlr` should be contained in the unit square\n");
	vector<int> idx = zsort(geom);
	tlr_aca_covM(geom, B, UV, kernel, m, idx, epslACA, allocSz);
	nugget = nugget / para[0] / para[0];
	if(nugget != 0.0) for(int i = 0; i < n; ++i) B[i/m](i%m, i%m) += nugget;
	if(para[0] != 1.0) {a.noalias() = a / para[0]; b.noalias() = b / para[0];}
	int lastBlkDim = n % m;
	if(lastBlkDim > 0)
	{
		VectorXd aCp = a;
		VectorXd bCp = b;
		vector<int> idxCp = idx;
		a.resize(n + m - lastBlkDim);
		b.resize(n + m - lastBlkDim);
		idx.resize(n + m - lastBlkDim);
		a.segment(0, n) = aCp;
		b.segment(0, n) = bCp;
		copy(idxCp.begin(), idxCp.end(), idx.begin());
		a.segment(n, m - lastBlkDim).setConstant(-10.0);
		b.segment(n, m - lastBlkDim).setConstant(10.0);
		iota(idx.begin()+n, idx.end(), n);
	}
	end = std::chrono::steady_clock::now();
	timeCovM = std::chrono::duration<double>(end - start).count();
	// mem alloc
	if(lastBlkDim > 0) n += m - lastBlkDim;
	int lworkDbl = max(19*m*m+16*m+4*n, (5*n + 4*m + 19)*N + ns + m);
	int lworkInt = 4*N + 2*n + 2*m + ns + 1;
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// recur blk reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *p = bcp + n;
	double *subworkDbl = p + n;
	int lsubworkDbl = 19*m*m + 16*m;
	if(subworkDbl + lsubworkDbl > workDbl + lworkDbl)
		Rcpp::stop("Memory overflow\n");
	for(int i = 0; i < n; i++)
	{
		acp[i] = a[idx[i]];
		bcp[i] = b[idx[i]];
	}
	fail = recur_blk_reorder(B, UV, acp, bcp, p, y, &(idx[0]), epsl, subworkDbl, 
		lsubworkDbl, workInt, lworkInt);
	if(fail)
		Rcpp::stop("TLR Cholesky failed. Either the original covariance "
			"matrix is not positive definite or the error in Cholesky "
			"factorization is significant\n");
	reorder(a.data(), &(idx[0]), n, workInt, lworkInt);
	reorder(b.data(), &(idx[0]), n, workInt, lworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the tlrmvn function
	start = std::chrono::steady_clock::now();
	tlrmvn(N, B, UV, a, b, v, e, ns, scaler, workDbl, lworkDbl, workInt, 
		lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// compute average rank
	int rkSum = 0;
	for(auto &tile : UV) rkSum += tile.crtColNum;
	int rkAvg = rkSum / UV.size();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}
}

// [[Rcpp::export]]
Rcpp::List mvt_internal(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd mu,
	double nu, Eigen::MatrixXd covM, bool useLog2, int N)
{
	// tmp var
	int n = covM.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeChol, timeInt;
	// scale a, b, covM
	{
		VectorXd diagVec = covM.diagonal();
		diagVec.noalias() = diagVec.unaryExpr([](double x){return 
			1.0/sqrt(x);});
		auto diagM = diagVec.asDiagonal();
		a.noalias() = diagM * a;
		b.noalias() = diagM * b;
		mu.noalias() = diagM * mu;
		covM.noalias() = diagM * covM;
		covM.noalias() = covM * diagM;
	}
	// mem alloc
	int lworkDbl = 9*n*N + n + ns + 17*N;
	int lworkInt = max(4*N + n + ns + 1, 2*n);
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// uni reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *subworkDbl = bcp + n;
	int lsubworkDbl = 6*n;
	int *idx = workInt;
	int *subworkInt = idx + n;
	int lsubworkInt = n;
	copy(a.data(), a.data()+n, acp);
	copy(b.data(), b.data()+n, bcp);
	transform(acp, acp+n, mu.data(), acp, [](double a, double mu){
		return a - mu;});
	transform(bcp, bcp+n, mu.data(), bcp, [](double b, double mu){
		return b - mu;});
	iota(idx, idx+n, 0);
	fail = uni_reorder(n, covM.data(), covM.rows(), acp, bcp, v, y, idx, 
		subworkDbl, lsubworkDbl);
	if(fail)
		Rcpp::stop("Cholesky failed. Please check the positive definiteness "
			"of the input covariance matrix\n");
	reorder(a.data(), idx, n, subworkInt, lsubworkInt);
	reorder(b.data(), idx, n, subworkInt, lsubworkInt);
	reorder(mu.data(), idx, n, subworkInt, lsubworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the mvt function
	start = std::chrono::steady_clock::now();
	mvt(N, nu, mu, covM, a, b, v, e, ns, scaler, workDbl, lworkDbl, workInt, 
		lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}
}

// [[Rcpp::export]]
Rcpp::List mvt_internal2(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd mu,
	double nu, Eigen::MatrixXd geom, int kernelType, 
	Eigen::VectorXd para, double nugget, bool useLog2, int N)
{
	// tmp var
	int n = geom.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeCovM, timeChol, timeInt;
	// build matern corr mat (btm half)
	start = std::chrono::steady_clock::now();
	function<double(double)> kernel;
	switch(kernelType){
		case 1:
			kernel = matern_kernel(para[1], para[2]);
			break;
		default:
			Rcpp::stop("Invalid kernel type\n");
			break;
	}
	MatrixXd covM = dense_covM(geom, kernel);
	nugget = nugget / para[0] / para[0];
	if(nugget != 0.0) for(int i = 0; i < n; i++) covM(i,i) += nugget;
	if(para[0] != 1.0) 
	{
		a.noalias() = a / para[0]; 
		b.noalias() = b / para[0];
		mu.noalias() = mu / para[0];
	}
	end = std::chrono::steady_clock::now();
	timeCovM = std::chrono::duration<double>(end - start).count();
	// mem alloc
	int lworkDbl = 9*n*N + n + ns + 17*N;
	int lworkInt = max(4*N + n + ns + 1, 2*n);
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// uni reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *subworkDbl = bcp + n;
	int lsubworkDbl = 6*n;
	int *idx = workInt;
	int *subworkInt = idx + n;
	int lsubworkInt = n;
	copy(a.data(), a.data()+n, acp);
	copy(b.data(), b.data()+n, bcp);
	transform(acp, acp+n, mu.data(), acp, [](double a, double mu){
		return a - mu;});
	transform(bcp, bcp+n, mu.data(), bcp, [](double b, double mu){
		return b - mu;});
	iota(idx, idx+n, 0);
	fail = uni_reorder(n, covM.data(), covM.rows(), acp, bcp, v, y, idx, 
		subworkDbl, lsubworkDbl);
	if(fail)
		Rcpp::stop("Cholesky failed. Please check the positive definiteness "
			"of the input covariance matrix\n");
	reorder(a.data(), idx, n, subworkInt, lsubworkInt);
	reorder(b.data(), idx, n, subworkInt, lsubworkInt);
	reorder(mu.data(), idx, n, subworkInt, lsubworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the mvt function
	start = std::chrono::steady_clock::now();
	mvt(N, nu, mu, covM, a, b, v, e, ns, scaler, workDbl, lworkDbl, workInt, 
		lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Building covariance matrix") = timeCovM,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Building covariance matrix") = timeCovM,
			Rcpp::Named("Univariate reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt);
	}
}

// [[Rcpp::export]]
Rcpp::List tlrmvt_internal(Eigen::VectorXd a, Eigen::VectorXd b, double nu, 
	Eigen::VectorXd mu,
	Eigen::MatrixXd covM, bool useLog2, int m, double epsl, int N)
{
	// tmp var
	int n = covM.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeCovM, timeChol, timeInt;
	// scale a, b, covM
	{
		VectorXd diagVec = covM.diagonal();
		diagVec.noalias() = diagVec.unaryExpr([](double x){return 
			1.0/sqrt(x);});
		auto diagM = diagVec.asDiagonal();
		a.noalias() = diagM * a;
		b.noalias() = diagM * b;
		mu.noalias() = diagM * mu;
		covM.noalias() = diagM * covM;
		covM.noalias() = covM * diagM;
	}
	// build tlr covM, extend a, b, and mu if necessary, resize covM after
	start = std::chrono::steady_clock::now();
	vector<MatrixXd> B;
	vector<TLRNode> UV;
	int allocSz = max(m >> 2, 20);
	double epslACA = epsl / (double) m;
	tlr_aca_covM(covM, B, UV, m, epslACA, allocSz);
	int lastBlkDim = n % m;
	if(lastBlkDim > 0)
	{
		VectorXd aCp = a;
		VectorXd bCp = b;
		VectorXd muCp = mu;
		a.resize(n + m - lastBlkDim);
		b.resize(n + m - lastBlkDim);
		mu.resize(n + m - lastBlkDim);
		a.segment(0, n) = aCp;
		b.segment(0, n) = bCp;
		mu.segment(0, n) = muCp;
		a.segment(n, m - lastBlkDim).setConstant(-20.0);
		b.segment(n, m - lastBlkDim).setConstant(20.0);
		mu.segment(n, m - lastBlkDim).setConstant(0.0);
	}
	covM.resize(0,0);
	end = std::chrono::steady_clock::now();
	timeCovM = std::chrono::duration<double>(end - start).count();
	// mem alloc
	if(lastBlkDim > 0) n += m - lastBlkDim;
	int lworkDbl = max(19*m*m+16*m+4*n, (5*n + 4*m + 19)*N + ns + m);
	int lworkInt = 4*N + 2*n + 2*m + ns + 1;
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// recur blk reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *p = bcp + n;
	double *subworkDbl = p + n;
	int lsubworkDbl = 19*m*m + 16*m;
	if(subworkDbl + lsubworkDbl > workDbl + lworkDbl)
		Rcpp::stop("Memory overflow\n");
	int *idx = workInt;
	int *subworkInt = idx + n;
	int lsubworkInt = max(2*m, n);
	if(subworkInt + lsubworkInt > workInt + lworkInt)
		Rcpp::stop("Memory overflow\n");
	copy(a.data(), a.data()+n, acp);
	copy(b.data(), b.data()+n, bcp);
	transform(acp, acp+n, mu.data(), acp, [](double a, double mu){
                return a - mu;});
        transform(bcp, bcp+n, mu.data(), bcp, [](double b, double mu){
                return b - mu;});
	iota(idx, idx+n, 0);
	fail = recur_blk_reorder(B, UV, acp, bcp, p, y, idx, epsl, subworkDbl, 
		lsubworkDbl, subworkInt, lsubworkInt);
	if(fail)
		Rcpp::stop("TLR Cholesky failed. Either the original covariance "
			"matrix is not positive definite or the error in Cholesky "
			"factorization is significant\n");
	reorder(a.data(), idx, n, subworkInt, lsubworkInt);
	reorder(b.data(), idx, n, subworkInt, lsubworkInt);
	reorder(mu.data(), idx, n, subworkInt, lsubworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the tlrmvt function
	start = std::chrono::steady_clock::now();
	tlrmvt(N, nu, mu, B, UV, a, b, v, e, ns, scaler, workDbl, lworkDbl, 
		workInt, lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// compute average rank
	int rkSum = 0;
	for(auto &tile : UV) rkSum += tile.crtColNum;
	int rkAvg = rkSum / UV.size();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}
}

// [[Rcpp::export]]
Rcpp::List tlrmvt_internal2(Eigen::VectorXd a, Eigen::VectorXd b, double nu, 
	Eigen::VectorXd mu,
	Eigen::MatrixXd geom, int kernelType, Eigen::VectorXd para, double nugget, 
	bool useLog2, int m, double epsl, int N)
{
	// tmp var
	int n = geom.rows();
	int ns = 10;
	int fail;
	double v, e;
	int scaler;
	TimeStamp start, end;
	double timeCovM, timeChol, timeInt;
	// build tlr covM, extend a and b if necessary, scale a and b
	start = std::chrono::steady_clock::now();
	vector<MatrixXd> B;
	vector<TLRNode> UV;
	int allocSz = max(m >> 2, 20);
	epsl = epsl * sqrt((double) m);
	double epslACA = epsl / (double) m;
	function<double(double)> kernel;
        switch(kernelType){
                case 1:
                        kernel = matern_kernel(para[1], para[2]);
                        break;
                default:
                        Rcpp::stop("Invalid kernel type\n");
                        break;
        }
	if(geom.cols() != 2) Rcpp::stop("The geometry for `mvn.tlr` should be 2D, "
		"each row represents one location\n");
	if(geom.maxCoeff() > 1.0 || geom.minCoeff() < 0.0) Rcpp::stop("The geometry "
		"for `mvn.tlr` should be contained in the unit square\n");
	vector<int> idx = zsort(geom);
	tlr_aca_covM(geom, B, UV, kernel, m, idx, epslACA, allocSz);
	nugget = nugget / para[0] / para[0];
	if(nugget != 0.0) for(int i = 0; i < n; ++i) B[i/m](i%m, i%m) += nugget;
	if(para[0] != 1.0) 
	{
                a.noalias() = a / para[0];
                b.noalias() = b / para[0];
                mu.noalias() = mu / para[0];
        }
	int lastBlkDim = n % m;
	if(lastBlkDim > 0)
	{
		VectorXd aCp = a;
		VectorXd bCp = b;
		VectorXd muCp = mu;
		vector<int> idxCp = idx;
		a.resize(n + m - lastBlkDim);
		b.resize(n + m - lastBlkDim);
		mu.resize(n + m - lastBlkDim);
		idx.resize(n + m - lastBlkDim);
		a.segment(0, n) = aCp;
		b.segment(0, n) = bCp;
		mu.segment(0, n) = muCp;
		copy(idxCp.begin(), idxCp.end(), idx.begin());
		a.segment(n, m - lastBlkDim).setConstant(-20.0);
		b.segment(n, m - lastBlkDim).setConstant(20.0);
		mu.segment(n, m - lastBlkDim).setConstant(0.0);
		iota(idx.begin()+n, idx.end(), n);
	}
	end = std::chrono::steady_clock::now();
	timeCovM = std::chrono::duration<double>(end - start).count();
	// mem alloc
	if(lastBlkDim > 0) n += m - lastBlkDim;
	int lworkDbl = max(19*m*m+16*m+4*n, (5*n + 4*m + 19)*N + ns + m);
	int lworkInt = 4*N + 2*n + 2*m + ns + 1;
	double *workDbl = new double[lworkDbl];
	int *workInt = new int[lworkInt];
	// recur blk reorder
	start = std::chrono::steady_clock::now();
	double *y = workDbl;
	double *acp = y + n;
	double *bcp = acp + n;
	double *p = bcp + n;
	double *subworkDbl = p + n;
	int lsubworkDbl = 19*m*m + 16*m;
	if(subworkDbl + lsubworkDbl > workDbl + lworkDbl)
		Rcpp::stop("Memory overflow\n");
	for(int i = 0; i < n; i++)
	{
		acp[i] = a[idx[i]] - mu[idx[i]];
		bcp[i] = b[idx[i]] - mu[idx[i]];
	}
	fail = recur_blk_reorder(B, UV, acp, bcp, p, y, &(idx[0]), epsl, subworkDbl, 
		lsubworkDbl, workInt, lworkInt);
	if(fail)
		Rcpp::stop("TLR Cholesky failed. Either the original covariance "
			"matrix is not positive definite or the error in Cholesky "
			"factorization is significant\n");
	reorder(a.data(), &(idx[0]), n, workInt, lworkInt);
	reorder(b.data(), &(idx[0]), n, workInt, lworkInt);
	reorder(mu.data(), &(idx[0]), n, workInt, lworkInt);
	end = std::chrono::steady_clock::now();
	timeChol = std::chrono::duration<double>(end - start).count();
	// call the tlrmvt function
	start = std::chrono::steady_clock::now();
	tlrmvt(N, nu, mu, B, UV, a, b, v, e, ns, scaler, workDbl, lworkDbl, 
		workInt, lworkInt);
	end = std::chrono::steady_clock::now();
	timeInt = std::chrono::duration<double>(end - start).count();
	// compute average rank
	int rkSum = 0;
	for(auto &tile : UV) rkSum += tile.crtColNum;
	int rkAvg = rkSum / UV.size();
	// mem release
	delete[] workDbl;
	delete[] workInt;
	// return a list
	if(useLog2)
	{
		v = log2(v) + (double) scaler * log2(FLT_RADIX);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}else
	{
		v = scalbn(v, scaler);
		e = scalbn(e, scaler);
		return Rcpp::List::create(Rcpp::Named("Estimation") = v, 
			Rcpp::Named("Error") = e,
			Rcpp::Named("Building TLR covariance matrix time") = timeCovM,
			Rcpp::Named("Recursive block reordering time") = timeChol,
			Rcpp::Named("Monte Carlo time") = timeInt,
			Rcpp::Named("Average rank") = rkAvg);
	}
}
