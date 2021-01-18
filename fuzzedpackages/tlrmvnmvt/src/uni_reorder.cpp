#include <RcppEigen.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include <cmath>
#include "mvphi.h"

using namespace std;
using namespace Eigen;

/*
	Univariate reordering
	The size of work should be at least 6n
	Corrected a bug so that the upper part of B is not referenced until the end
	2019/09/17
*/
int uni_reorder(int n, double *B, int ldB, double *a, double *b, double &p, double 
	*y, int *oldIdx, double *work, int lwork)
{
	// Caution: do not use variables before they are set
	if(lwork < 6*n)
		Rcpp::stop("Allocated memory for uni_reorder is insufficient\n");
	p = 1.0;
	double *sd = work;
	double *aScaled = work+n;
	double *bScaled = work+n*2;
	double *aScaledCDF = work+n*3;
	double *bScaledCDF = work+n*4;
	double *pAll = work+n*5;
	for(int i = 0; i < n; i++)
	{
		for(int j = i; j < n; j++) sd[j-i] = sqrt(B[ldB*j+j]); // sd i to n
		copy_n(a+i, n-i, aScaled); // cp a i to n
		copy_n(b+i, n-i, bScaled);
		transform(aScaled, aScaled+n-i, sd, aScaled, [](double aCoef, double 
			sdCoef){return aCoef/sdCoef;});
		transform(bScaled, bScaled+n-i, sd, bScaled, [](double bCoef, double 
			sdCoef){return bCoef/sdCoef;});
		lc_vdCdfNorm(n-i, aScaled, aScaledCDF);
		lc_vdCdfNorm(n-i, bScaled, bScaledCDF);
		transform(bScaledCDF, bScaledCDF+n-i, aScaledCDF, pAll, [](double p1,
			double p2){return p1 - p2;});
		double *minCoef = min_element(pAll, pAll+n-i);
		int minCoefIdx = minCoef - pAll + i;
		// swap B, a, b
		iter_swap(oldIdx+i, oldIdx+minCoefIdx);
		iter_swap(B+i+ldB*i, B+minCoefIdx+ldB*minCoefIdx);
		for(int j = 0; j < i; j++)
			iter_swap(B+i+ldB*j, B+minCoefIdx+ldB*j);
		for(int j = i+1; j < minCoefIdx; j++) 
			iter_swap(B+j+ldB*i, B+minCoefIdx+ldB*j);
		for(int j = minCoefIdx+1; j < n; j++)
			iter_swap(B+j+ldB*i, B+j+ldB*minCoefIdx);
		iter_swap(a+i, a+minCoefIdx);
		iter_swap(b+i, b+minCoefIdx);
		// Cholesky
		if(*(B+i+i*ldB) <= 0)
			return i;
		*(B+i+i*ldB) = sqrt(*(B+i+i*ldB));
		if(i < n-1)
		{
			for_each(B+i+1+i*ldB, B+n+i*ldB, [&quo = *(B+i+i*ldB)](
				double &Bcoef){Bcoef /= quo;});
			int nrowLeft = n-i-1;
			double alpha = -1.0;
			int step = 1;
			F77_CALL(dsyr)("L", &nrowLeft, &alpha, B+i+1+i*ldB, &step, 
				B+i+1+(i+1)*ldB, &ldB);
		} // i < n-1
		// compute y
		double aScaledCoef = a[i] / *(B+i+i*ldB);
		double bScaledCoef = b[i] / *(B+i+i*ldB);
		double aScaledCoefCDF;
		double bScaledCoefCDF;
		lc_vdCdfNorm(1, &aScaledCoef, &aScaledCoefCDF);
		lc_vdCdfNorm(1, &bScaledCoef, &bScaledCoefCDF);
		double pCoef = bScaledCoefCDF - aScaledCoefCDF;
		y[i] = ((exp(-aScaledCoef*aScaledCoef/2.0) - exp(-bScaledCoef*
			bScaledCoef/2.0)) / sqrt(2.0*M_PI)) / pCoef;
		p *= pCoef;
		// update a, b
		if(i < n-1)
		{
			transform(a+i+1, a+n, B+i+1+ldB*i, a+i+1, [&yCoef = y[i]](
				double aCoef, double BCoef){return aCoef - BCoef*
				yCoef;});
			transform(b+i+1, b+n, B+i+1+ldB*i, b+i+1, [&yCoef = y[i]](
				double bCoef, double BCoef){return bCoef - BCoef*
				yCoef;});
		}
	} // i = 0 : n
	for(int j = 1; j < n; j++) // col num
		for_each(B+j*ldB, B+j+j*ldB, [](double &BCoef){BCoef = 0.0;});
	return 0;
}

/*
	The old version of univariate reordering moved from the rtlrmvn project
	Used here for checking purpose
        2019/05/20
*/
//std::vector<int> uni_reorder_old(Eigen::MatrixXd &B, Eigen::VectorXd &a, 
//	Eigen::VectorXd &b, double &p, Eigen::VectorXd &y)
//{
//        int m = B.rows();
//        vector<int> idx(m);
//        iota(idx.begin(),idx.end(),0);
//
//        MatrixXf::Index minCoefPos;
//        int tmpInt;
//        double tmpDbl,tmpLowerb,tmpUpperb;
//        p = 1.0;
//        y = VectorXd::Zero(m);
//        VectorXd tmpVec(m),tmpVecRow(m);
//        VectorXd cdfATlowerb(m);
//        VectorXd cdfATupperb(m);
//        for(int j = 0; j < m; j++)
//        {
//                ArrayXd sd_vec = sqrt(B.diagonal().segment(j,m-j).array());
//                VectorXd lowerBounds = a.segment(j,m-j).array() / sd_vec;
//                VectorXd upperBounds = b.segment(j,m-j).array() / sd_vec;
//
//                lc_vdCdfNorm(m-j,lowerBounds.data(),cdfATlowerb.data());
//                lc_vdCdfNorm(m-j,upperBounds.data(),cdfATupperb.data());
//                VectorXd probVec = cdfATupperb.segment(0,m-j) -
//                        cdfATlowerb.segment(0,m-j);
//                probVec.minCoeff(&minCoefPos);
//                minCoefPos += j;
//                tmpInt = idx[j];
//                idx[j] = idx[minCoefPos];
//                idx[minCoefPos] = tmpInt;
//                tmpVecRow = B.row(j);
//                B.row(j) = B.row(minCoefPos);
//                B.row(minCoefPos) = tmpVecRow;
//                tmpVec = B.col(j);
//                B.col(j) = B.col(minCoefPos);
//                B.col(minCoefPos) = tmpVec;
//                tmpDbl = a(j);
//                a(j) = a(minCoefPos);
//                a(minCoefPos) = tmpDbl;
//                tmpDbl = b(j);
//                b(j) = b(minCoefPos);
//                b(minCoefPos) = tmpDbl;
//
//                assert(B(j,j) > 0);
//                B(j,j) = sqrt(B(j,j));
//                if(j < m-1)
//                {
//                        B.col(j).segment(j+1,m-1-j) = B.col(j).segment(j+1,m-1-j)
//                                / B(j,j);
//                        B.block(j+1,j+1,m-1-j,m-1-j) =
//                                B.block(j+1,j+1,m-1-j,m-1-j) -
//                                B.col(j).segment(j+1,m-1-j) *
//                                B.col(j).segment(j+1,m-1-j).transpose();
//                }
//                tmpLowerb = a(j) / B(j,j);
//                tmpUpperb = b(j) / B(j,j);
//                lc_vdCdfNorm(1,&tmpLowerb,cdfATlowerb.data());
//                lc_vdCdfNorm(1,&tmpUpperb,cdfATupperb.data());
//                tmpDbl = cdfATupperb(0) - cdfATlowerb(0);
//                y(j) = ((exp(-tmpLowerb * tmpLowerb / 2.0) -
//                exp(-tmpUpperb * tmpUpperb / 2.0)) / sqrt(2.0*M_PI)) / tmpDbl;
//                p *= tmpDbl;
//
//                if(j < m-1)
//                {
//                        a.segment(j+1,m-1-j) = a.segment(j+1,m-1-j) -
//                                B.col(j).segment(j+1,m-1-j)*y(j);
//                        b.segment(j+1,m-1-j) = b.segment(j+1,m-1-j) -
//                                B.col(j).segment(j+1,m-1-j)*y(j);
//                }
//        }
//        B.triangularView<StrictlyUpper>().setZero();
//        return idx;
//}
