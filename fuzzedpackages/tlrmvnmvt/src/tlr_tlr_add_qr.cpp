#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Print.h>
#include "cholesky.h"
#include "qr.h"
#include "svd.h"
#include "tlr.h"

using namespace std;
using namespace Eigen;

/*
	The leading dimension of U and V should be the same as chgNode.U.rows() and
		chgNode.V.rows()
	Total memory needed: 16*max(nrow,ncol)^2 + 12*max(nrow,ncol)
	Fix the bug that computing RU * RV^T, need to set the low tri of RU to 0 
		because RU is assumed to be a general matrix
	2019/09/17
	Add kNew = max(1, kNew)
*/
void tlr_tlr_add_qr(TLRNode &chgNode, const double *U, const double *V, int k2, double
	epsl, double *work, int lwork)
{
	int nrow = chgNode.U.rows();
	int ncol = chgNode.V.rows();
	int k1 = chgNode.crtColNum;
	int kCmb = k1 + k2;
	// tmp var
	double *UCmb = work;
	double *VCmb = UCmb + nrow*kCmb;
	double *RU = VCmb + ncol*kCmb;
	double *RV = RU + kCmb*kCmb;
	double *RUV = RV + kCmb*kCmb;
	double *tau = RUV + kCmb*kCmb;
	double *subwork = tau + kCmb; 
	int lsubwork = 10*max(nrow,ncol);
	int memTtl = (nrow+ncol)*kCmb + 3*kCmb*kCmb + kCmb + 10*max(nrow,ncol);
	if(lwork < memTtl)
		Rcpp::stop("Dimension of work is insufficient\n");
	int fail;
	// combine U and V
	copy_n(chgNode.U.data(), nrow*k1, UCmb);
	copy_n(U, nrow*k2, UCmb+nrow*k1);
	copy_n(chgNode.V.data(), ncol*k1, VCmb);
	copy_n(V, ncol*k2, VCmb+ncol*k1);
	// rk > dimension
	if(nrow < kCmb || ncol < kCmb)
	{
		if(nrow > kCmb)
		{
			// qr for Ucmb
			F77_CALL(dgeqrf)(&nrow, &kCmb, UCmb, &nrow, tau, subwork, 
				&lsubwork, &fail); // lsubwork >= kCmb * optimal bsz 
			if(fail)
				Rcpp::stop("`dgeqrf` failed code: %d\n", fail);
			for(int j = 0; j < kCmb; j++) 
				copy_n(UCmb+j*nrow, j+1, RU+j*kCmb);
			F77_CALL(dorgqr)(&nrow, &kCmb, &kCmb, UCmb, &nrow, tau,
				subwork, &lsubwork, &fail);
			if(fail)
				Rcpp::stop("`dorgqr` failed code: %d\n", fail);
			// VCmb = VCmb RU^T
			double alpha = 1.0;
			F77_CALL(dtrmm)("R", "U", "T", "N", &ncol, &kCmb, &alpha, RU, 
				&kCmb, VCmb, &ncol); 
			// svd(Vcmb)
			F77_CALL(dgesvd)("S", "S", &ncol, &kCmb, VCmb, &ncol, tau,
				RU, &kCmb, RV, &kCmb, subwork, &lsubwork, &fail);
				// RU stores the U factor
				// RV stores the VT factor
				// lsubwork >= 5 * kCmb
			if(fail)
				Rcpp::stop("`dgesvd` failed code: %d\n", fail);
			// count significant singular val
			int kNew = count_if(tau, tau+ncol, [&epsl](double &s){return
				s > epsl;});
			kNew = max(1, kNew);
			// check avail col in chgNode
			if(kNew > chgNode.maxColNum)
			{
//				Rprintf("Warning: the TLRNode needs resizing in TLR"
//					" addition. Consider enlarge their initial "
//					"allocations\n");
				chgNode.U.resize(nrow, kNew);
				chgNode.V.resize(ncol, kNew);
				chgNode.maxColNum = kNew;
			}
			// chgNode.U = UCmb * RV^T * sqrt(tau)
			// chgNode.V = RU * sqrt(tau)
			alpha = 1.0;
			double beta = 0.0;
			F77_CALL(dgemm)("N", "T", &nrow, &kNew, &kCmb, &alpha, UCmb, 
				&nrow, RV, &kCmb, &beta, chgNode.U.data(), &nrow);
			for(int j = 0; j < kNew; j++)
				copy_n(RU+j*kCmb, ncol, chgNode.V.data()+j*ncol);
			for(int j = 0; j < kNew; ++j)
			{
				double tauCoefSqrt = sqrt(tau[j]);
				for_each(chgNode.U.data()+j*nrow, chgNode.U.data()+
					j*nrow+nrow, [&tauCoefSqrt](double &UCoef){
					UCoef *= tauCoefSqrt;});
				for_each(chgNode.V.data()+j*ncol, chgNode.V.data()+
					j*ncol+ncol, [&tauCoefSqrt](double &VCoef){
					VCoef *= tauCoefSqrt;});
			}
			chgNode.crtColNum = kNew;
		} // nrow > kCmb && ncol <= kCmb
		else if(ncol > kCmb)
		{
			// qr for Vcmb
			F77_CALL(dgeqrf)(&ncol, &kCmb, VCmb, &ncol, tau, subwork, 
				&lsubwork, &fail); // lsubwork >= kCmb * optimal bsz 
			if(fail)
				Rcpp::stop("`dgeqrf` failed code: %d\n", fail);
			for(int j = 0; j < kCmb; j++) 
				copy_n(VCmb+j*ncol, j+1, RV+j*kCmb);
			F77_CALL(dorgqr)(&ncol, &kCmb, &kCmb, VCmb, &ncol, tau,
				subwork, &lsubwork, &fail);
			if(fail)
				Rcpp::stop("`dorgqr` failed code: %d\n", fail);
			// UCmb = UCmb RV^T
			double alpha = 1.0;
			F77_CALL(dtrmm)("R", "U", "T", "N", &nrow, &kCmb, &alpha, RV, 
				&kCmb, UCmb, &nrow); 
			// svd(Ucmb)
			F77_CALL(dgesvd)("S", "S", &nrow, &kCmb, UCmb, &nrow, tau,
				RU, &kCmb, RV, &kCmb, subwork, &lsubwork, &fail);
				// RU stores the U factor
				// RV stores the VT factor
				// lsubwork >= 5 * kCmb
			if(fail)
				Rcpp::stop("`dgesvd` failed code: %d\n", fail);
			// count significant singular val
			int kNew = count_if(tau, tau+nrow, [&epsl](double &s){return
				s > epsl;});
			kNew = max(1, kNew);
			// check avail col in chgNode
			if(kNew > chgNode.maxColNum)
			{
//				Rprintf("Warning: the TLRNode needs resizing in TLR"
//					" addition. Consider enlarge their initial "
//					"allocations\n");
				chgNode.U.resize(nrow, kNew);
				chgNode.V.resize(ncol, kNew);
				chgNode.maxColNum = kNew;
			}
			// chgNode.U = RU * sqrt(tau)
			// chgNode.V = VCmb * RV^T * sqrt(tau)
			alpha = 1.0;
			double beta = 0.0;
			F77_CALL(dgemm)("N", "T", &ncol, &kNew, &kCmb, &alpha, VCmb, 
				&ncol, RV, &kCmb, &beta, chgNode.V.data(), &ncol);
			for(int j = 0; j < kNew; j++)
				copy_n(RU+j*kCmb, nrow, chgNode.U.data()+j*nrow);
			for(int j = 0; j < kNew; ++j)
			{
				double tauCoefSqrt = sqrt(tau[j]);
				for_each(chgNode.U.data()+j*nrow, chgNode.U.data()+
					j*nrow+nrow, [&tauCoefSqrt](double &UCoef){
					UCoef *= tauCoefSqrt;});
				for_each(chgNode.V.data()+j*ncol, chgNode.V.data()+
					j*ncol+ncol, [&tauCoefSqrt](double &VCoef){
					VCoef *= tauCoefSqrt;});
			}	
			chgNode.crtColNum = kNew;
		} // ncol > kCmb
		else
		{
			// RUV = UCmb * VCmb^T
			double alpha = 1.0;
			double beta = 0.0;
			F77_CALL(dgemm)("N", "T", &nrow, &ncol, &kCmb, &alpha, UCmb,
			&nrow, VCmb, &ncol, &beta, RUV, &kCmb);
			// svd(RUV)
			F77_CALL(dgesvd)("S", "S", &nrow, &ncol, RUV, &kCmb, tau, 
				RU, &kCmb, RV, &kCmb, subwork, &lsubwork, &fail);
				// RU stores the U factor 
				// RV stores the VT factor
				// lsubwork >= 5 * max(nrow, ncol)
			// count significant singular val
			int kNew = count_if(tau, tau+min(nrow,ncol), [&epsl](double &
				s){return s > epsl;});
			kNew = max(1 , kNew);
			// check avail col in chgNode
			if(kNew > chgNode.maxColNum)
			{
//				Rprintf("Warning: the TLRNode needs resizing in TLR"
//					" addition. Consider enlarge their initial "
//					"allocations\n");
				chgNode.U.resize(nrow, kNew);
				chgNode.V.resize(ncol, kNew);
				chgNode.maxColNum = kNew;
			}
			// chgNode.U = RU * sqrt(tau)
			// chgNode.V = RV^T * sqrt(tau)
			for(int j = 0; j < kNew; j++)
			{
				copy_n(RU+j*kCmb, nrow, chgNode.U.data()+j*nrow);
				int incy = 1;
				F77_CALL(dcopy)(&ncol, RV+j, &kCmb, chgNode.V.data()+
					j*ncol, &incy);
			}
			for(int j = 0; j < kNew; ++j)
                        {
                                double tauCoefSqrt = sqrt(tau[j]);
                                for_each(chgNode.U.data()+j*nrow, chgNode.U.data()+
                                        j*nrow+nrow, [&tauCoefSqrt](double &UCoef){
                                        UCoef *= tauCoefSqrt;});
                                for_each(chgNode.V.data()+j*ncol, chgNode.V.data()+
                                        j*ncol+ncol, [&tauCoefSqrt](double &VCoef){
                                        VCoef *= tauCoefSqrt;});
                        }
			chgNode.crtColNum = kNew;
		} // ncol <= kCmb && nrow <= kCmb
	} // nrow < kCmb || ncol < kCmb
	else
	{
		// qr(UCmb)
		F77_CALL(dgeqrf)(&nrow, &kCmb, UCmb, &nrow, tau, subwork, &lsubwork, 
			&fail); // lsubwork >= kCmb * optimal bsz
		if(fail)
			Rcpp::stop("'dgeqrf' failed code: %d\n", fail);
		for(int j = 0; j < kCmb; j++)
		{
                        copy_n(UCmb+j*nrow, j+1, RU+j*kCmb);
			fill(RU+j*kCmb+j+1, RU+(j+1)*kCmb, 0.0);
		}
                F77_CALL(dorgqr)(&nrow, &kCmb, &kCmb, UCmb, &nrow, tau,
                        subwork, &lsubwork, &fail); // lsubwork >= kCmb * optimal bsz
                if(fail)
                        Rcpp::stop("`dorgqr` failed code: %d\n", fail);
		// qr(VCmb)
		F77_CALL(dgeqrf)(&ncol, &kCmb, VCmb, &ncol, tau, subwork, &lsubwork, 
			&fail); // lsubwork >= kCmb * optimal bsz
		if(fail)
			Rcpp::stop("'dgeqrf' failed code: %d\n", fail);
		for(int j = 0; j < kCmb; j++)
                        copy_n(VCmb+j*ncol, j+1, RV+j*kCmb);
                F77_CALL(dorgqr)(&ncol, &kCmb, &kCmb, VCmb, &ncol, tau,
                        subwork, &lsubwork, &fail); // lsubwork >= kCmb * optimal bsz
                if(fail)
                        Rcpp::stop("`dorgqr` failed code: %d\n", fail);
		// RU = RU * RV^T
		double alpha = 1.0;
		F77_CALL(dtrmm)("R", "U", "T", "N", &kCmb, &kCmb, &alpha, RV, &kCmb, 
			RU, &kCmb); // Assume RU is a general matrix 
			// Therefore, the lower tri of RU needs to be set to 0
		// svd(RU)
		F77_CALL(dgesvd)("S", "S", &kCmb, &kCmb, RU, &kCmb, tau, 
			RUV, &kCmb, RV, &kCmb, subwork, &lsubwork, &fail);
			// RUV stores the U factor 
			// RV stores the VT factor
			// lsubwork >= 5 * kCmb
		if(fail)
                        Rcpp::stop("`dgesvd` failed code: %d\n", fail);
		// count
		int kNew = count_if(tau, tau+kCmb, [&epsl](double &
			s){return s > epsl;});
		kNew = max(1, kNew);
		// check avail col in chgNode
		if(kNew > chgNode.maxColNum)
		{
//			Rprintf("Warning: the TLRNode needs resizing in TLR"
//				" addition. Consider enlarge their initial "
//				"allocations\n");
			chgNode.U.resize(nrow, kNew);
			chgNode.V.resize(ncol, kNew);
			chgNode.maxColNum = kNew;
		}
		// chgNode.U = UCmb * RUV * sqrt(tau)
		// chgNode.V = VCmb * RV^T * sqrt(tau)
		alpha = 1.0;
                double beta = 0.0;
                F77_CALL(dgemm)("N", "N", &nrow, &kNew, &kCmb, &alpha, UCmb,
                        &nrow, RUV, &kCmb, &beta, chgNode.U.data(), &nrow);
                F77_CALL(dgemm)("N", "T", &ncol, &kNew, &kCmb, &alpha, VCmb,
                        &ncol, RV, &kCmb, &beta, chgNode.V.data(), &ncol);
                for(int j = 0; j < kNew; ++j)
                {
                        double tauCoefSqrt = sqrt(tau[j]);
                        for_each(chgNode.U.data()+j*nrow, chgNode.U.data()+
                                j*nrow+nrow, [&tauCoefSqrt](double &UCoef){
                                UCoef *= tauCoefSqrt;});
                        for_each(chgNode.V.data()+j*ncol, chgNode.V.data()+
                                j*ncol+ncol, [&tauCoefSqrt](double &VCoef){
                                VCoef *= tauCoefSqrt;});
                }
		chgNode.crtColNum = kNew;
	} // nrow >= kCmb && ncol >= kCmb
}
		
