#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <algorithm>
#include "tlr.h"
#include "cholesky.h"
#include "inverse.h"
#include "uni_reorder.h"
#include "tlr_tlr_add_qr.h"
#include "misc.h"

using namespace std;
using namespace Eigen;

/*
	Perform TLR cholesky at the same time
	All block sizes should be the same
	The size of workInt should be at least 2*m
	The size of workDbl should be at least 19m^2 + 16m

	2019/09/09
*/
int recur_blk_reorder(std::vector<Eigen::MatrixXd> &B, std::vector<TLRNode> &UV, 
	double *a, double *b, double *p, double *y, int *oldIdx, double epsl, 
	double *workDbl, int lworkDbl, int *workInt, int lworkInt)
{
	int nb = B.size();
	int m = B[0].rows();
	int n = nb * m;
	// tmp var
	double *Bj = workDbl;
	double *aj = Bj + m*m;
	double *bj = aj + m; 
	double *Bi = bj + m;
	double *BSchur = Bi + m*m; 
	double *ySchur = BSchur + m*m; 
	double *ySchur2 = ySchur + m; 
	double *subworkDbl = ySchur2 + m; 
	int lsubworkDbl = 16*m*m + 12*m;
	int memTtlDbl = 3*m*m + 4*m + 16*m*m + 12*m;
	if(lworkDbl < memTtlDbl)
		Rcpp::stop("Dimension of workDbl is insufficient\n");
	int *idxj = workInt;
	int *subworkInt = workInt + m; 
	int lsubworkInt = m;
	int memTtlInt = m + m;
	if(lworkInt < memTtlInt)
		Rcpp::stop("Dimension of workInt is insufficient\n");
	// all diag blk
	for (int i = 0; i < nb; i++)
	{
		// remaining diag blk
		for (int j = i; j < nb; j++)
		{
			int offset = j*m;
			copy_n(B[j].data(), m*m, Bj);
			copy_n(a + offset, m, aj);
			copy_n(b + offset, m, bj);
			copy_n(oldIdx + offset, m, idxj);
			int fail = uni_reorder(m, Bj, m, aj, bj, p[j], y+offset, 
				idxj, subworkDbl, lsubworkDbl);
			if(fail)
				return j;
		} // j = i : nb
		int minIdx = distance(p, min_element(p+i, p+nb));
		if(minIdx != i) //switch blk order in B, UV, a, b, and oldIdx
		{
			int offset = i*m;
			swap_ranges(oldIdx + offset, oldIdx + offset + m, oldIdx +
				minIdx*m);
			swap_ranges(a + offset, a + offset + m, a + minIdx*m);
			swap_ranges(b + offset, b + offset + m, b + minIdx*m);
			for(int j = 0; j < nb; j++)
			{
				if(j < i)
				{
					int idx1 = (2*nb-1-j)*j/2 + (i-1-j);
					int idx2 = (2*nb-1-j)*j/2 + (minIdx-1-j);
					swap(UV[idx1], UV[idx2]);
				}
				else if(j == i)
				{
					int idx1 = (2*nb-1-i)*i/2 + (minIdx-1-i);
					swap(UV[idx1].U, UV[idx1].V);
				}
				else if(j < minIdx)
				{
					int idx1 = (2*nb-1-j)*j/2 + (minIdx-1-j);
					int idx2 = (2*nb-1-i)*i/2 + (j-1-i);
					swap(UV[idx1], UV[idx2]);
					swap(UV[idx1].U, UV[idx1].V);
					swap(UV[idx2].U, UV[idx2].V);
				}
				else if(j == minIdx)
				{
					swap(B[i], B[minIdx]);
				}else
				{
					int idx1 = (2*nb-1-i)*i/2 + (j-1-i);
					int idx2 = (2*nb-1-minIdx)*minIdx/2 + (j-1-
						minIdx);
					swap(UV[idx1], UV[idx2]);
				}
			} // j = 0 : nb
		} // minIdx != i
		// uni_reorder again
		int offset = i*m;
		iota(idxj, idxj+m, 0);
		int fail = uni_reorder(m, B[i].data(), m, a+offset, b+offset, p[i], 
			y+offset, idxj, subworkDbl, lsubworkDbl);
		if(fail)
			Rcpp::stop("uni_reorder failed after reordering diag blks "
				"called by recur_blk_reorder\n");
		reorder(oldIdx + offset, idxj, m, subworkInt, lsubworkInt);
		// reorder in i-th blk
		for(int j = 0; j < nb; j++)
		{
			if(j < i)
			{
				int idx1 = (2*nb-1-j)*j/2 + (i-1-j);	
				reorder_row(UV[idx1].U, idxj, subworkInt, lsubworkInt,
					UV[idx1].crtColNum);
			}
			else if(j > i)
			{
				int idx1 = (2*nb-1-i)*i/2 + (j-1-i);
				reorder_row(UV[idx1].V, idxj, subworkInt, lsubworkInt,
					UV[idx1].crtColNum);
			}
		}
		// Cholesky
		int success_code;
		copy_n(B[i].data(), m*m, Bi);
		F77_CALL(dtrtri)("L", "N", &m, Bi, &m, &success_code); // upper tri is
			// not referenced
		if(success_code)
			Rcpp::stop("Matrix inverse failed in recursive block "
				"reordering\n");
		for(int ii = i+1; ii < nb ; ii++)
		{
			int idx1 = (2*nb-1-i)*i/2 + (ii-1-i);
			// trailing blk in the same col
			double alpha = 1.0;
			F77_CALL(dtrmm)("L", "L", "N", "N", &m, &UV[idx1].crtColNum, 
				&alpha, Bi, &m, UV[idx1].V.data(), &m);
			// diag blk in Schur complement
			if(UV[idx1].crtColNum > m)
				Rcpp::stop("UV[%d]'s column number is greater than "
					"m\n", idx1);
			///////////////////////////
			double beta = 0.0;
			// Bj = V[idx1]^T V[idx1]
			F77_CALL(dsyrk)("L", "T", &UV[idx1].crtColNum, &m, &alpha, 
				UV[idx1].V.data(), &m, &beta, Bj, &m);
			// BSchur = U[idx1] Bj
			F77_CALL(dsymm)("R", "L", &m, &UV[idx1].crtColNum, &alpha, Bj,
				&m, UV[idx1].U.data(), &m, &beta, BSchur, &m);
			// B[ii] = B[ii] - BSchur * U[idx1].trans
			alpha = -1.0;
			beta = 1.0;
			F77_CALL(dgemm)("N", "T", &m, &m, &UV[idx1].crtColNum, &alpha,
				BSchur, &m, UV[idx1].U.data(), &m, &beta, B[ii].
				data(), &m);
			// integration limits
//			alpha = 1.0;
//			beta = 0.0;
//			int incx = 1;
//			int incy = 1;
//			F77_CALL(dgemv)("T", &m, &UV[idx1].crtColNum, &alpha, UV[
//				idx1].V.data(), &m, y+offset, &incx, &beta, ySchur,
//				&incy);
//			F77_CALL(dgemv)("N", &m, &UV[idx1].crtColNum, &alpha, UV[
//				idx1].U.data(), &m, ySchur, &incx, &beta, ySchur2,
//				&incy);
//			transform(a+ii*m, a+ii*m+m, ySchur2, a+ii*m, [](double &aCoef,
//				double &yCoef){return aCoef - yCoef;});
//			transform(b+ii*m, b+ii*m+m, ySchur2, b+ii*m, [](double &bCoef,
//				double &yCoef){return bCoef - yCoef;});
			///////////////////////////////////
			// update B[ii]
//			B[ ii ] -= UV[idx1].U.block(0,0,m,UV[idx1].crtColNum) * 
//				(UV[idx1].V.block(0,0,m,UV[idx1].crtColNum).
//				transpose() * 
//				UV[idx1].V.block(0,0,m,UV[idx1].crtColNum)) * 
//				UV[idx1].U.block(0,0,m,UV[idx1].crtColNum).
//				transpose();
                        // integration limits
			Map<VectorXd> yMin(y+offset,m);
                        VectorXd tmpVec = UV[idx1].U.block(0,0,m,UV[idx1].crtColNum)*
				(UV[idx1].V.block(0,0,m,UV[idx1].crtColNum).
				transpose()*yMin);
			transform(a+ii*m, a+ii*m+m, tmpVec.data(), a+ii*m, 
				[](double &aCoef, double &yCoef)
				{return aCoef - yCoef;});
			transform(b+ii*m, b+ii*m+m, tmpVec.data(), b+ii*m, 
				[](double &bCoef, double &yCoef)
				{return bCoef - yCoef;});
		} // ii = i+1 : nb
		// update Schur complement
		int dimRemain = nb - 1 - i;
		int nbRemain = dimRemain * ( dimRemain - 1 ) / 2;
		int idxOffset = (2*nb - 1 - (i+1)) * (i+1) / 2;
		for(int l = 0; l < nbRemain; l++)
		{
			int iRemain = l % dimRemain;
			int jRemain = l / dimRemain;
			if(iRemain <= jRemain)
			{
				jRemain = dimRemain - jRemain - 2;
				iRemain = dimRemain - iRemain - 1;
			}
			int idx = idxOffset + (2*dimRemain-1-jRemain)*jRemain/2 +
				(iRemain-1-jRemain);
			int idx1 = (2*nb-1-i)*i/2 + iRemain;
			int idx2 = (2*nb-1-i)*i/2 + jRemain;
			if(UV[idx1].crtColNum > m)
				Rcpp::stop("UV[%d]'s column number is greater than "
					"m\n", idx1);
			if(UV[idx2].crtColNum > m)
				Rcpp::stop("UV[%d]'s column number is greater than "
					"m\n", idx2);
			double *Ucmb = Bi;
			double *Utmp = Bj;
			double alpha = 1.0;
			double beta = 0.0;
			// Utmp = V[idx1]^T V[idx2]
			F77_CALL(dgemm)("T", "N", &UV[idx1].crtColNum, &UV[idx2].
				crtColNum, &m, &alpha, UV[idx1].V.data(), &m, 
				UV[idx2].V.data(), &m, &beta, Utmp, &m);
			alpha = -1.0;
			// Ucmb = - U[idx1] Utmp
			F77_CALL(dgemm)("N", "N", &m, &UV[idx2].crtColNum, &UV[idx1].
				crtColNum, &alpha, UV[idx1].U.data(), &m, Utmp, &m,
				&beta, Ucmb, &m);

			tlr_tlr_add_qr(UV[idx], Ucmb, UV[idx2].U.data(), UV[idx2].
				crtColNum, epsl, subworkDbl, lsubworkDbl);
		}
	} // i = 0 : nb
	return 0;
}

/*
        The old version of recursive block reordering moved from the rtlrmvn project
        Used here for checking purpose
	lwork should be no smaller than 16*max(nrowBlk, ncolBlk)^2 + 12*max(nrowBlk,
		ncolBlk)
        2019/09/17
*/
//std::vector<int> recur_blk_reorder_old(std::vector<Eigen::MatrixXd> &B, 
//	std::vector<Eigen::MatrixXd> &U, std::vector<Eigen::MatrixXd> &V, 
//	const Eigen::VectorXd &a, const Eigen::VectorXd &b, double epsl, double 
//	*work, int lwork)
//{
//        int num_blk = B.size();
//        int blk_sz = B[0].rows();
//        VectorXd a_cp = a;
//        VectorXd b_cp = b;
//        vector<int> odr(num_blk * blk_sz);
//        iota(odr.begin(),odr.end(),0);
//        for(int i = 0 ; i < num_blk ; i++)
//        {
//                int idx_min = -1;
//                double p_min = 2.0;
//                VectorXd y_min;
//                MatrixXd B_min;
//                vector<int> odr_min;
//                // find min
//                for(int j = i ; j < num_blk ; j++)
//                {
//                        MatrixXd B_j = B[j];
//                        VectorXd a_j = a_cp.segment(j*blk_sz,blk_sz);
//                        VectorXd b_j = b_cp.segment(j*blk_sz,blk_sz);
//                        double p_j;
//                        VectorXd y_j;
//                        vector<int> odr_j = uni_reorder_old(B_j , a_j , b_j , p_j , 
//				y_j);
//                        if(p_j < p_min)
//                        {
//                                idx_min = j;
//                                p_min = p_j;
//                                y_min = y_j;
//                                B_min = B_j;
//                                odr_min = odr_j;
//                        }
//                }
//                // in case all p_j are na
//                if(idx_min == -1)
//                {
//                        y_min = VectorXd::Zero(blk_sz);
//                }
//                // block-wise switch order
//                else if(idx_min != i)
//                {
//                        // odr switch
//                        for(int l = 0 ; l < blk_sz ; l++)
//                        {
//                                int tmp_int = odr[i*blk_sz+l];
//                                odr[i*blk_sz+l] = odr[idx_min*blk_sz+l];
//                                odr[idx_min*blk_sz+l] = tmp_int;
//                        }
//                        for(int j = 0 ; j < num_blk ; j++)
//                        {
//                                if(j < i)
//                                {
//                                        int idx1 = (2 * num_blk - 1 - j) * j / 2 +
//                                                ( i - 1 - j );
//                                        int idx2 = (2 * num_blk - 1 - j) * j / 2 +
//                                                ( idx_min - 1 - j );
//                                        MatrixXd tmp_mat = U[idx1];
//                                        U[idx1] = U[idx2];
//                                        U[idx2] = tmp_mat;
//                                        tmp_mat = V[idx1];
//                                        V[idx1] = V[idx2];
//                                        V[idx2] = tmp_mat;
//                                }
//                                if(j == i)
//                                {
//                                        int idx = (2 * num_blk - 1 - i) * i / 2 +
//                                                ( idx_min - 1 - i );
//                                        MatrixXd tmp_mat = U[idx];
//                                        U[idx] = V[idx];
//                                        V[idx] = tmp_mat;
//                                }
//                                if(i < j && j < idx_min)
//                                {
//                                        int idx1 = (2 * num_blk - 1 - j) * j / 2 +
//                                                ( idx_min - 1 - j );
//                                        int idx2 = (2 * num_blk - 1 - i) * i / 2 +
//                                                ( j - 1 - i );
//                                        MatrixXd tmp_mat = U[idx1];
//                                        U[idx1] = V[idx2];
//                                        V[idx2] = tmp_mat;
//                                        tmp_mat = V[idx1];
//                                        V[idx1] = U[idx2];
//                                        U[idx2] = tmp_mat;
//                                }
//                                if(j == idx_min)
//                                {
//                                        B[idx_min] = B[i];
//                                }
//                                if(j > idx_min)
//                                {
//                                        int idx1 = (2 * num_blk - 1 - i) * i / 2 +
//                                                ( j - 1 - i );
//                                        int idx2 = (2 * num_blk - 1 - idx_min) *
//                                                idx_min / 2 + ( j - 1 - idx_min );
//                                        MatrixXd tmp_mat = U[idx1];
//                                        U[idx1] = U[idx2];
//                                        U[idx2] = tmp_mat;
//                                        tmp_mat = V[idx1];
//                                        V[idx1] = V[idx2];
//                                        V[idx2] = tmp_mat;
//                                }
//                        }
//                        {
//                                VectorXd tmp_vec = a_cp.segment(i*blk_sz,blk_sz);
//                                a_cp.segment(i*blk_sz,blk_sz) =
//                                        a_cp.segment(idx_min*blk_sz,blk_sz);
//                                a_cp.segment(idx_min*blk_sz,blk_sz) =
//                                        tmp_vec;
//                                tmp_vec = b_cp.segment(i*blk_sz,blk_sz);
//                                b_cp.segment(i*blk_sz,blk_sz) =
//                                        b_cp.segment(idx_min*blk_sz,blk_sz);
//                                b_cp.segment(idx_min*blk_sz,blk_sz) =
//                                        tmp_vec;
//                        }
//                }
//
//
//                // within-blk reorder
//                if(idx_min != -1)
//                {
//                        vector<int> tmp_vec(odr.begin()+i*blk_sz,
//                                odr.begin()+i*blk_sz + blk_sz);
//                        for(int l = 0 ; l < blk_sz ; l++)
//                                odr[i*blk_sz+l] = tmp_vec[odr_min[l]];
//                }
//                if(idx_min != -1)
//                {
//                        for(int j = 0 ; j < num_blk ; j++)
//                        {
//                                if(j < i)
//                                {
//                                        int idx = (2 * num_blk - 1 - j) * j / 2 +
//                                                (i - 1 - j);
//                                        MatrixXd U_tmp = U[idx];
//                                        for(int l = 0 ; l < blk_sz ; l++)
//                                                U[idx].row(l) =
//                                                        U_tmp.row(odr_min[l]);
//                                }
//                                if(j == i)
//                                        B[i] = B_min;
//                                if(j > i)
//                                {
//                                        int idx = (2 * num_blk - 1 - i) * i / 2 +
//                                                (j - 1 - i);
//                                        MatrixXd V_tmp = V[idx];
//                                        for(int l = 0 ; l < blk_sz ; l++)
//                                                V[idx].row(l) =
//                                                        V_tmp.row(odr_min[l]);
//                                }
//                        }
//                }else
//                {
//                        int success_code = cholesky(B[i]);
//                        assert(success_code == 0);
//                        B_min = B[i];
//                }
//
//                // Cholesky routine
//		inverse(B_min);
//                B_min.triangularView<StrictlyUpper>().setZero();
//                for(int j = i + 1 ; j < num_blk ; j++)
//                {
//                        // update trailing tiles
//                        int idx = (2 * num_blk - 1 - i) * i / 2 + ( j - 1 - i );
//                        V[idx] = B_min * V[ idx ];
//                        // update diag tiles
//                        B[ j ] = B[ j ] - U[idx] * ( V[idx].transpose() *
//                        V[idx] ) * U[idx].transpose();
//                        // update integration limits
//                        VectorXd tmp_vec = U[idx]*(V[idx].transpose()*y_min);
//                        a_cp.segment(j*blk_sz,blk_sz) = a_cp.segment(
//                        j*blk_sz,blk_sz) - tmp_vec;
//                        b_cp.segment(j*blk_sz,blk_sz) = b_cp.segment(
//                        j*blk_sz,blk_sz) - tmp_vec;
//                }
//                // update other tiles
//                int dim_remain = num_blk - 1 - i;
//                int num_blk_remain = dim_remain * ( dim_remain - 1 ) / 2;
//                int idx_offset = (2 * num_blk - 1 - (i+1)) * (i+1) / 2;
//                for( int l = 0 ; l < num_blk_remain ; l++ )
//                {
//                        int i_remain = l % dim_remain;
//                        int j_remain = l / dim_remain;
//                        if( i_remain <= j_remain )
//                        {
//                                j_remain = dim_remain - j_remain - 2;
//                                i_remain = dim_remain - i_remain - 1;
//                        }
//                        int idx = idx_offset +
//                                ( 2 * dim_remain - 1 - j_remain ) * j_remain / 2 +
//                                ( i_remain - 1 - j_remain );
//                        int idx1 = (2 * num_blk - 1 - i) * i / 2 + i_remain;
//                        int idx2 = (2 * num_blk - 1 - i) * i / 2 + j_remain;
//                        MatrixXd U2 = U[ idx1 ] *
//                                ( V[ idx1 ].transpose() * V[ idx2 ] );
//                        MatrixXd V2 = - U[ idx2 ];
//                        TLRNode nodeNew;
//			nodeNew.U = U[idx];
//			nodeNew.V = V[idx];
//			nodeNew.crtColNum = U[idx].cols();
//			nodeNew.maxColNum = U[idx].cols();
//                        tlr_tlr_add_qr(nodeNew, U2.data(), V2.data() , U2.cols(), 
//				epsl, work, lwork);
//                        U[idx] = nodeNew.U.block(0,0,U2.rows(),nodeNew.crtColNum);
//                        V[idx] = nodeNew.V.block(0,0,V2.rows(),nodeNew.crtColNum);
//                }
//        }
//        return odr;
//}
//
