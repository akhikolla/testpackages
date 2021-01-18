#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include "uni_reorder.h"
#include "mvphi.h"

using namespace std;
using namespace Eigen;

/*
	Block reordering
	All block sizes should be the same
	The size of workDbl should be at least 6m
	The size of workInt should be at least n+nb
	2019/09/09
*/
void blk_reorder(std::vector<Eigen::MatrixXd> &B, double *a, double *b, double *p,
	double *y, int *oldIdx, double *workDbl, int lworkDbl, int *workInt, int 
	lworkInt)
{
	int nb = B.size();
	int m = B[0].rows();
	int n = nb * m;
	int offset = 0;
	if(lworkDbl < 6*m)
		Rcpp::stop("Dimension of workDbl is insufficient\n");
	if(lworkInt < n + nb)
		Rcpp::stop("Dimension of workInt is insufficient\n");
	for (int i = 0; i < nb; ++i)
	{
		uni_reorder(m, B[i].data(), m, a+offset, b+offset, p[i], y+offset,
			oldIdx+offset, workDbl, lworkDbl);
		offset += m;
	}
	int *blkOrder = workInt;
	iota(blkOrder, blkOrder+nb, 0);
	sort(blkOrder, blkOrder+nb, [&p](int i1, int i2) { return p[i1] < p[i2]; });
	int *oldIdxCp = workInt + nb;
	copy_n(oldIdx, n, oldIdxCp);
	offset = 0;
	for (int i = 0; i < nb; ++i)
	{
		copy_n(oldIdxCp+blkOrder[i]*m, m, oldIdx+offset);
		offset += m;
	}
}

/*
        The old version of block reordering moved from the rtlrmvn project
        Used here for checking purpose
        2019/05/20
*/
//std::vector<int> blk_reorder_old(std::vector<Eigen::MatrixXd> &B, const 
//	Eigen::VectorXd &a, const Eigen::VectorXd &b, int bsz)
//{
//    int bnum = B.size();
//    VectorXd probs(bnum);
//    vector<int> ind(bnum * B[0].rows());
//    int i1 = 0;
//    for(int i = 0; i < bnum; i++)
//    {
//        VectorXd a_ = a.segment(i1,bsz);
//        VectorXd b_ = b.segment(i1,bsz);
//        vector<int> subInd(bsz);
//        iota(subInd.begin(),subInd.end(),i1);
//
//        MatrixXf::Index minCoefPos;
//        int tempInt;
//        double tempDbl,tempLowerb,tempUpperb;
//        double prob = 1.0;
//        VectorXd y = VectorXd::Zero(bsz);
//        VectorXd tempVec(bsz),tempVecRow(bsz);
//        for(int j = 0; j < bsz; j++)
//        {
//            VectorXd s(bsz-j),t(bsz-j);
//            if(j > 0)
//            {
//                y(j-1) = ((exp(-tempLowerb * tempLowerb / 2.0) -
//                 exp(-tempUpperb * tempUpperb / 2.0)) / sqrt(2.0*M_PI)) / tempDbl;
//                s = B[i].block(j,0,bsz-j,j) * y.segment(0,j);
//                t = B[i].block(j,0,bsz-j,j).cwiseAbs2().rowwise().sum();
//            }else
//            {
//                s = VectorXd::Zero(bsz-j);
//                t = VectorXd::Zero(bsz-j);
//            }
//            VectorXd lowerBounds = ((a_.segment(j,bsz-j) - s).array() /
//             sqrt((B[i].diagonal().segment(j,bsz-j) - t).array())).matrix();
//            VectorXd upperBounds = ((b_.segment(j,bsz-j) - s).array() /
//             sqrt((B[i].diagonal().segment(j,bsz-j) - t).array())).matrix();
//            VectorXd cdfATlowerb(bsz-j);
//            VectorXd cdfATupperb(bsz-j);
//            lc_vdCdfNorm(bsz-j,lowerBounds.data(),cdfATlowerb.data());
//            lc_vdCdfNorm(bsz-j,upperBounds.data(),cdfATupperb.data());
//            VectorXd probVec = cdfATupperb - cdfATlowerb;
//            probVec.minCoeff(&minCoefPos);
//            minCoefPos += j;
//            tempInt = subInd[j];
//            subInd[j] = subInd[minCoefPos];
//            subInd[minCoefPos] = tempInt;
//            tempVecRow = B[i].row(j);
//            B[i].row(j) = B[i].row(minCoefPos);
//            B[i].row(minCoefPos) = tempVecRow;
//            tempVec = B[i].col(j);
//            B[i].col(j) = B[i].col(minCoefPos);
//            B[i].col(minCoefPos) = tempVec;
//            tempDbl = a_(j);
//            a_(j) = a_(minCoefPos);
//            a_(minCoefPos) = tempDbl;
//            tempDbl = b_(j);
//            b_(j) = b_(minCoefPos);
//            b_(minCoefPos) = tempDbl;
//
//            if(j > 0)
//            {
//                tempDbl = sqrt(B[i](j,j) - B[i].row(j).segment(0,j).cwiseAbs2().sum());
//                B[i].col(j).segment(j,bsz-j) = (B[i].col(j).segment(j,bsz-j) -
//                 B[i].block(j,0,bsz-j,j)*B[i].row(j).segment(0,j).transpose()) / tempDbl;
//            }else
//            {
//                tempDbl = sqrt(B[i](j,j));
//                B[i].col(j).segment(j,bsz-j) = B[i].col(j).segment(j,bsz-j) / tempDbl;
//            }
//
//            if(j > 0)
//            {
//                tempLowerb = (a_(j) - B[i].row(j).segment(0,j) * y.segment(0,j)) / B[i](j,j);
//                tempUpperb = (b_(j) - B[i].row(j).segment(0,j) * y.segment(0,j)) / B[i](j,j);
//            }else
//            {
//                tempLowerb = a_(j) / B[i](j,j);
//                tempUpperb = b_(j) / B[i](j,j);
//            }
//            lc_vdCdfNorm(1,&tempLowerb,cdfATlowerb.data());
//            lc_vdCdfNorm(1,&tempUpperb,cdfATupperb.data());
//            tempDbl = cdfATupperb(0) - cdfATlowerb(0);
//            prob *= tempDbl;
//        }
//        copy_n(subInd.begin(),bsz,ind.begin()+i1);
//        probs(i) = prob;
//        i1 += bsz;
//    }
//
//    vector<int> blkInd(bnum);
//    iota(blkInd.begin(),blkInd.end(),0);
//    sort(blkInd.begin(),blkInd.end(), [&probs](int i1, int i2) {return probs[i1] < probs[i2];});
//    vector<int> indCopy = ind;
//    for(int i = 0; i < bnum; i++)
//    {
//        int offset_from = bsz * blkInd[i];
//        int offset_to = bsz * i;
//        copy_n(indCopy.begin()+offset_from,bsz,ind.begin()+offset_to);
//    }
//
//    return ind;
//}
//
//
