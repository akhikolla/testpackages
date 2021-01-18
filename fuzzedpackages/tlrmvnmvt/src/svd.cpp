#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/Print.h>
#include "svd.h"

using namespace std;
using namespace Eigen;

/*
        svd decomposition
        Decomposition is performed to a leading block of A of size
                ARowNum X AColNum
        work should be of dimension larger than 5 X min(ARowNum,AColNum)
        lwork is the dimension of work
*/
int svd(Eigen::MatrixXd &A , int ARowNum , int AColNum , Eigen::MatrixXd &U ,
        Eigen::MatrixXd &VT , Eigen::VectorXd &S , double *work , int lwork)
{
        int ldA = A.rows();
        int rk = min(ARowNum , AColNum);
        if(U.rows() < ARowNum || U.cols() < rk)
        {
                Rprintf("Warning: the U factor in svd decomposition is resized. "
                        "Increasing the allocation for U will improve performance\n");
                U.resize(ARowNum , rk);
        }
        if(VT.rows() < rk || VT.cols() < AColNum)
        {
                Rprintf("Warning: the VT factor in svd decomposition is resized. "
                        "Increasing the allocation for VT will improve "
                        "performance\n");
                VT.resize(rk , AColNum);
        }
        if(S.size() < rk)
        {
                Rprintf("Warning: the S factor in svd decomposition is resized. "
                        "Increasing the allocation for S will improve performance\n");
                S.resize(rk);
        }
        int ldU = U.rows();
        int ldVT = VT.rows();
        int fail;
        F77_CALL(dgesvd)("S" , "S" , &ARowNum , &AColNum , A.data() , &ldA ,
                S.data() , U.data() , &ldU , VT.data() , &ldVT , work , &lwork ,
                &fail);
        return fail;
}
