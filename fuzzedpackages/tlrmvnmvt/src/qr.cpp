#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include <R_ext/Print.h>
#include "qr.h"

using namespace std;
using namespace Eigen;

/*
        qr factorization
        Factorization is performed to a leading block of A of size
                ARowNum X AColNum
        Resize R is R.rows() < AColNum or R.cols() < AColNum
        tau should be a double array of size at least AColNum
        work should be a double array of size at least AColNum
        lwork is the dimension of work
        2019/08/26
*/
int qr(Eigen::MatrixXd &A , int ARowNum , int AColNum , Eigen::MatrixXd &R ,
        Eigen::VectorXd &tau , double *work , int lwork)
{
        if(ARowNum < AColNum)
                Rcpp::stop("Number of rows should be no smaller than "
                        "the number of columns for the QR decomp");
        if(R.rows() < AColNum || R.cols() < AColNum)
        {
                Rprintf("Warning: R needs resizing in qr factorization. "
                        "Consider enlarge their initial allocations\n");
                R.resize(AColNum , AColNum);
        }
        if(tau.size() < AColNum)
        {
                Rprintf("Warning: tau needs resizing in qr factorization. "
                        "Consider enlarge their initial allocations\n");
                tau.resize(AColNum);
        }
        int fail; 
        int ldA = A.rows();
        F77_CALL(dgeqrf)(&ARowNum , &AColNum  , A.data() , &ldA , tau.data() , work ,
                &lwork , &fail);
        if(fail)
                return fail;
        for(int i = 0 ; i < AColNum*AColNum ; i++)
        {
                unsigned int row_num = i%AColNum;
                unsigned int col_num = i/AColNum;
                R(row_num , col_num) =
                        col_num < row_num ? 0 : A(row_num , col_num);
        }
        F77_CALL(dorgqr)(&ARowNum , &AColNum , &AColNum , A.data() , &ARowNum , tau.
                data() , work , &lwork , &fail);
        if(fail)
                return fail;
        return 0;
}
