// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


#ifdef _OPENMP
#include <omp.h>
//   [[Rcpp::plugins(openmp)]]
#endif







// internal function to remove a row from a dynamic matrix
void removeRow(Eigen::MatrixXd& matrix, unsigned long rowToRemove)
{
    unsigned long numRows = matrix.rows()-1;
    unsigned long numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}


