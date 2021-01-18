#ifndef RNEWMAT_H_
#define RNEWMAT_H_

#include <R.h> // in /usr/lib/R/include. For R
#include <Rinternals.h>
#include <Rmath.h>

#include <newmat.h>  // For Newmat
#include <newmatap.h> // for cholesky decomposition

#include <set>   


ReturnMatrix
getMultipleCols(const Matrix&, const std::set<int>&); // get different concatenated columns of matrix
 
#endif /*RNEWMAT_H_*/
