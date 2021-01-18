#ifndef _AdMatAlgFoo_h
#define _AdMatAlgFoo_h

#include "RcppArmadillo.h"
using namespace arma;

const double MINEGVAL = std::numeric_limits<double>::min();   //  minimum value of smallest eigenvalue for a square symmetric matrix to be considered numerically positive-definite 
const double MAXLNK2 = log(1e12);                              //  maximum (l2 norm) condition number for a matrix to be considered numerically non-singular 
const double MINLNDET = -500;                                 //  minimum log-determinant for a matrix to be considered numerically non-singular 


bool pdsolve(const mat& M, mat& MInv, double* logDet);
bool pdsolve(const mat& M, mat& rhs, mat& res, double* logDet);  // Try replacing pointer to logDet by a constant passed by reference
void SetZero(vec& v,const int n,bool cheksize);
void SetZero(mat& M,const int m,const int n,bool cheksize);
bool chcksing(const mat& M, double& logDet, double& viol, double& minegv, double& maxegv,
  const double minlndet=MINLNDET, const double maxlnk2=MAXLNK2, const bool scale=false);
bool safepdsolve(const mat& M, mat& MInv, double& logDet, double& viol, double& minegv, double& maxegv, 
  const double minlndet=MINLNDET, const double maxlnk2=MAXLNK2, const bool scale=false); 
bool safepdsolve(const mat& M, const mat& rhs, mat& res, double& logDet, double& viol, double& minegv, double& maxegv,  
  const double minlndet=MINLNDET, const double maxlnk2=MAXLNK2, const bool scale=false);

#endif




