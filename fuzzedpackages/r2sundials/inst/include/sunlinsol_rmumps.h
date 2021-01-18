// cf sunlinsol_rmumps.cpp for legal information

#ifndef _SUNLINSOL_RMUMPS_H
#define _SUNLINSOL_RMUMPS_H

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h> // access to serial N_Vector
#include <sunmatrix/sunmatrix_sparse.h>
#include <rmumps.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#define RMUMPS_CONTENT(S)     ( (SUNLinearSolverContent_RMUMPS)(S->content) )
#define LASTFLAG(S)        ( RMUMPS_CONTENT(S)->last_flag )
#define RMU(S)          ( *(RMUMPS_CONTENT(S)->rmu) )

/* Default RMUMPS solver parameters */
#define SUNRMUMPS_ORDERING_DEFAULT  "auto"

/* Interfaces to match 'sunindextype' with the correct RMUMPS types/functions */
/*
#if defined(SUNDIALS_INT64_T)
#error Incompatible SUNDIALS int type for rmumps
#elif !defined(SUNDIALS_INT32_T)
#error  Incompatible SUNDIALS int type for rmumps
#endif
*/

#if defined(SUNDIALS_DOUBLE_PRECISION)
#else
#error  Incompatible realtype for RMUMPS
#endif
 
struct _SUNLinearSolverContent_RMUMPS {
  long int last_flag;
  XPtr<Rmumps> *rmu;
  Col<MUMPS_INT> *irp;
  Col<MUMPS_INT> *jcp;
};

#undef SUNDIALS_EXPORT
#define SUNDIALS_EXPORT
typedef struct _SUNLinearSolverContent_RMUMPS *SUNLinearSolverContent_RMUMPS;
SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_RMUMPS(N_Vector y, SUNMatrix A, int permutation);
SUNDIALS_EXPORT int SUNLinSol_RMUMPSSetOrdering(SUNLinearSolver S, std::string ordering_choice);
  
/*
 * -----------------------------------------------------------------
 * RMUMPS implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_RMUMPS(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_RMUMPS(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_RMUMPS(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_RMUMPS(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_RMUMPS(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolFree_RMUMPS(SUNLinearSolver S);
#ifdef __cplusplus
}
#endif

#endif
