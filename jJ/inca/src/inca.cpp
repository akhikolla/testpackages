#if defined _OPENMP
  #if (_OPENMP > 200800)
    #include <omp.h>
    #define __VOPENMP 1
  #else
    #define __VOPENMP 0
  #endif
#else
  #define __VOPENMP 0
#endif

#include <RcppArmadillo.h>
#include <Rcpp.h>

void R_init_inca(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

using namespace arma;
using namespace Rcpp;

#ifndef _DEBUG // Debug modality (usaully passed via gpp -D_DEBUG=1)
    #define _DEBUG 0
/* DEFINITION OF "_DEBUG" preprocessor expression
 *   0 : no log is printed
 *   1 : the optimization-log is printed
 *   2 : a colorful log is printed (with colors based on "ANSI escape code")
 */
#endif
#if _DEBUG
    #ifndef _POSCHK // Position Check
        #define _POSCHK !0
/* DEFINITION OF "_POSCHK" preprocessor expression
 *   !0 : it prints (only when _DEBUG != 0) the positions that changed (df, Grad, Wts)
 *    0 : it does not print the positions that changed
 */
    #endif // _POSCHK
#endif // _DEBUG

#include "bestRound.h"
#include "intCalib.h"
//#include "qininCalib.h"

colvec de_round(const mat& A, const colvec& y, colvec& w,
             const mat& Bnds, const vec& scale,
             const std::string lossType) {
    // Optimal Rounding for Dense Matrices
    return bestRound(A, y, w, Bnds, scale, lossType);
}

RcppExport SEXP dense_round(SEXP ASEXP, SEXP ySEXP, SEXP wSEXP, SEXP BndsSEXP, SEXP scaleSEXP, SEXP lossTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Bnds(BndsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const std::string >::type lossType(lossTypeSEXP);
    __result = Rcpp::wrap(de_round(A, y, w, Bnds, scale, lossType));
    return __result;
END_RCPP
}

colvec sp_round(const sp_mat& A, const colvec& y, colvec& w,
             const mat& Bnds, const vec& scale,
             const std::string lossType) {
    // Optimal Rounding for Sparse Matrices
    return bestRound(A, y, w, Bnds, scale, lossType);
}

RcppExport SEXP sparse_round(SEXP ASEXP, SEXP ySEXP, SEXP wSEXP, SEXP BndsSEXP, SEXP scaleSEXP, SEXP lossTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Bnds(BndsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const std::string >::type lossType(lossTypeSEXP);
    __result = Rcpp::wrap(sp_round(A, y, w, Bnds, scale, lossType));
    return __result;
END_RCPP
}

colvec sp_ipc(const sp_mat& A, const colvec& y, colvec& w,
              const vec lower, const vec upper, const mat& Bnds,
              const vec& scale, const std::string lossType) {
    // Integer Programming Calibration for Sparse Matrices
    int lt = 0;
    if(lossType == "L1") lt = LOSS_L1;
    if(lossType == "aL1") lt = LOSS_aL1;
    if(lossType == "rL1") lt = LOSS_rL1;
    if(lossType == "LB1") lt = LOSS_LB1;
    if(lossType == "rB1") lt = LOSS_rB1;
    if(lossType == "rbLasso1") lt = LOSS_RBLASSO1;
    if(lossType == "L2") lt = LOSS_L2;
    if(lossType == "aL2") lt = LOSS_aL2;
    if(lossType == "rL2") lt = LOSS_rL2;
    if(lossType == "LB2") lt = LOSS_LB2;
    if(lossType == "rB2") lt = LOSS_rB2;
    if(lossType == "rbLasso2") lt = LOSS_RBLASSO2;
    colvec dse(IS_RBLASSO(lt) ? w.size() : 1);
    if (IS_RBLASSO(lt)) dse = w + 0.0;
    w = bestRound(A, y, w, Bnds, scale, lossType);
    return IntProgCalib(A, y, w, dse, lower, upper, Bnds, scale, lossType);
}

RcppExport SEXP sparse_ipc(SEXP ASEXP, SEXP ySEXP, SEXP wSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP BndsSEXP, SEXP scaleSEXP, SEXP lossTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const vec >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Bnds(BndsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const std::string >::type lossType(lossTypeSEXP);
    __result = Rcpp::wrap(sp_ipc(A, y, w, lower, upper, Bnds, scale, lossType));
    return __result;
END_RCPP
}

colvec de_ipc(const mat& A, const colvec& y, colvec& w,
              const vec lower, const vec upper, const mat& Bnds,
              const vec& scale, const std::string lossType) {
    // Integer Programming Calibration for Sparse Matrices
    int lt = 0;
    if(lossType == "L1") lt = LOSS_L1;
    if(lossType == "aL1") lt = LOSS_aL1;
    if(lossType == "rL1") lt = LOSS_rL1;
    if(lossType == "LB1") lt = LOSS_LB1;
    if(lossType == "rB1") lt = LOSS_rB1;
    if(lossType == "rbLasso1") lt = LOSS_RBLASSO1;
    if(lossType == "L2") lt = LOSS_L2;
    if(lossType == "aL2") lt = LOSS_aL2;
    if(lossType == "rL2") lt = LOSS_rL2;
    if(lossType == "LB2") lt = LOSS_LB2;
    if(lossType == "rB2") lt = LOSS_rB2;
    if(lossType == "rbLasso2") lt = LOSS_RBLASSO2;
    colvec dse(IS_RBLASSO(lt) ? w.size() : 1);
    if (IS_RBLASSO(lt)) dse = w + 0.0;
    w = bestRound(A, y, w, Bnds, scale, lossType);
    return IntProgCalib(A, y, w, dse, lower, upper, Bnds, scale, lossType);
}

RcppExport SEXP dense_ipc(SEXP ASEXP, SEXP ySEXP, SEXP wSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP BndsSEXP, SEXP scaleSEXP, SEXP lossTypeSEXP) {
    BEGIN_RCPP
        Rcpp::RObject __result;
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const mat& >::type A(ASEXP);
        Rcpp::traits::input_parameter< const colvec& >::type y(ySEXP);
        Rcpp::traits::input_parameter< colvec& >::type w(wSEXP);
        Rcpp::traits::input_parameter< const vec >::type lower(lowerSEXP);
        Rcpp::traits::input_parameter< const vec >::type upper(upperSEXP);
        Rcpp::traits::input_parameter< const mat& >::type Bnds(BndsSEXP);
        Rcpp::traits::input_parameter< const vec& >::type scale(scaleSEXP);
        Rcpp::traits::input_parameter< const std::string >::type lossType(lossTypeSEXP);
        __result = Rcpp::wrap(de_ipc(A, y, w, lower, upper, Bnds, scale, lossType));
        return __result;
    END_RCPP
}
/*
colvec sp_qiipc(const sp_mat& A, const colvec& y, colvec& w,
                const vec lower, const vec upper, const mat& Bnds,
                const vec& scale, const size_t nsim, const double memrate,
                const std::string lossType) {
    // Quantum Inspired Integer Programming Calibration for Sparse Matrices
    int lt = 0;
    if(lossType == "L1") lt = LOSS_L1;
    if(lossType == "aL1") lt = LOSS_aL1;
    if(lossType == "rL1") lt = LOSS_rL1;
    if(lossType == "LB1") lt = LOSS_LB1;
    if(lossType == "rB1") lt = LOSS_rB1;
    if(lossType == "rbLasso1") lt = LOSS_RBLASSO1;
    if(lossType == "L2") lt = LOSS_L2;
    if(lossType == "aL2") lt = LOSS_aL2;
    if(lossType == "rL2") lt = LOSS_rL2;
    if(lossType == "LB2") lt = LOSS_LB2;
    if(lossType == "rB2") lt = LOSS_rB2;
    if(lossType == "rbLasso2") lt = LOSS_RBLASSO2;
    colvec dse(IS_RBLASSO(lt) ? w.size() : 1);
    if (IS_RBLASSO(lt)) dse = w + 0.0;
    w = bestRound(A, y, w, Bnds, scale, lossType);
    w = IntProgCalib(A, y, w, dse, lower, upper, Bnds, scale, lossType);
    return QInInPrCalib(A, y, w, dse, lower, upper, Bnds, scale, nsim, memrate, lossType);
}

RcppExport SEXP sparse_qiipc(SEXP ASEXP, SEXP ySEXP, SEXP wSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP BndsSEXP, SEXP scaleSEXP, SEXP nsimSEXP, SEXP memrateSEXP, SEXP lossTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const sp_mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const vec >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Bnds(BndsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< const double >::type memrate(memrateSEXP);
    Rcpp::traits::input_parameter< const std::string >::type lossType(lossTypeSEXP);
    __result = Rcpp::wrap(sp_qiipc(A, y, w, lower, upper, Bnds, scale, nsim, memrate, lossType));
    return __result;
END_RCPP
}

colvec de_qiipc(const mat& A, const colvec& y, colvec& w,
                const vec lower, const vec upper, const mat& Bnds,
                const vec& scale, const size_t nsim, const double memrate,
                const std::string lossType) {
    // Quantum Inspired Integer Programming Calibration for Sparse Matrices
    int lt = 0;
    if(lossType == "L1") lt = LOSS_L1;
    if(lossType == "aL1") lt = LOSS_aL1;
    if(lossType == "rL1") lt = LOSS_rL1;
    if(lossType == "LB1") lt = LOSS_LB1;
    if(lossType == "rB1") lt = LOSS_rB1;
    if(lossType == "rbLasso1") lt = LOSS_RBLASSO1;
    if(lossType == "L2") lt = LOSS_L2;
    if(lossType == "aL2") lt = LOSS_aL2;
    if(lossType == "rL2") lt = LOSS_rL2;
    if(lossType == "LB2") lt = LOSS_LB2;
    if(lossType == "rB2") lt = LOSS_rB2;
    if(lossType == "rbLasso2") lt = LOSS_RBLASSO2;
    colvec dse(IS_RBLASSO(lt) ? w.size() : 1);
    if (IS_RBLASSO(lt)) dse = w + 0.0;
    w = bestRound(A, y, w, Bnds, scale, lossType);
    w = IntProgCalib(A, y, w, dse, lower, upper, Bnds, scale, lossType);
    return QInInPrCalib(A, y, w, dse, lower, upper, Bnds, scale, nsim, memrate, lossType);
}

RcppExport SEXP dense_qiipc(SEXP ASEXP, SEXP ySEXP, SEXP wSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP BndsSEXP, SEXP scaleSEXP, SEXP nsimSEXP, SEXP memrateSEXP, SEXP lossTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< colvec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const vec >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const vec >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Bnds(BndsSEXP);
    Rcpp::traits::input_parameter< const vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< const double >::type memrate(memrateSEXP);
    Rcpp::traits::input_parameter< const std::string >::type lossType(lossTypeSEXP);
    __result = Rcpp::wrap(de_qiipc(A, y, w, lower, upper, Bnds, scale, nsim, memrate, lossType));
    return __result;
END_RCPP
}*/
