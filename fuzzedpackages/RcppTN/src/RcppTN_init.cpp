// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include "etn1.h"
#include "rtn1.h"
#include "vtn1.h"
#include "dtn1.h"
#include "enttn1.h"

RcppExport void R_init_RcppTN (DllInfo *info) {
    R_RegisterCCallable("RcppTN", "RcppTN_etn1", (DL_FUNC) etn1) ;
    R_RegisterCCallable("RcppTN", "RcppTN_rtn1", (DL_FUNC) rtn1) ;
    R_RegisterCCallable("RcppTN", "RcppTN_vtn1", (DL_FUNC) vtn1) ;
    R_RegisterCCallable("RcppTN", "RcppTN_dtn1", (DL_FUNC) dtn1) ;
    R_RegisterCCallable("RcppTN", "RcppTN_enttn1", (DL_FUNC) enttn1) ;
    R_registerRoutines(info, NULL, NULL, NULL, NULL) ;
    R_useDynamicSymbols(info, TRUE) ;

}
