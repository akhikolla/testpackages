#ifndef EMGAUSSH
#define EMGAUSSH

RcppExport
SEXP CEMGauss(SEXP X_s, SEXP k_s, SEXP Config_s, SEXP Homoc_s, 
    SEXP maxiter_s, SEXP tautol_s, SEXP convtol_s, SEXP k2max_s, SEXP MaxSctEgvlRt_s,
//  SEXP InitSolz_s, SEXP InitSoltau_s, SEXP InitSolmuk_s, SEXP InitSolSigma_s, SEXP InitSolSigmak_s, SEXP InitSolLnLik_s, SEXP startwithM_s);       
  SEXP InitSolz_s, SEXP InitSoltau_s, SEXP InitSolmuk_s, SEXP InitSolSigma_s, SEXP InitSolSigmak_s,
  SEXP InitSolLnLik_s, SEXP startwithM_s, SEXP SctEgvCnstr_s, SEXP MaxVarGRt_s);       

#endif
