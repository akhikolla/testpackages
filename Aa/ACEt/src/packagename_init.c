#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void CWrapper_mcmc_atctet(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void CWrapper_mcmc_atctet_2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP gr_AtCtEt_epsp_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gr_AtCtEt_epsp_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gr_AtCtEt_esp_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP hessian_AtCtEt_esp_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP loglik_AtCtEt_epsp_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP loglik_AtCtEt_epsp_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP loglik_AtCtEt_esp_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"CWrapper_mcmc_atctet",   (DL_FUNC) &CWrapper_mcmc_atctet,   32},
    {"CWrapper_mcmc_atctet_2", (DL_FUNC) &CWrapper_mcmc_atctet_2, 23},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"gr_AtCtEt_epsp_c",       (DL_FUNC) &gr_AtCtEt_epsp_c,       17},
    {"gr_AtCtEt_epsp_g_c",     (DL_FUNC) &gr_AtCtEt_epsp_g_c,     17},
    {"gr_AtCtEt_esp_c",        (DL_FUNC) &gr_AtCtEt_esp_c,        11},
    {"hessian_AtCtEt_esp_c",   (DL_FUNC) &hessian_AtCtEt_esp_c,   11},
    {"loglik_AtCtEt_epsp_c",   (DL_FUNC) &loglik_AtCtEt_epsp_c,   17},
    {"loglik_AtCtEt_epsp_g_c", (DL_FUNC) &loglik_AtCtEt_epsp_g_c, 17},
    {"loglik_AtCtEt_esp_c",    (DL_FUNC) &loglik_AtCtEt_esp_c,    11},
    {NULL, NULL, 0}
};

void R_init_ACEt(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

