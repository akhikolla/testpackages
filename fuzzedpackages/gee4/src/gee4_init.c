#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expresssions for .NAME have been omitted

    geefun
    gee_jmcm_method("new")
    ipwfun
    ipw_method("new")

  Most likely possible values need to be added below.
*/

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP gee4_geerfit_ar1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gee4_geerfit_cs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gee4_geerfit_id(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gee4_gees_estimation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gee4_ipw_estimation(SEXP, SEXP, SEXP, SEXP);
extern SEXP gee_jmcm__new(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gee_jmcm__get_m(SEXP, SEXP);
extern SEXP gee_jmcm__get_Y(SEXP, SEXP);
extern SEXP gee_jmcm__get_X(SEXP, SEXP);
extern SEXP gee_jmcm__get_Z(SEXP, SEXP);
extern SEXP gee_jmcm__get_W(SEXP, SEXP);
extern SEXP gee_jmcm__get_D(SEXP, SEXP, SEXP);
extern SEXP gee_jmcm__get_T(SEXP, SEXP, SEXP);
extern SEXP gee_jmcm__get_mu(SEXP, SEXP, SEXP);
extern SEXP gee_jmcm__get_Sigma(SEXP, SEXP, SEXP);
extern SEXP gee_jmcm__get_fim(SEXP, SEXP);
extern SEXP gee_jmcm__get_sd(SEXP, SEXP);
extern SEXP ipw__new(SEXP, SEXP, SEXP);
extern SEXP ipw__get_p(SEXP, SEXP);
extern SEXP ipw__get_Pi(SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"gee4_geerfit_ar1",     (DL_FUNC) &gee4_geerfit_ar1,     10},
    {"gee4_geerfit_cs",      (DL_FUNC) &gee4_geerfit_cs,      10},
    {"gee4_geerfit_id",      (DL_FUNC) &gee4_geerfit_id,      10},
    {"gee4_gees_estimation", (DL_FUNC) &gee4_gees_estimation, 13},
    {"gee4_ipw_estimation",  (DL_FUNC) &gee4_ipw_estimation,   4},
    {"gee_jmcm__new",        (DL_FUNC) &gee_jmcm__new,         7},
    {"gee_jmcm__get_m",      (DL_FUNC) &gee_jmcm__get_m,       2},
    {"gee_jmcm__get_Y",      (DL_FUNC) &gee_jmcm__get_Y,       2},
    {"gee_jmcm__get_X",      (DL_FUNC) &gee_jmcm__get_X,       2},
    {"gee_jmcm__get_Z",      (DL_FUNC) &gee_jmcm__get_Z,       2},
    {"gee_jmcm__get_W",      (DL_FUNC) &gee_jmcm__get_W,       2},
    {"gee_jmcm__get_D",      (DL_FUNC) &gee_jmcm__get_D,       3},
    {"gee_jmcm__get_T",      (DL_FUNC) &gee_jmcm__get_T,       3},
    {"gee_jmcm__get_mu",     (DL_FUNC) &gee_jmcm__get_mu,      3},
    {"gee_jmcm__get_Sigma",  (DL_FUNC) &gee_jmcm__get_Sigma,   3},
    {"gee_jmcm__get_fim",    (DL_FUNC) &gee_jmcm__get_fim,     2},
    {"gee_jmcm__get_sd",     (DL_FUNC) &gee_jmcm__get_sd,      2},
    {"ipw__new",             (DL_FUNC) &ipw__new,              3},
    {"ipw__get_p",           (DL_FUNC) &ipw__get_p,            2},
    {"ipw__get_Pi",          (DL_FUNC) &ipw__get_Pi,           2},
    {NULL, NULL, 0}
};

void R_init_gee4(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
