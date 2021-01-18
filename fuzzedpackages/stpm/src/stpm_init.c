#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Call calls */
  extern SEXP complik_gen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP complikGenNonGenetic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP complikMD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP complikMD_linear(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP complikMD_quadratic_linear(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_a0_abqff1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_a2_abqff1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_b0_bqff1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_b2_bqff1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_f0_ff1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_f1_f1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_f2_ff1mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_m0_mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_m2_mt_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_q0_qff1mtb_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_q2_qff1mtb_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP d_f_i1_t_t_g_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP devlik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP devlik_g_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mloglik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mloglik_g_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simCont(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simdata_ND(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simGenCont(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"complik_gen",                (DL_FUNC) &complik_gen,                23},
  {"complikGenNonGenetic",       (DL_FUNC) &complikGenNonGenetic,       21},
  {"complikMD",                  (DL_FUNC) &complikMD,                  14},
  {"complikMD_linear",           (DL_FUNC) &complikMD_linear,           13},
  {"complikMD_quadratic_linear", (DL_FUNC) &complikMD_quadratic_linear, 14},
  {"d_f_i1_a0_abqff1mt_g_c",     (DL_FUNC) &d_f_i1_a0_abqff1mt_g_c,     21},
  {"d_f_i1_a2_abqff1mt_g_c",     (DL_FUNC) &d_f_i1_a2_abqff1mt_g_c,     21},
  {"d_f_i1_b0_bqff1mt_g_c",      (DL_FUNC) &d_f_i1_b0_bqff1mt_g_c,      21},
  {"d_f_i1_b2_bqff1mt_g_c",      (DL_FUNC) &d_f_i1_b2_bqff1mt_g_c,      21},
  {"d_f_i1_f0_ff1mt_g_c",        (DL_FUNC) &d_f_i1_f0_ff1mt_g_c,        21},
  {"d_f_i1_f1_f1mt_g_c",         (DL_FUNC) &d_f_i1_f1_f1mt_g_c,         21},
  {"d_f_i1_f2_ff1mt_g_c",        (DL_FUNC) &d_f_i1_f2_ff1mt_g_c,        21},
  {"d_f_i1_m0_mt_g_c",           (DL_FUNC) &d_f_i1_m0_mt_g_c,           21},
  {"d_f_i1_m2_mt_g_c",           (DL_FUNC) &d_f_i1_m2_mt_g_c,           21},
  {"d_f_i1_q0_qff1mtb_g_c",      (DL_FUNC) &d_f_i1_q0_qff1mtb_g_c,      21},
  {"d_f_i1_q2_qff1mtb_g_c",      (DL_FUNC) &d_f_i1_q2_qff1mtb_g_c,      21},
  {"d_f_i1_t_t_g_c",             (DL_FUNC) &d_f_i1_t_t_g_c,             21},
  {"devlik",                     (DL_FUNC) &devlik,                      9},
  {"devlik_g_2",                 (DL_FUNC) &devlik_g_2,                 15},
  {"mloglik",                    (DL_FUNC) &mloglik,                     9},
  {"mloglik_g_2",                (DL_FUNC) &mloglik_g_2,                15},
  {"simCont",                    (DL_FUNC) &simCont,                    16},
  {"simdata_ND",                 (DL_FUNC) &simdata_ND,                 14},
  {"simGenCont",                 (DL_FUNC) &simGenCont,                 25},
  {NULL, NULL, 0}
};

void R_init_stpm(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
