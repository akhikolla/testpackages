#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP JSM_calc_bi_st(SEXP, SEXP, SEXP);
extern SEXP JSM_calc_expM2(SEXP);
extern SEXP JSM_calc_M_v(SEXP, SEXP);
extern SEXP JSM_calc_M1_a_M2_Hadamard(SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_M1_M2_Hadamard(SEXP, SEXP);
extern SEXP JSM_calc_M1_M2_Hadamard_a(SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_M1_M2_M3_Hadamard(SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_M1timesM2v(SEXP, SEXP, SEXP);
extern SEXP JSM_calc_muB(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_muBMult(SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_mult_rowsum1(SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_mult_rowsum2(SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_mult_rowsum3(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP JSM_calc_MVND(SEXP, SEXP, SEXP);
extern SEXP JSM_calc_rowsum(SEXP, SEXP);
extern SEXP JSM_calc_rowsum_mult(SEXP, SEXP, SEXP);
extern SEXP JSM_calc_tapply_vect_sum(SEXP, SEXP);
extern SEXP JSM_calc_v_a(SEXP, SEXP);
extern SEXP JSM_calc_VB(SEXP, SEXP, SEXP);
extern SEXP JSM_calc_VY(SEXP, SEXP, SEXP);
extern SEXP JSM_fast_lapply_length(SEXP, SEXP, SEXP);
extern SEXP JSM_fast_rbind_lapply_outerprod(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"JSM_calc_bi_st",                  (DL_FUNC) &JSM_calc_bi_st,                  3},
    {"JSM_calc_expM2",                  (DL_FUNC) &JSM_calc_expM2,                  1},
    {"JSM_calc_M_v",                    (DL_FUNC) &JSM_calc_M_v,                    2},
    {"JSM_calc_M1_a_M2_Hadamard",       (DL_FUNC) &JSM_calc_M1_a_M2_Hadamard,       4},
    {"JSM_calc_M1_M2_Hadamard",         (DL_FUNC) &JSM_calc_M1_M2_Hadamard,         2},
    {"JSM_calc_M1_M2_Hadamard_a",       (DL_FUNC) &JSM_calc_M1_M2_Hadamard_a,       4},
    {"JSM_calc_M1_M2_M3_Hadamard",      (DL_FUNC) &JSM_calc_M1_M2_M3_Hadamard,      4},
    {"JSM_calc_M1timesM2v",             (DL_FUNC) &JSM_calc_M1timesM2v,             3},
    {"JSM_calc_muB",                    (DL_FUNC) &JSM_calc_muB,                    6},
    {"JSM_calc_muBMult",                (DL_FUNC) &JSM_calc_muBMult,                4},
    {"JSM_calc_mult_rowsum1",           (DL_FUNC) &JSM_calc_mult_rowsum1,           4},
    {"JSM_calc_mult_rowsum2",           (DL_FUNC) &JSM_calc_mult_rowsum2,           4},
    {"JSM_calc_mult_rowsum3",           (DL_FUNC) &JSM_calc_mult_rowsum3,           5},
    {"JSM_calc_MVND",                   (DL_FUNC) &JSM_calc_MVND,                   3},
    {"JSM_calc_rowsum",                 (DL_FUNC) &JSM_calc_rowsum,                 2},
    {"JSM_calc_rowsum_mult",            (DL_FUNC) &JSM_calc_rowsum_mult,            3},
    {"JSM_calc_tapply_vect_sum",        (DL_FUNC) &JSM_calc_tapply_vect_sum,        2},
    {"JSM_calc_v_a",                    (DL_FUNC) &JSM_calc_v_a,                    2},
    {"JSM_calc_VB",                     (DL_FUNC) &JSM_calc_VB,                     3},
    {"JSM_calc_VY",                     (DL_FUNC) &JSM_calc_VY,                     3},
    {"JSM_fast_lapply_length",          (DL_FUNC) &JSM_fast_lapply_length,          3},
    {"JSM_fast_rbind_lapply_outerprod", (DL_FUNC) &JSM_fast_rbind_lapply_outerprod, 1},
    {NULL, NULL, 0}
};

void R_init_JSM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
