#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Cggmfit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _gRim_C_ghk2pms(SEXP);
extern SEXP _gRim_C_pms2ghk(SEXP);
extern SEXP _gRim_fit2way_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRim_parm_ghk2pms_(SEXP);
extern SEXP _gRim_parm_normalize_ghk_(SEXP);
extern SEXP _gRim_parm_pms2ghk_(SEXP);
extern SEXP _gRim_parm_update_ghk_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRim_updateA(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_ghk2pms(SEXP);
extern SEXP C_pms2ghk(SEXP);

static const R_CMethodDef CEntries[] = {
    {"Cggmfit", (DL_FUNC) &Cggmfit, 14},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_gRim_C_ghk2pms",           (DL_FUNC) &_gRim_C_ghk2pms,           1},
    {"_gRim_C_pms2ghk",           (DL_FUNC) &_gRim_C_pms2ghk,           1},
    {"_gRim_fit2way_",            (DL_FUNC) &_gRim_fit2way_,            4},
    {"_gRim_parm_ghk2pms_",       (DL_FUNC) &_gRim_parm_ghk2pms_,       1},
    {"_gRim_parm_normalize_ghk_", (DL_FUNC) &_gRim_parm_normalize_ghk_, 1},
    {"_gRim_parm_pms2ghk_",       (DL_FUNC) &_gRim_parm_pms2ghk_,       1},
    {"_gRim_parm_update_ghk_",    (DL_FUNC) &_gRim_parm_update_ghk_,    9},
    {"_gRim_updateA",             (DL_FUNC) &_gRim_updateA,             4},
    {"C_ghk2pms",                 (DL_FUNC) &C_ghk2pms,                 1},
    {"C_pms2ghk",                 (DL_FUNC) &C_pms2ghk,                 1},
    {NULL, NULL, 0}
};

void R_init_gRim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
