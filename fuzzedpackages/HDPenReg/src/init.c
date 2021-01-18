// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP cvEMfusedLasso1D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMfusedLasso2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlogisticFusedLasso1D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlogisticFusedLasso2D(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlogisticLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvlars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMfusedLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMlassoC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMlogisticFusedLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMlogisticLasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fusion(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP lars(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP EMlassoMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMlogisticLassoMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMfusedLassoMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP EMlogisticFusedLassoMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlassoMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMfusedLasso1DMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMfusedLasso2DMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlogisticLassoMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlogisticFusedLasso1DMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cvEMlogisticFusedLasso2DMain(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
  {"cvEMfusedLasso1D",         (DL_FUNC) &cvEMfusedLasso1D,         12},
  {"cvEMfusedLasso2D",         (DL_FUNC) &cvEMfusedLasso2D,         11},
  {"cvEMlasso",                (DL_FUNC) &cvEMlasso,                10},
  {"cvEMlogisticFusedLasso1D", (DL_FUNC) &cvEMlogisticFusedLasso1D, 12},
  {"cvEMlogisticFusedLasso2D", (DL_FUNC) &cvEMlogisticFusedLasso2D, 11},
  {"cvEMlogisticLasso",        (DL_FUNC) &cvEMlogisticLasso,        10},
  {"cvlars",                   (DL_FUNC) &cvlars,                   11},
  {"EMfusedLasso",             (DL_FUNC) &EMfusedLasso,             10},
  {"EMlassoC",                  (DL_FUNC) &EMlassoC,                   9},
  {"EMlogisticFusedLasso",     (DL_FUNC) &EMlogisticFusedLasso,     10},
  {"EMlogisticLasso",          (DL_FUNC) &EMlogisticLasso,           9},
  {"fusion",                   (DL_FUNC) &fusion,                    7},
  {"lars",                     (DL_FUNC) &lars,                      7},
  {"EMlassoMain",              (DL_FUNC) &EMlassoMain,               9},
  {"EMlogisticLassoMain",      (DL_FUNC) &EMlogisticLassoMain,       9},
  {"EMfusedLassoMain",         (DL_FUNC) &EMfusedLassoMain,         10},
  {"EMlogisticFusedLassoMain", (DL_FUNC) &EMlogisticFusedLassoMain, 10},
  {"cvEMlassoMain",            (DL_FUNC) &cvEMlassoMain,            10},
  {"cvEMfusedLasso1DMain",     (DL_FUNC) &cvEMfusedLasso1DMain,     12},
  {"cvEMfusedLasso2DMain",     (DL_FUNC) &cvEMfusedLasso2DMain,     11},
  {"cvEMlogisticLassoMain",    (DL_FUNC) &cvEMlogisticLassoMain,    10},
  {"cvEMlogisticFusedLasso1DMain", (DL_FUNC) &cvEMlogisticFusedLasso1DMain, 12},
  {"cvEMlogisticFusedLasso2DMain", (DL_FUNC) &cvEMlogisticFusedLasso2DMain, 11},
  {NULL, NULL, 0}
};

void R_init_HDPenReg(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
