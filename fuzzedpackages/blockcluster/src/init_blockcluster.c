// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "Rcoclustmain.h"

// declare RcppExport SEXP CoClustmain(SEXP robj)
//
static const R_CallMethodDef callMethods[]  =
{
  {"CoClustmain", (DL_FUNC) &CoClustmain, 2},
  {NULL}
};

void R_init_myRoutines(DllInfo *info)
{
	/* Register the .Call routines.
	No .C  .Fortran() or .External() routines,
	so pass those arrays as NULL.
	*/
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}
