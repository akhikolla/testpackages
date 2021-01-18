#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP Rlabkey_listToMatrix(SEXP, SEXP);

R_CallMethodDef callMethods[] = {
	{ "Rlabkey_listToMatrix", (DL_FUNC) &Rlabkey_listToMatrix, 2},
	{NULL, NULL, 0}
};

void
R_init_Rlabkey(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
}
