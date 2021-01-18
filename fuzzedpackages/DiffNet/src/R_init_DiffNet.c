#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP DiffNet_GHD_Fast(SEXP ASEXP, SEXP BSEXP);

SEXP DiffNet_MU_Fast(SEXP ASEXP, SEXP BSEXP);

SEXP DiffNet_STD_Fast(SEXP ASEXP, SEXP BSEXP);


static R_CallMethodDef callMethods[] = {
  {"DiffNet_GHD_Fast", (DL_FUNC) &DiffNet_GHD_Fast, 2},
  {"DiffNet_MU_Fast", (DL_FUNC) &DiffNet_MU_Fast, 2},
  {"DiffNet_STD_Fast",(DL_FUNC) &DiffNet_STD_Fast, 2},
  {NULL, NULL, 0}
};

void R_init_DiffNet(DllInfo *info)
{
  /* Register routines,
   allocate resources. */
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_DiffNet(DllInfo *info)
{
  /* Release resources. */
}
