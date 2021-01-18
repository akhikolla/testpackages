#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void R_FastOGK(void *, void *, void *, void *, void *);
extern void R_FastR(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_inQn(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"R_FastOGK",  (DL_FUNC) &R_FastOGK,   5},
    {"R_FastR",    (DL_FUNC) &R_FastR,    16},
    {"R_inQn",     (DL_FUNC) &R_inQn,      3},
    {NULL, NULL, 0}
};

void R_init_DetMCD(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
