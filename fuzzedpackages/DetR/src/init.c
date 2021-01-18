#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void R_extCstep(void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_FastOGK(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void R_inFM(void *, void *, void *);
extern void R_inQn(void *, void *, void *);
extern void R_unimcd(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"R_extCstep", (DL_FUNC) &R_extCstep,  8},
    {"R_FastOGK",  (DL_FUNC) &R_FastOGK,  14},
    {"R_inFM",     (DL_FUNC) &R_inFM,      3},
    {"R_inQn",     (DL_FUNC) &R_inQn,      3},
    {"R_unimcd",   (DL_FUNC) &R_unimcd,    7},
    {NULL, NULL, 0}
};

void R_init_DetR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
