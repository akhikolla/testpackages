#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(bivrecur)(int *n, double *gtime, double *ctime, double *count,
                     int *mc, int *m,  double *cen,  double *ucen, int *nd,
                     double *udt,  int *tot,  double *gap,  double *event,
                     double *r,  double *d,  double *sest,
                     double *var,  double *markvar1,  double *markvar2,  double *mark1,
                     double *mark2, double *u1,  double *u2,  double *Fest,
                     int *tmpindex,  double *prob,  double *std);
extern void F77_NAME(mprovar)(int *n,  int *nparams,
                     double *xmati, double *ymati,
                     double *gmatx,  double *gmaty,
                     double *l1, double *l2,
                     double *expAx, double *expAy, double *subsumx,    double *subsumy,
                     double *dx, double *dy, double *mstar, int *mc);
extern void F77_NAME(onesamp)(int *n, double *gtime, double *ctime, double *count,
                     int *mc, int *m,  double *cen,  double *ucen, int *nd,
                     double *udt,  int *tot,  double *gap,  double *event,
                     double *r,  double *d,  double *sest, double *std);
extern void F77_NAME(xmproee)(int *n,  int *nparams,
                     double *di, double *xmati, double *gmati,
                     double *L, double *expA, double *subsum,
                     int *kcount);
extern void F77_NAME(ymproee)(int *n,  int *nparams,
                     double *di, double *xmati, double *ymati, double *gmati,
                     double *L, double *expA, double *subsum,
                     int *kcount);

static const R_FortranMethodDef FortranEntries[] = {
    {"bivrecur",  (DL_FUNC) &F77_NAME(bivrecur),  27},
    {"mprovar",   (DL_FUNC) &F77_NAME(mprovar),   16},
    {"onesamp",   (DL_FUNC) &F77_NAME(onesamp),   17},
    {"xmproee",   (DL_FUNC) &F77_NAME(xmproee),    9},
    {"ymproee",   (DL_FUNC) &F77_NAME(ymproee),   10},
    {NULL, NULL, 0}
};

void R_init_BivRec(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
