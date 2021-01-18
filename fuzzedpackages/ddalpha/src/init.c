#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void AlphaLearn(double *, int *, int *, int *, int *, int *, double *);
extern void AlphaClassify(void *, void *, void *, void *, void *, void *);
extern void AlphaLearnCV(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void DKnnClassify(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void DKnnLearnCv(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HDepth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void HDepthEx(void *, void *, void *, void *, void *, void *, void *);
extern void HDepthSpaceEx(void *, void *, void *, void *, void *, void *, void *, void *);
extern void HDSpace(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void IsInConvexes(void *, void *, void *, void *, void *, void *, void *, void *);
extern void KnnAffInvClassify(void *, void *, void *, void *, void *, void *, void *);
extern void KnnAffInvLearnJK(void *, void *, void *, void *, void *);
extern void KnnClassify(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void KnnLearnJK(void *, void *, void *, void *, void *, void *, void *);
extern void MahalanobisDepth(void *, void *, void *, void *, void *, void *, void *);
extern void OjaDepth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PolynomialLearnCV(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void PotentialDepthsCount(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ProjectionDepth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void SimplicialDepth(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ZDepth(void *, void *, void *, void *, void *, void *, void *);
extern void BetaSkeletonDepth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void SimplicialBandDepthF(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(adjc)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(adjlp)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(bd)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cvkernsm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(diffd)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dpth1)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dpth2)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fund1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(fund2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(funmd)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(funrpd1)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(funrpd2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(hrd)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(metrc)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(metrl2)(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"AlphaLearn",           (DL_FUNC) &AlphaLearn,            7},
    {"AlphaClassify",        (DL_FUNC) &AlphaClassify,         6},
    {"AlphaLearnCV",         (DL_FUNC) &AlphaLearnCV,          9},
    {"DKnnClassify",         (DL_FUNC) &DKnnClassify,         10},
    {"DKnnLearnCv",          (DL_FUNC) &DKnnLearnCv,           9},
    {"HDepth",               (DL_FUNC) &HDepth,               12},
    {"HDepthEx",             (DL_FUNC) &HDepthEx,              7},
    {"HDepthSpaceEx",        (DL_FUNC) &HDepthSpaceEx,         8},
    {"HDSpace",              (DL_FUNC) &HDSpace,              10},
    {"IsInConvexes",         (DL_FUNC) &IsInConvexes,          8},
    {"KnnAffInvClassify",    (DL_FUNC) &KnnAffInvClassify,     7},
    {"KnnAffInvLearnJK",     (DL_FUNC) &KnnAffInvLearnJK,      5},
    {"KnnClassify",          (DL_FUNC) &KnnClassify,           9},
    {"KnnLearnJK",           (DL_FUNC) &KnnLearnJK,            7},
    {"MahalanobisDepth",     (DL_FUNC) &MahalanobisDepth,      7},
    {"OjaDepth",             (DL_FUNC) &OjaDepth,             11},
    {"PolynomialLearnCV",    (DL_FUNC) &PolynomialLearnCV,    10},
    {"PotentialDepthsCount", (DL_FUNC) &PotentialDepthsCount, 11},
    {"ProjectionDepth",      (DL_FUNC) &ProjectionDepth,      12},
    {"SimplicialDepth",      (DL_FUNC) &SimplicialDepth,       9},
    {"ZDepth",               (DL_FUNC) &ZDepth,                7},
    {"BetaSkeletonDepth",    (DL_FUNC) &BetaSkeletonDepth,    10},
    {"SimplicialBandDepthF", (DL_FUNC) &SimplicialBandDepthF, 10},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"adjc",     (DL_FUNC) &F77_NAME(adjc),      8},
    {"adjlp",    (DL_FUNC) &F77_NAME(adjlp),     8},
    {"bd",       (DL_FUNC) &F77_NAME(bd),        6},
    {"cvkernsm", (DL_FUNC) &F77_NAME(cvkernsm), 15},
    {"diffd",    (DL_FUNC) &F77_NAME(diffd),    15},
    {"dpth1",    (DL_FUNC) &F77_NAME(dpth1),     6},
    {"dpth2",    (DL_FUNC) &F77_NAME(dpth2),     8},
    {"fund1",    (DL_FUNC) &F77_NAME(fund1),    11},
    {"fund2",    (DL_FUNC) &F77_NAME(fund2),    13},
    {"funmd",    (DL_FUNC) &F77_NAME(funmd),     7},
    {"funrpd1",  (DL_FUNC) &F77_NAME(funrpd1),  11},
    {"funrpd2",  (DL_FUNC) &F77_NAME(funrpd2),  15},
    {"hrd",      (DL_FUNC) &F77_NAME(hrd),       6},
    {"metrc",    (DL_FUNC) &F77_NAME(metrc),     6},
    {"metrl2",   (DL_FUNC) &F77_NAME(metrl2),    6},
    {NULL, NULL, 0}
};

void R_init_ddalpha(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
