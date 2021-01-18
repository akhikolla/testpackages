#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP Acpp(SEXP, SEXP);
extern SEXP Anlmcpp(SEXP, SEXP);
extern SEXP beetlecpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Dcpp(SEXP, SEXP);
extern SEXP distcpp(SEXP);
extern SEXP Dnlmcpp(SEXP, SEXP);
extern SEXP Ecpp(SEXP, SEXP);
extern SEXP Enlmcpp(SEXP, SEXP);
extern SEXP GPpredcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP HLRAcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP HLRDcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LMcpp(SEXP);
extern SEXP LRAcpp(SEXP, SEXP);
extern SEXP LRDcpp(SEXP, SEXP);
extern SEXP LRLAP2cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LRNSELLAP2cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nselhlrcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP nsellrcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP NSELnlmcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nselprcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PRLAP2cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PRNSELLAP2cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pvalcpp(SEXP, SEXP);
extern SEXP rowSumscpp(SEXP);
extern SEXP sighlrcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP siglrcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP SIGnlmcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sigprcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP utilcomp15badcpp(SEXP, SEXP);
extern SEXP utilcomp15sigcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP utilcomp18badcpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Acpp",             (DL_FUNC) &Acpp,             2},
    {"Anlmcpp",          (DL_FUNC) &Anlmcpp,          2},
    {"beetlecpp",        (DL_FUNC) &beetlecpp,        5},
    {"Dcpp",             (DL_FUNC) &Dcpp,             2},
    {"distcpp",          (DL_FUNC) &distcpp,          1},
    {"Dnlmcpp",          (DL_FUNC) &Dnlmcpp,          2},
    {"Ecpp",             (DL_FUNC) &Ecpp,             2},
    {"Enlmcpp",          (DL_FUNC) &Enlmcpp,          2},
    {"GPpredcpp",        (DL_FUNC) &GPpredcpp,        5},
    {"HLRAcpp",          (DL_FUNC) &HLRAcpp,          5},
    {"HLRDcpp",          (DL_FUNC) &HLRDcpp,          5},
    {"LMcpp",            (DL_FUNC) &LMcpp,            1},
    {"LRAcpp",           (DL_FUNC) &LRAcpp,           2},
    {"LRDcpp",           (DL_FUNC) &LRDcpp,           2},
    {"LRLAP2cpp",        (DL_FUNC) &LRLAP2cpp,        5},
    {"LRNSELLAP2cpp",    (DL_FUNC) &LRNSELLAP2cpp,    5},
    {"nselhlrcpp",       (DL_FUNC) &nselhlrcpp,       4},
    {"nsellrcpp",        (DL_FUNC) &nsellrcpp,        4},
    {"NSELnlmcpp",       (DL_FUNC) &NSELnlmcpp,       5},
    {"nselprcpp",        (DL_FUNC) &nselprcpp,        5},
    {"PRLAP2cpp",        (DL_FUNC) &PRLAP2cpp,        5},
    {"PRNSELLAP2cpp",    (DL_FUNC) &PRNSELLAP2cpp,    5},
    {"pvalcpp",          (DL_FUNC) &pvalcpp,          2},
    {"rowSumscpp",       (DL_FUNC) &rowSumscpp,       1},
    {"sighlrcpp",        (DL_FUNC) &sighlrcpp,        7},
    {"siglrcpp",         (DL_FUNC) &siglrcpp,         4},
    {"SIGnlmcpp",        (DL_FUNC) &SIGnlmcpp,        5},
    {"sigprcpp",         (DL_FUNC) &sigprcpp,         5},
    {"utilcomp15badcpp", (DL_FUNC) &utilcomp15badcpp, 2},
    {"utilcomp15sigcpp", (DL_FUNC) &utilcomp15sigcpp, 4},
    {"utilcomp18badcpp", (DL_FUNC) &utilcomp18badcpp, 2},
    {NULL, NULL, 0}
};

void R_init_acebayes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
