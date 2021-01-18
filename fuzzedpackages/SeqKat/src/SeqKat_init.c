#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SeqKat_cget_nucleotide_chunk_counts(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SeqKat_cget_trinucleotide_counts(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SeqKat_cpp_get_context(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SeqKat_cget_nucleotide_chunk_counts", (DL_FUNC) &_SeqKat_cget_nucleotide_chunk_counts, 6},
    {"_SeqKat_cget_trinucleotide_counts",    (DL_FUNC) &_SeqKat_cget_trinucleotide_counts,    5},
    {"_SeqKat_cpp_get_context",              (DL_FUNC) &_SeqKat_cpp_get_context,              3},
    {NULL, NULL, 0}
};

void R_init_SeqKat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
