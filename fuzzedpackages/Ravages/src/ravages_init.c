// with the help of package_native_routine_registration_skeleton(".", "src/ravages_init.c")

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP block_diag(SEXP);
extern SEXP gg_geno_stats_snps(SEXP, SEXP, SEXP);
extern SEXP label_multiple_genes(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP moments(SEXP, SEXP);
extern SEXP oz_burden2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP oz_jaccard(SEXP, SEXP, SEXP);
extern SEXP oz_new_bed_matrix(SEXP, SEXP);
extern SEXP oz_random_filling_bed_matrix_noHW(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP oz_sum_fst(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP oz_sum_fst_max_perm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP oz_sum_fst1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP oz_sum_fst1_max_perm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rbm_haplos_freqs(SEXP, SEXP, SEXP, SEXP);
extern SEXP rbm_haplos_thresholds_filling(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP skat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP skat_bootstrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"block_diag",                        (DL_FUNC) &block_diag,                         1},
    {"gg_geno_stats_snps",                (DL_FUNC) &gg_geno_stats_snps,                 3},
    {"label_multiple_genes",              (DL_FUNC) &label_multiple_genes,               5},
    {"moments",                           (DL_FUNC) &moments,                            2},
    {"oz_burden2",                        (DL_FUNC) &oz_burden2,                         6},
    {"oz_jaccard",                        (DL_FUNC) &oz_jaccard,                         3},
    {"oz_new_bed_matrix",                 (DL_FUNC) &oz_new_bed_matrix,                  2},
    {"oz_random_filling_bed_matrix_noHW", (DL_FUNC) &oz_random_filling_bed_matrix_noHW,  5},
    {"oz_sum_fst",                        (DL_FUNC) &oz_sum_fst,                         6},
    {"oz_sum_fst_max_perm",               (DL_FUNC) &oz_sum_fst_max_perm,                6},
    {"oz_sum_fst1",                       (DL_FUNC) &oz_sum_fst1,                        6},
    {"oz_sum_fst1_max_perm",              (DL_FUNC) &oz_sum_fst1_max_perm,               6},
    {"rbm_haplos_freqs",                  (DL_FUNC) &rbm_haplos_freqs,                   4},
    {"rbm_haplos_thresholds_filling",     (DL_FUNC) &rbm_haplos_thresholds_filling,      9},
    {"skat",                              (DL_FUNC) &skat,                               9},
    {"skat_bootstrap",                    (DL_FUNC) &skat_bootstrap,                    10},
    {NULL, NULL, 0}
};

void R_init_Ravages(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
