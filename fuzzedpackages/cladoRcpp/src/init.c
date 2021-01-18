#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*

NOTE (2018-09-26, NJM)

This init.c file was built by package_native_routine_registration_skeleton().

It was put in, copying this rexpokit update:

It was a change requested by Kurt in rexpokit updates:

Registering the C++ and FORTRAN calls:
https://stat.ethz.ch/pipermail/r-devel/2017-February/073755.html

library(tools)
package_native_routine_registration_skeleton(dir="/GitHub/rexpokit")

*/



/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP convolve3cpp(SEXP a, SEXP b);

extern SEXP cpp_areas_list_to_states_list(SEXP R_areas_indices, SEXP R_maxareas, SEXP R_include_null_range);

extern SEXP cpp_calc_anclikes_sp(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP states_indices, SEXP s,
SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP maxent01y, SEXP
max_minsize_as_function_of_ancsize, SEXP Rsp_rowsums);

extern SEXP cpp_calc_anclikes_sp_COOprobs(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP states_indices,
SEXP s, SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP maxent01y, SEXP
max_minsize_as_function_of_ancsize);

extern SEXP cpp_calc_anclikes_sp_COOweights_faster(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP
states_indices, SEXP s, SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP
maxent01y, SEXP max_minsize_as_function_of_ancsize);

extern SEXP cpp_calc_anclikes_sp_rowsums(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP states_indices,
SEXP s, SEXP v, SEXP j, SEXP y, SEXP dmat, SEXP maxent01s, SEXP maxent01v, SEXP maxent01j, SEXP maxent01y, SEXP
max_minsize_as_function_of_ancsize);

extern SEXP cpp_calc_anclikes_sp_using_COOprobs(SEXP Rprintmat, SEXP leftprobs, SEXP rightprobs, SEXP
RCOO_left_i_list, SEXP RCOO_right_j_list, SEXP RCOO_probs_list, SEXP Rsp_rowsums);

extern SEXP cpp_calc_rowsums_for_COOweights_columnar(SEXP RCOO_weights_columnar_anc_i_list, SEXP RCOO_probs_list, SEXP Rnumstates);

extern SEXP cpp_calc_splitlikes_using_COOweights_columnar(SEXP leftprobs, SEXP rightprobs, SEXP RCOO_weights_columnar_anc_i_list, SEXP RCOO_left_i_list, SEXP RCOO_right_j_list, SEXP RCOO_probs_list, SEXP Rsp_rowsums);

extern SEXP cpp_combn_zerostart(SEXP R_n, SEXP R_m, SEXP R_maxval);

extern SEXP cpp_states_list_to_DEmat(SEXP R_areas_indices, SEXP R_states_indices, SEXP R_dmat, SEXP R_elist, 
SEXP R_amat, SEXP R_normalize_TF);

extern SEXP cpp_states_list_to_DEmat_COO(SEXP R_areas_indices, SEXP R_states_indices, SEXP R_dmat, SEXP R_elist, SEXP R_amat, SEXP R_normalize_TF, SEXP R_min_precision);

extern SEXP mult2probvect(SEXP leftprobs, SEXP rightprobs);

static const R_CallMethodDef CallEntries[] = {
    {"convolve3cpp",                                  (DL_FUNC) &convolve3cpp,                                   2},
    {"cpp_areas_list_to_states_list",                 (DL_FUNC) &cpp_areas_list_to_states_list,                  3},
    {"cpp_calc_anclikes_sp",                          (DL_FUNC) &cpp_calc_anclikes_sp,                          15},
    {"cpp_calc_anclikes_sp_COOprobs",                 (DL_FUNC) &cpp_calc_anclikes_sp_COOprobs,                 14},
    {"cpp_calc_anclikes_sp_COOweights_faster",        (DL_FUNC) &cpp_calc_anclikes_sp_COOweights_faster,        14},
    {"cpp_calc_anclikes_sp_rowsums",                  (DL_FUNC) &cpp_calc_anclikes_sp_rowsums,                  14},
    {"cpp_calc_anclikes_sp_using_COOprobs",           (DL_FUNC) &cpp_calc_anclikes_sp_using_COOprobs,            7},
    {"cpp_calc_rowsums_for_COOweights_columnar",      (DL_FUNC) &cpp_calc_rowsums_for_COOweights_columnar,       3},
    {"cpp_calc_splitlikes_using_COOweights_columnar", (DL_FUNC) &cpp_calc_splitlikes_using_COOweights_columnar,  7},
    {"cpp_combn_zerostart",                           (DL_FUNC) &cpp_combn_zerostart,                            3},
    {"cpp_states_list_to_DEmat",                      (DL_FUNC) &cpp_states_list_to_DEmat,                       6},
    {"cpp_states_list_to_DEmat_COO",                  (DL_FUNC) &cpp_states_list_to_DEmat_COO,                   7},
    {"mult2probvect",                                 (DL_FUNC) &mult2probvect,                                  2},
    {NULL, NULL, 0}
};

void R_init_cladoRcpp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
