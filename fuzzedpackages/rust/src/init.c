#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rust_all_pos(SEXP);
extern SEXP _rust_any_col_nonneg(SEXP);
extern SEXP _rust_any_col_nonpos(SEXP);
extern SEXP _rust_any_naC(SEXP);
extern SEXP _rust_any_neg(SEXP);
extern SEXP _rust_any_nonpos(SEXP);
extern SEXP _rust_any_pos(SEXP);
extern SEXP _rust_bc_log_j(SEXP, SEXP);
extern SEXP _rust_bc_no_trans(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_bc_phi_to_theta(SEXP, SEXP);
extern SEXP _rust_cpp_a_obj(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_a_obj_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_logf(SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_logf_rho(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_logf_rho_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_logf_rho_3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_logf_rho_4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_logf_scaled(SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_lower_box(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_lower_box_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_psi_to_phi(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_psi_to_phi_0(SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_rho_to_psi(SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_upper_box(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_cpp_upper_box_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_create_log_j_xptr(SEXP);
extern SEXP _rust_create_log_jac_xptr(SEXP);
extern SEXP _rust_create_phi_to_theta_xptr(SEXP);
extern SEXP _rust_create_psi_to_phi_xptr(SEXP);
extern SEXP _rust_create_trans_xptr(SEXP);
extern SEXP _rust_create_xptr(SEXP);
extern SEXP _rust_exptrans(SEXP, SEXP);
extern SEXP _rust_gp_phi_to_theta(SEXP, SEXP);
extern SEXP _rust_log_none_jac(SEXP, SEXP);
extern SEXP _rust_logcauchy(SEXP, SEXP);
extern SEXP _rust_logdgamma(SEXP, SEXP);
extern SEXP _rust_logdlnorm(SEXP, SEXP);
extern SEXP _rust_logdmvnorm(SEXP, SEXP);
extern SEXP _rust_logdN01(SEXP, SEXP);
extern SEXP _rust_logdnorm2(SEXP, SEXP);
extern SEXP _rust_loggp(SEXP, SEXP);
extern SEXP _rust_loghalfcauchy(SEXP, SEXP);
extern SEXP _rust_lognormalmix(SEXP, SEXP);
extern SEXP _rust_lognormt(SEXP, SEXP);
extern SEXP _rust_neglog(SEXP, SEXP);
extern SEXP _rust_no_naC(SEXP);
extern SEXP _rust_no_trans(SEXP, SEXP);
extern SEXP _rust_null_phi_to_theta_xptr(SEXP);
extern SEXP _rust_rcpp_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_RcppExport_registerCCallable();
extern SEXP _rust_ru_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_ru_cpp_2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_ru_cpp_3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_ru_cpp_4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rust_vecpow(SEXP, SEXP);
extern SEXP _rust_vecpower(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rust_all_pos",                      (DL_FUNC) &_rust_all_pos,                       1},
    {"_rust_any_col_nonneg",               (DL_FUNC) &_rust_any_col_nonneg,                1},
    {"_rust_any_col_nonpos",               (DL_FUNC) &_rust_any_col_nonpos,                1},
    {"_rust_any_naC",                      (DL_FUNC) &_rust_any_naC,                       1},
    {"_rust_any_neg",                      (DL_FUNC) &_rust_any_neg,                       1},
    {"_rust_any_nonpos",                   (DL_FUNC) &_rust_any_nonpos,                    1},
    {"_rust_any_pos",                      (DL_FUNC) &_rust_any_pos,                       1},
    {"_rust_bc_log_j",                     (DL_FUNC) &_rust_bc_log_j,                      2},
    {"_rust_bc_no_trans",                  (DL_FUNC) &_rust_bc_no_trans,                   4},
    {"_rust_bc_phi_to_theta",              (DL_FUNC) &_rust_bc_phi_to_theta,               2},
    {"_rust_cpp_a_obj",                    (DL_FUNC) &_rust_cpp_a_obj,                     9},
    {"_rust_cpp_a_obj_2",                  (DL_FUNC) &_rust_cpp_a_obj_2,                  15},
    {"_rust_cpp_logf",                     (DL_FUNC) &_rust_cpp_logf,                      3},
    {"_rust_cpp_logf_rho",                 (DL_FUNC) &_rust_cpp_logf_rho,                  6},
    {"_rust_cpp_logf_rho_2",               (DL_FUNC) &_rust_cpp_logf_rho_2,               11},
    {"_rust_cpp_logf_rho_3",               (DL_FUNC) &_rust_cpp_logf_rho_3,               11},
    {"_rust_cpp_logf_rho_4",               (DL_FUNC) &_rust_cpp_logf_rho_4,               11},
    {"_rust_cpp_logf_scaled",              (DL_FUNC) &_rust_cpp_logf_scaled,               3},
    {"_rust_cpp_lower_box",                (DL_FUNC) &_rust_cpp_lower_box,                10},
    {"_rust_cpp_lower_box_2",              (DL_FUNC) &_rust_cpp_lower_box_2,              16},
    {"_rust_cpp_psi_to_phi",               (DL_FUNC) &_rust_cpp_psi_to_phi,                4},
    {"_rust_cpp_psi_to_phi_0",             (DL_FUNC) &_rust_cpp_psi_to_phi_0,              4},
    {"_rust_cpp_rho_to_psi",               (DL_FUNC) &_rust_cpp_rho_to_psi,                3},
    {"_rust_cpp_upper_box",                (DL_FUNC) &_rust_cpp_upper_box,                10},
    {"_rust_cpp_upper_box_2",              (DL_FUNC) &_rust_cpp_upper_box_2,              16},
    {"_rust_create_log_j_xptr",            (DL_FUNC) &_rust_create_log_j_xptr,             1},
    {"_rust_create_log_jac_xptr",          (DL_FUNC) &_rust_create_log_jac_xptr,           1},
    {"_rust_create_phi_to_theta_xptr",     (DL_FUNC) &_rust_create_phi_to_theta_xptr,      1},
    {"_rust_create_psi_to_phi_xptr",       (DL_FUNC) &_rust_create_psi_to_phi_xptr,        1},
    {"_rust_create_trans_xptr",            (DL_FUNC) &_rust_create_trans_xptr,             1},
    {"_rust_create_xptr",                  (DL_FUNC) &_rust_create_xptr,                   1},
    {"_rust_exptrans",                     (DL_FUNC) &_rust_exptrans,                      2},
    {"_rust_gp_phi_to_theta",              (DL_FUNC) &_rust_gp_phi_to_theta,               2},
    {"_rust_log_none_jac",                 (DL_FUNC) &_rust_log_none_jac,                  2},
    {"_rust_logcauchy",                    (DL_FUNC) &_rust_logcauchy,                     2},
    {"_rust_logdgamma",                    (DL_FUNC) &_rust_logdgamma,                     2},
    {"_rust_logdlnorm",                    (DL_FUNC) &_rust_logdlnorm,                     2},
    {"_rust_logdmvnorm",                   (DL_FUNC) &_rust_logdmvnorm,                    2},
    {"_rust_logdN01",                      (DL_FUNC) &_rust_logdN01,                       2},
    {"_rust_logdnorm2",                    (DL_FUNC) &_rust_logdnorm2,                     2},
    {"_rust_loggp",                        (DL_FUNC) &_rust_loggp,                         2},
    {"_rust_loghalfcauchy",                (DL_FUNC) &_rust_loghalfcauchy,                 2},
    {"_rust_lognormalmix",                 (DL_FUNC) &_rust_lognormalmix,                  2},
    {"_rust_lognormt",                     (DL_FUNC) &_rust_lognormt,                      2},
    {"_rust_neglog",                       (DL_FUNC) &_rust_neglog,                        2},
    {"_rust_no_naC",                       (DL_FUNC) &_rust_no_naC,                        1},
    {"_rust_no_trans",                     (DL_FUNC) &_rust_no_trans,                      2},
    {"_rust_null_phi_to_theta_xptr",       (DL_FUNC) &_rust_null_phi_to_theta_xptr,        1},
    {"_rust_rcpp_apply",                   (DL_FUNC) &_rust_rcpp_apply,                   11},
    {"_rust_RcppExport_registerCCallable", (DL_FUNC) &_rust_RcppExport_registerCCallable,  0},
    {"_rust_ru_cpp",                       (DL_FUNC) &_rust_ru_cpp,                       11},
    {"_rust_ru_cpp_2",                     (DL_FUNC) &_rust_ru_cpp_2,                     16},
    {"_rust_ru_cpp_3",                     (DL_FUNC) &_rust_ru_cpp_3,                     16},
    {"_rust_ru_cpp_4",                     (DL_FUNC) &_rust_ru_cpp_4,                     16},
    {"_rust_vecpow",                       (DL_FUNC) &_rust_vecpow,                        2},
    {"_rust_vecpower",                     (DL_FUNC) &_rust_vecpower,                      2},
    {NULL, NULL, 0}
};

void R_init_rust(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
