#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _revdbayes_any_nonpos(SEXP);
extern SEXP _revdbayes_cpp_gev_beta(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_flat(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_flatflat(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_loglik(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_loglognorm(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_mdi(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_norm(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_prob(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gev_quant(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_beta(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_flat(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_flatflat(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_jeffreys(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_loglik(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_mdi(SEXP, SEXP);
extern SEXP _revdbayes_cpp_gp_norm(SEXP, SEXP);
extern SEXP _revdbayes_cpp_os_loglik(SEXP, SEXP);
extern SEXP _revdbayes_cpp_pp_loglik(SEXP, SEXP);
extern SEXP _revdbayes_create_prior_xptr(SEXP);
extern SEXP _revdbayes_gev_beta_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_beta_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_flat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_flat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_flatflat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_flatflat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_loglognorm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_loglognorm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_logpost_phi_xptr(SEXP);
extern SEXP _revdbayes_gev_logpost_xptr(SEXP);
extern SEXP _revdbayes_gev_mdi_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_mdi_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_norm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_norm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_phi_to_theta(SEXP, SEXP);
extern SEXP _revdbayes_gev_prob_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_prob_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_quant_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_quant_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gev_user_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gev_user_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_beta_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_beta_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_flat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_flat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_flatflat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_flatflat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_jeffreys_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_jeffreys_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_logpost_phi_xptr(SEXP);
extern SEXP _revdbayes_gp_logpost_xptr(SEXP);
extern SEXP _revdbayes_gp_mdi_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_mdi_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_norm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_norm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_gp_phi_to_theta(SEXP, SEXP);
extern SEXP _revdbayes_gp_user_logpost(SEXP, SEXP);
extern SEXP _revdbayes_gp_user_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_kgaps_log_j(SEXP, SEXP);
extern SEXP _revdbayes_kgaps_logpost(SEXP, SEXP);
extern SEXP _revdbayes_kgaps_logpost_xptr(SEXP);
extern SEXP _revdbayes_kgaps_phi_to_theta(SEXP, SEXP);
extern SEXP _revdbayes_lgdgev_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _revdbayes_log_j_xptr(SEXP);
extern SEXP _revdbayes_os_beta_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_beta_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_flat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_flat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_flatflat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_flatflat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_loglognorm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_loglognorm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_logpost_phi_xptr(SEXP);
extern SEXP _revdbayes_os_logpost_xptr(SEXP);
extern SEXP _revdbayes_os_mdi_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_mdi_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_norm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_norm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_prob_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_prob_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_quant_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_quant_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_os_user_logpost(SEXP, SEXP);
extern SEXP _revdbayes_os_user_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pgev_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _revdbayes_phi_to_theta_xptr(SEXP);
extern SEXP _revdbayes_pp_beta_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_beta_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_flat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_flat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_flatflat_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_flatflat_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_loglognorm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_loglognorm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_logpost_phi_xptr(SEXP);
extern SEXP _revdbayes_pp_logpost_xptr(SEXP);
extern SEXP _revdbayes_pp_mdi_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_mdi_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_norm_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_norm_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_phi_to_theta(SEXP, SEXP);
extern SEXP _revdbayes_pp_prob_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_prob_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_quant_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_quant_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_pp_user_logpost(SEXP, SEXP);
extern SEXP _revdbayes_pp_user_logpost_phi(SEXP, SEXP);
extern SEXP _revdbayes_qgev_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _revdbayes_RcppExport_registerCCallable();
extern SEXP _revdbayes_user_gev_flat(SEXP, SEXP);
extern SEXP _revdbayes_user_gev_norm(SEXP, SEXP);
extern SEXP _revdbayes_user_gp_flat(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_revdbayes_any_nonpos",                   (DL_FUNC) &_revdbayes_any_nonpos,                   1},
    {"_revdbayes_cpp_gev_beta",                 (DL_FUNC) &_revdbayes_cpp_gev_beta,                 2},
    {"_revdbayes_cpp_gev_flat",                 (DL_FUNC) &_revdbayes_cpp_gev_flat,                 2},
    {"_revdbayes_cpp_gev_flatflat",             (DL_FUNC) &_revdbayes_cpp_gev_flatflat,             2},
    {"_revdbayes_cpp_gev_loglik",               (DL_FUNC) &_revdbayes_cpp_gev_loglik,               2},
    {"_revdbayes_cpp_gev_loglognorm",           (DL_FUNC) &_revdbayes_cpp_gev_loglognorm,           2},
    {"_revdbayes_cpp_gev_mdi",                  (DL_FUNC) &_revdbayes_cpp_gev_mdi,                  2},
    {"_revdbayes_cpp_gev_norm",                 (DL_FUNC) &_revdbayes_cpp_gev_norm,                 2},
    {"_revdbayes_cpp_gev_prob",                 (DL_FUNC) &_revdbayes_cpp_gev_prob,                 2},
    {"_revdbayes_cpp_gev_quant",                (DL_FUNC) &_revdbayes_cpp_gev_quant,                2},
    {"_revdbayes_cpp_gp_beta",                  (DL_FUNC) &_revdbayes_cpp_gp_beta,                  2},
    {"_revdbayes_cpp_gp_flat",                  (DL_FUNC) &_revdbayes_cpp_gp_flat,                  2},
    {"_revdbayes_cpp_gp_flatflat",              (DL_FUNC) &_revdbayes_cpp_gp_flatflat,              2},
    {"_revdbayes_cpp_gp_jeffreys",              (DL_FUNC) &_revdbayes_cpp_gp_jeffreys,              2},
    {"_revdbayes_cpp_gp_loglik",                (DL_FUNC) &_revdbayes_cpp_gp_loglik,                2},
    {"_revdbayes_cpp_gp_mdi",                   (DL_FUNC) &_revdbayes_cpp_gp_mdi,                   2},
    {"_revdbayes_cpp_gp_norm",                  (DL_FUNC) &_revdbayes_cpp_gp_norm,                  2},
    {"_revdbayes_cpp_os_loglik",                (DL_FUNC) &_revdbayes_cpp_os_loglik,                2},
    {"_revdbayes_cpp_pp_loglik",                (DL_FUNC) &_revdbayes_cpp_pp_loglik,                2},
    {"_revdbayes_create_prior_xptr",            (DL_FUNC) &_revdbayes_create_prior_xptr,            1},
    {"_revdbayes_gev_beta_logpost",             (DL_FUNC) &_revdbayes_gev_beta_logpost,             2},
    {"_revdbayes_gev_beta_logpost_phi",         (DL_FUNC) &_revdbayes_gev_beta_logpost_phi,         2},
    {"_revdbayes_gev_flat_logpost",             (DL_FUNC) &_revdbayes_gev_flat_logpost,             2},
    {"_revdbayes_gev_flat_logpost_phi",         (DL_FUNC) &_revdbayes_gev_flat_logpost_phi,         2},
    {"_revdbayes_gev_flatflat_logpost",         (DL_FUNC) &_revdbayes_gev_flatflat_logpost,         2},
    {"_revdbayes_gev_flatflat_logpost_phi",     (DL_FUNC) &_revdbayes_gev_flatflat_logpost_phi,     2},
    {"_revdbayes_gev_loglognorm_logpost",       (DL_FUNC) &_revdbayes_gev_loglognorm_logpost,       2},
    {"_revdbayes_gev_loglognorm_logpost_phi",   (DL_FUNC) &_revdbayes_gev_loglognorm_logpost_phi,   2},
    {"_revdbayes_gev_logpost_phi_xptr",         (DL_FUNC) &_revdbayes_gev_logpost_phi_xptr,         1},
    {"_revdbayes_gev_logpost_xptr",             (DL_FUNC) &_revdbayes_gev_logpost_xptr,             1},
    {"_revdbayes_gev_mdi_logpost",              (DL_FUNC) &_revdbayes_gev_mdi_logpost,              2},
    {"_revdbayes_gev_mdi_logpost_phi",          (DL_FUNC) &_revdbayes_gev_mdi_logpost_phi,          2},
    {"_revdbayes_gev_norm_logpost",             (DL_FUNC) &_revdbayes_gev_norm_logpost,             2},
    {"_revdbayes_gev_norm_logpost_phi",         (DL_FUNC) &_revdbayes_gev_norm_logpost_phi,         2},
    {"_revdbayes_gev_phi_to_theta",             (DL_FUNC) &_revdbayes_gev_phi_to_theta,             2},
    {"_revdbayes_gev_prob_logpost",             (DL_FUNC) &_revdbayes_gev_prob_logpost,             2},
    {"_revdbayes_gev_prob_logpost_phi",         (DL_FUNC) &_revdbayes_gev_prob_logpost_phi,         2},
    {"_revdbayes_gev_quant_logpost",            (DL_FUNC) &_revdbayes_gev_quant_logpost,            2},
    {"_revdbayes_gev_quant_logpost_phi",        (DL_FUNC) &_revdbayes_gev_quant_logpost_phi,        2},
    {"_revdbayes_gev_user_logpost",             (DL_FUNC) &_revdbayes_gev_user_logpost,             2},
    {"_revdbayes_gev_user_logpost_phi",         (DL_FUNC) &_revdbayes_gev_user_logpost_phi,         2},
    {"_revdbayes_gp_beta_logpost",              (DL_FUNC) &_revdbayes_gp_beta_logpost,              2},
    {"_revdbayes_gp_beta_logpost_phi",          (DL_FUNC) &_revdbayes_gp_beta_logpost_phi,          2},
    {"_revdbayes_gp_flat_logpost",              (DL_FUNC) &_revdbayes_gp_flat_logpost,              2},
    {"_revdbayes_gp_flat_logpost_phi",          (DL_FUNC) &_revdbayes_gp_flat_logpost_phi,          2},
    {"_revdbayes_gp_flatflat_logpost",          (DL_FUNC) &_revdbayes_gp_flatflat_logpost,          2},
    {"_revdbayes_gp_flatflat_logpost_phi",      (DL_FUNC) &_revdbayes_gp_flatflat_logpost_phi,      2},
    {"_revdbayes_gp_jeffreys_logpost",          (DL_FUNC) &_revdbayes_gp_jeffreys_logpost,          2},
    {"_revdbayes_gp_jeffreys_logpost_phi",      (DL_FUNC) &_revdbayes_gp_jeffreys_logpost_phi,      2},
    {"_revdbayes_gp_logpost_phi_xptr",          (DL_FUNC) &_revdbayes_gp_logpost_phi_xptr,          1},
    {"_revdbayes_gp_logpost_xptr",              (DL_FUNC) &_revdbayes_gp_logpost_xptr,              1},
    {"_revdbayes_gp_mdi_logpost",               (DL_FUNC) &_revdbayes_gp_mdi_logpost,               2},
    {"_revdbayes_gp_mdi_logpost_phi",           (DL_FUNC) &_revdbayes_gp_mdi_logpost_phi,           2},
    {"_revdbayes_gp_norm_logpost",              (DL_FUNC) &_revdbayes_gp_norm_logpost,              2},
    {"_revdbayes_gp_norm_logpost_phi",          (DL_FUNC) &_revdbayes_gp_norm_logpost_phi,          2},
    {"_revdbayes_gp_phi_to_theta",              (DL_FUNC) &_revdbayes_gp_phi_to_theta,              2},
    {"_revdbayes_gp_user_logpost",              (DL_FUNC) &_revdbayes_gp_user_logpost,              2},
    {"_revdbayes_gp_user_logpost_phi",          (DL_FUNC) &_revdbayes_gp_user_logpost_phi,          2},
    {"_revdbayes_kgaps_log_j",                  (DL_FUNC) &_revdbayes_kgaps_log_j,                  2},
    {"_revdbayes_kgaps_logpost",                (DL_FUNC) &_revdbayes_kgaps_logpost,                2},
    {"_revdbayes_kgaps_logpost_xptr",           (DL_FUNC) &_revdbayes_kgaps_logpost_xptr,           1},
    {"_revdbayes_kgaps_phi_to_theta",           (DL_FUNC) &_revdbayes_kgaps_phi_to_theta,           2},
    {"_revdbayes_lgdgev_cpp",                   (DL_FUNC) &_revdbayes_lgdgev_cpp,                   4},
    {"_revdbayes_log_j_xptr",                   (DL_FUNC) &_revdbayes_log_j_xptr,                   1},
    {"_revdbayes_os_beta_logpost",              (DL_FUNC) &_revdbayes_os_beta_logpost,              2},
    {"_revdbayes_os_beta_logpost_phi",          (DL_FUNC) &_revdbayes_os_beta_logpost_phi,          2},
    {"_revdbayes_os_flat_logpost",              (DL_FUNC) &_revdbayes_os_flat_logpost,              2},
    {"_revdbayes_os_flat_logpost_phi",          (DL_FUNC) &_revdbayes_os_flat_logpost_phi,          2},
    {"_revdbayes_os_flatflat_logpost",          (DL_FUNC) &_revdbayes_os_flatflat_logpost,          2},
    {"_revdbayes_os_flatflat_logpost_phi",      (DL_FUNC) &_revdbayes_os_flatflat_logpost_phi,      2},
    {"_revdbayes_os_loglognorm_logpost",        (DL_FUNC) &_revdbayes_os_loglognorm_logpost,        2},
    {"_revdbayes_os_loglognorm_logpost_phi",    (DL_FUNC) &_revdbayes_os_loglognorm_logpost_phi,    2},
    {"_revdbayes_os_logpost_phi_xptr",          (DL_FUNC) &_revdbayes_os_logpost_phi_xptr,          1},
    {"_revdbayes_os_logpost_xptr",              (DL_FUNC) &_revdbayes_os_logpost_xptr,              1},
    {"_revdbayes_os_mdi_logpost",               (DL_FUNC) &_revdbayes_os_mdi_logpost,               2},
    {"_revdbayes_os_mdi_logpost_phi",           (DL_FUNC) &_revdbayes_os_mdi_logpost_phi,           2},
    {"_revdbayes_os_norm_logpost",              (DL_FUNC) &_revdbayes_os_norm_logpost,              2},
    {"_revdbayes_os_norm_logpost_phi",          (DL_FUNC) &_revdbayes_os_norm_logpost_phi,          2},
    {"_revdbayes_os_prob_logpost",              (DL_FUNC) &_revdbayes_os_prob_logpost,              2},
    {"_revdbayes_os_prob_logpost_phi",          (DL_FUNC) &_revdbayes_os_prob_logpost_phi,          2},
    {"_revdbayes_os_quant_logpost",             (DL_FUNC) &_revdbayes_os_quant_logpost,             2},
    {"_revdbayes_os_quant_logpost_phi",         (DL_FUNC) &_revdbayes_os_quant_logpost_phi,         2},
    {"_revdbayes_os_user_logpost",              (DL_FUNC) &_revdbayes_os_user_logpost,              2},
    {"_revdbayes_os_user_logpost_phi",          (DL_FUNC) &_revdbayes_os_user_logpost_phi,          2},
    {"_revdbayes_pgev_cpp",                     (DL_FUNC) &_revdbayes_pgev_cpp,                     4},
    {"_revdbayes_phi_to_theta_xptr",            (DL_FUNC) &_revdbayes_phi_to_theta_xptr,            1},
    {"_revdbayes_pp_beta_logpost",              (DL_FUNC) &_revdbayes_pp_beta_logpost,              2},
    {"_revdbayes_pp_beta_logpost_phi",          (DL_FUNC) &_revdbayes_pp_beta_logpost_phi,          2},
    {"_revdbayes_pp_flat_logpost",              (DL_FUNC) &_revdbayes_pp_flat_logpost,              2},
    {"_revdbayes_pp_flat_logpost_phi",          (DL_FUNC) &_revdbayes_pp_flat_logpost_phi,          2},
    {"_revdbayes_pp_flatflat_logpost",          (DL_FUNC) &_revdbayes_pp_flatflat_logpost,          2},
    {"_revdbayes_pp_flatflat_logpost_phi",      (DL_FUNC) &_revdbayes_pp_flatflat_logpost_phi,      2},
    {"_revdbayes_pp_loglognorm_logpost",        (DL_FUNC) &_revdbayes_pp_loglognorm_logpost,        2},
    {"_revdbayes_pp_loglognorm_logpost_phi",    (DL_FUNC) &_revdbayes_pp_loglognorm_logpost_phi,    2},
    {"_revdbayes_pp_logpost_phi_xptr",          (DL_FUNC) &_revdbayes_pp_logpost_phi_xptr,          1},
    {"_revdbayes_pp_logpost_xptr",              (DL_FUNC) &_revdbayes_pp_logpost_xptr,              1},
    {"_revdbayes_pp_mdi_logpost",               (DL_FUNC) &_revdbayes_pp_mdi_logpost,               2},
    {"_revdbayes_pp_mdi_logpost_phi",           (DL_FUNC) &_revdbayes_pp_mdi_logpost_phi,           2},
    {"_revdbayes_pp_norm_logpost",              (DL_FUNC) &_revdbayes_pp_norm_logpost,              2},
    {"_revdbayes_pp_norm_logpost_phi",          (DL_FUNC) &_revdbayes_pp_norm_logpost_phi,          2},
    {"_revdbayes_pp_phi_to_theta",              (DL_FUNC) &_revdbayes_pp_phi_to_theta,              2},
    {"_revdbayes_pp_prob_logpost",              (DL_FUNC) &_revdbayes_pp_prob_logpost,              2},
    {"_revdbayes_pp_prob_logpost_phi",          (DL_FUNC) &_revdbayes_pp_prob_logpost_phi,          2},
    {"_revdbayes_pp_quant_logpost",             (DL_FUNC) &_revdbayes_pp_quant_logpost,             2},
    {"_revdbayes_pp_quant_logpost_phi",         (DL_FUNC) &_revdbayes_pp_quant_logpost_phi,         2},
    {"_revdbayes_pp_user_logpost",              (DL_FUNC) &_revdbayes_pp_user_logpost,              2},
    {"_revdbayes_pp_user_logpost_phi",          (DL_FUNC) &_revdbayes_pp_user_logpost_phi,          2},
    {"_revdbayes_qgev_cpp",                     (DL_FUNC) &_revdbayes_qgev_cpp,                     4},
    {"_revdbayes_RcppExport_registerCCallable", (DL_FUNC) &_revdbayes_RcppExport_registerCCallable, 0},
    {"_revdbayes_user_gev_flat",                (DL_FUNC) &_revdbayes_user_gev_flat,                2},
    {"_revdbayes_user_gev_norm",                (DL_FUNC) &_revdbayes_user_gev_norm,                2},
    {"_revdbayes_user_gp_flat",                 (DL_FUNC) &_revdbayes_user_gp_flat,                 2},
    {NULL, NULL, 0}
};

void R_init_revdbayes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
