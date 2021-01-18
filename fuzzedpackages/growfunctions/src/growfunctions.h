/* 
 * File:   growfunctions.h
 * Author: savitsky
 *
 * Created on January 28, 2014, 11:45 AM
 */

#ifndef   GROWFUNCTIONS_H
#define   GROWFUNCTIONS_H

#include <RcppArmadillo.h>
#include <time.h>

RcppExport SEXP GPDPMIX(SEXP Ymat, SEXP Otrend, SEXP Oseas, SEXP o_gp_mod, SEXP o_jitter,
           SEXP o_a, SEXP o_b, SEXP o_atau, SEXP o_btau,SEXP _lower, SEXP _upper,
           SEXP o_w_star, SEXP o_w, SEXP o_n_slice_iter, 
           SEXP o_y_index, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
           SEXP ntuneInt, SEXP Minit, SEXP shapealph, SEXP ratebeta, SEXP o_progress, SEXP o_ipr);
RcppExport SEXP GPDPMIXMIS(SEXP Ymat, SEXP Otrend, SEXP Oseas, SEXP o_gp_mod, SEXP o_jitter,
           SEXP o_b_move, SEXP o_a, SEXP o_b, SEXP o_atau, SEXP o_btau, SEXP o_lower, SEXP o_upper,
           SEXP o_w_star, SEXP o_w, SEXP o_n_slice_iter, 
           SEXP o_y_index, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
           SEXP ntuneInt, SEXP Minit, SEXP shapealph, SEXP ratebeta, SEXP o_progress, SEXP o_ipr);
RcppExport SEXP GP(SEXP Ymat, SEXP Otrend, SEXP Oseas, SEXP o_gp_mod, SEXP o_jitter,
           SEXP o_a, SEXP o_b, SEXP o_atau, SEXP o_btau, SEXP o_lower, SEXP o_upper,
           SEXP o_w, SEXP o_n_slice_iter, 
           SEXP o_y_index, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
           SEXP ntuneInt, SEXP o_progress, SEXP o_ipr);
RcppExport SEXP GPFIX(SEXP Ymat, SEXP Otrend, SEXP Oseas, SEXP o_gp_mod, SEXP o_jitter,
           SEXP o_a, SEXP o_b, SEXP o_atau, SEXP o_btau, SEXP o_lower, SEXP o_upper,
           SEXP o_w, SEXP o_n_slice_iter, 
           SEXP o_y_index, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
           SEXP ntuneInt, SEXP o_progress, SEXP o_s, SEXP o_ipr);
RcppExport SEXP GPBFIX(SEXP Ymat, SEXP Otrend, SEXP Oseas, SEXP o_gp_mod, SEXP o_jitter,
           SEXP o_a, SEXP o_b, SEXP o_atau, SEXP o_btau, SEXP o_lower, SEXP o_upper,
           SEXP o_w, SEXP o_n_slice_iter, 
           SEXP o_y_index, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
           SEXP ntuneInt, SEXP o_progress, SEXP o_s, SEXP o_ipr);
RcppExport SEXP IGMRFDPMIX(SEXP Ymat, SEXP o_ksi, SEXP o_C, SEXP o_D, SEXP o_order,  
           SEXP niterInt, SEXP nburnInt, SEXP nthinInt,
           SEXP Minit, SEXP o_w_star, SEXP o_a, SEXP o_b, SEXP o_a_tau, SEXP o_b_tau,
           SEXP shapealph, SEXP ratebeta, SEXP o_nu,
           SEXP o_progress, SEXP o_jitter, SEXP o_kappa_fast, SEXP o_ipr);
RcppExport SEXP IGMRFDPMIXCOUNT(SEXP Ymat, SEXP o_E, SEXP o_ksi, SEXP o_C, SEXP o_D, 
           SEXP o_order, SEXP niterInt, SEXP nburnInt, SEXP nthinInt,
           SEXP Minit, SEXP o_w_star, SEXP o_a, SEXP o_b, SEXP o_a_tau, SEXP o_b_tau,
           SEXP shapealph, SEXP ratebeta, SEXP o_nu, SEXP o_Rep, 
           SEXP o_progress, SEXP o_jitter, SEXP o_kappa_fast, SEXP o_stable_launch,
           SEXP o_ipr);
RcppExport SEXP predict_bb(SEXP res, SEXP o_Omegas_tetr, SEXP o_Omegas_tete,
                    SEXP o_Omegat_tetr, SEXP o_Omegat_tete, SEXP o_J);
RcppExport SEXP predict_gmrf_bb(SEXP res, SEXP o_R, SEXP o_J);
/* utility functions */
SEXP wishrnd(arma::mat& L, const arma::mat& V, double nu);
arma::colvec rdirich(arma::colvec& shape);
double logmatrixdens(const arma::mat& B_i, const arma::mat& P, const arma::mat& Q);
double logmvndens(const arma::colvec& b_i, const arma::colvec& m, const arma::mat& Q);
double loggmrfdens_full(const arma::colvec& b_i, const arma::colvec& m, const arma::mat& Q,
                    const arma::vec& eigraw, double kappa);
double loggmrfdens(const arma::colvec& b_i, const arma::colvec& m, const arma::mat& Q,
                    int df, double kappa);
double log_dnorm_vec(const arma::rowvec& y, const arma::rowvec& b, double tau_e);
double dev(const arma::colvec& resid, double taue);
SEXP dmarg(const arma::colvec& resid, double taue, arma::rowvec& devmarg);
SEXP cpo(const arma::mat& Devmarg, arma::rowvec& logcpo, double& lpml);
SEXP dic3comp(const arma::colvec& Deviance, const arma::mat& Devmarg, arma::colvec& devres);
SEXP rmvnsample(const arma::mat& phi, const arma::colvec& h, arma::colvec& b);
SEXP rmvncov(const arma::mat& phi_inv, const arma::colvec& h, arma::colvec& b);
SEXP rmvnchol(const arma::mat& U, const arma::colvec& h, arma::colvec& b);
SEXP rmvnbasic(const arma::mat& phi, const arma::colvec& e, arma::colvec& b);
unsigned long rdrawone(const arma::colvec& pr, unsigned long k);
/* gpdpmix_moves */
SEXP gen_P(arma::uvec& gp_mod, arma::uvec& P_vec);
arma::mat gen_C(const arma::colvec& thetastar_m, double tau_e, const arma::mat& Omega_t, 
          const arma::cube& Omega_s, double jitter, arma::uvec& gp_mod, 
          arma::uvec& n_parms, arma::uvec& pos_s, int noise);
arma::mat gen_Cterm(const arma::colvec& thetastar_m, double tau_e, const arma::mat& Omega_t, 
               const arma::mat& Omega_s, double jitter, 
               int gp_mod_term, int noise);
arma::mat gen_Casym(const arma::colvec& theta_i, const arma::mat& Omega_t, const arma::cube& Omega_s,
               arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s);
double logy_like(int i, const arma::mat& y,  arma::ucolvec& s, 
                    const arma::cube& U_last);
arma::mat compute_Um(const arma::colvec& thetastar_m, double tau_e, double jitter,  
                    arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s, 
                    const arma::mat& Omega_t, const arma::cube& Omega_s, int noise);
arma::mat compute_Upm(double thetastar_pm, const arma::mat& theta_star, double tau_e, double jitter, int p, 
                    int m, const arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s, const arma::mat& Omega_t, 
                    const arma::cube& Omega_s, int noise);
SEXP compute_U(const arma::mat& theta_star, double tau_e, double jitter, arma::uvec& gp_mod,
                    arma::uvec& n_parms, arma::uvec& pos_s, const arma::mat& Omega_t, 
                    const arma::cube& Omega_s, int noise, arma::cube& U);
double logFm_like(const arma::mat& y, int m, const arma::ucolvec& s, 
                    const arma::mat& U_m, const arma::vec& ipr);
double logFpm_post(double thetastar_pm, const arma::mat& theta_star, double tau_e, double jitter, int p, 
                    int m, arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s, 
                    const arma::mat& Omega_t, const arma::cube& Omega_s, int noise, 
                    const arma::mat& y, const arma::ucolvec& s, double a, double b, const arma::vec& ipr);
double logFtau_post(const arma::mat& theta_star, double tau_e, double jitter, arma::uvec& gp_mod,
                    arma::uvec& n_parms, arma::uvec& pos_s, const arma::mat& Omega_t, 
                    const arma::cube& Omega_s, int noise, const arma::mat& y, 
                    const arma::ucolvec& s, double a, double b, const arma::vec& ipr);
double logFtau_like(const arma::mat& y, const arma::ucolvec& s, const arma::cube& U_last,
                    const arma::vec& ipr);
double log_prior(double theta_star_pm, double a, double b);
SEXP uni_slice_pm(arma::colvec& theta_updown, arma::mat& Pi_n, int pos_ud, int dist,
          const arma::mat& theta_star, double tau_e, int p, int m, int& slice_evals_pm,
          double lower, double upper, int n_slice_iter, const arma::mat& w_tot,
          const arma::ucolvec& s, const arma::mat& Omegat_n, const arma::cube& Omegas_n,
          arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s, int noise, double jitter,
          const arma::mat& y_n, double a, double b,  int transition, const arma::vec& ipr);
arma::ucolvec temper_dist_alt(int n);
int temper_dist_compute(int i, int n);
double temper_logpmove_compute(arma::mat& Pi_n);
double temper_logpmove_alt(arma::mat& Pi_n, const arma::ucolvec& dist);
SEXP temper(arma::cube& U_last, arma::mat& theta_star, double& tau_e, arma::uvec& gp_mod, 
          arma::uvec& n_parms, arma::uvec& pos_s, double& slice_levals_theta, double& slice_levals_tau, 
          double& n_slice_theta, double& n_slice_tau, double& accept_temper, 
          double& n_temper_evals, double lower, double upper, 
          int n_slice_iter, const arma::mat& w_tot, const arma::ucolvec& s, 
          const arma::field<arma::mat>& Omegat_ns, 
          const arma::field<arma::cube>& Omegas_ns, const arma::mat& Omega_t, const arma::cube& Omega_s, 
          const arma::field<arma::mat>& y_ns, const arma::mat& y, int noise, double jitter,
          double a, double b, double atau, double btau, const arma::vec& ipr);
SEXP temper_b(arma::cube& U_last, arma::mat& theta_star, double tau_e, arma::uvec& gp_mod, 
          arma::uvec& n_parms, arma::uvec& pos_s, double& slice_levals_theta, 
          double& n_slice_theta, double& accept_temper, 
          double& n_temper_evals, double lower, double upper, 
          int n_slice_iter, const arma::mat& w_tot, const arma::ucolvec& s, 
          const arma::field<arma::mat>& Omegat_ns, 
          const arma::field<arma::cube>& Omegas_ns, const arma::mat& Omega_t, const arma::cube& Omega_s, 
          const arma::field<arma::mat>& bb_ns, const arma::mat& bb, double jitter,
          double a, double b, const arma::vec& ipr);
SEXP auxclusterstep(arma::mat& theta_star, arma::mat& wpm, arma::cube& U_last, const arma::mat& Omega_t, 
            const arma::cube& Omega_s, const arma::mat& y, double tau_e, int noise, double jitter, 
            arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s,
            arma::ucolvec& s, arma::ucolvec& num, unsigned int& M, const int& w_star, double& conc,
            double a, double b, const arma::vec& ipr, arma::colvec& Num);
SEXP concstep(double& conc, int M, int N,  double a6, double b6);
SEXP wp_tune(arma::cube& Theta_tune, arma::colvec& wtune);
SEXP wpm_tune(arma::cube& Theta_tune, arma::mat& wtune);
SEXP wtau_tune(arma::colvec& Taue_tune, double& wtune);
arma::mat update_w(const arma::mat& wpm, double wtau);
SEXP wpm_aux(arma::mat& wpm_aux_h, const arma::mat& theta_aux_h);
SEXP gen_bb(const arma::mat& y, const arma::mat& Tau_e, arma::umat& S, 
               const arma::field<arma::mat>& Theta_star, const arma::mat& Omega_t, 
               const arma::cube& Omega_s, arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s,
               arma::uvec& P_vec, arma::mat& bb, arma::field<arma::mat>& f, 
               arma::field<arma::cube>& invG_star, double jitter);
SEXP gen_f(const arma::mat& BB, const arma::mat& Tau_e, arma::umat& S, 
               const arma::field<arma::mat>& Theta_star, const arma::mat& Omega_t, 
               const arma::cube& Omega_s, arma::uvec& gp_mod, arma::uvec& n_parms, arma::uvec& pos_s,
               arma::uvec& P_vec, arma::field<arma::mat>& f, 
               const arma::field<arma::cube>& invG_star, double jitter);  
SEXP move_b(arma::mat& bb, arma::cube& invG_star, const arma::cube& U_last, const arma::ucolvec& s, 
               const arma::mat& y, double tau_e);
SEXP move_bslice(arma::mat& bb, const arma::cube& U_last, const arma::ucolvec& s, 
               const arma::mat& y, double tau_e, int R);
SEXP gen_bb_ns(const arma::mat& bb, const arma::field<arma::uvec>& index, 
               arma::field<arma::mat>& bb_ns);
SEXP lsqcluster(arma::umat& S, arma::ucolvec& ordscore, arma::mat& phat, 
                arma::field<arma::ucolvec>& bigS);
SEXP pop_Num(const arma::ucolvec& s, const arma::vec& ipr, arma::colvec& Num);
/* dpmix_moves */
SEXP clusterstep(const arma::cube& B, arma::mat& kappa_star, arma::mat& B1, 
            const arma::uvec& o, const arma::field<arma::mat>& c, 
            //const arma::field<arma::sp_mat>& c,
            const arma::mat& D, arma::ucolvec& s, 
            arma::ucolvec& num, unsigned int& M, double& conc, int a, int b, const arma::vec& ipr,
            arma::colvec& Num);
SEXP clusterstep_alt(const arma::cube& B, arma::mat& kappa_star, arma::mat& B1, const arma::uvec& o,
            const arma::cube& Q, arma::ucolvec& s, 
            arma::ucolvec& num, unsigned int& M, double& conc, int a, int b, const arma::vec& ipr,
            arma::colvec& Num);
SEXP move_kappastar(arma::mat& kappa_star, const arma::mat& B1, const arma::ucolvec& s, 
                    arma::uvec& o, int T, int a, int b, const arma::vec& ipr);
SEXP move_kappastar_alt(arma::mat& kappa_star, const arma::cube& B, const arma::cube& Q, 
                    const arma::ucolvec& s, 
                    arma::uvec& o, int T, int a, int b, const arma::vec& ipr);
SEXP move_B(const arma::mat& y, arma::cube& B, const arma::mat& kappa_star, 
               //const arma::field<arma::sp_mat>& C, 
               const arma::field<arma::mat>& C,
               arma::mat& gamma, const arma::mat& D, 
               const arma::ucolvec& s, double tau_e);
SEXP move_B_alt(const arma::mat& y, arma::cube& B, const arma::mat& kappa_star, 
               //const arma::field<arma::sp_mat>& C, 
               const arma::field<arma::mat>& C,
               arma::mat& gamma, const arma::mat& D, 
               const arma::ucolvec& s, double tau_e);
SEXP move_taue(const arma::mat& y, const arma::mat& gamma, double& tau_e, double a, double b,
               const arma::vec& ipr);
SEXP move_taue_jitter(const arma::mat& y, const arma::mat& gamma, double& tau_e, double a, double b,
                         double jitter, const arma::vec& ipr);
SEXP miss_ystep(arma::mat& y_rep, const arma::mat& y, const arma::mat& gamma, double tau_e);
/* sppmdpmix_moves */
SEXP auxclusterstep_gmrf(const arma::cube& B, const arma::mat& ksi, 
                    arma::mat& kappa_star, const arma::uvec& o,
                    const arma::field<arma::mat>& C, const arma::mat& D, 
                    arma::mat& u_star, arma::cube& Lambda_star,
                    arma::mat& as_star, double nu,
                    const arma::colvec& u_bar, const arma::mat& P_bar, 
                    arma::ucolvec& s, 
                    arma::ucolvec& num, unsigned int& M, const int& w_star,
                    double& conc, int a, int b,
                    const arma::vec& ipr, arma::colvec& Num);
SEXP move_ustar(arma::mat& u_star, const arma::mat& ksi,
                const arma::cube& Lambda_star, const arma::colvec& u_bar,
                const arma::mat& P_bar, const arma::ucolvec& s);
SEXP move_Lambdastar(arma::cube& Lambda_star, const arma::mat& as_star, 
                     const arma::mat& Ksi, const arma::mat& u_star, 
                     const arma::ucolvec& s, double nu);
SEXP move_ubar(arma::colvec& u_bar, const arma::mat& u_star,
               const arma::mat& P_bar);
SEXP move_Pbar(arma::mat& P_bar, const arma::colvec& u_bar, const arma::mat& u_star);
SEXP move_as(const arma::mat& P_8, arma::vec& as, double nu, double b);
double loglike_psi_i(const arma::colvec& psi_i, const arma::mat& Y, 
                     const arma::mat& E,
                     int i);
SEXP move_Psi_i(arma::mat& Psi, const arma::mat& Y, const arma::mat& E,  
                const arma::mat& gamma, double tau_e, int R);
SEXP miss_ycount(arma::mat& Y_rep, const arma::mat& Y, 
                 const arma::mat& E, const arma::mat& Psi);
SEXP dmarg_count(const arma::colvec& y_vec, const arma::colvec& mu_vec, 
                 arma::rowvec& devmarg);

#endif     /* GROWFUNCTIONS_H */

