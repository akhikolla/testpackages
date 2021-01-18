#ifndef dblpm_H
#define dblpm_H

#include <RcppArmadillo.h>

class dblpm
{
public:
  
  dblpm (unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, // GLOBAL PARS
         arma::mat, // DATA
         arma::mat, arma::cube, arma::vec, arma::vec, // LIKELIHOOD PARS
         double, double, double, double, double, double, // OTHER PARS
         double, double, double, double, double, double, double, double, // HYPERPARS
         unsigned int, unsigned int, unsigned int, // MCMC SETTINGS
         double, double, double, double); // PROPOSAL PARS
  
  // GLOBAL PARAMETERS
  unsigned int T;
  unsigned int N;
  unsigned int M;
  unsigned int L;
  unsigned int D;
  
  // DATA
  arma::mat edgelist;
  
  // LIKELIHOOD PARAMETERS
  arma::mat x;
  arma::cube w;
  arma::vec gamma;
  arma::vec beta;
  
  // OTHER PARAMETERS
  double tauw;
  double tauw0;
  double taugamma;
  double taugamma0;
  double taubeta;
  double taubeta0;
  
  // HYPERPARAMETERS
  double taux;
  double delta;
  double aw;
  double bw;
  double agamma;
  double bgamma;
  double abeta;
  double bbeta;
  
  // STATISTICS
  arma::cube y;
  arma::mat out_degrees;
  arma::vec out_tot_degrees;
  arma::mat in_degrees;
  arma::vec in_tot_degrees;
  arma::vec j_first_activity;
  arma::vec j_last_activity;
  arma::mat j_activity_table;
  arma::field < arma::vec > i_activity_list;
  arma::mat i_activity_table;
  double N_active;
  arma::vec i_active;
  double M_active;
  arma::vec j_active;
  double w0_ss;
  double w_innovation_ss;
  double gamma_innovation_ss;
  double beta_innovation_ss;
  
  // INFERENCE
  double likelihood_value;
  double posterior_value;
  
  // BASIC METHODS
  void Print ();
  void EvaluateSumOfSquares ();
  void SetNoMissingData ();
  void FillActivity ();
  void FillY ();
  void Likelihood ();
  void Posterior ();
  
  // MCMC METHODS
  void UpdateX (unsigned int, double);
  void UpdateW (unsigned int, unsigned int, double);
  void UpdateGamma (unsigned int, double);
  void UpdateBeta (unsigned int, double);
  void UpdateTauw ();
  void UpdateTauw0 ();
  void UpdateTaugamma ();
  void UpdateTaugamma0 ();
  void UpdateTaubeta ();
  void UpdateTaubeta0 ();
  
  // MCMC SAMPLES
  arma::field < arma::mat > x_store;
  arma::field < arma::cube > w_store;
  arma::mat gamma_store;
  arma::mat beta_store;
  arma::vec tauw_store;
  arma::vec tauw0_store;
  arma::vec taugamma_store;
  arma::vec taugamma0_store;
  arma::vec taubeta_store;
  arma::vec taubeta0_store;
  arma::vec posterior_store;
  
  // MCMC SETTINGS
  unsigned int niter;
  unsigned int burnin;
  unsigned int thin;
  unsigned int total_niter;
  
  // PROPOSAL PARAMETERS
  double x_var;
  double w_var;
  double gamma_var;
  double beta_var;
  
  void MCMC(bool);
  
  bool debug = false;
  
};

#endif

