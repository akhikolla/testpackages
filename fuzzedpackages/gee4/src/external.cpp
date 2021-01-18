#include <RcppArmadillo.h>
#include <string>

#include "newton.h"
#include "ipw.h"
#include "gee_jmcm.h"

//'@title Fit Inverse Probability Weights Model
//'@description Fit inverse probability weights model.
//'@param m an integer vector of number of measurements for each subject.
//'@param Y a vector of responses for all subjects.
//'@param order the order for MAR remaining model.
//'@param trace whether or not the optimization iteration should be printed.
//'@export
// [[Rcpp::export]]
Rcpp::List ipw_estimation(arma::uvec m, arma::vec Y, arma::uword order, bool trace = false) {
  int debug = 0;
  if (debug) Rcpp::Rcout << "ipw_estimation()" << std::endl;

  gee::ipw weights(m, Y, order);
  dragonwell::Newton<gee::ipw> newt(weights);


  //trace = true;
  arma::vec x = arma::zeros<arma::vec>(order + 1);
  newt.Optimize(x, 1.0e-6, trace);
  arma::vec result = weights.CalWeights(x);
  //result.print("Weight = ");

  return Rcpp::List::create(Rcpp::Named("alpha") = x,
                            Rcpp::Named("pij") = weights.get_p(),
                            Rcpp::Named("cpij") = weights.get_Pi(),
                            Rcpp::Named("weights") = result);

}

//'@title Fit (Weighted) Generalized Estimating Equations based on MCD
//'@description Fit (weighted) generalized estimating equations based on MCD.
//'@param m an integer vector of number of measurements for each subject.
//'@param Y a vector of responses for all subjects.
//'@param X model matrix for mean structure model.
//'@param Z model matrix for the diagonal matrix.
//'@param W model matrix for the lower triangular matrix.
//'@param H a vector of weights used in WGEE-MCD.
//'@param method choose 'gee-mcd' (Ye and Pan, 2006) or 'wgee-mcd' (Pan et al. 2012).
//'@param corrStruct choose 'id' (independent), 'cs' (compound symmetry) or ar1' (AR(1)).
//'@param rho a parameter used in the 'working' covariance structure.
//'@param start starting values for the parameters in the model.
//'@param trace the values of the objective function and the parameters are
//'       printed for all the trace'th iterations.
//'@param profile whether parameters should be estimated sequentially using the
//'       idea of profile likelihood or not.
//'@param errorMsg whether or not the error message should be print.
//'@export
// [[Rcpp::export]]
Rcpp::List gees_estimation(arma::uvec m, arma::vec Y, arma::mat X, arma::mat Z, arma::mat W,
                           arma::vec H,
                           std::string method, std::string corrStruct, double rho, arma::vec start,
                           bool trace = false, bool profile = true, bool errorMsg = false)
{
  int debug = 0;
  int debug_test = 0;

  if (debug) Rcpp::Rcout << "gees_estimation(): " << std::endl;

  int n_bta = X.n_cols;
  int n_lmd = Z.n_cols;
  int n_gma = W.n_cols;

  if (debug) Rcpp::Rcout << "gees_estimation(): setting corr_mode..." << std::endl;
  gee_corr_mode corr_mode(0);
  if (corrStruct == "id") corr_mode.setid(1);
  else if (corrStruct == "cs") corr_mode.setid(2);
  else if (corrStruct == "ar1") corr_mode.setid(3);

  if (debug) Rcpp::Rcout << "gees_estimation(): creating gees object..." << std::endl;
  gee::gee_jmcm gees(m, Y, X, Z, W, rho, identity_link, corr_mode);
  if (method == "wgee-mcd") {
    //Rcpp::Rcout << "length(H) = " << H.n_rows << std::endl;
    //Rcpp::Rcout << "length(Y) = " << Y.n_elem << std::endl;
    gees.set_weights(H);
  }

  dragonwell::Newton<gee::gee_jmcm> newt(gees);
  dragonwell::LineSearch<dragonwell::NRfmin<gee::gee_jmcm>> linesearch;
  linesearch.set_message(errorMsg);

  arma::vec x = start;
  //double f_min = 0.0;
  int n_iters  = 0;

  if (profile) {
    if(debug) {
      Rcpp::Rcout << "gees_estimation(): Start profile optimization ... " << std::endl;
      x.t().print("start value: ");
    }
    gees.set_params(x);

    bool check = true;                  // check is true if the routine has converged to
                                        // a local minimum of the function fmin
    const arma::uword kMaxIters = 200;  // maximum number of iterations
    const double kTolF = 1.0e-8;        // convergence criterion on the function value
    const double kTolMin = 1.0e-12;     // criterion for deciding whether spurious
                                        // convergence to a minimum of fmin has occurred
    const double kStpMax = 100.0;       // scaled maximum step length allowed in line searches
    // const double kTolX = std::numeric_limits<double>::epsilon(); // convergence criterion on
    // delta x
    const double kTolX = 1e-10;

    const arma::uword n = x.n_elem;
    dragonwell::NRfmin<gee::gee_jmcm> fmin(gees); // Set up NRfmin object
    dragonwell::NRfdjac<gee::gee_jmcm> fdjac(gees); // Set up NRfdjac object
    double f = fmin(x);                             // 0.5 * fvec.t() * fvec
    arma::vec fvec = gees(x);                       // gradient vector
    arma::mat fjac = arma::zeros<arma::mat>(n, n);  // jacobian matrix

    // Test for initial guess being a root
    double test = 0.0;
    for (arma::uword i = 0; i != n; ++i) {
      if (std::abs(fvec(i)) > test) test = std::abs(fvec(i));
    }
    if (test < 0.01 * kTolF) {
      check = false;
      if (debug_test)
        Rcpp::Rcout << "Testing initial guess... " << std::endl
                    << "test = " << test << std::endl;
    }

    double sum = arma::as_scalar(x.t() * x);
    double stpmax = kStpMax * std::max(std::sqrt(sum), (double)n);
    for (arma::uword iter = 0; iter < kMaxIters; ++iter) {
      //if (!check) break;


      fjac = fdjac(x, fvec);
      arma::vec g = fjac.t() * fvec; // Compute delta f for the line search

      // arma::vec g = fvec;
      arma::vec xold = x;            // store x
      // double fold = f;               // store f

      if (debug) Rcpp::Rcout << "Update beta..." << std::endl;
      gees.UpdateBeta();
      if (debug) gees.get_beta().t().print("beta = ");
      if (debug) Rcpp::Rcout << "Update lambda..." << std::endl;
      gees.UpdateLambda();
      if (debug) gees.get_lambda().t().print("lambda = ");
      if (debug) Rcpp::Rcout << "Update gamma..." << std::endl;
      gees.UpdateGamma();
      if (debug) gees.get_gamma().t().print("gamma = ");
      if (debug) Rcpp::Rcout << "Update theta..." << std::endl;
      x = gees.get_theta();
      if (debug) x.t().print("xnew = ");
      arma::vec p = x - xold;

      check = linesearch.GetStep(fmin, &f, &x, g, p, stpmax);

      f = fmin(x);
      fvec = gees(x);

      if (trace) {
        Rcpp::Rcout << iter << ": " << f << ": " << x.t() << std::endl;
      }

      // Test for convergence on function values
      test = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        if (std::abs(fvec(i)) > test) test = std::abs(fvec(i));
      }

      if (test < kTolF) {
        check = false;
        if (debug_test)
          Rcpp::Rcout << "Testing convergence on function values... " << std::endl
                      << "test = " << test << std::endl;

        break;
      }

      // Check for gradient of f zero, i.e., spurious convergence
      if (check) {
        test = 0.0;
        double den = std::max(f, 0.5*n);
        for(arma::uword i = 0; i < n; ++i) {
          double temp = std::abs(g(i)) * std::max(std::abs(x(i)),1.0)/den;
          if(temp > test) test = temp;
        }
        check = (test < kTolMin);

        if (debug_test)
          Rcpp::Rcout << "Testing convergence on gradient of f zero... " << std::endl
                      << "test = " << test << std::endl;

        // break;  // return check;
      }

      // Test for convergence on delta x
      test = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        double temp = std::abs(x(i) - xold(i)) / std::max(std::abs(x(i)), 1.0);
        if (temp > test) test = temp;
      }
      if (debug_test)
        Rcpp::Rcout << "kTolX = " << kTolX << std::endl
                    << "test = " << test << std::endl;
      if (test < kTolX) {
        if (debug_test)
          Rcpp::Rcout << "Testing convergence on delta x... " << std::endl
                      << "test = " << test << std::endl;
        break;
      }
      if (debug_test)
        Rcpp::Rcout << "kTolX = " << kTolX << std::endl;
    }
  } else {
    if (debug) Rcpp::Rcout << "gees_estimation(): starting non-profile estimation..." << std::endl;
    newt.Optimize(x, 1.0e-6, trace);
    //f_min = newt.f_min();
    n_iters = newt.n_iters();
  }

  arma::vec beta   = x.rows(0, n_bta-1);
  arma::vec lambda = x.rows(n_bta, n_bta+n_lmd-1);
  arma::vec gamma  = x.rows(n_bta+n_lmd, n_bta+n_lmd+n_gma-1);

  if (debug) {
    Rcpp::Rcout << "Estimation finished" << std::endl;
    beta.t().print("beta = ");
    lambda.t().print("lambda = ");
    gamma.t().print("gamma = ");
  }

  int n_par = n_bta + n_lmd + n_gma;
  int n_sub = m.n_rows;

  double quasilik = gees.get_quasi_likelihood(x);

  return Rcpp::List::create(Rcpp::Named("par") = x,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("quasilik") = quasilik,
                            Rcpp::Named("QIC") = -2 * quasilik / n_sub + 2 * n_par / n_sub,
                            Rcpp::Named("iter") = n_iters);

}

// [[Rcpp::export]]
Rcpp::List geerfit_id(arma::uvec m,
                      arma::vec Y,
                      arma::mat X,
                      arma::mat Z,
                      arma::mat W,
                      double rho,
                      arma::vec start,
                      bool trace = false,
                      bool profile = true,
                      bool errorMsg = false) {
  //int debug = 1;

  int n_bta = X.n_cols;
  int n_lmd = Z.n_cols;
  int n_gma = W.n_cols;

  gee::gee_jmcm gees(m, Y, X, Z, W, rho, identity_link, Identity_corr);
  dragonwell::Newton<gee::gee_jmcm> newt(gees);
  arma::vec x = start;

  int n_iters  = 0;
  //double f_min = 0.0;
  newt.Optimize(x, 1.0e-6, trace);
  //f_min = newt.f_min();
  n_iters = newt.n_iters();
  //else {
  //   gees.learn(m, Y, X, Z, W, identity_link, Identity_corr, rho, x, 1000, trace);
  // }

  arma::vec beta   = x.rows(0, n_bta-1);
  arma::vec lambda = x.rows(n_bta, n_bta+n_lmd-1);
  arma::vec gamma  = x.rows(n_bta+n_lmd, n_bta+n_lmd+n_gma-1);

  int n_par = n_bta + n_lmd + n_gma;
  int n_sub = m.n_rows;

  double quasilik = gees.get_quasi_likelihood(x);

  return Rcpp::List::create(Rcpp::Named("par") = x,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("quasilik") = quasilik,
                            Rcpp::Named("QIC") = -2 * quasilik / n_sub + 2 * n_par / n_sub,
                            Rcpp::Named("iter") = n_iters);

}

// [[Rcpp::export]]
Rcpp::List geerfit_cs(arma::uvec m,
                      arma::vec Y,
                      arma::mat X,
                      arma::mat Z,
                      arma::mat W,
                      double rho,
                      arma::vec start,
                      bool trace = false,
                      bool profile = true,
                      bool errorMsg = false) {
  int n_bta = X.n_cols;
  int n_lmd = Z.n_cols;
  int n_gma = W.n_cols;

  gee::gee_jmcm gees(m, Y, X, Z, W, rho, identity_link, CompSymm_corr);
  dragonwell::Newton<gee::gee_jmcm> newt(gees);
  arma::vec x = start;

  //double f_min = 0.0;
  int n_iters  = 0;

  newt.Optimize(x, 1.0e-6, trace);
  //f_min = newt.f_min();
  n_iters = newt.n_iters();

  arma::vec beta   = x.rows(0, n_bta-1);
  arma::vec lambda = x.rows(n_bta, n_bta+n_lmd-1);
  arma::vec gamma  = x.rows(n_bta+n_lmd, n_bta+n_lmd+n_gma-1);

  int n_par = n_bta + n_lmd + n_gma;
  int n_sub = m.n_rows;

  double quasilik = gees.get_quasi_likelihood(x);

  return Rcpp::List::create(Rcpp::Named("par") = x,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("quasilik") = quasilik,
                            Rcpp::Named("QIC") = -2 * quasilik / n_sub + 2 * n_par / n_sub,
                            Rcpp::Named("iter") = n_iters);

}

// [[Rcpp::export]]
Rcpp::List geerfit_ar1(arma::uvec m,
                       arma::vec Y,
                       arma::mat X,
                       arma::mat Z,
                       arma::mat W,
                       double rho,
                       arma::vec start,
                       bool trace = false,
                       bool profile = true,
                       bool errorMsg = false) {
  int n_bta = X.n_cols;
  int n_lmd = Z.n_cols;
  int n_gma = W.n_cols;

  gee::gee_jmcm gees(m, Y, X, Z, W, rho, identity_link, AR1_corr);
  dragonwell::Newton<gee::gee_jmcm> newt(gees);
  arma::vec x = start;

  //double f_min = 0.0;
  int n_iters  = 0;

  newt.Optimize(x, 1.0e-6, trace);
  //f_min = newt.f_min();
  n_iters = newt.n_iters();

  arma::vec beta   = x.rows(0, n_bta-1);
  arma::vec lambda = x.rows(n_bta, n_bta+n_lmd-1);
  arma::vec gamma  = x.rows(n_bta+n_lmd, n_bta+n_lmd+n_gma-1);

  int n_par = n_bta + n_lmd + n_gma;
  int n_sub = m.n_rows;

  double quasilik = gees.get_quasi_likelihood(x);

  return Rcpp::List::create(Rcpp::Named("par") = x,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("gamma") = gamma,
                            Rcpp::Named("quasilik") = quasilik,
                            Rcpp::Named("QIC") = -2 * quasilik / n_sub + 2 * n_par / n_sub,
                            Rcpp::Named("iter") = n_iters);

}

// RCPP_EXPOSED_CLASS(gee_link_mode);
// RCPP_EXPOSED_CLASS(gee_corr_mode);
// RCPP_MODULE(mod_gee_jmcm) {
//   using namespace Rcpp;
//   using namespace arma;
//   using namespace gee;

//   class_<gee_jmcm>("gee_jmcm")

//     .constructor<uvec,vec,mat,mat,mat,double,gee_link_mode,gee_corr_mode>()

//     .property("m_", &gee_jmcm::get_m)
//     .property("Y_", &gee_jmcm::get_Y)
//     .property("X_", &gee_jmcm::get_X)
//     ;
// }

RcppExport SEXP gee_jmcm__new(SEXP m_, SEXP Y_, SEXP X_, SEXP Z_, SEXP W_,
                              SEXP corrStruct_, SEXP rho_) {
  arma::uvec m = Rcpp::as<arma::uvec>(m_);
  arma::vec Y = Rcpp::as<arma::vec>(Y_);
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Z = Rcpp::as<arma::mat>(Z_);
  arma::mat W = Rcpp::as<arma::mat>(W_);

  std::string corrStruct = Rcpp::as<std::string>(corrStruct_);
  double rho = Rcpp::as<double>(rho_);

  gee_corr_mode corr_mode(0);
  if (corrStruct == "id") corr_mode.setid(1);
  else if (corrStruct == "cs") corr_mode.setid(2);
  else if (corrStruct == "ar1") corr_mode.setid(3);

  Rcpp::XPtr<gee::gee_jmcm>
    ptr(new gee::gee_jmcm(m, Y, X, Z, W, rho, identity_link, corr_mode), true);

  return ptr;
}

RcppExport SEXP gee_jmcm__get_m(SEXP xp, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);
  int i = Rcpp::as<int>(i_) - 1;

  return Rcpp::wrap(ptr->get_m(i));
}

RcppExport SEXP gee_jmcm__get_Y(SEXP xp, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);
  int i = Rcpp::as<int>(i_) - 1;

  return Rcpp::wrap(ptr->get_Y(i));
}

RcppExport SEXP gee_jmcm__get_X(SEXP xp, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);
  int i = Rcpp::as<int>(i_) - 1;

  return Rcpp::wrap(ptr->get_X(i));
}

RcppExport SEXP gee_jmcm__get_Z(SEXP xp, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);
  int i = Rcpp::as<int>(i_) - 1;

  return Rcpp::wrap(ptr->get_Z(i));
}

RcppExport SEXP gee_jmcm__get_W(SEXP xp, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);
  int i = Rcpp::as<int>(i_) - 1;

  return Rcpp::wrap(ptr->get_W(i));
}

RcppExport SEXP gee_jmcm__get_D(SEXP xp, SEXP x_, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);

  arma::vec x = Rcpp::as<arma::vec>(x_);
  int i = Rcpp::as<int>(i_) - 1;

  ptr->UpdateGEES(x);

  return Rcpp::wrap(ptr->get_D(i));
}

RcppExport SEXP gee_jmcm__get_T(SEXP xp, SEXP x_, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);

  arma::vec x = Rcpp::as<arma::vec>(x_);
  int i = Rcpp::as<int>(i_) - 1;

  ptr->UpdateGEES(x);

  return Rcpp::wrap(ptr->get_T(i));
}

RcppExport SEXP gee_jmcm__get_mu(SEXP xp, SEXP x_, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);

  arma::vec x = Rcpp::as<arma::vec>(x_);
  int i = Rcpp::as<int>(i_) - 1;

  ptr->UpdateGEES(x);

  return Rcpp::wrap(ptr->get_mu(i));
}

RcppExport SEXP gee_jmcm__get_Sigma(SEXP xp, SEXP x_, SEXP i_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);

  arma::vec x = Rcpp::as<arma::vec>(x_);
  int i = Rcpp::as<int>(i_) - 1;

  arma::mat Sigmai;

  ptr->UpdateGEES(x);

  return Rcpp::wrap(ptr->get_Sigma(i));
}

RcppExport SEXP gee_jmcm__get_fim(SEXP xp, SEXP x_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);

  arma::vec x = Rcpp::as<arma::vec>(x_);

  ptr->UpdateGEES(x);

  return Rcpp::wrap(ptr->get_fim());
}

RcppExport SEXP gee_jmcm__get_sd(SEXP xp, SEXP x_) {
  Rcpp::XPtr<gee::gee_jmcm> ptr(xp);

  arma::vec x = Rcpp::as<arma::vec>(x_);

  ptr->UpdateGEES(x);

  return Rcpp::wrap(ptr->get_sd());
}

RcppExport SEXP ipw__new(SEXP m_, SEXP Y_, SEXP order_) {
  arma::uvec m = Rcpp::as<arma::uvec>(m_);
  arma::vec Y = Rcpp::as<arma::vec>(Y_);
  arma::uword order = Rcpp::as<arma::uword>(order_);

  Rcpp::XPtr<gee::ipw>
    ptr(new gee::ipw(m, Y, order), true);

  return ptr;
}

RcppExport SEXP ipw__get_p(SEXP xp, SEXP alpha_) {
  Rcpp::XPtr<gee::ipw> ptr(xp);

  arma::vec alpha = Rcpp::as<arma::vec>(alpha_);
  arma::vec weights = ptr->CalWeights(alpha);

  return Rcpp::wrap(ptr->get_p());
}

RcppExport SEXP ipw__get_Pi(SEXP xp, SEXP alpha_) {
  Rcpp::XPtr<gee::ipw> ptr(xp);

  arma::vec alpha = Rcpp::as<arma::vec>(alpha_);
  arma::vec weights = ptr->CalWeights(alpha);

  return Rcpp::wrap(ptr->get_Pi());
}
