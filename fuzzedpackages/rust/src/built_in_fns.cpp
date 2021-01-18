// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Miscellaneous functions.

// [[Rcpp::export]]
bool any_naC(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(Rcpp::is_na(x)));
}

// [[Rcpp::export]]
bool no_naC(const Rcpp::NumericVector& x) {
  return Rcpp::is_false(Rcpp::any(Rcpp::is_na(x)));
}

// [[Rcpp::export]]
bool any_nonpos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x <= 0));
}

// [[Rcpp::export]]
bool all_pos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::all(x > 0));
}

// [[Rcpp::export]]
bool any_neg(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x < 0));
}

// [[Rcpp::export]]
bool any_col_nonneg(const Rcpp::NumericMatrix& x) {
  int ncols = x.ncol();
  for(int i=0; i < ncols; i++) {
    if (!any_neg(x(_,i))) {
      return true ;
    }
  }
  return false ;
}

// [[Rcpp::export]]
bool any_pos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x > 0));
}

// [[Rcpp::export]]
bool any_col_nonpos(const Rcpp::NumericMatrix& x) {
  int ncols = x.ncol();
  for(int i=0; i < ncols; i++) {
    if (!any_pos(x(_,i))) {
      return true ;
    }
  }
  return false ;
}

// [[Rcpp::export]]
Rcpp::NumericVector vecpow(const Rcpp::NumericVector& base,
                           const Rcpp::NumericVector& exp) {
  int n = base.size() ;
  Rcpp::NumericVector res(n) ;
  for(int i=0; i < n; i++) {
    res[i] = std::pow(base[i], exp[i]) ;
  }
  return res;
}

// Functions to perform (inverse) transformations of variable.

// Rotation and mode relocation.

// [[Rcpp::export]]
arma::vec cpp_rho_to_psi(const arma::vec& rho, const arma::vec& psi_mode,
                         const arma::mat& rot_mat) {
  arma::vec y = rot_mat * rho + psi_mode ;
  return Rcpp::as<Rcpp::NumericVector>(wrap(y)) ;
}

// Inverse Box-Cox transformation.

// No lambda = 0.

// [[Rcpp::export]]
Rcpp::NumericVector cpp_psi_to_phi(const Rcpp::NumericVector& psi,
                                   const Rcpp::NumericVector& lambda,
                                   const Rcpp::NumericVector& gm,
                                   const Rcpp::NumericVector& con) {
  return Rcpp::wrap(vecpow(psi * con + 1.0, 1.0 / lambda)) ;
}

// At least one lambda = 0.

// [[Rcpp::export]]
Rcpp::NumericVector cpp_psi_to_phi_0(const Rcpp::NumericVector& psi,
                                     const Rcpp::NumericVector& lambda,
                                     const Rcpp::NumericVector& gm,
                                     const Rcpp::NumericVector& con) {
  return Rcpp::wrap( ifelse(lambda == 0, exp(psi / gm),
                            vecpow(psi * con + 1.0, 1.0 / lambda) )) ;
}

// Functions to calculate target log-density.

// Original logf (no transformation).

// [[Rcpp::export]]
double cpp_logf(const Rcpp::NumericVector& theta, const SEXP& logf,
                const Rcpp::List& pars) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  double val ;
  val = fun(theta, pars) ;
  return val ;
}

// Original logf (no transformation) but scaled by the maximum value
// of the logged target density (perhaps after transformation).
// Used only in plot.ru() when d = 2 (to avoid over/under-flow).

// [[Rcpp::export]]
double cpp_logf_scaled(const Rcpp::NumericVector& theta, const SEXP& logf,
                       const Rcpp::List& pars) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  double val ;
  double hscale = pars["hscale"] ;
  val = fun(theta, pars) - hscale ;
  return val ;
}

// Case 1: rotation and relocation only.

// [[Rcpp::export]]
double cpp_logf_rho(const arma::vec& rho, const arma::vec& psi_mode,
                    const arma::mat& rot_mat, const double& hscale,
                    const SEXP& logf, const Rcpp::List& pars) {
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  Rcpp::NumericVector theta ;
  double val ;
  theta = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
  val = fun(theta, pars) - hscale ;
  return val ;
}

// Case 2: Box-Cox, rotation and relocation.

// [[Rcpp::export]]
double cpp_logf_rho_2(const arma::vec& rho, const arma::vec& psi_mode,
                      const arma::mat& rot_mat, const double& hscale,
                      const SEXP& logf, const Rcpp::List& pars,
                      const Rcpp::List& tpars, const SEXP& ptpfun,
                      const SEXP& phi_to_theta, const SEXP& log_j,
                      const Rcpp::List& user_args) {
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  // Unwrap pointer to psi_to_phi transformation function.
  typedef Rcpp::NumericVector (*ptpPtr)(const Rcpp::NumericVector& psi,
                               const Rcpp::NumericVector& lambda,
                               const Rcpp::NumericVector& gm,
                               const Rcpp::NumericVector& con) ;
  Rcpp::XPtr<ptpPtr> xptpfun(ptpfun) ;
  ptpPtr psi_to_phi_fun = *xptpfun ;
  Rcpp::NumericVector lambda = tpars["lambda"] ;
  Rcpp::NumericVector gm = tpars["gm"] ;
  Rcpp::NumericVector con = tpars["con"] ;
  Rcpp::IntegerVector which_lam = tpars["which_lam"] ;
  Rcpp::NumericVector phi, psi, phi2, temp, temp2 ;
  double val, log_bc_jac ;
  psi = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
  temp = psi * con + 1.0 ;
  temp = temp[which_lam] ;
  if (any_nonpos(temp)) {
    return R_NegInf ;
  }
  phi = psi_to_phi_fun(psi, lambda, gm, con) ;
  phi2 = phi[which_lam] ;
  temp = Rcpp::log(phi2) ;
  temp2 = lambda[which_lam] ;
  log_bc_jac = sum((temp2 - 1.0) * temp) ;
  val = fun(phi, pars) - log_bc_jac - hscale ;
  return val ;
}

// Case 3: Transformation to positivity, Box-Cox, rotation and relocation.

// [[Rcpp::export]]
double cpp_logf_rho_3(const arma::vec& rho, const arma::vec& psi_mode,
                      const arma::mat& rot_mat, const double& hscale,
                      const SEXP& logf, const Rcpp::List& pars,
                      const Rcpp::List& tpars, const SEXP& ptpfun,
                      const SEXP& phi_to_theta, const SEXP& log_j,
                      const Rcpp::List& user_args) {
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  // Unwrap pointer to psi_to_phi transformation function.
  typedef Rcpp::NumericVector (*ptpPtr)(const Rcpp::NumericVector& psi,
                               const Rcpp::NumericVector& lambda,
                               const Rcpp::NumericVector& gm,
                               const Rcpp::NumericVector& con) ;
  Rcpp::XPtr<ptpPtr> xptpfun(ptpfun) ;
  ptpPtr psi_to_phi_fun = *xptpfun ;
  // Unwrap pointer to phi_to_theta transformation function.
  typedef Rcpp::NumericVector (*pttPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  Rcpp::XPtr<pttPtr> xpttfun(phi_to_theta) ;
  pttPtr phi_to_theta_fun = *xpttfun ;
  // Unwrap pointer to log_j function.
  typedef double (*logjacPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<logjacPtr> xlogjfun(log_j) ;
  logjacPtr log_j_fun = *xlogjfun ;
  Rcpp::NumericVector lambda = tpars["lambda"] ;
  Rcpp::NumericVector gm = tpars["gm"] ;
  Rcpp::NumericVector con = tpars["con"] ;
  Rcpp::IntegerVector which_lam = tpars["which_lam"] ;
  Rcpp::NumericVector theta, phi, psi, phi2, temp, temp2 ;
  double val, log_bc_jac, logj ;
  psi = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
  temp = psi * con + 1.0 ;
  temp = temp[which_lam] ;
  if (any_nonpos(temp)) {
    return R_NegInf ;
  }
  phi = psi_to_phi_fun(psi, lambda, gm, con) ;
  theta = phi_to_theta_fun(phi, user_args) ;
  if (any_naC(theta)) {
    return R_NegInf ;
  }
  logj = log_j_fun(theta, user_args) ;
  phi2 = phi[which_lam] ;
  temp = Rcpp::log(phi2) ;
  temp2 = lambda[which_lam] ;
  log_bc_jac = sum((temp2 - 1.0) * temp) ;
  val = fun(theta, pars) - log_bc_jac - logj - hscale ;
  return val ;
}

// Case 4: User-supplied transformation, rotation and relocation.

// [[Rcpp::export]]
double cpp_logf_rho_4(const arma::vec& rho, const arma::vec& psi_mode,
                      const arma::mat& rot_mat, const double& hscale,
                      const SEXP& logf, const Rcpp::List& pars,
                      const Rcpp::List& tpars, const SEXP& ptpfun,
                      const SEXP& phi_to_theta, const SEXP& log_j,
                      const Rcpp::List& user_args) {
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  // Unwrap pointer to phi_to_theta transformation function.
  typedef Rcpp::NumericVector (*pttPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  Rcpp::XPtr<pttPtr> xpttfun(phi_to_theta) ;
  pttPtr phi_to_theta_fun = *xpttfun ;
  // Unwrap pointer to log_j function.
  typedef double (*logjacPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<logjacPtr> xlogjfun(log_j) ;
  logjacPtr log_j_fun = *xlogjfun ;
  Rcpp::NumericVector theta, phi ;
  double val, logj ;
  phi = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
  theta = phi_to_theta_fun(phi, user_args) ;
  if (any_naC(theta)) {
    return R_NegInf ;
  }
  logj = log_j_fun(theta, user_args) ;
  val = fun(theta, pars) - logj - hscale ;
  return val ;
}

// Function to vectorize cpp_logf_rho_4() for use in find_lambda_rcpp()
// and find_lambda_one_d_rcpp().

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_apply(const Rcpp::NumericMatrix& x,
                               const arma::vec& psi_mode,
                               const arma::mat& rot_mat, const double& hscale,
                               const SEXP& logf, const Rcpp::List& pars,
                               const Rcpp::List& tpars, const SEXP& ptpfun,
                               const SEXP& phi_to_theta, const SEXP& log_j,
                               const Rcpp::List& user_args) {
  int nRows = x.nrow() ;
  Rcpp::NumericVector out = no_init(nRows) ;
  for(int i=0; i < nRows; i++) {
    arma::vec rho = x(i, _) ;
    out[i] = cpp_logf_rho_4(rho, psi_mode, rot_mat, hscale, logf, pars, tpars,
                            ptpfun, phi_to_theta, log_j, user_args) ;
  }
  return out ;
}

// [[Rcpp::export]]
double cpp_a_obj(const arma::vec& psi, const arma::vec& psi_mode,
                 const arma::mat& rot_mat, const double& hscale,
                 const SEXP& logf, const int& d, const double& r,
                 const double& big_val, const Rcpp::List& pars) {
  if (any_naC(Rcpp::as<Rcpp::NumericVector>(wrap(psi)))) {
    return big_val ;
  }
  double val ;
  val = cpp_logf_rho(psi, psi_mode, rot_mat, hscale, logf, pars) ;
  if (val == R_NegInf) {
    return big_val ;
  }
  return -val / (d * r + 1) ;
}

// [[Rcpp::export]]
double cpp_a_obj_2(const arma::vec& psi, const arma::vec& psi_mode,
                   const arma::mat& rot_mat, const double& hscale,
                   const int& d, const double& r, const double& big_val,
                   const SEXP& tfun, const Rcpp::List& tpars,
                   const SEXP& logf, const Rcpp::List& pars,
                   const SEXP& ptpfun,
                   const SEXP& phi_to_theta, const SEXP& log_j,
                   const Rcpp::List& user_args) {
  if (any_naC(Rcpp::as<Rcpp::NumericVector>(wrap(psi)))) {
    return big_val ;
  }
  // Unwrap pointer to transformation function.
  typedef double (*transPtr)(const arma::vec& rho, const arma::vec& psi_mode,
                  const arma::mat& rot_mat, const double& hscale,
                  const SEXP& logf, const Rcpp::List& pars,
                  const Rcpp::List& tpars, const SEXP& ptpfun,
                  const SEXP& phi_to_theta, const SEXP& log_j,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<transPtr> xtfun(tfun) ;
  transPtr transfun = *xtfun ;
  double val ;
  val = transfun(psi, psi_mode, rot_mat, hscale, logf, pars, tpars, ptpfun,
                 phi_to_theta, log_j, user_args) ;
  if (val == R_NegInf) {
    return big_val ;
  }
  return -val / (d * r + 1) ;
}

// Why did I have if (check == 0L | is.infinite(check)) { in lower and upper box?
// i.e. why the is.infinite bit?
// ... because the density could be infinite

// [[Rcpp::export]]
double cpp_lower_box(const arma::vec& rho, int j, const arma::vec& psi_mode,
                     const arma::mat& rot_mat, const double& hscale,
                     const SEXP& logf, const int& d, const double& r,
                     const double& big_val, const Rcpp::List& pars) {
  if (rho(j) > 0) {
    return(big_val) ;
  }
  if (any_naC(Rcpp::as<Rcpp::NumericVector>(wrap(rho)))) {
    return(big_val) ;
  }
  double val ;
  val = cpp_logf_rho(rho, psi_mode, rot_mat, hscale, logf, pars) ;
  if (val == R_NegInf) {
    return(big_val) ;
  }
  return rho(j) * pow(exp(val), (r / (d * r + 1))) ;
}

// [[Rcpp::export]]
double cpp_lower_box_2(const arma::vec& rho, int j, const arma::vec& psi_mode,
                       const arma::mat& rot_mat, const double& hscale,
                       const SEXP& tfun, const Rcpp::List& tpars,
                       const SEXP& logf, const Rcpp::List& pars,
                       const int& d, const double& r,
                       const double& big_val, const SEXP& ptpfun,
                       const SEXP& phi_to_theta, const SEXP& log_j,
                       const Rcpp::List& user_args) {
  if (rho(j) > 0) {
    return(big_val) ;
  }
  if (any_naC(Rcpp::as<Rcpp::NumericVector>(wrap(rho)))) {
    return(big_val) ;
  }
  // Unwrap pointer to transformation function.
  typedef double (*transPtr)(const arma::vec& rho, const arma::vec& psi_mode,
                  const arma::mat& rot_mat, const double& hscale,
                  const SEXP& logf, const Rcpp::List& pars,
                  const Rcpp::List& tpars, const SEXP& ptpfun,
                  const SEXP& phi_to_theta, const SEXP& log_j,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<transPtr> xtfun(tfun) ;
  transPtr transfun = *xtfun ;
  double val ;
  val = transfun(rho, psi_mode, rot_mat, hscale, logf, pars, tpars, ptpfun,
                 phi_to_theta, log_j, user_args) ;
  if (val == R_NegInf) {
    return(big_val) ;
  }
  return rho(j) * pow(exp(val), (r / (d * r + 1))) ;
}

// [[Rcpp::export]]
double cpp_upper_box(const arma::vec& rho, int j, const arma::vec& psi_mode,
                     const arma::mat& rot_mat, const double& hscale,
                     const SEXP& logf, const int& d, const double& r,
                     const double& big_val, const Rcpp::List& pars) {
  if (rho(j) < 0) {
    return big_val ;
  }
  if (any_naC(Rcpp::as<Rcpp::NumericVector>(wrap(rho)))) {
    return big_val ;
  }
  double val ;
  val = cpp_logf_rho(rho, psi_mode, rot_mat, hscale, logf, pars) ;
  if (val == R_NegInf) {
    return big_val ;
  }
  return -rho(j) * pow(exp(val), (r / (d * r + 1))) ;
}

// [[Rcpp::export]]
double cpp_upper_box_2(const arma::vec& rho, int j, const arma::vec& psi_mode,
                       const arma::mat& rot_mat, const double& hscale,
                       const SEXP& tfun, const Rcpp::List& tpars,
                       const SEXP& logf, const Rcpp::List& pars,
                       const int& d, const double& r,
                       const double& big_val, const SEXP& ptpfun,
                       const SEXP& phi_to_theta, const SEXP& log_j,
                       const Rcpp::List& user_args) {
  if (rho(j) < 0) {
    return(big_val) ;
  }
  if (any_naC(Rcpp::as<Rcpp::NumericVector>(wrap(rho)))) {
    return(big_val) ;
  }
  // Unwrap pointer to transformation function.
  typedef double (*transPtr)(const arma::vec& rho, const arma::vec& psi_mode,
                  const arma::mat& rot_mat, const double& hscale,
                  const SEXP& logf, const Rcpp::List& pars,
                  const Rcpp::List& tpars, const SEXP& ptpfun,
                  const SEXP& phi_to_theta, const SEXP& log_j,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<transPtr> xtfun(tfun) ;
  transPtr transfun = *xtfun ;
  double val ;
  val = transfun(rho, psi_mode, rot_mat, hscale, logf, pars, tpars, ptpfun,
                 phi_to_theta, log_j, user_args) ;
  if (val == R_NegInf) {
    return(big_val) ;
  }
  return -rho(j) * pow(exp(val), (r / (d * r + 1))) ;
}

// [[Rcpp::export]]
Rcpp::List ru_cpp(const int& n, const int& d, const double& r,
                  const double& a_box, const Rcpp::NumericVector& l_box,
                  const Rcpp::NumericVector& u_box, const SEXP& logf,
                  const arma::vec& psi_mode, const arma::mat& rot_mat,
                  const double& hscale, const Rcpp::List& pars) {
  RNGScope scope; // ensure RNG gets set/reset
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  Rcpp::NumericMatrix sim_vals_rho(n, d), sim_vals(n, d) ;
  int ntry = 0, nacc = 0 ;
  double u, rhs, d_r;
  Rcpp::NumericVector d_box, vs, rho, theta ;
  d_r = d * r + 1 ;
  d_box = u_box - l_box ;
  while (nacc < n) {
    if (ntry % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
    u = runif(1, 0, a_box)[0] ;
    vs = d_box * Rcpp::runif(d) + l_box ;
    rho = vs / pow(u, r) ;
    theta = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
    rhs = fun(theta, pars) - hscale ;
    ntry++ ;
    if (d_r * log(u) < rhs) {
      sim_vals_rho(nacc,_) = rho ;
      sim_vals(nacc,_) = theta ;
      nacc++ ;
    }
  }
  return List::create(Named("sim_vals_rho") = sim_vals_rho,
                      Named("sim_vals") = sim_vals,
                      Named("ntry") = ntry) ;
}

// [[Rcpp::export]]
Rcpp::List ru_cpp_2(const int& n, const int& d, const double& r,
                    const double& a_box, const Rcpp::NumericVector& l_box,
                    const Rcpp::NumericVector& u_box,
                    const arma::vec& psi_mode, const arma::mat& rot_mat,
                    const double& hscale, const SEXP& logf,
                    const Rcpp::List& pars, const Rcpp::List& tpars,
                    const SEXP& ptpfun, const SEXP& phi_to_theta,
                    const SEXP& log_j, const Rcpp::List& user_args) {
  RNGScope scope; // ensure RNG gets set/reset
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  // Unwrap pointer to psi_to_phi transformation function.
  typedef Rcpp::NumericVector (*ptpPtr)(const Rcpp::NumericVector& psi,
                               const Rcpp::NumericVector& lambda,
                               const Rcpp::NumericVector& gm,
                               const Rcpp::NumericVector& con) ;
  Rcpp::XPtr<ptpPtr> xptpfun(ptpfun) ;
  ptpPtr psi_to_phi_fun = *xptpfun ;
  Rcpp::NumericVector lambda = tpars["lambda"] ;
  Rcpp::NumericVector gm = tpars["gm"] ;
  Rcpp::NumericVector con = tpars["con"] ;
  Rcpp::IntegerVector which_lam = tpars["which_lam"] ;
  Rcpp::NumericMatrix sim_vals_rho(n, d), sim_vals(n, d) ;
  int ntry = 0, nacc = 0 ;
  double u, rhs, d_r, log_bc_jac;
  Rcpp::NumericVector d_box, vs, rho, phi, psi, temp, temp2, phi2 ;
  d_r = d * r + 1 ;
  d_box = u_box - l_box ;
  while (nacc < n) {
    if (ntry % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
    u = runif(1, 0, a_box)[0] ;
    vs = d_box * Rcpp::runif(d) + l_box ;
    rho = vs / pow(u, r) ;
    psi = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
    temp = psi * con + 1.0 ;
    temp = temp[which_lam] ;
    ntry++ ;
    if (all_pos(temp)) {
      phi = psi_to_phi_fun(psi, lambda, gm, con) ;
      phi2 = phi[which_lam] ;
      temp = Rcpp::log(phi2) ;
      temp2 = lambda[which_lam] ;
      log_bc_jac = sum((temp2 - 1.0) * temp) ;
      rhs = fun(phi, pars) - log_bc_jac - hscale ;
      if (d_r * log(u) < rhs) {
        sim_vals_rho(nacc,_) = rho ;
        sim_vals(nacc,_) = phi ;
        nacc++ ;
      }
    }
  }
  return List::create(Named("sim_vals_rho") = sim_vals_rho,
                      Named("sim_vals") = sim_vals,
                      Named("ntry") = ntry) ;
}

// [[Rcpp::export]]
Rcpp::List ru_cpp_3(const int& n, const int& d, const double& r,
                    const double& a_box, const Rcpp::NumericVector& l_box,
                    const Rcpp::NumericVector& u_box,
                    const arma::vec& psi_mode, const arma::mat& rot_mat,
                    const double& hscale, const SEXP& logf,
                    const Rcpp::List& pars, const Rcpp::List& tpars,
                    const SEXP& ptpfun, const SEXP& phi_to_theta,
                    const SEXP& log_j, const Rcpp::List& user_args) {
  RNGScope scope; // ensure RNG gets set/reset
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  // Unwrap pointer to psi_to_phi transformation function.
  typedef Rcpp::NumericVector (*ptpPtr)(const Rcpp::NumericVector& psi,
                               const Rcpp::NumericVector& lambda,
                               const Rcpp::NumericVector& gm,
                               const Rcpp::NumericVector& con) ;
  Rcpp::XPtr<ptpPtr> xptpfun(ptpfun) ;
  ptpPtr psi_to_phi_fun = *xptpfun ;
  // Unwrap pointer to phi_to_theta transformation function.
  typedef Rcpp::NumericVector (*pttPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  Rcpp::XPtr<pttPtr> xpttfun(phi_to_theta) ;
  pttPtr phi_to_theta_fun = *xpttfun ;
  // Unwrap pointer to log_j function.
  typedef double (*logjacPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<logjacPtr> xlogjfun(log_j) ;
  logjacPtr log_j_fun = *xlogjfun ;
  Rcpp::NumericVector lambda = tpars["lambda"] ;
  Rcpp::NumericVector gm = tpars["gm"] ;
  Rcpp::NumericVector con = tpars["con"] ;
  Rcpp::IntegerVector which_lam = tpars["which_lam"] ;
  Rcpp::NumericMatrix sim_vals_rho(n, d), sim_vals(n, d) ;
  int ntry = 0, nacc = 0 ;
  double u, rhs, d_r, log_bc_jac, logj;
  Rcpp::NumericVector d_box, vs, rho, theta, phi, psi, temp, temp2, phi2 ;
  d_r = d * r + 1 ;
  d_box = u_box - l_box ;
  while (nacc < n) {
    if (ntry % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
    u = runif(1, 0, a_box)[0] ;
    vs = d_box * Rcpp::runif(d) + l_box ;
    rho = vs / pow(u, r) ;
    psi = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
    temp = psi * con + 1.0 ;
    temp = temp[which_lam] ;
    ntry++ ;
    if (all_pos(temp)) {
      phi = psi_to_phi_fun(psi, lambda, gm, con) ;
      theta = phi_to_theta_fun(phi, user_args) ;
      if (no_naC(theta)) {
        logj = log_j_fun(theta, user_args) ;
        phi2 = phi[which_lam] ;
        temp = Rcpp::log(phi2) ;
        temp2 = lambda[which_lam] ;
        log_bc_jac = sum((temp2 - 1.0) * temp) ;
        rhs = fun(theta, pars) - log_bc_jac - logj - hscale ;
        if (d_r * log(u) < rhs) {
          sim_vals_rho(nacc,_) = rho ;
          sim_vals(nacc,_) = theta ;
          nacc++ ;
        }
      }
    }
  }
  return List::create(Named("sim_vals_rho") = sim_vals_rho,
                      Named("sim_vals") = sim_vals,
                      Named("ntry") = ntry) ;
}

// [[Rcpp::export]]
Rcpp::List ru_cpp_4(const int& n, const int& d, const double& r,
                    const double& a_box, const Rcpp::NumericVector& l_box,
                    const Rcpp::NumericVector& u_box,
                    const arma::vec& psi_mode, const arma::mat& rot_mat,
                    const double& hscale, const SEXP& logf,
                    const Rcpp::List& pars, const Rcpp::List& tpars,
                    const SEXP& ptpfun, const SEXP& phi_to_theta,
                    const SEXP& log_j, const Rcpp::List& user_args) {
  RNGScope scope; // ensure RNG gets set/reset
  // Unwrap pointer to untransformed target log-density.
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  Rcpp::XPtr<funcPtr> xpfun(logf) ;
  funcPtr fun = *xpfun ;
  // Unwrap pointer to phi_to_theta transformation function.
  typedef Rcpp::NumericVector (*pttPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  Rcpp::XPtr<pttPtr> xpttfun(phi_to_theta) ;
  pttPtr phi_to_theta_fun = *xpttfun ;
  // Unwrap pointer to log_j function.
  typedef double (*logjacPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  Rcpp::XPtr<logjacPtr> xlogjfun(log_j) ;
  logjacPtr log_j_fun = *xlogjfun ;
  Rcpp::NumericMatrix sim_vals_rho(n, d), sim_vals(n, d) ;
  int ntry = 0, nacc = 0 ;
  double u, rhs, d_r, logj;
  Rcpp::NumericVector d_box, vs, rho, theta, phi, psi, temp, temp2, phi2 ;
  d_r = d * r + 1 ;
  d_box = u_box - l_box ;
  while (nacc < n) {
    if (ntry % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
    u = runif(1, 0, a_box)[0] ;
    vs = d_box * Rcpp::runif(d) + l_box ;
    rho = vs / pow(u, r) ;
    phi = cpp_rho_to_psi(rho, psi_mode, rot_mat) ;
    ntry++ ;
    theta = phi_to_theta_fun(phi, user_args) ;
    if (no_naC(theta)) {
      logj = log_j_fun(theta, user_args) ;
      rhs = fun(theta, pars) - logj - hscale ;
      if (d_r * log(u) < rhs) {
        sim_vals_rho(nacc,_) = rho ;
        sim_vals(nacc,_) = theta ;
        nacc++ ;
      }
    }
  }
  return List::create(Named("sim_vals_rho") = sim_vals_rho,
                      Named("sim_vals") = sim_vals,
                      Named("ntry") = ntry) ;
}

// [[Rcpp::export]]
SEXP create_trans_xptr(std::string fstr) {
  typedef double (*transPtr)(const arma::vec& rho, const arma::vec& psi_mode,
                  const arma::mat& rot_mat, const double& hscale,
                  const SEXP& logf, const Rcpp::List& pars,
                  const Rcpp::List& tpars, const SEXP& ptpfun,
                  const SEXP& phi_to_theta, const SEXP& log_j,
                  const Rcpp::List& user_args) ;
  if (fstr == "case_2")
    return(Rcpp::XPtr<transPtr>(new transPtr(&cpp_logf_rho_2))) ;
  else if (fstr == "case_3")
    return(Rcpp::XPtr<transPtr>(new transPtr(&cpp_logf_rho_3))) ;
  else if (fstr == "case_4")
    return(Rcpp::XPtr<transPtr>(new transPtr(&cpp_logf_rho_4))) ;
  else
    return(Rcpp::XPtr<transPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
Rcpp::NumericVector bc_no_trans(const Rcpp::NumericVector& psi,
                                const Rcpp::NumericVector& lambda,
                                const Rcpp::NumericVector& gm,
                                const Rcpp::NumericVector& con) {
  return psi ;
}

// [[Rcpp::export]]
SEXP create_psi_to_phi_xptr(std::string fstr) {
  typedef Rcpp::NumericVector (*ptpPtr)(const Rcpp::NumericVector& psi,
                               const Rcpp::NumericVector& lambda,
                               const Rcpp::NumericVector& gm,
                               const Rcpp::NumericVector& con) ;
  if (fstr == "no_zero")
    return(Rcpp::XPtr<ptpPtr>(new ptpPtr(&cpp_psi_to_phi))) ;
  else if (fstr == "has_zero")
    return(Rcpp::XPtr<ptpPtr>(new ptpPtr(&cpp_psi_to_phi_0))) ;
  else if (fstr == "no_trans")
    return(Rcpp::XPtr<ptpPtr>(new ptpPtr(&bc_no_trans))) ;
  else
    return(Rcpp::XPtr<ptpPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
double log_none_jac(const Rcpp::NumericVector& theta,
                    const Rcpp::List& user_args) {
  return 0.0 ;
}

// [[Rcpp::export]]
SEXP create_log_jac_xptr(std::string fstr) {
  typedef double (*logjacPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  return(Rcpp::XPtr<logjacPtr>(new logjacPtr(&log_none_jac))) ;
}

// [[Rcpp::export]]
Rcpp::NumericVector no_trans(const Rcpp::NumericVector& theta,
                             const Rcpp::List& user_args) {
  return theta ;
}

// [[Rcpp::export]]
SEXP null_phi_to_theta_xptr(std::string fstr) {
  typedef Rcpp::NumericVector (*p2tPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&no_trans))) ;
}
