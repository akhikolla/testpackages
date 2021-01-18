#include "psqn.h"
#include "psqn-reporter.h"
#include <stdexcept>

using namespace Rcpp;

/**
 simple wrapper for an R function which takes three argument.

 Caution: it is much faster but is not gaurded against errors:
   https://stackoverflow.com/a/37846827/5861244
 */
class simple_R_func3 {
  SEXP fn, env;

public:
  simple_R_func3(SEXP fn, SEXP env):
  fn(fn), env(env) { }

  SEXP operator()(SEXP a1, SEXP a2, SEXP a3) const {
    SEXP R_fcall, out;
    PROTECT(R_fcall = Rf_lang4(fn, a1, a2, a3));
    PROTECT(out = Rf_eval(R_fcall, env));
    UNPROTECT(2);
    return out;
  }
};

/** same as simple_R_func3 but with only one argument */
class simple_R_func1 {
  SEXP fn, env;

public:
  simple_R_func1(SEXP fn, SEXP env):
  fn(fn), env(env) { }

  SEXP operator()(SEXP a1) const {
    SEXP R_fcall, out;
    PROTECT(R_fcall = Rf_lang2(fn, a1));
    PROTECT(out = Rf_eval(R_fcall, env));
    UNPROTECT(2);
    return out;
  }
};

class r_worker {
  simple_R_func3 f;
  IntegerVector f_idx;
  LogicalVector mutable scomp_grad = LogicalVector(1L);
  size_t const g_dim, p_dim,
               n_ele = g_dim + p_dim;

  NumericVector mutable par = NumericVector(g_dim + p_dim);

public:
  r_worker(SEXP func, int iarg, SEXP rho):
  f(func, rho),
  f_idx(([&]() -> IntegerVector {
    IntegerVector out(1L);
    out[0] = iarg + 1L;
    return out;
  })()),
  g_dim(([&]() -> size_t {
    NumericVector dum(0);
    scomp_grad[0] = false;
    SEXP res = PROTECT(f(f_idx, dum, scomp_grad));
    if(!Rf_isInteger(res) or !Rf_isVector(res) or Rf_xlength(res) != 2L){
      UNPROTECT(1);
      throw std::invalid_argument(
          "fn returns invalid lengths with zero length par");
    }
    int const out = *INTEGER(res);
    UNPROTECT(1);
    return out;
  })()),
  p_dim(([&]() -> size_t {
    NumericVector dum(0);
    scomp_grad[0] = false;
    SEXP res =  PROTECT(f(f_idx, dum, scomp_grad));
    if(!Rf_isInteger(res) or !Rf_isVector(res) or Rf_xlength(res) != 2L){
      UNPROTECT(1);
      throw std::invalid_argument(
          "fn returns invalid lengths with zero length par");
    }
    int const out = INTEGER(res)[1L];
    UNPROTECT(1);
    return out;
  })())
  { };

  size_t global_dim() const {
    return g_dim;
  };
  size_t private_dim() const {
    return p_dim;
  }

  double func(double const *point) const {
    double *p = &par[0];
    for(size_t j = 0; j < n_ele; ++j, ++point, ++p)
      *p = *point;
    scomp_grad[0] = false;
    SEXP res =  PROTECT(f(f_idx, par, scomp_grad));
    if(!Rf_isReal(res) or !Rf_isVector(res) or Rf_xlength(res) != 1L){
      UNPROTECT(1);
      throw std::invalid_argument(
          "fn returns invalid output with comp_grad = FALSE");
    }
    double const out = *REAL(res);
    UNPROTECT(1);
    return out;
  }

  double grad
    (double const * __restrict__ point, double * __restrict__ gr) const {
    for(size_t j = 0; j < n_ele; ++j, ++point)
      par[j] = *point;
    scomp_grad[0] = true;
    SEXP res =  PROTECT(f(f_idx, par, scomp_grad));
    CharacterVector what("grad");
    SEXP gr_val = PROTECT(Rf_getAttrib(res, what));

    if(!Rf_isReal(res) or !Rf_isVector(res) or Rf_xlength(res) != 1L or
         Rf_isNull(gr_val) or !Rf_isReal(gr_val) or
         static_cast<size_t>(Rf_xlength(gr_val)) != n_ele){
      UNPROTECT(2);
      throw std::invalid_argument(
          "fn returns invalid output with comp_grad = TRUE");
    }

    lp::copy(gr, REAL(gr_val), n_ele);
    double const out = *REAL(res);
    UNPROTECT(2);
    return out;
  };

  virtual bool thread_safe() const {
    return false;
  };
};

List wrap_optim_info(NumericVector par_res, PSQN::optim_info res){
  NumericVector counts = NumericVector::create(
    res.n_eval, res.n_grad,  res.n_cg);
  counts.names() = CharacterVector::create("function", "gradient", "n_cg");

  int const info = static_cast<int>(res.info);
  return List::create(
    _["par"] = par_res, _["value"] = res.value, _["info"] = info,
      _["counts"] = counts, _["convergence"] = info >= 0L);
}

//' Partially Separable Function Optimization
//'
//' @description
//' Optimization method for specially structured partially separable
//' functions.
//'
//' @param par Initial values for the parameters. It is a concatenated
//' vector of the global parameters and all the private parameters.
//' @param fn Function to compute the element functions and their
//' derivatives. Each call computes on element function. See the examples
//' section.
//' @param n_ele_func Number of element functions.
//' @param rel_eps Relative convergence threshold.
//' @param n_threads Number of threads to use.
//' @param max_it Maximum number of iterations.
//' @param c1,c2 Thresholds for the Wolfe condition.
//' @param use_bfgs Logical for whether to use BFGS updates or SR1 updates.
//' @param trace Integer where larger values gives more information during the
//' optimization.
//' @param cg_tol Threshold for the conjugate gradient method.
//' @param strong_wolfe \code{TRUE} if the strong Wolfe condition should be used.
//' @param env Environment to evaluate \code{fn} in. \code{NULL} yields the
//' global environment.
//' @param max_cg maximum number of conjugate gradient iterations in each
//' iteration. Use zero if there should not be a limit.
//' @param pre_method preconditioning method in conjugate gradient method.
//' zero yields no preconditioning and one yields diagonal preconditioning.
//'
//' @details
//' The function follows the method described by Nocedal and Wright (2006)
//' and mainly what is described in Section 7.4. Details are provided
//' in the psqn vignette. See \code{vignette("psqn", package = "psqn")}.
//'
//' The partially separable function we consider are special in that the
//' function to be minimized is a sum of so-called element functions which
//' only depend on few shared (global) parameters and some
//' private parameters which are particular to each element function.
//'
//' The optimization function is also available in C++ as a header-only
//' library. Using C++ may reduce the computation time substantially. See
//' the vignette in the package for examples.
//'
//' @return
//' An object with the following elements:
//' \item{par}{the estimated global and private parameters.}
//' \item{value}{function value at \code{par}.}
//' \item{info}{information code. 0 implies convergence.
//' -1 implies that the maximum number iterations is reached.
//' -2 implies that the conjugate gradient method failed.
//' -3 implies that the line search failed.
//' -4 implies that the user interrupted the optimization.}
//' \item{counts}{An integer vector with the number of function evaluations,
//' gradient evaluations, and the number of conjugate gradient iterations.}
//' \item{convergence}{\code{TRUE} if \code{info == 0}.}
//'
//' @references
//' Nocedal, J. and Wright, S. J. (2006). \emph{Numerical Optimization}
//' (2nd ed.). Springer.
//'
//' @examples
//' # example with inner problem in a Taylor approximation for a mixed GLMM as
//' # in the vignette
//'
//' # assign model parameters, number of random effects, and fixed effects
//' q <- 2 # number of private parameters per cluster
//' p <- 1 # number of global parameters
//' beta <- sqrt((1:p) / sum(1:p))
//' Sigma <- diag(q)
//'
//' # simulate a data set
//' set.seed(66608927)
//' n_clusters <- 20L # number of clusters
//' sim_dat <- replicate(n_clusters, {
//'   n_members <- sample.int(8L, 1L) + 2L
//'   X <- matrix(runif(p * n_members, -sqrt(6 / 2), sqrt(6 / 2)),
//'               p)
//'   u <- drop(rnorm(q) %*% chol(Sigma))
//'   Z <- matrix(runif(q * n_members, -sqrt(6 / 2 / q), sqrt(6 / 2 / q)),
//'               q)
//'   eta <- drop(beta %*% X + u %*% Z)
//'   y <- as.numeric((1 + exp(-eta))^(-1) > runif(n_members))
//'
//'   list(X = X, Z = Z, y = y, u = u, Sigma_inv = solve(Sigma))
//' }, simplify = FALSE)
//'
//' # evalutes the negative log integrand.
//' #
//' # Args:
//' #   i cluster/element function index.
//' #   par the global and private parameter for this cluster. It has length
//' #       zero if the number of parameters is requested. That is, a 2D integer
//' #       vector the number of global parameters as the first element and the
//' #       number of private parameters as the second element.
//' #   comp_grad logical for whether to compute the gradient.
//' r_func <- function(i, par, comp_grad){
//'   dat <- sim_dat[[i]]
//'   X <- dat$X
//'   Z <- dat$Z
//'
//'   if(length(par) < 1)
//'     # requested the dimension of the parameter
//'     return(c(global_dim = NROW(dat$X), private_dim = NROW(dat$Z)))
//'
//'   y <- dat$y
//'   Sigma_inv <- dat$Sigma_inv
//'
//'   beta <- par[1:p]
//'   uhat <- par[1:q + p]
//'   eta <- drop(beta %*% X + uhat %*% Z)
//'   exp_eta <- exp(eta)
//'
//'   out <- -sum(y * eta) + sum(log(1 + exp_eta)) +
//'     sum(uhat * (Sigma_inv %*% uhat)) / 2
//'   if(comp_grad){
//'     d_eta <- -y + exp_eta / (1 + exp_eta)
//'     grad <- c(X %*% d_eta,
//'               Z %*% d_eta + dat$Sigma_inv %*% uhat)
//'     attr(out, "grad") <- grad
//'   }
//'
//'   out
//' }
//'
//' # optimize the log integrand
//' res <- psqn(par = rep(0, p + q * n_clusters), fn = r_func,
//'             n_ele_func = n_clusters)
//' head(res$par, p)              # the estimated global parameters
//' tail(res$par, n_clusters * q) # the estimated private parameters
//'
//' # compare with
//' beta
//' c(sapply(sim_dat, "[[", "u"))
//' @export
// [[Rcpp::export]]
List psqn
  (NumericVector par, SEXP fn, unsigned const n_ele_func,
   double const rel_eps = .00000001,
   unsigned const max_it = 100L, unsigned const n_threads = 1L,
   double const c1 = .0001, double const c2 = .9,
   bool const use_bfgs = true, int const trace = 0L,
   double const cg_tol = .5, bool const strong_wolfe = true,
   SEXP env = R_NilValue, int const max_cg = 0L,
   int const pre_method = 1L){
  if(n_ele_func < 1L)
    throw std::invalid_argument("psqn: n_ele_func < 1L");

  if(Rf_isNull(env))
    env = Environment::global_env();
  if(!Rf_isEnvironment(env))
    throw std::invalid_argument("psqn: env is not an environment");
  if(!Rf_isFunction(fn))
    throw std::invalid_argument("psqn: fn is not a function");
  if(pre_method < 0L or pre_method > 1L)
    throw std::invalid_argument("psqn: invalid pre_method");

  std::vector<r_worker> funcs;
  funcs.reserve(n_ele_func);
  for(size_t i = 0; i < n_ele_func; ++i)
    funcs.emplace_back(fn, i, env);

  PSQN::optimizer<r_worker, PSQN::R_reporter,
                  PSQN::R_interrupter> optim(funcs, n_threads);

  // check that we pass a parameter value of the right length
  if(optim.n_par != static_cast<size_t>(par.size()))
    throw std::invalid_argument("psqn: invalid parameter size");

  NumericVector par_arg = clone(par);
  optim.set_n_threads(n_threads);
  auto res = optim.optim(&par_arg[0], rel_eps, max_it, c1, c2,
                         use_bfgs, trace, cg_tol, strong_wolfe, max_cg,
                         static_cast<PSQN::precondition>(pre_method));

  return wrap_optim_info(par_arg, res);
}


class r_worker_bfgs : public PSQN::problem {
  simple_R_func1 f, g;
  size_t const n_ele;
  NumericVector par = NumericVector(n_ele);

public:
  r_worker_bfgs(SEXP f, SEXP g, size_t const n_ele, SEXP env):
  f(f, env), g(g, env), n_ele(n_ele) { }

  size_t size() const {
    return n_ele;
  }

  double func(double const *val){
    lp::copy(&par[0], val, n_ele);
    SEXP res = PROTECT(f(par));
    if(!Rf_isReal(res) or !Rf_isVector(res) or Rf_xlength(res) != 1L){
      UNPROTECT(1);
      throw std::invalid_argument("fn returns invalid output");
    }
    double const out = *REAL(res);
    UNPROTECT(1);
    return out;
  }

  double grad(double const * __restrict__ val,
              double       * __restrict__ gr){
    lp::copy(&par[0], val, n_ele);

    SEXP res = PROTECT(g(par));
    CharacterVector what("value");
    SEXP func_val = PROTECT(Rf_getAttrib(res, what));

    if(!Rf_isReal(res) or !Rf_isVector(res) or
         static_cast<size_t>(Rf_xlength(res)) != n_ele or
         Rf_isNull(func_val) or !Rf_isReal(func_val) or
         Rf_xlength(func_val) != 1L){
      UNPROTECT(2);
      throw std::invalid_argument("gr returns invalid output");
    }

    double const out = *REAL(func_val);
    lp::copy(gr, REAL(res), n_ele);
    UNPROTECT(2);
    return out;
  }
};

//' BFGS Implementation Used Internally in the psqn Package
//'
//' @description
//' The method seems to mainly differ from \code{\link{optim}} by the line search
//' method. This version uses the interpolation method with a zoom phase
//' using cubic interpolation as described by Nocedal and Wright (2006).
//'
//' @references
//' Nocedal, J. and Wright, S. J. (2006). \emph{Numerical Optimization}
//' (2nd ed.). Springer.
//'
//' @return
//' An object like the object returned by \code{\link{psqn}}.
//'
//' @inheritParams psqn
//' @param par Initial values for the parameters.
//' @param fn Function to evaluate the function to be minimized.
//' @param gr Gradient of \code{fn}. Should return the function value as an
//' attribute called \code{"value"}.
//' @param env Environment to evaluate \code{fn} and \code{gr} in.
//' \code{NULL} yields the global environment.
//' @export
//'
//' @examples
//' # declare function and gradient from the example from help(optim)
//' fn <- function(x) {
//'   x1 <- x[1]
//'   x2 <- x[2]
//'   100 * (x2 - x1 * x1)^2 + (1 - x1)^2
//' }
//' gr <- function(x) {
//'   x1 <- x[1]
//'   x2 <- x[2]
//'   c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
//'      200 *      (x2 - x1 * x1))
//' }
//'
//' # we need a different function for the method in this package
//' gr_psqn <- function(x) {
//'   x1 <- x[1]
//'   x2 <- x[2]
//'   out <- c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
//'             200 *      (x2 - x1 * x1))
//'   attr(out, "value") <- 100 * (x2 - x1 * x1)^2 + (1 - x1)^2
//'   out
//' }
//'
//' # we get the same
//' optim    (c(-1.2, 1), fn, gr, method = "BFGS")
//' psqn_bfgs(c(-1.2, 1), fn, gr_psqn)
//'
//' # compare the computation time
//' system.time(replicate(1000,
//'                       optim    (c(-1.2, 1), fn, gr, method = "BFGS")))
//' system.time(replicate(1000,
//'                       psqn_bfgs(c(-1.2, 1), fn, gr_psqn)))
// [[Rcpp::export]]
List psqn_bfgs
  (NumericVector par, SEXP fn, SEXP gr,
   double const rel_eps = .00000001, size_t const max_it = 100,
   double const c1 = .0001, double const c2 = .9, int const trace = 0L,
   SEXP env = R_NilValue){
  if(Rf_isNull(env))
    env = Environment::global_env();
  if(!Rf_isEnvironment(env))
    throw std::invalid_argument("psqn_bfgs: env is not an environment");
  if(!Rf_isFunction(fn))
    throw std::invalid_argument("psqn_bfgs: fn is not a function");
  if(!Rf_isFunction(gr))
    throw std::invalid_argument("psqn_bfgs: gr is not a function");

  r_worker_bfgs problem(fn, gr, par.size(), env);

  NumericVector par_res = clone(par);
  auto const out = PSQN::bfgs<PSQN::R_reporter, PSQN::R_interrupter>
    (problem, &par_res[0], rel_eps, max_it, c1, c2, trace);

  return wrap_optim_info(par_res, out);
}
