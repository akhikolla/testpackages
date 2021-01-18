// we will use openMP to perform the comptutation in parallel
// [[Rcpp::plugins(openmp, cpp11)]]

// we use RcppArmadillo to simplify the code
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(psqn)]]
#include "psqn.h"
#include "psqn-reporter.h"
using namespace Rcpp;

/// simple function to avoid copying a vector. You can ignore this
inline arma::vec vec_no_cp(double const * x, size_t const n_ele){
  return arma::vec(const_cast<double *>(x), n_ele, false);
}

/***
 implements the element function for a given cluster. The class must provide
 the member functions which we provide here.

 We do not need to inherit from the element_function class but we can do it
 to ensure that we have implemented all the member functions.
 */
class m_logit_func final : public PSQN::element_function {
  /// design matrices
  arma::mat const X, Z;
  /// outcomes
  arma::vec const y;
  /// inverse covariance matrix
  arma::mat const Sigma_inv;

public:
  m_logit_func(List data):
  X        (as<arma::mat>(data["X"        ])),
  Z        (as<arma::mat>(data["Z"        ])),
  y        (as<arma::vec>(data["y"        ])),
  Sigma_inv(as<arma::mat>(data["Sigma_inv"])) { }

  /// dimension of the global parameters
  size_t global_dim() const {
    return X.n_rows;
  }
  /// dimension of the private parameters
  size_t private_dim() const {
    return Z.n_rows;
  }

  /***
   computes the element function.
   @param point point to compute function at.
   */
  double func(double const *point) const {
    arma::vec const beta = vec_no_cp(point           , X.n_rows),
                       u = vec_no_cp(point + X.n_rows, Z.n_rows);

    double out(0);
    for(size_t i = 0; i < y.n_elem; ++i){
      double const eta =
        arma::dot(beta, X.col(i)) + arma::dot(u, Z.col(i));
      out -= y[i] * eta - log(1 + exp(eta));
    }

    out += arma::dot(u, Sigma_inv * u) * .5;

    return out;
  }

  /***
   computes the element function and its gradient.
   @param point point to compute function at.
   @param gr gradient vector with respect to global and private parameters.
   */
  double grad
    (double const * point, double * gr) const {
    arma::vec const beta = vec_no_cp(point           , X.n_rows),
                       u = vec_no_cp(point + X.n_rows, Z.n_rows);

    // create objects to write to for the gradient
    std::fill(gr, gr + beta.n_elem + u.n_elem, 0.);
    arma::vec dbeta(gr              , beta.n_elem, false),
              du   (gr + beta.n_elem, u.n_elem   , false);

    double out(0);
    for(size_t i = 0; i < y.n_elem; ++i){
      arma::vec const xi = X.unsafe_col(i),
                      zi = Z.unsafe_col(i);
      double const eta = arma::dot(beta, xi) + arma::dot(u, zi),
               exp_eta = exp(eta),
               d_eta   = y[i] - exp_eta / (1 + exp_eta);
      out -= y[i] * eta - log(1 + exp_eta);
      dbeta -= d_eta * xi;
      du    -= d_eta * zi;
    }

    arma::vec u_scaled = Sigma_inv * u;
    out += arma::dot(u, u_scaled) * .5;
    du += u_scaled;

    return out;
  }

  /***
   returns true if the member functions are thread-safe.
   */
  bool thread_safe() const {
    return true;
  }
};

using mlogit_topim = PSQN::optimizer<m_logit_func, PSQN::R_reporter,
                                     PSQN::R_interrupter>;

/***
 creates a pointer to an object which is needed in the optim_mlogit
 function.
 @param data list with data for each element function.
 @param max_threads maximum number of threads to use.
 */
// [[Rcpp::export]]
SEXP get_mlogit_optimizer(List data, unsigned const max_threads){
  size_t const n_elem_funcs = data.size();
  std::vector<m_logit_func> funcs;
  funcs.reserve(n_elem_funcs);
  for(auto dat : data)
    funcs.emplace_back(List(dat));

  // create an XPtr to the object we will need
  XPtr<mlogit_topim> ptr(new mlogit_topim(funcs, max_threads));

  // return the pointer to be used later
  return ptr;
}

/***
 performs the optimization.
 @param val vector with starting value for the global and private
 parameters.
 @param ptr returned object from get_mlogit_optimizer.
 @param rel_eps relative convergence threshold.
 @param max_it maximum number iterations.
 @param n_threads number of threads to use.
 @param c1,c2 thresholds for Wolfe condition.
 @param use_bfgs boolean for whether to use SR1 or BFGS updates.
 @param trace integer where larger values gives more information during the
 optimization.
 @param cg_tol threshold for conjugate gradient method.
 @param strong_wolfe true if the strong Wolfe condition should be used.
 @param max_cg maximum number of conjugate gradient iterations in each
 iteration. Use zero if there should not be a limit.
 @param pre_method preconditioning method in conjugate gradient method.
 zero yields no preconditioning and one yields diagonal preconditioning.
 */
// [[Rcpp::export]]
List optim_mlogit
  (NumericVector val, SEXP ptr, double const rel_eps, unsigned const max_it,
   unsigned const n_threads, double const c1,
   double const c2, bool const use_bfgs = true, int const trace = 0L,
   double const cg_tol = .5, bool const strong_wolfe = true,
   size_t const max_cg = 0L, int const pre_method = 1L){
  XPtr<mlogit_topim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("optim_mlogit: invalid parameter size");

  NumericVector par = clone(val);
  optim->set_n_threads(n_threads);
  auto res = optim->optim(&par[0], rel_eps, max_it, c1, c2,
                          use_bfgs, trace, cg_tol, strong_wolfe, max_cg,
                          static_cast<PSQN::precondition>(pre_method));
  NumericVector counts = NumericVector::create(
    res.n_eval, res.n_grad,  res.n_cg);
  counts.names() = CharacterVector::create("function", "gradient", "n_cg");

  int const info = static_cast<int>(res.info);
  return List::create(
    _["par"] = par, _["value"] = res.value, _["info"] = info,
    _["counts"] = counts,
    _["convergence"] =  res.info == PSQN::info_code::converged );
}

/***
 performs the optimization but only for the private parameters.
 @param val vector with starting value for the global and private
 parameters.
 @param ptr returned object from get_mlogit_optimizer.
 @param rel_eps relative convergence threshold.
 @param max_it maximum number iterations.
 @param n_threads number of threads to use.
 @param c1,c2 thresholds for Wolfe condition.
 */
// [[Rcpp::export]]
NumericVector optim_mlogit_private
  (NumericVector val, SEXP ptr, double const rel_eps, unsigned const max_it,
   unsigned const n_threads, double const c1, double const c2){
  XPtr<mlogit_topim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("optim_mlogit_private: invalid parameter size");

  NumericVector par = clone(val);
  optim->set_n_threads(n_threads);
  double const res = optim->optim_priv(&par[0], rel_eps, max_it, c1, c2);
  par.attr("value") = res;
  return par;
}

/***
 evaluates the partially separable function.
 @param val vector with global and private parameters to evaluate the
 function at.
 @param ptr returned object from get_mlogit_optimizer.
 @param n_threads number of threads to use.
 */
// [[Rcpp::export]]
double eval_mlogit(NumericVector val, SEXP ptr, unsigned const n_threads){
  XPtr<mlogit_topim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("eval_mlogit: invalid parameter size");

  optim->set_n_threads(n_threads);
  return optim->eval(&val[0], nullptr, false);
}

/***
 evaluates the gradient of a partially separable function.
 @param val vector with global and private parameters to evaluate the
 function at.
 @param ptr returned object from get_mlogit_optimizer.
 @param n_threads number of threads to use.
 */
// [[Rcpp::export]]
NumericVector grad_mlogit(NumericVector val, SEXP ptr,
                          unsigned const n_threads){
  XPtr<mlogit_topim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("grad_mlogit: invalid parameter size");

  NumericVector grad(val.size());
  optim->set_n_threads(n_threads);
  grad.attr("value") = optim->eval(&val[0], &grad[0], true);

  return grad;
}

/***
 returns the current Hessian approximation.
 @param ptr returned object from get_mlogit_optimizer.
 */
// [[Rcpp::export]]
NumericMatrix get_Hess_approx_mlogit(SEXP ptr){
  XPtr<mlogit_topim> optim(ptr);

  NumericMatrix out(optim->n_par, optim->n_par);
  optim->get_hess(&out[0]);

  return out;
}
