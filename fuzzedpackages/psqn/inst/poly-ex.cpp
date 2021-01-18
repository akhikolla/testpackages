// see `mlogit-ex.cpp` for an example with more comments

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

class poly_func final : public PSQN::element_function {
  /// global parameter indices
  arma::uvec const g_idx;
  /// centroid vector
  arma::vec const mu_cluster;
  /// matrix used to transform subset of global parameters
  arma::mat const Psi;
  /// number of global parameters
  size_t const n_global;
  /// global parameter centroid vector
  arma::vec const mu_global;
  /**
   true if this element function should compute the terms from the global
   paramaters */
  bool const comp_global;

public:
  poly_func(List data, arma::vec const &mu_g, bool const comp_global):
  g_idx     (as<arma::uvec>(data["g_idx"    ]) - 1L),
  mu_cluster(as<arma::vec>(data["mu_cluster"])     ),
  Psi       (as<arma::mat>(data["Psi"       ])     ),
  n_global(mu_g.n_elem),
  mu_global(comp_global ? mu_g : arma::vec() ),
  comp_global(comp_global)
  { }

  size_t global_dim() const {
    return n_global;
  }
  size_t private_dim() const {
    return mu_cluster.n_elem;
  }

  double func(double const *point) const {
    arma::vec const x_glob = vec_no_cp(point           , n_global),
                    x_priv = vec_no_cp(point + n_global, mu_cluster.n_elem),
                     delta = x_priv - Psi * x_glob(g_idx) - mu_cluster;

    // compute the function
    double out(0);
    out += arma::dot(delta, delta);

    if(comp_global)
      out += arma::dot(x_glob - mu_global, x_glob - mu_global);

    return out;
  }

  double grad
  (double const * point, double * gr) const {
    arma::vec const x_glob = vec_no_cp(point           , n_global),
                    x_priv = vec_no_cp(point + n_global, mu_cluster.n_elem),
                     delta = x_priv - Psi * x_glob(g_idx) - mu_cluster;

    // create objects to write to for the gradient
    std::fill(gr, gr + x_glob.n_elem, 0.);
    arma::vec d_glob(gr                , x_glob.n_elem, false),
              d_priv(gr + x_glob.n_elem, x_priv.n_elem, false);

    // compute the function and the gradient
    double out(0);
    out += arma::dot(delta, delta);
    d_glob(g_idx) -= 2 * Psi.t() * delta;
    d_priv         = 2 * delta;

    if(comp_global){
      out += arma::dot(x_glob - mu_global, x_glob - mu_global);
      d_glob += 2. * x_glob;
      d_glob -= 2 * mu_global;
    }

    return out;
  }

  bool thread_safe() const {
    return true;
  }
};

using poly_optim = PSQN::optimizer<poly_func, PSQN::R_reporter,
                                   PSQN::R_interrupter>;

// [[Rcpp::export]]
SEXP get_poly_optimizer(List data, arma::vec const &mu_global,
                        unsigned const max_threads){
  size_t const n_elem_funcs = data.size();
  std::vector<poly_func> funcs;
  funcs.reserve(n_elem_funcs);
  bool comp_global(true);
  for(auto dat : data){
    funcs.emplace_back(List(dat), mu_global, comp_global);
    comp_global = false;
  }

  // create an XPtr to the object we will need
  XPtr<poly_optim>ptr(new poly_optim(funcs, max_threads));

  // return the pointer to be used later
  return ptr;
}

// [[Rcpp::export]]
List optim_poly
  (NumericVector val, SEXP ptr, double const rel_eps, unsigned const max_it,
   unsigned const n_threads, double const c1,
   double const c2, bool const use_bfgs = true, int const trace = 0L,
   double const cg_tol = .5, bool const strong_wolfe = true,
   size_t const max_cg = 0L, int const pre_method = 1L){
  XPtr<poly_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("optim_poly: invalid parameter size");

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
      _["convergence"] =  res.info == PSQN::info_code::converged);
}

// [[Rcpp::export]]
double eval_poly(NumericVector val, SEXP ptr, unsigned const n_threads){
  XPtr<poly_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("eval_poly: invalid parameter size");

  optim->set_n_threads(n_threads);
  return optim->eval(&val[0], nullptr, false);
}

// [[Rcpp::export]]
NumericVector grad_poly(NumericVector val, SEXP ptr,
                        unsigned const n_threads){
  XPtr<poly_optim> optim(ptr);

  // check that we pass a parameter value of the right length
  if(optim->n_par != static_cast<size_t>(val.size()))
    throw std::invalid_argument("grad_poly: invalid parameter size");

  NumericVector grad(val.size());
  optim->set_n_threads(n_threads);
  grad.attr("value") = optim->eval(&val[0], &grad[0], true);

  return grad;
}
