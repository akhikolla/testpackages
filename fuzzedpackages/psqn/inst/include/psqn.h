#ifndef PSQN_H
#define PSQN_H
#include <vector>
#include <array>
#include <memory>
#include "lp.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include "constant.h"
#include <cmath>
#include "intrapolate.h"
#include "psqn-misc.h"
#include "psqn-bfgs.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace PSQN {
using std::abs;
using std::sqrt;

/***
 virtual base class which computes an element function and its gradient. The
 virtual class is mainly used as a check to ensure that all the member
 fucntions are implemented.
 */
class element_function {
public:
  /// dimension of the global parameters
  virtual size_t global_dim() const = 0;
  /// dimension of the private parameters
  virtual size_t private_dim() const = 0;

  /***
   computes the element function.
   @param point point to compute function at.
   */
  virtual double func(double const *point) const = 0;

  /***
   computes the element function and its gradient.
   @param point point to compute function at.
   @param gr gradient vector with respect to global and private parameters.
   */
  virtual double grad
    (double const * PSQN_RESTRICT point, double * PSQN_RESTRICT gr)
    const = 0;

  /***
   returns true if the member functions are thread-safe.
   */
  virtual bool thread_safe() const = 0;

  virtual ~element_function() = default;
};

/***
 the reporter class can be used to report results during the optimization.
 */
template<class EFunc, class Reporter = dummy_reporter,
         class interrupter = dummy_interrupter>
class optimizer {
  /***
   worker class to hold an element function and the element function''s
   Hessian approximation.
   */
  class worker {
    /// logical for whether it the first call
    bool first_call = true;

  public:
    /// the element function for this worker
    EFunc const func;
    /// number of elements
    size_t const n_ele = func.global_dim() + func.private_dim();
    /// memory for the Hessian approximation
    double * const PSQN_RESTRICT B;
    /// memory for the gradient
    double * const PSQN_RESTRICT gr = B + (n_ele * (n_ele + 1)) / 2L;
    /// memory for the old gradient
    double * const PSQN_RESTRICT gr_old = gr + n_ele;
    /// memory for the old value
    double * const PSQN_RESTRICT x_old = gr_old + n_ele;
    /// memory for the current value
    double * const PSQN_RESTRICT x_new = x_old + n_ele;
    /// indices of first set of private parameters
    size_t const par_start;
    /// bool for whether to use BFGS or SR1 updates
    bool use_bfgs = true;

    /***
     save the current parameter values and gradient in order to do update
     the Hessian approximation.
     */
    void record() noexcept {
      lp::copy(x_old , static_cast<double const*>(x_new), n_ele);
      lp::copy(gr_old, static_cast<double const*>(gr   ), n_ele);
    }

    /***
     resets the Hessian approximation.
     */
    void reset() noexcept {
      first_call = true;

      std::fill(B, B + n_ele * n_ele, 0.);
      // set diagonal entries to one
      double *b = B;
      for(size_t i = 0; i < n_ele; ++i, b += i + 1)
        *b = 1.;
    }

    worker(EFunc &&func, double * mem, size_t const par_start):
      func(func), B(mem), par_start(par_start) {
      reset();
    }

    /***
     computes the element function and possibly its gradient.

     @param global values for the global parameters
     @param vprivate values for private parameters.
     @param comp_grad logical for whether to compute the gradient
     */
    double operator()
      (double const * PSQN_RESTRICT global,
       double const * PSQN_RESTRICT vprivate, bool const comp_grad){
      // copy values
      size_t const d_global  = func.global_dim(),
                   d_private = func.private_dim();

      lp::copy(x_new           , global  , d_global);
      lp::copy(x_new + d_global, vprivate, d_private);

      if(!comp_grad)
        return func.func(static_cast<double const *>(x_new));

      double const out =  func.grad(
        static_cast<double const *>(x_new), gr);

      return out;
    }

    /***
     updates the Hessian approximation. Assumes that the () operator have
     been called at-least twice at different values.
     @param wmem working memory to use.
     */
    void update_Hes(double * const wmem){
      // differences in parameters and gradient
      double * const PSQN_RESTRICT s   = wmem,
             * const PSQN_RESTRICT y   = s + n_ele,
             * const PSQN_RESTRICT wrk = y + n_ele;

      lp::vec_diff(x_new, x_old , s, n_ele);

      bool all_unchanged(true);
      for(size_t i = 0; i < n_ele; ++i)
        if(abs(s[i]) > abs(x_new[i]) *
           std::numeric_limits<double>::epsilon() * 100){
          all_unchanged = false;
          break;
        }

      if(!all_unchanged){
        lp::vec_diff(gr, gr_old, y, n_ele);

        if(use_bfgs){
          double const s_y = lp::vec_dot(y, s, n_ele);
          if(first_call){
            first_call = false;
            // make update on page 143
            double const scal = lp::vec_dot(y, n_ele) / s_y;
            double *b = B;
            for(size_t i = 0; i < n_ele; ++i, b += i + 1)
              *b = scal;
          }

          // perform BFGS update
          std::fill(wrk, wrk + n_ele, 0.);
          lp::mat_vec_dot(B, s, wrk, n_ele);
          double const s_B_s = lp::vec_dot(s, wrk, n_ele);

          lp::rank_one_update(B, wrk, -1. / s_B_s, n_ele);

          if(s_y < .2 * s_B_s){
            // damped BFGS
            double const theta = .8 * s_B_s / (s_B_s - s_y);
            double *yi = y,
                   *wi = wrk;
            for(size_t i = 0; i < n_ele; ++i, ++yi, ++wi)
              *yi = theta * *yi + (1 - theta) * *wi;
            double const s_r = lp::vec_dot(y, s, n_ele);
            lp::rank_one_update(B, y, 1. / s_r, n_ele);

          } else
            // regular BFGS
            lp::rank_one_update(B, y, 1. / s_y, n_ele);

        } else {
          if(first_call){
            first_call = false;
            // make update on page 143
            double const scal =
              lp::vec_dot(y, n_ele) / lp::vec_dot(y, s, n_ele);
            double *b = B;
            for(size_t i = 0; i < n_ele; ++i, b += i + 1)
              *b = scal;
          }

          /// maybe perform SR1
          std::fill(wrk, wrk + n_ele, 0.);
          lp::mat_vec_dot(B, s, wrk, n_ele);
          for(size_t i = 0; i < n_ele; ++i){
            wrk[i] *= -1;
            wrk[i] += y[i];
          }
          double const s_w = lp::vec_dot(s, wrk, n_ele),
                    s_norm = sqrt(abs(lp::vec_dot(s, n_ele))),
                  wrk_norm = sqrt(abs(lp::vec_dot(wrk, n_ele)));
          constexpr double const r = 1e-8;
          if(abs(s_w) > r * s_norm * wrk_norm)
            lp::rank_one_update(B, wrk, 1. / s_w, n_ele);

        }
      } else
        // essentially no change in the input. Reset the Hessian
        // approximation
        reset();

      record();
    }
  };

  /**
   class to optimize one set of private parameters given the global
   parameters. */
  class sub_problem final : public problem {
    worker &w;
    double const * g_val;
    size_t const p_dim = w.func.private_dim(),
                 g_dim = w.func.global_dim();

  public:
    sub_problem(worker &w, double const *g_val): w(w), g_val(g_val) { }

    size_t size() const {
      return p_dim;
    }
    double func(double const *val){
      return w(g_val, val, false);
    }
    double grad(double const * PSQN_RESTRICT val,
                double       * PSQN_RESTRICT gr){
      double const out = w(g_val, val, true);
      for(size_t i = 0; i < p_dim; ++i)
        gr[i] = w.gr[i + g_dim];
      return out;
    }
  };

public:
  /// dimension of the global parameters
  size_t const global_dim;
  /// true if the element functions are thread-safe
  bool const is_ele_func_thread_safe;
  /// total number of parameters
  size_t const n_par;

private:
  /***
   size of the allocated working memory. The first element is needed for
   the workers. The second element is needed during the computation for the
   master thread. The third element is number required per thread.
  */
  std::array<size_t, 3L> const n_mem;
  /// maximum number of threads to use
  size_t const max_threads;
  /// working memory
  std::unique_ptr<double[]> mem =
    std::unique_ptr<double[]>(
        new double[n_mem[0] + n_mem[1] + max_threads * n_mem[2]]);
  /// pointer to temporary memory to use on the master thread
  double * const temp_mem = mem.get() + n_mem[0];
  /// pointer to temporray memory to be used by the threads
  double * const temp_thread_mem = temp_mem + n_mem[1];
  /// element functions
  std::vector<worker> funcs;
  /// number of function evaluations
  size_t n_eval = 0L;
  /// number of gradient evaluations
  size_t n_grad = 0L;
  /// number of iterations of conjugate gradient
  size_t n_cg = 0L;

  /***
   reset the counters for the number of evaluations
   */
  void reset_counters() {
    n_eval = 0L;
    n_grad = 0L;
    n_cg = 0L;
  }

  /// returns the thread number.
  int get_thread_num() const noexcept {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0L;
#endif
  }

  /// returns working memory for this thread
  double * get_thread_mem(int const thread_num) const noexcept {
    return temp_thread_mem + thread_num * n_mem[2];
  }

  double * get_thread_mem() const noexcept {
    return get_thread_mem(get_thread_num());
  }

  /// number of threads to use
  int n_threads = 1L;

public:
  /// set the number of threads to use
  void set_n_threads(size_t const n_threads_new) noexcept {
#ifdef _OPENMP
    n_threads = std::max(
      static_cast<size_t>(1L), std::min(n_threads_new, max_threads));
    omp_set_num_threads(n_threads);
    omp_set_dynamic(0L);
#endif
  }

  /***
   takes in a vector with element functions and constructs the optimizer.
   The members are moved out of the vector.
   @param funcs_in vector with element functions.
   @param max_threads maximum number of threads to use.
   */
  optimizer(std::vector<EFunc> &funcs_in, size_t const max_threads):
  global_dim(([&]() -> size_t {
    if(funcs_in.size() < 1L)
      throw std::invalid_argument(
          "optimizer<EFunc>::optimizer: no functions supplied");
    return funcs_in[0].global_dim();
  })()),
  is_ele_func_thread_safe(funcs_in[0].thread_safe()),
  n_par(([&]() -> size_t {
    size_t out(global_dim);
    for(auto &f : funcs_in)
      out += f.private_dim();
    return out;
  })()),
  n_mem(([&]() -> std::array<size_t, 3L> {
    size_t out(0L),
           max_priv(0L);
    for(auto &f : funcs_in){
      if(f.global_dim() != global_dim)
        throw std::invalid_argument(
            "optimizer<EFunc>::optimizer: global_dim differs");
      if(f.thread_safe() != is_ele_func_thread_safe)
        throw std::invalid_argument(
            "optimizer<EFunc>::optimizer: thread_safe differs");
      size_t const private_dim = f.private_dim(),
                   n_ele       = private_dim + global_dim;
      if(max_priv < private_dim)
        max_priv = private_dim;

      out += n_ele * 4L + (n_ele * (n_ele + 1L)) / 2L;
    }

    constexpr size_t const mult = cacheline_size() / sizeof(double),
                       min_size = 2L * mult;

    size_t thread_mem = std::max(3 * (global_dim + max_priv), min_size);
    thread_mem = (thread_mem + mult - 1L) / mult;
    thread_mem *= mult;

    std::array<size_t, 3L> ret = {
      out,
      5L * n_par + (global_dim * (global_dim + 1L)) / 2L,
      thread_mem };
    return ret;
  })()),
  max_threads(max_threads > 0 ? max_threads : 1L),
  funcs(([&]() -> std::vector<worker> {
    std::vector<worker> out;
    size_t const n_ele(funcs_in.size());
    out.reserve(funcs_in.size());

    double * mem_ptr = mem.get();
    size_t i_start(global_dim);
    for(size_t i = 0; i < n_ele; ++i){
      worker new_func(std::move(funcs_in[i]), mem_ptr, i_start);
      out.emplace_back(std::move(new_func));
      size_t const n_ele = out.back().n_ele;
      mem_ptr += n_ele * 4L + (n_ele * (n_ele + 1L)) / 2L;
      i_start += out.back().func.private_dim();
    }

    return out;
  })()) { }

  /***
   evalutes the partially separable and also the gradient if requested.
   @param val pointer to the value to evaluate the function at.
   @param gr pointer to store gradient in.
   @param comp_grad boolean for whether to compute the gradient.
   */
  double eval(double const * val, double * PSQN_RESTRICT gr,
              bool const comp_grad){
    if(comp_grad)
      n_grad++;
    else
      n_eval++;

    size_t const n_funcs = funcs.size();
    auto serial_version = [&]() -> double {
      double out(0.);
      for(size_t i = 0; i < n_funcs; ++i){
        auto &f = funcs[i];
        out += f(val, val + f.par_start, comp_grad);
      }

      if(comp_grad){
        std::fill(gr, gr + global_dim, 0.);
        for(size_t i = 0; i < n_funcs; ++i){
          auto const &f = funcs[i];
          for(size_t j = 0; j < global_dim; ++j)
            gr[j] += f.gr[j];

          size_t const iprivate = f.func.private_dim();
          lp::copy(gr + f.par_start, f.gr + global_dim, iprivate);
        }
      }

      return out;
    };

    if(n_threads < 2L or !is_ele_func_thread_safe)
      return serial_version();

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
    {
    double * r_mem = get_thread_mem(),
           * v_mem =
               r_mem + global_dim + 1L /* leave 1 ele for func value*/;
    lp::copy(v_mem, val, global_dim);
    if(comp_grad)
      std::fill(r_mem, r_mem + global_dim, 0.);

    double &thread_terms = *(r_mem + global_dim);
    thread_terms = 0;
#pragma omp for schedule(static)
    for(size_t i = 0; i < n_funcs; ++i){
      auto &f = funcs[i];
      thread_terms += f(v_mem, val + f.par_start, comp_grad);

      if(comp_grad){
        // update global
        double *lhs = r_mem;
        double const *rhs = f.gr;
        for(size_t j = 0; j < global_dim; ++j, ++lhs, ++rhs)
          *lhs += *rhs;

        // update private
        lp::copy(gr + f.par_start, f.gr + global_dim, f.func.private_dim());
      }
    }
    }

    if(comp_grad)
      std::fill(gr, gr + global_dim, 0.);

    // add to global parameters
    double out(0.);
    for (int t = 0; t < n_threads; t++){
      double const *r_mem = get_thread_mem(t);
      if(comp_grad)
        for(size_t i = 0; i < global_dim; ++i)
          gr[i] += r_mem[i];
      out += r_mem[global_dim];
    }

    return out;
#else
    return serial_version();
#endif
  }

  /***
   computes y <- y + B.x where B is the current Hessian approximation.
   @param val vector on the right-hand side.
   @param res output vector on the left-hand side.
   @param B_start memory with the [global_dim] x [global_dim] part of B.
   @param comp_B_start true if B_start should be computed.
   ***/
  void B_vec(double const * const PSQN_RESTRICT val,
             double * const PSQN_RESTRICT res,
             double * const PSQN_RESTRICT B_start,
             bool const comp_B_start) const noexcept {
    size_t const n_funcs = funcs.size();

    // aggregate the first part of B if needed
    if(comp_B_start){
      size_t const B_sub_ele = (global_dim * (global_dim + 1)) / 2;
      std::fill(B_start, B_start + B_sub_ele, 0);
      for(size_t i = 0; i < n_funcs; ++i){
        auto &f = funcs[i];
        double * b = B_start;
        double const *b_inc = f.B;
        for(size_t j = 0; j < B_sub_ele; ++j)
          *b++ += *b_inc++;
      }
    }

    // compute the first part
    lp::mat_vec_dot(B_start, val, res, global_dim);

    // the serial version
    auto serial_version = [&]() -> void {
      for(size_t i = 0; i < n_funcs; ++i){
        auto &f = funcs[i];
        size_t const iprivate = f.func.private_dim(),
               private_offset = f.par_start;

        lp::mat_vec_dot_excl_first(f.B, val, val + private_offset, res,
                                   res + private_offset, global_dim,
                                   iprivate);
      }
    };

    if(n_threads < 2L){
      serial_version();
      return;
    }

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
    {
    double * r_mem = get_thread_mem(),
           * v_mem = r_mem + global_dim;
    lp::copy(v_mem, val, global_dim);
    std::fill(r_mem, r_mem + global_dim, 0.);

#pragma omp for schedule(static)
    for(size_t i = 0; i < n_funcs; ++i){
      auto &f = funcs[i];
      size_t const iprivate = f.func.private_dim(),
             private_offset = f.par_start;

      lp::mat_vec_dot_excl_first(f.B, v_mem, val + private_offset, r_mem,
                                 res + private_offset, global_dim,
                                 iprivate);
    }
    }

    // add to global parameters
    for (int t = 0; t < n_threads; t++){
      double const *r_mem = get_thread_mem(t);
      for(size_t i = 0; i < global_dim; ++i)
        res[i] += r_mem[i];
    }
#else
    serial_version();
#endif
  }

  /**
    computes the diagonl of diag(B).
   */
  void get_diag(double * x){
    std::fill(x, x + global_dim, 0.);
    double * PSQN_RESTRICT x_priv = x + global_dim;

    for(size_t i = 0; i < funcs.size(); ++i){
      auto &f = funcs[i];
      size_t const iprivate = f.func.private_dim();

      // add to the global parameters
      double const * b_diag = f.B;
      size_t j = 0L;
      for(; j <            global_dim; ++j, b_diag += j + 1)
        x[j] += *b_diag;

      for(; j < iprivate + global_dim; ++j, b_diag += j + 1)
        *x_priv++ = *b_diag;
    }
  }

  /***
    conjugate gradient method with diagonal preconditioning. Solves B.y = x
    where B is the Hessian approximation.
    @param tol convergence threshold.
    @param max_cg maximum number of conjugate gradient iterations.
    @param trace controls the amount of tracing information.
    @param pre_method preconditioning method.
   */
  bool conj_grad(double const * PSQN_RESTRICT x, double * PSQN_RESTRICT y,
                 double const tol, size_t const max_cg,
                 int const trace, precondition const pre_method){
    double * PSQN_RESTRICT r       = temp_mem,
           * PSQN_RESTRICT p       = r      + n_par,
           * PSQN_RESTRICT B_p     = p      + n_par,
           * PSQN_RESTRICT v       = B_p    + n_par,
           * PSQN_RESTRICT B_diag  = v      + n_par,
           * PSQN_RESTRICT B_start = B_diag + n_par;
    bool const do_pre = pre_method == 1L;

    // setup before first iteration
    if(do_pre){
      get_diag(B_diag);
      double *b = B_diag;
      for(size_t i = 0; i < n_par; ++i, ++b)
        *b = 1. / *b; // want to use multiplication rather than division
    }

    auto diag_solve = [&](double       * PSQN_RESTRICT vy,
                          double const * PSQN_RESTRICT vx) -> void {
      double * di = B_diag;
      for(size_t i = 0; i < n_par; ++i)
        *vy++ = *vx++ * *di++;
    };

    std::fill(y, y + n_par, 0.);
    for(size_t i = 0; i < n_par; ++i){
      r[i] = -x[i];
      if(!do_pre)
        p[i] = x[i];
    }

    if(do_pre){
      diag_solve(v, r);
      for(size_t i = 0; i < n_par; ++i)
        p[i] = -v[i];
    }

    auto get_r_v_dot = [&]() -> double {
      return do_pre ? lp::vec_dot(r, v, n_par) : lp::vec_dot(r, n_par);
    };

    double old_r_v_dot = get_r_v_dot();
    for(size_t i = 0; i < max_cg; ++i){
      ++n_cg;
      std::fill(B_p, B_p + n_par, 0.);
      B_vec(p, B_p, B_start, i == 0);

      double const p_B_p = lp::vec_dot(p, B_p, n_par);
      if(p_B_p <= 0){
        // negative curvature. Thus, exit
        if(i < 1L)
          // set output to be the gradient
          for(size_t j = 0; j < n_par; ++j)
            y[j] = x[j];

        break;
      }
      double const alpha = old_r_v_dot / p_B_p;

      for(size_t j = 0; j < n_par; ++j){
        y[j] += alpha *   p[j];
        r[j] += alpha * B_p[j];
      }

      if(do_pre)
        diag_solve(v, r);
      double const r_v_dot = get_r_v_dot(),
                   t_val   = do_pre ? sqrt(abs(lp::vec_dot(r, n_par))) :
                                      sqrt(abs(r_v_dot));
      Reporter::cg_it(trace, i, max_cg, t_val, tol);
      if(t_val < tol)
        break;

      double const beta = r_v_dot / old_r_v_dot;
      old_r_v_dot = r_v_dot;
      for(size_t j = 0; j < n_par; ++j){
        p[j] *= beta;
        p[j] -= do_pre ? v[j] : r[j];
      }
    }

    return true;
  }

  /***
   performs line search to satisfy the Wolfe condition.
   @param f0 value of the functions at the current value.
   @param x0 value the function is evaluted.
   @param gr0 value of the current gradient.
   @param dir direction to search in.
   @param fnew the function value at the found solution.
   @param c1,c2 thresholds for Wolfe condition.
   @param strong_wolfe true if the strong Wolfe condition should be used.
   @param trace controls the amount of tracing information.

   x0 and gr0 contains the new value and gradient on return. The method
   returns false if the line search fails.
   */
  bool line_search(
      double const f0, double * PSQN_RESTRICT x0, double * PSQN_RESTRICT gr0,
      double * PSQN_RESTRICT dir, double &fnew, double const c1,
      double const c2, bool const strong_wolfe, int const trace){
    double * const x_mem = temp_mem;

    // declare 1D functions
    auto psi = [&](double const alpha) -> double {
      for(size_t i = 0; i < n_par; ++i)
        x_mem[i] = x0[i] + alpha * dir[i];

      return eval(x_mem, nullptr, false);
    };

    // returns the function value and the gradient
    auto dpsi = [&](double const alpha) -> double {
      for(size_t i = 0; i < n_par; ++i)
        x_mem[i] = x0[i] + alpha * dir[i];

      fnew = eval(x_mem, gr0, true);
      return lp::vec_dot(gr0, dir, n_par);
    };

    // the above at alpha = 0
    double const dpsi_zero = lp::vec_dot(gr0, dir, n_par);
    if(dpsi_zero > 0)
      // not a descent direction
      return false;

    constexpr size_t const max_it = 20L;
    static double const NaNv = std::numeric_limits<double>::quiet_NaN();
    auto zoom =
      [&](double a_low, double a_high, intrapolate &inter) -> bool {
        double f_low = psi(a_low);
        for(size_t i = 0; i < max_it; ++i){
          double const ai = inter.get_value(a_low, a_high),
                       fi = psi(ai);
          if(!std::isfinite(fi)){
            // this is very bad! We do not know in which direction to go.
            // We try to go in a lower direction
            if(a_low < a_high)
              a_high = ai;
            else
              a_low = ai;
            continue;
          }

          inter.update(ai, fi);
          Reporter::line_search_inner(trace, a_low, ai, fi, true,
                                      NaNv, a_high);

          if(fi > f0 + c1 * ai * dpsi_zero or fi >= f_low){
            a_high = ai;
            continue;
          }

          double const dpsi_i = dpsi(ai);
          Reporter::line_search_inner(trace, a_low, ai, fi, true,
                                      dpsi_i, a_high);
          double const test_val = strong_wolfe ? abs(dpsi_i) : -dpsi_i;
          if(test_val <= - c2 * dpsi_zero)
            return true;

          if(dpsi_i * (a_high - a_low) >= 0.)
            a_high = a_low;

          a_low = ai;
          f_low = fi;
        }

        return false;
      };

    double fold(f0),
         a_prev(0),
             ai(.5);
    bool found_ok_prev = false,
         failed_once   = false;
    double mult = 2;
    for(size_t i = 0; i < max_it; ++i){
      ai *= mult;
      double fi = psi(ai);
      Reporter::line_search_inner(trace, a_prev, ai, fi, false,
                                  NaNv, NaNv);
      if(!std::isfinite(fi)){
        // handle inf/nan case
        failed_once = true;
        mult = .5;

        if(!found_ok_prev)
          // no valid previous value yet to use
          continue;
        else {
          // the previous value was ok. Use that one
          fi = fold;
          ai = a_prev;
        }
      }

      if(fi > f0 + c1 * ai * dpsi_zero or (found_ok_prev and fi > fold)){
        intrapolate inter(f0, dpsi_zero, ai, fi);
        bool const out = zoom(a_prev, ai, inter);
        lp::copy(x0, x_mem, n_par);
        return out;
      }

      double const dpsi_i = dpsi(ai);
      Reporter::line_search_inner(trace, a_prev, ai, fi, false,
                                  dpsi_i, NaNv);

      double const test_val = strong_wolfe ? abs(dpsi_i) : -dpsi_i;
      if(test_val <= - c2 * dpsi_zero){
        lp::copy(x0, x_mem, n_par);
        return true;
      }

      if(failed_once and fi < f0){
        // effectively just line search
        lp::copy(x0, x_mem, n_par);
        return false;
      }

      if(dpsi_i >= 0){
        intrapolate inter = ([&]() -> intrapolate {
          if(found_ok_prev){
            // we have two values that we can use
            intrapolate out(f0, dpsi_zero, a_prev, fold);
            out.update(ai, fi);
            return out;
          }

          return intrapolate(f0, dpsi_zero, ai, fi);
        })();
        bool const out = zoom(ai, a_prev, inter);
        lp::copy(x0, x_mem, n_par);
        return out;
      }

      found_ok_prev = true;
      a_prev = ai;
      fold = fi;
    }

    return false;
  }

  /***
   optimizes the partially separable function.
   @param val pointer to starting value. Set to the final estimate at the
   end.
   @param rel_eps relative convergence threshold.
   @param max_it maximum number of iterations.
   @param c1,c2 thresholds for Wolfe condition.
   @param use_bfgs bool for whether to use BFGS updates or SR1 updates.
   @param trace integer with info level passed to reporter.
   @param cg_tol threshold for conjugate gradient method.
   @param strong_wolfe true if the strong Wolfe condition should be used.
   @param max_cg maximum number of conjugate gradient iterations in each
   iteration. Use zero if there should not be a limit.
   @param pre_method preconditioning method.
   */
  optim_info optim
    (double * val, double const rel_eps, size_t const max_it,
     double const c1, double const c2,
     bool const use_bfgs = true, int const trace = 0,
     double const cg_tol = .5, bool const strong_wolfe = true,
     size_t const max_cg = 0,
     precondition const pre_method = precondition::diag){
    reset_counters();
    for(auto &f : funcs){
      f.reset();
      f.use_bfgs = use_bfgs;
    }

    std::unique_ptr<double[]> gr(new double[n_par]),
                             dir(new double[n_par]);

    // evaluate the gradient at the current value
    double fval = eval(val, gr.get(), true);
    for(auto &f : funcs)
      f.record();

    info_code info = info_code::max_it_reached;
    int n_line_search_fail = 0;
    for(size_t i = 0; i < max_it; ++i){
      if(i % 10 == 0)
        if(interrupter::check_interrupt()){
          info = info_code::user_interrupt;
          break;
        }

      double const fval_old = fval,
                     gr_nom = sqrt(abs(lp::vec_dot(gr.get(), n_par))),
                 cg_tol_use = std::min(cg_tol, sqrt(gr_nom)) * gr_nom;
      if(!conj_grad(gr.get(), dir.get(), cg_tol_use,
                    max_cg < 1 ? n_par : max_cg, trace,
                    pre_method)){
        info = info_code::conjugate_gradient_failed;
        Reporter::cg(trace, i, n_cg, false);
        break;
      } else
        Reporter::cg(trace, i, n_cg, true);

      for(double * d = dir.get(); d != dir.get() + n_par; ++d)
        *d *= -1;

      double const x1 = *val;
      if(!line_search(fval_old, val, gr.get(), dir.get(), fval, c1, c2,
                      strong_wolfe, trace)){
        info = info_code::line_search_failed;
        Reporter::line_search
          (trace, i, n_eval, n_grad, fval_old, fval, false,
           std::numeric_limits<double>::quiet_NaN(),
           const_cast<double const *>(val), global_dim);
        if(++n_line_search_fail > 2)
          break;
      } else{
        n_line_search_fail = 0;
        Reporter::line_search
          (trace, i, n_eval, n_grad, fval_old, fval, true,
           (*val - x1) / *dir.get(), const_cast<double const *>(val),
           global_dim);

      }

      bool const has_converged =
        abs(fval - fval_old) < rel_eps * (abs(fval_old) + rel_eps);
      if(has_converged){
        info = info_code::converged;
        break;
      }

      // update the Hessian and take another iteation
      if(n_line_search_fail < 1){
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
        {
          double * const tmp_mem_use = get_thread_mem();
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
          for(size_t i = 0; i < funcs.size(); ++i)
            funcs[i].update_Hes(tmp_mem_use);
        }
      } else
        for(size_t i = 0; i < funcs.size(); ++i){
          funcs[i].reset();
          funcs[i].record();
        }
    }

    return { fval, info, n_eval, n_grad, n_cg };
  }

  /***
   returns the current Hessian approximation.
   */
  void get_hess(double * const PSQN_RESTRICT hess) const {
    // TODO: make an implementation which returns a sparse Hessian
    std::fill(hess, hess + n_par * n_par, 0.);

    size_t private_offset(global_dim);
    for(auto &f : funcs){
      size_t const iprivate = f.func.private_dim();

      auto get_i = [&](size_t const i, size_t const j) -> size_t {
        size_t const ii = std::min(i, j),
                     jj = std::max(j, i);

        return ii + (jj * (jj + 1L)) / 2L;
      };

      double const * const b = f.B;
      {
        double *h1 = hess,
               *h2 = hess + private_offset;
        for(size_t j = 0; j < global_dim;
            ++j, h1 += n_par, h2 += n_par){
          for(size_t i = 0; i < global_dim; ++i)
            h1[i] += b[get_i(i             , j)];
          for(size_t i = 0; i < iprivate; ++i)
            h2[i] += b[get_i(i + global_dim, j)];
        }
      }

      double *h1 = hess + private_offset * n_par,
             *h2 = h1 + private_offset;
      for(size_t j = 0; j < iprivate;
          ++j, h1 += n_par, h2 += n_par){
        for(size_t i = 0; i < global_dim; ++i)
          h1[i] += b[get_i(i             , j + global_dim)];
        for(size_t i = 0; i < iprivate; ++i)
          h2[i] += b[get_i(i + global_dim, j + global_dim)];
      }

      private_offset += iprivate;
    }
  }

  /***
   optimizes the private parameters given the global parameters.
   @param val pointer to starting value. Set to the final estimate at the
   end.
   @param rel_eps relative convergence threshold.
   @param max_it maximum number of iterations.
   @param c1,c2 thresholds for Wolfe condition.
   */
  double optim_priv
  (double * val, double const rel_eps, size_t const max_it,
   double const c1, double const c2){
    double out(0.);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads) reduction(+:out) if(is_ele_func_thread_safe)
#endif
    for(size_t i = 0; i < funcs.size(); ++i){
      auto &f = funcs[i];
      sub_problem prob(f, val);
      double * const p_val = val + f.par_start;

      auto const opt_out = bfgs(prob, p_val, rel_eps, max_it, c1, c2, 0L);
      out += opt_out.value;
    }

    return out;
  }
};

template<class EFunc>
class optimizer<EFunc, dummy_reporter, dummy_interrupter>;

} // namespace PSQN

#endif
