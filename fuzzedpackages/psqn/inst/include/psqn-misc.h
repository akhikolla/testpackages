#ifndef PSQN_MISC_H
#define PSQN_MISC_H

namespace PSQN {
enum info_code : int {
  max_it_reached = -1L,
  conjugate_gradient_failed = -2,
  line_search_failed = -3,
  user_interrupt = -4,
  converged = 0
};

enum precondition : int {
   non = 0,
   diag = 1
};

struct optim_info {
  double value;
  info_code info;
  size_t n_eval, n_grad, n_cg;
};

struct dummy_reporter {
  /**
   reporting during the conjugate gradient method.

   @param trace info level passed by the user.
   @param iteration itteration number starting at zero.
   @param maxit maximum number of iterations.
   @param r_norm norm of the residual vector.
   @param threshold convergence threshold.
   */
  static void cg_it(int const trace, size_t const iteration,
                    size_t const maxit, double const r_norm,
                    double const threshold) { }

  /**
   reporting after conjugate gradient method.

   @param trace info level passed by the user.
   @param iteration itteration number starting at zero.
   @param n_cg number of conjugate gradient iterations.
   @param successful true of conjugate gradient method succeeded.
   */
  static void cg(int const trace, size_t const iteration,
                 size_t const n_cg, bool const successful) { }

  /**
   reporting during line search.

   @param trace info level passed by the user.
   @param a_old lower value in zoom or previous value if not zoom.
   @param a_new new value to test.
   @param f_new function value at a_new.
   @param is_zoom true if the call is from the zoom function.
   @param d_new derivative at a_new.
   @param a_high upper value in zoom.
   */
  static void line_search_inner
  (int const trace, double const a_old, double const a_new,
   double const f_new, bool const is_zoom, double const d_new,
   double const a_high) { }

  /**
   reporting after line search.

   @param trace info level passed by the user.
   @param iteration itteration number starting at zero.
   @param n_eval number of function evaluations.
   @param n_grad number of gradient evaluations.
   @param fval_old old function value.
   @param fval new function value.
   @param successful true of conjugate gradient method succeeded.
   @param step_size found step size.
   @param new_x new parameter value.
   @param n_global number of global parameters.
   */
  static void line_search
  (int const trace, size_t const iteration, size_t const n_eval,
   size_t const n_grad, double const fval_old,
   double const fval, bool const successful, double const step_size,
   double const *new_x, size_t const n_global) { }
};

class dummy_interrupter {
public:
   /**
    use during the algorithm to interup the computation. Should be
    thread-safe and return true if the optimization should stop. */
   static bool check_interrupt() {
      return false;
   }
};

} // namespace PSQN

#endif
