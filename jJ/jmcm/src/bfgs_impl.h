// Copyright (C) 2015-2016 The University of Manchester
//
// Written by Yi Pan - ypan1988@gmail.com

/**
 * Constructor
 */
template <typename T>
BFGS<T>::BFGS() : LineSearch<T>() {}

/**
 * Destructor
 */
template <typename T>
BFGS<T>::~BFGS() {}

/**
 * Implement the quasi-Newton method
 *
 * @param func Instance of function to be optimized
 * @param x Parameters
 * @param grad_tol The convergence requirement on zeroing the gradient
 */
template <typename T>
void BFGS<T>::Optimize(T &func, arma::vec &x, const double grad_tol) {
  int debug = 0;

  // Maximum number of iterations
  const int kIterMax = 200;

  // Machine precision
  const double kEpsilon = std::numeric_limits<double>::epsilon();

  // Convergence criterion on x values
  const double kTolX = 4 * kEpsilon;

  // Scaled maximum step length allowed in line searches
  const double kScaStepMax = 100;

  const int n_pars = x.n_rows;  // number of parameters

  // Calculate starting function value and gradient
  if (debug) {
    Rcpp::Rcout << "Initializing the function value" << std::endl;
  }

  double f = func(x);

  if (debug) {
    Rcpp::Rcout << "Initializing the gradient" << std::endl;
  }

  arma::vec grad;
  func.Gradient(x, grad);

  // Initialize the inverse Hessian to a unit matrix
  arma::mat hess_inv = arma::eye<arma::mat>(n_pars, n_pars);

  // Initialize Newton Step
  arma::vec p = -hess_inv * grad;

  // Calculate the maximum step length
  double sum = sqrt(arma::dot(x, x));
  const double kStepMax = kScaStepMax * std::max(sum, double(n_pars));

  // Main loop over the iterations
  for (int iter = 0; iter != kIterMax; ++iter) {
    n_iters_ = iter;

    arma::vec x2 = x;  // Save the old point

    this->GetStep(func, x, p, kStepMax);

    f = func(x);  // Update function value
    p = x - x2;   // Update line direction
    x2 = x;       // Update the current point
    f_min_ = f;

    if (trace_) {
      Rcpp::Rcout << std::setw(5) << iter << ": " << std::setw(10) << f << ": ";
      x.t().print();
    }

    // Test for convergence on Delta x
    double test = 0.0;
    for (int i = 0; i != n_pars; ++i) {
      double temp = std::abs(p(i)) / std::max(std::abs(x(i)), 1.0);
      if (temp > test) test = temp;
    }
    if (test < kTolX) return;

    arma::vec grad2 = grad;  // Save the old gradient
    func.Gradient(x, grad);  // Get the new gradient

    // Test for convergence on zero gradient
    test = 0.0;
    double den = std::max(f, 1.0);
    for (int i = 0; i != n_pars; ++i) {
      double temp = std::abs(grad(i)) * std::max(std::abs(x(i)), 1.0) / den;
      if (temp > test) test = temp;
    }
    if (test < grad_tol) return;

    // Compute difference of gradients
    arma::vec dg = grad - grad2;
    // Compute difference times current matrix
    arma::vec hdg = hess_inv * dg;

    // Calculate dot products for the denominators
    double fac = dot(dg, p);
    double fad = 0.0;
    double fae = dot(dg, hdg);
    double sumdg = dot(dg, dg);
    double sump = dot(p, p);

    // Skip update if fac not sufficiently positive
    if (fac > sqrt(kEpsilon * sumdg * sump)) {
      fac = 1.0 / fac;
      fad = 1.0 / fae;

      // The vector that makes BFGS different from DFP
      arma::vec u = fac * p - fad * hdg;

      hess_inv += fac * p * p.t() - fad * hdg * hdg.t() + fae * u * u.t();
    }

    // Calculate the next direction to go
    p = -hess_inv * grad;
  }
  if (this->message_) {
    Rcpp::Rcerr << "too many iterations in bfgs" << std::endl;
  }
}

template <typename T>
int BFGS<T>::n_iters() const {
  return n_iters_;
}

template <typename T>
double BFGS<T>::f_min() const {
  return f_min_;
}
