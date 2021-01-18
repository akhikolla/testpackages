// Copyright (C) 2015-2016 The University of Manchester
//
// Written by Yi Pan - ypan1988@gmail.com

template <typename T>
LineSearch<T>::LineSearch() {}

template <typename T>
LineSearch<T>::~LineSearch() {}

/**
 * Line Searches and backtracking
 *
 * @param func Instance of function to be optimized
 * @param x Parameters
 * @param p Newton step
 * @param stepmax Maximum step length
 */
template <typename T>
void LineSearch<T>::GetStep(T &func, arma::vec &x, arma::vec &p,
                            const double stepmax) {
  int debug = 0;

  // Maximum number of iterations
  const int kIterMax = 200;

  // Ensure sufficient decrease in function value
  const double kAlpha = 1.0e-4;

  // The convergence criterion on Delta X
  const double kTolX = std::numeric_limits<double>::epsilon();

  const int n_pars = x.n_rows;  // number of parameters

  const arma::vec xold = x;
  const double fold = func(xold);
  arma::vec grad;
  func.Gradient(xold, grad);

  // Scale if attempted step is too big
  double sum = sqrt(arma::dot(p, p));
  if (sum > stepmax) p *= stepmax / sum;

  double slope = dot(grad, p);
  if (slope >= 0.0 && message_)
    Rcpp::Rcerr << "Roundoff problem in linesearch." << std::endl;

  // Calculate the minimum step length
  double test = 0.0;
  for (int i = 0; i != n_pars; ++i) {
    double temp = std::abs(p(i)) / std::max(std::abs(xold(i)), 1.0);
    if (temp > test) test = temp;
  }
  double stepmin = kTolX / test;

  double lambda, lambda2, lambda_tmp, f, f2;
  lambda = 1.0;  // Always try full Newton step first
  lambda2 = lambda_tmp = f = f2 = 0.0;
  for (int iter = 0; iter != kIterMax; ++iter) {
    // Start of iteration loop
    x = xold + lambda * p;
    f = func(x);

    if (debug) {
      Rcpp::Rcout << "iter " << iter << ": "
                  << "\tlambda = " << lambda << "\tf = " << f << std::endl;
    }

    if (lambda < stepmin) {
      // x is too close to xold, ignored
      x = xold;
      return;
    } else if (f <= fold + kAlpha * lambda * slope) {
      // Sufficient function decrease
      return;
    } else if (IsInfOrNaN(f)) {
      // f is INF or NAN
      while (!IsInfOrNaN(lambda) && IsInfOrNaN(f)) {
        lambda *= 0.5;
        x = xold + lambda * p;
        f = func(x);
      }

      if (debug) {
        Rcpp::Rcout << "iter " << iter << ":(INF) "
                    << "\tlambda = " << lambda << "\tf = " << f << std::endl;
      }

      lambda_tmp = 0.5 * lambda;
    } else {
      // Backtrack
      if (lambda == 1.0) {
        lambda_tmp = -slope / (2.0 * (f - fold - slope));
      } else {
        double rhs1 = f - fold - lambda * slope;
        double rhs2 = f2 - fold - lambda2 * slope;
        double a = rhs1 / (lambda * lambda) / (lambda - lambda2) -
                   rhs2 / (lambda2 * lambda2) / (lambda - lambda2);
        double b = -lambda2 * rhs1 / (lambda * lambda) / (lambda - lambda2) +
                   lambda * rhs2 / (lambda2 * lambda2) / (lambda - lambda2);
        if (IsInfOrNaN(a) || IsInfOrNaN(b)) {
          lambda_tmp = 0.5 * lambda;
        } else if (a == 0.0) {
          lambda_tmp = -slope / (2.0 * b);
        } else {
          double disc = b * b - 3.0 * a * slope;
          if (disc < 0.0) {
            lambda_tmp = 0.5 * lambda;
          } else if (b <= 0.0) {
            lambda_tmp = (-b + sqrt(disc)) / (3.0 * a);
          } else {
            lambda_tmp = -slope / (b + sqrt(disc));
          }
        }
        if (lambda_tmp > 0.5 * lambda || IsInfOrNaN(lambda_tmp)) {
          lambda_tmp = 0.5 * lambda;
        }
        if (debug) {
          Rcpp::Rcout << "a = " << a << "\tb = " << b
                      << "\tlambda_tmp = " << lambda_tmp << std::endl;
        }
      }
    }
    lambda2 = lambda;
    f2 = f;
    lambda = std::max(lambda_tmp, 0.1 * lambda);

    if (debug) {
      Rcpp::Rcout << "iter " << iter << ": "
                  << "\tlambda = " << lambda << "\tf = " << f << std::endl;
    }
  }
}

template <typename T>
bool LineSearch<T>::IsInfOrNaN(double x) {
  return (x == std::numeric_limits<double>::infinity() ||
          x == -std::numeric_limits<double>::infinity() || x != x);
}
