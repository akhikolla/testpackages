/* Copyright 2016 The University of Manchester.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Written by Yi Pan - ypan1988@gmail.com
   ==============================================================================*/

/**
 * Line Searches and backtracking
 *
 * @param func Instance of function to be optimized
 * @param x Parameters
 * @param p Newton step
 * @param stepmax Maximum step length
 */
template <typename T>
bool LineSearch<T>::GetStep(T &func, double *f, arma::vec *x,
                            const arma::vec &g, const arma::vec &p,
                            const double stepmax) {
  bool check = false;

  const double fold = *f;
  const arma::vec xold = *x;
  const arma::vec gradient = g;
  arma::vec direction = p;

  const int kIterMax = 200;      // Maximum number of iterations
  const double kAlpha = 1.0e-4;  // Ensure sufficient decrease in function value
  const double kTolX =
    std::numeric_limits<double>::epsilon();  // The convergence criterion on
  // Delta X
  const arma::uword n_param = (*x).n_rows;     // number of parameters

  // Scale if attempted step is too big
  double sum = std::sqrt(arma::dot(direction, direction));
  if (sum > stepmax) direction *= stepmax / sum;

  double slope = dot(gradient, direction);
  if (slope >= 0.0 && message_)
    Rcpp::Rcerr << "LineSearch<T>::GetStep(): Roundoff problem." << std::endl;

  // Calculate the minimum step length
  double test = 0.0;
  for (arma::uword i = 0; i < n_param; ++i) {
    double temp = std::abs(direction(i)) / std::max(std::abs(xold(i)), 1.0);
    if (temp > test) test = temp;
  }
  double stepmin = kTolX / test;

  double f2 = 0.0;
  double lambda2, lambda_tmp;
  double lambda = 1.0;  // Always try full Newton step first
  lambda2 = lambda_tmp = 0.0;
  for (int iter = 0; iter < kIterMax; ++iter) {
    // Start of iteration loop
    *x = xold + lambda * direction;
    *f = func(*x);

    if (lambda < stepmin) {  // x is too close to xold, ignored
      *x = xold;
      check = true;
      return check;
    } else if (*f <= fold + kAlpha * lambda * slope)
      return check;
    else if (IsInfOrNaN(*f)) {
      while (!IsInfOrNaN(lambda) && IsInfOrNaN(*f)) {
        lambda *= 0.5;
        *x = xold + lambda * direction;
        *f = func(*x);
      }

      lambda_tmp = 0.5 * lambda;
    } else {                // Backtrack
      if (lambda == 1.0) {  // First time
        lambda_tmp = -slope / (2.0 * ((*f) - fold - slope));
      } else {  // Subsequent backtracks
        double rhs1 = *f - fold - lambda * slope;
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
          if (disc < 0.0)
            lambda_tmp = 0.5 * lambda;
          else if (b <= 0.0)
            lambda_tmp = (-b + sqrt(disc)) / (3.0 * a);
          else
            lambda_tmp = -slope / (b + sqrt(disc));
        }
        if (lambda_tmp > 0.5 * lambda || IsInfOrNaN(lambda_tmp)) {
          lambda_tmp = 0.5 * lambda;
        }
      }
    }
    lambda2 = lambda;
    f2 = *f;
    lambda = std::max(lambda_tmp, 0.1 * lambda);
  }  // for
  return check;
}

template <typename T>
bool LineSearch<T>::IsInfOrNaN(double x) {
  return (x == std::numeric_limits<double>::infinity() ||
          x == -std::numeric_limits<double>::infinity() || x != x);
}
