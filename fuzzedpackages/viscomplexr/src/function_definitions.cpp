// R package viscomplexr - phase portraits of functions in the
// complex number plane
// Copyright (C) 2020 Peter Biber
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>



#include <Rcpp.h>
using namespace Rcpp;


//' Mandelbrot iteration with a given number of steps
//'
//' This function is provided as a basis for visualizing the Mandelbrot set with
//' \code{\link{phasePortrait}}. While usual visualizations color the points
//' \emph{outside} the Mandelbrot set dependent on the velocity of divergence,
//' this function produces the information required for coloring the Mandelbrot
//' set itself. For numbers that can be identified as not being elements of the
//' Mandelbrot set, we obtain a \code{NaN+NaNi} value; for all other numbers,
//' the function gives back the value after a user-defined number of iterations.
//' The function has been implemented in C++; it runs fairly fast.
//'
//' The Mandelbrot set comprises all complex numbers \code{z} for which the
//' sequence \code{a[n+1] = a[n]^2 + z} starting with \code{a[0] = 0} remains
//' bounded for all \code{n > 0}. This condition is certainly not true, if, at
//' any time, \code{abs(a[]) >= 2}. The function \code{mandelbrot} performs the
//' iteration for \code{n = 0, ..., itDepth - 1} and permanently checks for
//' \code{abs(a[n+1]) >= 2}. If this is the case, it stops the iteration and
//' returns \code{NaN+NaNi}. In all other cases, it returns \code{a[itDepth]}.
//'
//' @param z Complex number; the point in the complex plane to which the output
//'   of the function is mapped
//'
//' @param itDepth An integer which defines the depth of the iteration, i.e. the
//'   maximum number of iteration (default: \code{itDepth =  500})
//'
//' @return Either \code{NaN+NaNi} or the complex number obtained after
//'   \code{itDepth} iterations
//'
//' @family fractals
//' @family maths
//'
//' @examples
//' # This code shows the famous Mandelbrot figure in total, just in the
//' # opposite way as usual: the Mandelbrot set itself is colored, while the
//' # points outside are uniformly black.
//' # Adjust xlim and ylim to zoom in wherever you like.
//' \donttest{
//' phasePortrait(mandelbrot,
//'   xlim = c(-2.3, 0.7),
//'   ylim = c(-1.2, 1.2),
//'   hsvNaN = c(0, 0, 0),
//'   nCores = 1)          # Max. two cores on CRAN, not a limit for your use
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//' @export
// [[Rcpp::export]]
std::complex<double> mandelbrot(std::complex<double> z,
                                int itDepth = 500) {
  std::complex<double> zz = 0;
  for(int i = 0; i < itDepth; ++i) {
    zz = pow(zz, 2) + z;
    if((pow(real(zz), 2) + pow(imag(zz), 2)) >= 4) {
      zz = std::complex<double>(NAN, NAN);
      break;
    }
  }
  return zz;
}


//' Julia iteration with a given number of steps
//'
//' This function is designed as the basis for visualizing normal Julia sets
//' with \code{\link{phasePortrait}}. In contrast to usual visualizations of
//' Julia sets, this requires coloring the actual member points of the set and
//' not the points outside. Therefore, for numbers that can be identified as not
//' being parts of the Julia set, this function returns \code{NaN+NaNi}. All
//' other numbers are mapped to the complex value obtained after a user-defined
//' number of iterations. This function has been implemented in C++; therefore
//' it is fairly fast.
//'
//' Normal Julia sets are closely related to the Mandelbrot set. A normal Julia
//' set comprises all complex numbers \code{z} for which the following sequence
//' is bounded for all \code{n > 0}: \code{a[n+1] = a[n]^2 + c}, starting with
//' \code{a[0] = z}. The parameter \code{c} is a complex number, and the
//' sequence is certainly unbounded if \code{abs(a[]) >= R} with \code{R} being
//' an escape Radius which matches the inequality \code{R^2 - R >= abs(c)}. As
//' the visualization with this package gives interesting pictures (i.e. other
//' than a blank screen) only for \code{c} which are elements of the Mandelbrot
//' set, \code{R = 2} is a good choice. For the author's taste, the Julia
//' visualizations become most interesting for \code{c} located in the border
//' zone of the Mandelbrot set.
//'
//' @param z Complex number; the point in the complex plane to which the output
//'   of the function is mapped
//'
//' @param c Complex number; a parameter whose choice has an enormous effect on
//'   the shape of the Julia set. For obtaining useful results with
//'   \code{\link{phasePortrait}}, \code{c} must be an element of the Mandelbrot
//'   set.
//'
//' @param R_esc Real number; the espace radius. If the absolute value of a
//'   number obtained during iteration attains or excels the value of
//'   \code{R_esc}, \code{juliaNormal} will return \code{NaN+NaNi}. \code{R_esc
//'   = 2} is a good choice for \code{c} being an element of the Mandelbrot set.
//'   See Details for more information.
//'
//' @param itDepth An integer which defines the depth of the iteration, i.e. the
//'   maximum number of iteration (default: \code{itDepth =  500})
//'
//' @return Either \code{NaN+NaNi} or the complex number obtained after
//'   \code{itDepth} iterations
//'
//' @family fractals
//' @family maths
//'
//' @examples
//' # This code visualizes a Julia set with some appeal (for the author's
//' # taste). Zoom in as you like by adjusting xlim and ylim.
//' \donttest{
//' phasePortrait(juliaNormal,
//'   moreArgs = list(c = -0.09 - 0.649i, R_esc = 2),
//'   xlim = c(-2, 2),
//'   ylim = c(-1.3, 1.3),
//'   hsvNaN = c(0, 0, 0),
//'   nCores = 1)          # Max. two cores on CRAN, not a limit for your use
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//'
//' @export
// [[Rcpp::export]]
std::complex<double> juliaNormal(std::complex<double> z,
                                 std::complex<double> c,
                                 double R_esc,
                                 int itDepth = 500) {

  // squared escape radius
  // double Rq     = pow(1/2 * (1 + sqrt(1 + 4 * std::abs(c))), 2);
  double Rq     = pow(R_esc, 2);

  std::complex<double> zz = z;
  for(int i = 0; i < itDepth; ++i) {
    zz = pow(zz, 2) + c;
    if((pow(real(zz), 2) + pow(imag(zz), 2)) >= Rq) {
      zz = std::complex<double>(NAN, NAN);
      break;
    }
  }
  return zz;
}


//' Calculate Blaschke products
//'
//' This function calculates Blaschke products
//' (\url{https://en.wikipedia.org/wiki/Blaschke_product}) for a complex number
//' \code{z} given a sequence \code{a} of complex numbers inside the unit disk,
//' which are the zeroes of the Blaschke product.
//'
//' A sequence of points \code{a[n]} located inside the unit disk satisfies the
//' Blaschke condition, if \code{sum[1:n] (1 - abs(a[n])) < Inf}. For each
//' element \code{a != 0} of such a sequence, \code{B(a, z) = abs(a)/a * (a -
//' z)/(1 - conj(a) * z)} can be calculated. For \code{a = 0}, \code{B(a, z) =
//' z}. The Blaschke product \code{B(z)} results as \code{B(z) = prod[1:n]
//' (B(a[n], z))}.
//'
//' @param z Complex number; the point in the complex plane to which the output
//'   of the function is mapped
//'
//' @param a Vector of complex numbers located inside the unit disk. At each
//'   \code{a}, the Blaschke product will have a zero.
//'
//' @return The value of the Blaschke product at \code{z}.
//'
//' @family maths
//'
//' @examples
//' # Generate random vector of 17 zeroes inside the unit disk
//' n <- 17
//' a <- complex(modulus = runif(n, 0, 1), argument = runif(n, 0, 2*pi))
//' \donttest{
//' # Portrait the Blaschke product
//' phasePortrait(blaschkeProd, moreArgs = list(a = a),
//'   xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
//'   nCores = 1) # Max. two cores on CRAN, not a limit for your use
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//' @export
// [[Rcpp::export]]
std::complex<double> blaschkeProd(std::complex<double> z,
                                  std::vector<std::complex<double>> a) {
  int n = a.size();
  std::complex<double> zz = 1;
  std::complex<double> fact;

  for(int i = 0; i < n; ++i) {
    if(std::abs(a[i]) != 0) {
      fact = std::abs(a[i]) / a[i] * (a[i] - z) /
        (std::complex<double>(1, 0) - z * std::conj(a[i]));
    }
    else {
      fact = z;
    }
    zz = zz * fact;
  }
  return zz;
}


//' Jacobi theta function
//'
//' Approximation of "the" Jacobi theta function using the first \code{nn}
//' factors in its triple product version
//'
//' This function approximates the Jacobi theta function theta(z; tau) which is
//' the sum of exp(pi*i*n^2*tau + 2*pi*i*n*z) for n in -Inf, Inf. It uses,
//' however, the function's triple product representation. See
//' \url{https://en.wikipedia.org/wiki/Theta_function} for details. This function
//' has been implemented in C++, but it is only slightly faster than well-crafted
//' R versions, because the calculation can be nicely vectorized in R.
//'
//' @param z Complex number; the point in the complex plane to which the output
//'   of the function is mapped
//'
//' @param tau Complex number; the so-called half-period ratio, must have a
//'   positive imaginary part
//'
//' @param nn Integer; number of factors to be used when approximating the
//'   triple product (default = 30)
//'
//' @return The value of the function for \code{z} and \code{tau}.
//'
//' @family maths
//'
//'
//' @examples
//' \donttest{
//' phasePortrait(jacobiTheta, moreArgs = list(tau = 1i/2-1/4),
//' pType = "p", xlim = c(-2, 2), ylim = c(-2, 2),
//' nCores = 1) # Max. two cores on CRAN, not a limit for your use
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//' \donttest{
//' phasePortrait(jacobiTheta, moreArgs = list(tau = 1i/2-1/2),
//' pType = "p", xlim = c(-2, 2), ylim = c(-2, 2),
//' nCores = 1)
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//' \donttest{
//' phasePortrait(jacobiTheta, moreArgs = list(tau = 1i/3+1/3),
//' pType = "p", xlim = c(-2, 2), ylim = c(-2, 2),
//' nCores = 1)
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//' \donttest{
//' phasePortrait(jacobiTheta, moreArgs = list(tau = 1i/4+1/2),
//' pType = "p", xlim = c(-2, 2), ylim = c(-2, 2),
//' nCores = 1)
//'   \dontshow{
//'   # R CMD check: make sure any open connections are closed afterward
//'   foreach::registerDoSEQ()
//'   doParallel::stopImplicitCluster()
//'   }
//' }
//'
//'
//' @export
// [[Rcpp::export]]
std::complex<double> jacobiTheta(std::complex<double> z,
                                 std::complex<double> tau,
                                 int nn = 30) {

  std::complex<double> one   = std::complex<double>(1, 0);
  std::complex<double> I     = std::complex<double>(0, 1);
  std::complex<double> pi    = std::complex<double>(M_PI, 0);
  std::complex<double> two   = std::complex<double>(2, 0);
  std::complex<double> theta = std::complex<double>(1, 0);

  for(int k = 0; k < nn; k++) {
    std::complex<double> n = std::complex<double>(k + 1, 0);
    theta *= (one - exp(two*pi*I*n*tau)) *
      (one + exp(pi*I*((two*n - one)*tau + two*z))) *
      (one + exp(pi*I*((two*n - one)*tau - two*z)));
  }

  return theta;
}








