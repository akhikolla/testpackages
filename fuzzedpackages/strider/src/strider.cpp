// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
#include <algorithm>
using std::for_each;
using std::transform;
using std::accumulate;
using std::begin;
using std::plus;
using std::end;
#include <strider.h>
using strider::make_strided;
using strider::make_strided_range;

//' Fast row sums
//' 
//' Demonstration of fast row and columns sums in C++
//' 
//' @param x a numeric matrix
//' 
//' @details A very efficient row summing algorithm that demonstrates
//' the use of the strided pointer concept. The \code{row_sum} algorithm is
//' roughly twice as fast as \code{\link{rowSums}}. The \code{col_sum} algorithm 
//' matches \code{\link{colSums}} for speed.
//'
//' @examples
//' row_sums(matrix(1:9, 3))
//' col_sums(matrix(1:9, 3))
//' 
//' @seealso \code{\link{rowSums}}, \code{\link{colSums}}
//' @rdname sums
//' @export
// [[Rcpp::export]]
NumericVector row_sums(const NumericMatrix& x)
{
  auto nr = x.nrow();
  NumericVector res(nr, 0.0);
  for_each(make_strided(begin(x), nr), make_strided(end(x)), [&](const double& y) {
    transform(&y, &y + nr, begin(res), begin(res), plus<double>());               });
  return res;
}

//' @rdname sums
//' @export
// [[Rcpp::export]]
NumericVector col_sums(const NumericMatrix& x)
{
  auto nr = x.nrow();
  NumericVector res(x.ncol());
  transform(make_strided(begin(x), nr), make_strided(end(x), nr), begin(res),
            [&](const double& v){ return accumulate(&v, &v + nr, 0.0); });
  return res;
}

//' Convolve Matrices
//' 
//' Demonstration of fast matrix convolution C++
//' 
//' @param a a numeric matrix
//' @param b a numeric matrix
//' 
//' @details A very efficient matric convolution implementation that demonstrates
//' the use of the strided pointer and strided range concepts. Performance will be
//' improved if the smaller matrix is given as the second argument.
//'
//' @examples
//' a = matrix(c(1, 2, 1,
//'              1, 1, 1), 2, 3, byrow = TRUE)
//' b = matrix(c(0, 0, 0,
//'              0, 0, 0,
//'              0, 1, 0,
//'              0, 0, 0), 4, 3, byrow = TRUE)
//' convolve2(a, b)
//' 
//' @seealso \code{\link{convolve}}
//' @rdname convolve
//' @export
// [[Rcpp::export]]
NumericMatrix
convolve2(const NumericMatrix& a,
          const NumericMatrix& b)
{
  int
    nra = a.nrow(), nca = a.ncol(),
    nrb = b.nrow(), ncb = b.ncol(),
    nrab = nra + nrb - 1,
    ncab = nca + ncb - 1;
  NumericMatrix ab(nrab, ncab);
  transform(make_strided(begin(a), nra), make_strided(end(a)),
            make_strided(begin(ab), nrab), make_strided(begin(ab), nrab),
            [&](const double& t, double& u) {
              transform(make_strided(begin(b), nrb), make_strided(end(b)),
                        make_strided(&u, nrab), make_strided(&u, nrab),
                        [&](const double& v, double& w) {
                          transform(&t, &t + nra, &w, &w,
                                    [&](const double x, double& y) {
                                      transform(&v, &v + nrb, &y, &y,
                                                [&](const double z, const double zz) {
                                                  return zz + x * z; });
                                      return y; });
                          return w; });
              return u; });
  return ab;
}

