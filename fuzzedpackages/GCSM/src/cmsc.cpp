#include <RcppArmadillo.h>
using namespace arma;

//' @rdname gcsm
//' @export
// [[Rcpp::export]]
double cmsc(arma::vec x, arma::vec y, bool rescale = false,
            double xmin = NA_REAL, double xmax = NA_REAL,
            double ymin = NA_REAL, double ymax = NA_REAL,
            std::string comp = "si") {
  if (comp != "si" && comp != "s1" && comp != "s2" && comp != "s3")
    Rcpp::stop("comp should be 'si' or 's1', 's2', 's3'!");
  if (x.has_nan() || y.has_nan()) {
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
    if (find_finite(x).is_empty()) Rcpp::stop("x and y have no valid overlap!");
  }
  if (!std::isfinite(xmin)) xmin = x.min();
  if (!std::isfinite(xmax)) xmax = x.max();
  if (!std::isfinite(ymin)) ymin = y.min();
  if (!std::isfinite(ymax)) ymax = y.max();
  if (xmin > xmax) Rcpp::stop("xmin > xmax, please reset them!");
  if (ymin > ymax) Rcpp::stop("ymin > ymax, please reset them!");
  if (xmax < x.min() || xmin > x.max()) Rcpp::stop("[xmin, xmax] is beyond the range of x!");
  if (ymax < y.min() || ymin > y.max()) Rcpp::stop("[ymin, ymax] is beyond the range of x!");
  double xymin = (xmin < ymin) ? xmin : ymin;
  double xymax = (xmax > ymax) ? xmax : ymax;
  if (xymin == xymax) Rcpp::stop("global min equals to global max!");
  double maxmmin = xymax - xymin;
  if (rescale) {
    // avoid zero-denominator
    if (xmin == xmax) x.elem(find_finite(x)).ones(); else x = (x-xmin) / (xmax-xmin);
    if (ymin == ymax) y.elem(find_finite(y)).ones(); else y = (y-ymin) / (ymax-ymin);
    maxmmin = 1;
  }

  uvec inx = find_finite(x);
  arma::vec xsub = vectorise(x.elem(inx));
  arma::vec ysub = vectorise(y.elem(inx));
  double d2 = pow((stddev(xsub) - stddev(ysub)) / (maxmmin / 2), 2);
  // special cases e.g., stddev(x) is near to 0 and stddev(y) > maxmmin/2
  // because the denominator of the standard deviation is n-1.
  if (d2 > 1) d2 = 1;
  double cc = as_scalar(cor(xsub, ysub));
  // special cases when standard deviation of x or y is zero.
  if (!std::isfinite(cc)) cc = (d2 == 0) ? 1 : 0;
  double d1 = pow((mean(xsub) - mean(ysub)) / maxmmin, 2);
  // special cases when the specified rescale range is narrow
  // because the raw data would not be clamped
  if (d1 > 1) d1 = 1;
  double s1 = 1 - d1;
  double s2 = 1 - d2;
  double s3 = cc;
  if (comp == "si") return s1 * s2 * s3;
  if (comp == "s1") return s1;
  if (comp == "s2") return s2;
  if (comp == "s3") return s3;
  // unnecessary
  Rcpp::stop("comp should be 'si' or 's1', 's2', 's3'!");
}


//' @rdname gcsm_sw
//' @export
// [[Rcpp::export]]
arma::mat cmsc_sw(arma::mat x, arma::mat y, bool rescale = false,
                  double xmin = NA_REAL, double xmax = NA_REAL,
                  double ymin = NA_REAL, double ymax = NA_REAL,
                  double ksize = 9, bool globe = false,
                  std::string comp = "si") {
  if (comp != "si" && comp != "s1" && comp != "s2" && comp != "s3")
    Rcpp::stop("comp should be 'si' or 's1', 's2', 's3'!");
  if (x.has_nan() || y.has_nan()) {
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
    if (find_finite(x).is_empty()) Rcpp::stop("x and y have no valid overlap!");
  }
  if (!std::isfinite(xmin)) xmin = x.min();
  if (!std::isfinite(xmax)) xmax = x.max();
  if (!std::isfinite(ymin)) ymin = y.min();
  if (!std::isfinite(ymax)) ymax = y.max();
  if (xmin > xmax) Rcpp::stop("xmin > xmax, please reset them!");
  if (ymin > ymax) Rcpp::stop("ymin > ymax, please reset them!");
  if (xmax < x.min() || xmin > x.max()) Rcpp::stop("[xmin, xmax] is beyond the range of x!");
  if (ymax < y.min() || ymin > y.max()) Rcpp::stop("[ymin, ymax] is beyond the range of x!");
  double xymin = (xmin < ymin) ? xmin : ymin;
  double xymax = (xmax > ymax) ? xmax : ymax;
  if (xymin == xymax) Rcpp::stop("global min equals to global max!");
  double maxmmin = xymax - xymin;
  if (rescale) {
    // avoid zero-denominator
    if (xmin == xmax) x.elem(find_finite(x)).ones(); else x = (x-xmin) / (xmax-xmin);
    if (ymin == ymax) y.elem(find_finite(y)).ones(); else y = (y-ymin) / (ymax-ymin);
    maxmmin = 1;
  }

  arma::mat s1 = arma::mat(size(x)).fill(datum::nan);
  arma::mat s2 = arma::mat(size(x)).fill(datum::nan);
  arma::mat s3 = arma::mat(size(x)).fill(datum::nan);
  int rad = floor(0.5 * ksize);
  // pad cols
  arma::mat xheadcols = x.head_cols(rad);
  x.insert_cols(0, x.tail_cols(rad));
  x.insert_cols(x.n_cols, xheadcols);
  arma::mat yheadcols = y.head_cols(rad);
  y.insert_cols(0, y.tail_cols(rad));
  y.insert_cols(y.n_cols, yheadcols);
  if (!globe) {
    x.head_cols(rad).fill(datum::nan);
    x.tail_cols(rad).fill(datum::nan);
    y.head_cols(rad).fill(datum::nan);
    y.tail_cols(rad).fill(datum::nan);
  }
  // pad rows
  arma::mat nanrows = arma::mat(rad, x.n_cols).fill(datum::nan);
  x.insert_rows(0, nanrows);
  x.insert_rows(x.n_rows, nanrows);
  y.insert_rows(0, nanrows);
  y.insert_rows(y.n_rows, nanrows);
  uvec locs = find_finite(x);
  locs = locs.elem(find(locs >= rad * x.n_rows && locs <= x.n_elem - 1 - rad * x.n_rows));
  int nlocs = locs.n_elem;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int loc = 0; loc < nlocs; ++loc) {
    // decalre inside as private variables for parallel
    // "%" here is modulo operator, not elewise multiplication
    int col = locs(loc) / x.n_rows;
    int row = locs(loc) % x.n_rows;
    arma::mat xsubmat = x.submat(row-rad, col-rad, row+rad, col+rad);
    arma::mat ysubmat = y.submat(row-rad, col-rad, row+rad, col+rad);
    uvec inx = find_finite(xsubmat);
    arma::vec xsub = vectorise(xsubmat.elem(inx));
    arma::vec ysub = vectorise(ysubmat.elem(inx));
    double d2 = pow((stddev(xsub) - stddev(ysub)) / (maxmmin / 2), 2);
    // special cases e.g., stddev(x) is near to 0 and stddev(y) > maxmmin/2
    // because the denominator of the standard deviation is n-1.
    if (d2 > 1) d2 = 1;
    double cc = as_scalar(cor(xsub, ysub));
    // special cases when standard deviation of x or y is zero.
    if (!std::isfinite(cc)) cc = (d2 == 0) ? 1 : 0;
    double d1 = pow((mean(xsub) - mean(ysub)) / maxmmin, 2);
    // special cases when the specified rescale range is narrow
    // because the raw data would not be clamped
    if (d1 > 1) d1 = 1;
    s1(row-rad, col-rad) = (1 - d1);
    s2(row-rad, col-rad) = (1 - d2);
    s3(row-rad, col-rad) = cc;
  }
  if (comp == "si") return s1 % s2 % s3;
  if (comp == "s1") return s1;
  if (comp == "s2") return s2;
  if (comp == "s3") return s3;
  // unnecessary
  Rcpp::stop("comp should be 'si' or 's1', 's2', 's3'!");
}


//' @rdname gcsm_tw
//' @export
// [[Rcpp::export]]
arma::mat cmsc_tw(arma::cube xxx, arma::cube yyy, bool rescale = false,
                  double xmin = NA_REAL, double xmax = NA_REAL,
                  double ymin = NA_REAL, double ymax = NA_REAL,
                  std::string comp = "si") {
  if (comp != "si" && comp != "s1" && comp != "s2" && comp != "s3")
    Rcpp::stop("comp should be 'si' or 's1', 's2', 's3'!");
  // make sure x and y are not altered, use const make nonsense since they are SEXP, i.e. pointer
  // see [Rcpp FAQ 5.1](https://cloud.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf)
  // and discussion in https://github.com/RcppCore/RcppArmadillo/pull/109#issuecomment-255484817
  cube x = xxx;
  cube y = yyy;
  if (x.has_nan() || y.has_nan()) {
    x.elem(find_nonfinite(y)).fill(datum::nan);
    y.elem(find_nonfinite(x)).fill(datum::nan);
    if (find_finite(x).is_empty())Rcpp::stop("x and y have no valid overlap!");
  }
  if (!std::isfinite(xmin)) xmin = x.min();
  if (!std::isfinite(xmax)) xmax = x.max();
  if (!std::isfinite(ymin)) ymin = y.min();
  if (!std::isfinite(ymax)) ymax = y.max();
  if (xmin > xmax) Rcpp::stop("xmin > xmax, please reset them!");
  if (ymin > ymax) Rcpp::stop("ymin > ymax, please reset them!");
  if (xmax < x.min() || xmin > x.max()) Rcpp::stop("[xmin, xmax] is beyond the range of x!");
  if (ymax < y.min() || ymin > y.max()) Rcpp::stop("[ymin, ymax] is beyond the range of x!");
  double xymin = (xmin < ymin) ? xmin : ymin;
  double xymax = (xmax > ymax) ? xmax : ymax;
  if (xymin == xymax) Rcpp::stop("global min equals to global max!");
  double maxmmin = xymax - xymin;
  if (rescale) {
    // avoid zero-denominator
    if (xmin == xmax) x.elem(find_finite(x)).ones(); else x = (x-xmin) / (xmax-xmin);
    if (ymin == ymax) y.elem(find_finite(y)).ones(); else y = (y-ymin) / (ymax-ymin);
    maxmmin = 1;
  }

  arma::mat s1 = arma::mat(x.n_rows, x.n_cols).fill(datum::nan);
  arma::mat s2 = arma::mat(x.n_rows, x.n_cols).fill(datum::nan);
  arma::mat s3 = arma::mat(x.n_rows, x.n_cols).fill(datum::nan);
  int nlocs = x.slice(0).n_elem;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int loc = 0; loc < nlocs; ++loc) {
    // decalre inside as private variables for parallel
    // "%" here is modulo operator, not elewise multiplication
    int col = loc / x.n_rows;
    int row = loc % x.n_rows;
    arma::vec xtube = x.tube(row, col);
    arma::vec ytube = y.tube(row, col);
    uvec inx = find_finite(xtube);
    if (inx.is_empty()) continue;
    arma::vec xsub = vectorise(xtube.elem(inx));
    arma::vec ysub = vectorise(ytube.elem(inx));
    double d2 = pow((stddev(xsub) - stddev(ysub)) / (maxmmin / 2), 2);
    // special cases e.g., stddev(x) is near to 0 and stddev(y) > maxmmin/2
    // because the denominator of the standard deviation is n-1.
    if (d2 > 1) d2 = 1;
    double cc = as_scalar(cor(xsub, ysub));
    // special cases when standard deviation of x or y is zero.
    if (!std::isfinite(cc)) cc = (d2 == 0) ? 1 : 0;
    double d1 = pow((mean(xsub) - mean(ysub)) / maxmmin, 2);
    // special cases when the specified rescale range is narrow
    // because the raw data would not be clamped
    if (d1 > 1) d1 = 1;
    s1(row, col) = 1 - d1;
    s2(row, col) = 1 - d2;
    s3(row, col) = cc;
  }
  if (comp == "si") return s1 % s2 % s3;
  if (comp == "s1") return s1;
  if (comp == "s2") return s2;
  if (comp == "s3") return s3;
  // unnecessary
  Rcpp::stop("comp should be 'si' or 's1', 's2', 's3'!");
}
