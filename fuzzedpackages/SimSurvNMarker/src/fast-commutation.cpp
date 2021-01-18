#include <Rcpp.h>

void get_commutation_unequal
  (unsigned const n, unsigned const m, double * o) noexcept{
  unsigned const nm = n * m,
             nnm_p1 = n * nm + 1L,
              nm_pm = nm + m;

  for(unsigned i = 0; i < n; ++i, o += nm_pm){
    double *o1 = o;
    for(unsigned j = 0; j < m; ++j, o1 += nnm_p1)
      *o1 = 1.;
  }
}

void get_commutation_equal(unsigned const m, double * const o) noexcept {
  unsigned const mm = m * m,
                mmm = mm * m,
             mmm_p1 = mmm + 1L,
              mm_pm = mm + m;

  unsigned inc_i(0L);
  for(unsigned i = 0; i < m; ++i, inc_i += m){
    double *o1 = o + inc_i + i * mm,
           *o2 = o + i     + inc_i * mm;
    for(unsigned j = 0; j < i; ++j, o1 += mmm_p1, o2 += mm_pm){
      *o1 = 1.;
      *o2 = 1.;
    }
    *o1 += 1.;
  }
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_commutation(unsigned const n, unsigned const m){
  size_t const nm = n * m;
  Rcpp::NumericMatrix out(nm, nm);
  if(n == m)
    get_commutation_equal(n, &out[0]);
  else
    get_commutation_unequal(n, m, &out[0]);

  return out;
}

/*** R
options(digits = 3)

library(matrixcalc)
get_commutation_R <- function(m){
  out <- matrix(0., nrow = m * m, ncol = m * m)
  for(i in 1:m)
    for(j in 1:m)
      out[(i - 1L) * m + j, (j - 1L) * m + i] <- 1.

    return(out)
}

# equal
for(i in 2:10){
  stopifnot(all.equal(commutation.matrix(i), get_commutation_R(i)))
  stopifnot(all.equal(commutation.matrix(i), get_commutation  (i, i)))
}

# unequal
for(i in 3:10)
  for(j in 2:(i - 1L)){
    stopifnot(all.equal(commutation.matrix(i, j), get_commutation  (i, j)))
    stopifnot(all.equal(commutation.matrix(j, i), get_commutation  (j, i)))
  }

# benchmark: equal
library(microbenchmark)
microbenchmark(
  matrixcalc = commutation.matrix(4L),
  R          = get_commutation_R (4L),
  cpp        = get_commutation   (4L, 4L),
  times = 1000)
#R> Unit: nanoseconds
#R>        expr    min     lq   mean median     uq     max neval
#R>  matrixcalc 403333 422449 461087 433004 444748 3127575  1000
#R>           R   4205   4661   5330   5224   5624   54460  1000
#R>         cpp    638    848   1288   1264   1510   20620  1000

microbenchmark(
  matrixcalc = commutation.matrix(20L),
  R          = get_commutation_R (20L),
  cpp        = get_commutation   (20L, 20L),
  times = 25)
#R> Unit: microseconds
#R>        expr      min       lq   mean   median       uq    max neval
#R>  matrixcalc 531256.1 565798.3 571792 572429.9 579425.0 622050    25
#R>           R    159.4    172.4    184    179.4    183.1    266    25
#R>         cpp     39.1     44.3    103     48.4     52.9   1381    25

# benchmark: unequal
microbenchmark(
  matrixcalc = commutation.matrix(17L, 20L),
  cpp        = get_commutation   (17L, 20L),
  times = 25)
#R> Unit: microseconds
#R>        expr    min     lq   mean   median     uq    max neval
#R>  matrixcalc 318443 324061 338861 329943.3 358941 363992    25
#R>         cpp     29     37    140     42.1     46   1368    25
*/
