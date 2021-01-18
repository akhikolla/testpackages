#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List make_CDF( int lp,
                     arma::uvec col0_tot,
                     arma::uvec row1_tot,
                     arma::vec min_margin,
                     Rcpp::List ALL_PROBS,
                     bool compute_all_holm,
                     double min_p) {

  Rcpp::List pCDFlist(lp);

  for (int s=0; s<lp; s++) {
    int lo = col0_tot(s) - row1_tot(s);
    if (lo<0) lo=0;
    int hi = min_margin(s);
    IntegerVector support = seq(lo, hi);
    vec ap = Rcpp::as<vec>(ALL_PROBS(s));
    NumericVector all_probs = Rcpp::NumericVector(ap.begin(), ap.end());
    int lsupp = ap.n_elem;
    vec pvals(lsupp);

    uvec idx = Rcpp::as<uvec>(Rcpp::duplicated(all_probs));
    if (accu(ap.elem(find(idx==1))) < 1e-20) {
      pvals = cumsum(sort(ap));
    } else {
      // if (max(abs(ap.rows(0, lsupp/2 - 1)-reverse(ap.rows(lsupp-lsupp/2, lsupp - 1))))< 1e-7) {
      //   vec u = cumsum(ap.rows(0, lsupp/2-1));
      //   IntegerVector ii = seq(0, u.n_elem-1);
      //   uvec idx = vectorise(repmat(Rcpp::as<urowvec>(ii), 2, 1));
      //   if (lsupp%2==0) {
      //     pvals = 2*u.rows(idx);
      //   } else {
      //     pvals.rows(0, lsupp-2) = 2*u.rows(idx);
      //     pvals(pvals.n_elem-1) = 1;
      //   }
      // } else {
        double relErr = 1.0 + 1e-7;
        for (int i=0; i<lsupp; i++) {
          pvals(i) = accu(ap(find(ap <= ap(support(i)-lo)*relErr)));
        }
        pvals=sort(pvals);
      }
    // }

    vec mxmp = pvals(find(pvals<=min_p));
    double mp = 0;
    if (mxmp.n_elem>0) mp = max(mxmp);

    Rcpp::List tempList = Rcpp::List::create(
                              Rcpp::Named( "pvals" ) = pvals,
                              Rcpp::Named( "mp" ) = mp );
    pCDFlist(s) = tempList;
  }

return pCDFlist;
}
