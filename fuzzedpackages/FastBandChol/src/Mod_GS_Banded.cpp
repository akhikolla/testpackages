// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]

List BandCholcpp(arma::mat x,  int k){


  int p = x.n_cols;
  int s = 0;

  arma::mat Q = x;
  arma::mat R = zeros(p, p);


  for(int i=0; i<(p-2); i++){
    if(i+k > p-1){
    	s = p-1;
    } else {
    	s = i+k;
    }
  R(span(i,i), span(i,i)) = norm(Q(span::all,span(i,i)), 2);
  Q(span::all,span(i,i)) = Q(span::all,span(i,i))/as_scalar(R(span(i,i), span(i,i)));
  R(span(i,i), span(i+1,s)) = trans(Q(span::all,span(i,i)))*Q(span::all,span(i+1,s));
  Q(span::all, span(i+1,s)) = Q(span::all,span(i+1, s)) - Q(span::all,span(i,i))*R(span(i,i), span(i+1,s));
  }

  R(span(p-1,p-1),span(p-1,p-1)) = norm(Q(span::all,span(p-1, p-1)), 2);
  Q(span::all,span(p-1, p-1)) = Q(span::all,span(p-1,p-1))/as_scalar(R(span(p-1,p-1), span(p-1,p-1)));


  return List::create(Named("Q") = wrap(Q),
                      Named("R") = wrap(R) );
}









