#include <RcppArmadillo.h>
#include <Rmath.h>
#include <Rcpp.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double cont(arma::vec x,  arma::vec y, const arma::vec& w) {
  double sumw = arma::sum(w);
  double xb = arma::sum(w%x)/sumw;
  double yb = arma::sum(w%y)/sumw;
  const arma::vec temp1 = x-xb;
  const arma::vec temp2 = y-yb;
  double numerator = arma::sum(w%temp1%temp2);
  double denominator = pow(arma::sum(w%arma::square(temp1)) * sum(w%arma::square(temp2)), 0.5);
  return numerator/denominator;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec wrankFast(arma::vec x, const arma::vec& w) {
  const arma::uvec ord = arma::sort_index(x);
  int size = ord.size();
  arma::vec temp = arma::linspace(0, size-1, size);
  const arma::uvec ind = arma::sort_index(ord);
  
  temp = temp.elem(ind);
  arma::uvec rord = arma::conv_to<arma::uvec>::from(temp);
  const arma::vec xp = x.elem(ord);
  const arma::vec wp = w.elem(ord);
  arma::vec rnk(size);
  double t1 = 0;
  int i = 0;
  double t2 = 0;
  double n = 0 ;
  while(i < size-1) {
    
    t2 = t2 + wp[i];
    n = n + 1;
    if(xp[i+1] != xp[i]) {
      double rnki = t1 + (n+1)/(2*n)*t2;
      for(int ii =0; ii < n; ii++) {
        rnk[i-ii] = rnki;
        
      }
      t1 = t1 + t2; 
      t2 = 0;
      n = 0;
    }
    i = i + 1;
  }
  n = n + 1;
  t2 = t2 + wp[i];
  double rnki = t1 + (n+1)/(2*n)*t2; 
  for(int ii = 0; ii < n; ii++) {
      rnk[i-ii] = rnki;
    }
    rnk = rnk.elem(rord);
  return(rnk);
}

// arma::vec wrankFastOld(arma::vec x, const arma::vec& w) {
//   int size = x.size();
//   arma::vec ranks(size);
//   for(int i =0; i < size; i++ ) {
//     const arma::uvec ind = arma::find(x<x[i]);
//     double t1 = arma::sum(w.elem(ind));
//     const arma::uvec ind1 = arma::find(x==x[i]);
//     arma::vec t2 = w.elem(ind1);
//     ranks[i] = t1 + arma::mean(t2) + (arma::sum(t2) - arma::mean(t2))/2;
//   }
//   return(ranks);
// }
