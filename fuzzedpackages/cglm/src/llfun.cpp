#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]


double llfun(arma::vec b, arma::vec y, arma::mat design, arma::vec nidcumsum){
  
  design.each_row() %= b.t(); 
  mat designb = design;
  int m = nidcumsum.size()-1;
  double ll = 0;
  for(int i = 0; i<m; ++i){
    span ind = span(nidcumsum(i), nidcumsum(i+1)-1);
    mat designbi = designb(ind, span::all); 
    vec yi = y(ind); 
    double lli = 0;
    double numi = 0;
    int j = 0;  
    vector<int> perm(yi.size());
    for (size_t k=0; k<yi.size(); k++) perm[k] = k;
    uvec permij(yi.size());
    do{ 
      for(size_t k=0; k<yi.size(); k++) permij[k] = perm[k];
      vec yij = yi(permij);  
      mat deni = sum(yij.t()*designbi, 1);
      if(j==0){
        numi = deni(0, 0);
        j = j+1;
      }
      mat llimat = exp(deni-numi);
      lli += llimat(0, 0);
    }while(next_permutation(perm.begin(), perm.end()));  
    ll += log(lli);  
  }
  return ll; 
    
}



 

