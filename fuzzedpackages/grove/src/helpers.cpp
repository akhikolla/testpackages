#include <RcppArmadillo.h>
//#include <omp.h>

using namespace Rcpp;
using namespace arma;

const double log2pi = std::log(2.0 * M_PI);

double log_exp_x_plus_exp_y(double x, double y) 
{
  double result;
  if ( x - y >= 100 ) result = x;
  else if ( x - y <= -100 ) result = y;
  else {
    if (x > y) {
      result = y + log( 1 + exp(x-y) );
    }
    else result = x + log( 1 + exp(y-x) );
  }
  return result;
}

arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = true) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) 
{
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

int sampling(int Kf, arma::vec probs)
{  int v;
   Rcpp::NumericVector probsum(Kf);
   double x = R::runif(0.0, sum(probs));
   probsum(0)=probs(0);
   for(int k = 1; k < Kf; k++)
   {  probsum(k)=probs(k)+probsum(k-1);
   }
   if(x < probsum(0)){ v=0;}
   for (int k = 0; k < (Kf-1); k++)
   {  if(x > probsum(k)){ if(x< probsum(k+1)){v=k+1;} }
   }
   if(x > probsum(Kf-1)){v=Kf-1;}
   return v;
}


arma::uvec integerToBinary(int s, int l = 16)
{
  uvec output(l); 
  for(int i=0; i<l; i++)
    output(i) = (s >> i) & 1;
  return output;  
}
