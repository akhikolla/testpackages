#include <RcppArmadillo.h>
// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

List tvpcpp(mat x, vec y, mat xe, mat theta, mat E, double lambda, double V, Nullable<double> kappa) {

  int N = x.n_rows;
  int t;
  double yhat, ei, tv, Vu, temp;
  mat R(E.n_rows,E.n_cols,fill::zeros);
  vec pdensi(N,fill::zeros), ytvp(N,fill::zeros);
  mat thetas(theta.n_rows,N+1,fill::zeros);
  mat xx(1,xe.n_cols,fill::zeros);
  mat Rf(E.n_rows,E.n_cols,fill::zeros);
  thetas.col(0) <- theta;

  for (t=0;t!=N;t++)
    {
      xx = xe.row(t);
      xx = xx.t();
      yhat = as_scalar(xx.t() * theta);
      ei = as_scalar(y(t) - yhat);
      R = E / lambda;
      Rf = sqrtmat_sympd(R);
      tv = as_scalar(dot(xx.t() * Rf,xx.t() * Rf));
      Vu = V + tv;
      E = R - (R * xx) *  (xx.t() * R) / Vu;
      
      pdensi(t) = exp(-0.5 * ei * ei / Vu ) / sqrt(2*M_PI*Vu);
      
      theta = theta + ( R * xx ) * ei / Vu;
      
      if (!kappa.isNull())
        {
          V = V * as<double>(kappa) + ((double)1-as<double>(kappa)) * ei * ei;
        }
      else
        {
          temp = ( t * V + (ei * ei - tv) ) / (t+1);
          if (temp>0)
           {
             V = temp;
           }
        }

      ytvp(t) = yhat;
      thetas.col(t+1) = theta;
    }

  return List::create(thetas,ytvp,pdensi);

}
