#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List l0gen_c(arma::vec y, arma::mat D, int T0, int max_steps, arma::mat ddinv)
{
  arma::vec beta = y;
  arma::vec z = D*y;
  int n = y.n_elem;
  int p=z.n_elem;

  double rho = n*n;

  arma::uvec A(T0);
  A=A.fill(0);
  int m=z.n_elem;

  arma::vec u(m);
  u.zeros();

  int k=0;


  while(k < max_steps)
  {
    arma::vec bd=abs(z+u/rho);
    arma::uvec bd_order=sort_index(bd,"descend");
    arma::uvec A_new=bd_order.head(T0);
    arma::uvec I_new=bd_order.tail(p-T0);

    z=z.zeros(m);
    u=z.zeros(m);


    arma::mat di_inv=ddinv(I_new,I_new)-ddinv(I_new,A_new)*inv(ddinv(A_new,A_new))*ddinv(A_new,I_new);
    u(I_new)=di_inv*D.rows(I_new)*y;

    z(A_new)=D.rows(A_new)*(y-D.st()*u);

    A_new=sort(A_new);
    A=sort(A);
    arma::uvec if_true=A_new==A;
    if(sum(if_true)==A.n_elem) break;
    else{
      k=k+1;
      A=A_new;
    }

  }

  beta=y-D.st()*u;

   return List::create(Named("beta")=beta,Named("z")=z,Named("u")=u);

}
