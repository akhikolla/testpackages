#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List l0tf_c(arma::vec y,int k0,int T0,int max_steps)
{
  arma::vec beta=y;
  int n = y.n_elem;
  int m = n-k0-1; 
  arma::vec z = diff(y, k0+1);
  k0 = k0+1;
  arma::vec u(m);
  u.zeros();
  
  double rho = n*n;
  arma::uvec A(T0);
  A=A.fill(0);

  int k=0;
  int k_00=1;

  
  arma::mat D(n-k0,n,arma::fill::zeros);
  arma::mat D_ori(n,n,arma::fill::zeros);
  arma::vec D_temp=arma::ones<arma::vec>(n);
  D_temp=-D_temp;
  
  
  arma::uvec D_index = arma::linspace<arma::uvec>(1,n-1,n-1);
  arma::uvec D_index2 = D_index-1;
  arma::vec D_temp2 = arma::ones<arma::vec>(n-1);
  arma::mat D0(n-1,n-1,arma::fill::zeros);
  D0.diag()=D_temp2;
  D_ori(D_index2,D_index)=D0;
  D_ori.diag() = D_temp;
  
  arma::mat D_save = D_ori;
  while (k_00<k0)
  {
    D_save = D_ori*D_save;
    k_00=k_00+1;
    
  }

  arma::uvec ddindex = arma::linspace<arma::uvec>(1,n-k0,n-k0);
  arma::uvec ddindex2 = arma::linspace<arma::uvec>(1,n,n);
  D=D_save(ddindex-1,ddindex2-1);
  arma::mat dd = D*D.t();
  arma::mat ddinv = dd.i();
  

  while(k < max_steps)
  {
    arma::vec bd=abs(z+u/rho);
    arma::uvec bd_order=sort_index(bd,"descend");
    arma::uvec A_new=bd_order.head(T0);
    arma::uvec I_new=bd_order.tail(n-k0-T0);

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

