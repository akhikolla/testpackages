#include "RcppArmadillo.h"
// [[Rcpp::export]]

arma::vec CQRPMMCPP (arma::mat xr,arma::vec yr,arma::vec betar,arma::vec betaoldr,
double to,int m,arma::vec ta,double lamdar) {

double toler = (to),lamda=(lamdar);
int maxit=(m);
arma:: vec tau=(ta);
arma:: mat x=(xr),r,product,xt,denominator;
arma:: vec W,uv,v,y=(yr),delta,E;
arma:: vec betaold,beta=(betar),quantile,u,yh,qrbeta=(betaoldr);
arma::uvec order, index;


int n=x.n_rows,k=tau.n_elem;
int p=x.n_cols;
double error=10000,epsilon=0.9999,epsilonprime;
//toler=1e-2;
epsilonprime=toler*p/2;
int iteration=1;
u.zeros(k);
r.zeros(n,k);

//beta= solve(x.t()*x, x.t()*y);
product.ones(p,n);
xt=x.t();



while (iteration<=maxit&& error>toler)
{
  betaold=beta;
  yh=x*beta;
  uv=arma::sort(y-yh);

  // u is vec of the quantiles of given vector
  quantile=(n-1)*tau-floor((n-1)*tau);

	for (int i=0;i<k;i++)
  	{	
		u(i)=quantile(i)*uv(ceil((n-1)*tau(i)))+(1-quantile(i))*uv(floor((n-1)*tau(i))); 
	}

	for (int i=0;i<k;i++)
  	{       r.col(i)=y-u(i)-yh;    
 
	}
        denominator=1/(arma::abs(r)+epsilon);
	W=arma::sum(denominator,1);
        v=k-2*arma::sum(tau)-arma::sum(r%denominator,1);
	E=lamda*arma::abs(beta)/beta/(epsilonprime+arma::abs(beta))/qrbeta/qrbeta;
        for (int i=0;i<n;i++)
		{	product.col(i)=xt.col(i)*W(i);}

 	delta=arma::solve(product*x-diagmat(E),x.t()*v-E%beta); 
  
	beta=beta-delta;
error=arma::sum(arma::abs(delta));
iteration++;
}

return beta;


}
