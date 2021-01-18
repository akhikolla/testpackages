#include "RcppArmadillo.h"
// [[Rcpp::export]]

arma::vec QRPMMCPP (arma::mat xr,arma::vec yr,arma::vec betar,arma::vec betaoldr,
double to,int m,double ta,double lamdar){

double toler = (to);
int maxit=(m);
double tau=(ta),lamda=(lamdar);
arma:: mat x=(xr),product,xt;
arma:: vec W,newX,z,signw,v,r,y=(yr),E,delta;
arma:: vec betaold,beta=(betar),qrbeta=(betaoldr);
arma::uvec order, index;

int n=x.n_rows;
x.insert_cols( 0, arma::ones(n) );
int p=x.n_cols;
double error=10000,epsilon=0.9999,epsilonprime;
//toler=1e-3;
int iteration=1;
//beta= solve(x.t()*x, x.t()*y);
product.ones(p,n);
epsilonprime=toler*p/2;
xt=x.t();
qrbeta.insert_rows(0,1);

while (iteration<=maxit&& error>toler)
{
	betaold=beta;
 
	r=y-x*beta;
	v=1-2*tau-r/(arma::abs(r)+epsilon);
	E=lamda*arma::abs(beta)/beta/(epsilonprime+arma::abs(beta))/qrbeta/qrbeta;
        E(0)=0;
	W=1/(epsilon+abs(r));

	for (int i=0;i<n;i++)
		{	product.col(i)=xt.col(i)*W(i);}

 	delta=arma::solve(product*x-diagmat(E),xt*v-E%beta); 
	beta=beta-delta;


error=sum(abs(delta));
iteration++;
}
arma::vec betanew=beta;
return betanew;

}
