#include "RcppArmadillo.h"
// [[Rcpp::export]]
arma::vec QRMMCPP(arma::mat xr,arma:: vec yr,arma::vec betar,double to,int m,double ta) 
{

double toler =(to);
int maxit=(m);
double tau=(ta);
arma:: mat x=(xr),product,xt;
arma:: vec W,newX,z,signw,v,r,y=(yr);
arma:: vec betaold,beta=(betar),delta;
arma::uvec order, index;

int n=x.n_rows;
x.insert_cols( 0, arma::ones(n) );
int p=x.n_cols;
double error=10000,epsilon=0.9999;
//toler=1e-3;
int iteration=1;
//beta= solve(x.t()*x, x.t()*y);
product.ones(p,n);
xt=x.t();

while (iteration<=maxit&& error>toler)
{
	betaold=beta;
 
	r=y-x*beta;
	v=1-2*tau-r/(arma::abs(r)+epsilon);

	W=1/(epsilon+arma::abs(r));

	for (int i=0;i<n;i++)
		{	product.col(i)=xt.col(i)*W(i);}

 	delta=arma::solve(product*x,xt*v); 
	beta=beta-delta;


error=sum(abs(delta));
iteration++;
}
arma::vec betanew=beta;

return beta;

}

