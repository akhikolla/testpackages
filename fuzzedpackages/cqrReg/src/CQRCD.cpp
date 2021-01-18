#include "RcppArmadillo.h"
// [[Rcpp::export]]

arma::vec CQRCDCPP (arma::mat xr,arma:: vec yr,arma::vec betar,double to,int m,arma::vec ta) 
{

double toler = (to);
int maxit=(m);
arma:: vec tau=(ta);
arma:: mat x=(xr),r,signw,z,newX;
arma:: vec sortz,w,uv,y=(yr),vz,vnewX;
arma:: vec betaold,beta=(betar),quantile,u,yh;
arma::uvec order, index;


int place,n=x.n_rows,k=tau.n_elem;
int p=x.n_cols;
double error=10000;
//toler=1e-2;
int iteration=1;
u.zeros(k);
r.zeros(n,k);
signw.zeros(n,k);
z.zeros(n,k);
newX.zeros(n,k);

//beta= solve(x, y);;



while (iteration<=maxit&& error>toler)
{
  betaold=beta;
  uv=arma::sort(y-x*beta);

  // u is vec of the quantiles of given vector
  quantile=(n-1)*tau-floor((n-1)*tau);

	for (int i=0;i<k;i++)
  	{	
		u(i)=quantile(i)*uv(ceil((n-1)*tau(i)))+(1-quantile(i))*uv(floor((n-1)*tau(i))); 
	}
  yh=x*beta;

	for (int i=0;i<k;i++)
  	{       r.col(i)=y-u(i)-yh;    
  		signw.col(i)=(1-arma::sign(r.col(i)))/2*(1-tau(i))+(arma::sign(r.col(i))+1)*tau(i)/2;	 
	}

	for (int j=0;j<p;j++)
	{       
               arma::vec xbeta=beta(j)*x.col(j);
               for(int i=0;i<k;i++)
		{
			z.col(i)=(r.col(i)+xbeta)/x.col(j);	
			newX.col(i)=x.col(j)%signw.col(i);	
		}
		vz=arma::vectorise(z);	
		order=arma::sort_index(vz);
		sortz=vz(order);
                vnewX=arma::vectorise(newX)	;			
		w=arma::abs(vnewX(order));		
                index=arma::find(cumsum(w)>(sum(w)/2),1,"first");
                place=int(index(0));
		beta(j)=sortz(place);
		
	}
error=sum(abs(beta-betaold));
iteration++;
}

return beta;

}
