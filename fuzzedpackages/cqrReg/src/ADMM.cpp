#include "RcppArmadillo.h"
using namespace std;

arma::vec shrinkcpp(arma::vec u, arma::vec v)
{	
	arma::vec w=(1+sign(u-v))/2%(u-v)-(1+sign(-u-v))/2%(-u-v);
	return w;
}

// [[Rcpp::export]]
arma::vec QRADMMCPP(arma::mat xr,arma:: vec yr,arma::vec betar,double to,int mr,double ta,double rhor) 
{
double toler = (to);
int maxit=(mr);
double tau=(ta),rho=(rhor);
arma:: mat x=(xr);
arma:: vec signw,y=(yr),comparev=arma::zeros(3);
arma:: vec betaold,beta=(betar);

int n=x.n_rows;
int p=x.n_cols;
x.insert_cols( 0, arma::ones(n) );
double ABSTOL = 1e-4,RELTOL = 1e-2,alpha=1.4;
double rnorm,epspri,snorm,epsdual;
int iteration=1;
//toler=0.01;


//im=inverse matrix (x^T x)
arma::mat m= (x.t()*x);
arma::mat im=m.i()*x.t();
arma::vec betah=beta,gammah,uh,betai=betah,gammaold;
arma::vec u=arma::zeros(n),z=arma::zeros(n),zold,xbetah;
arma::vec yh=x*betah,uold,newr;
z=y-yh;


while (iteration<=maxit)
{
	betaold=betai;

//update beta
	betai = betah + im*(u/rho-z);
	beta = betai.subvec(1,p);

// updata z
	zold = z;
	xbetah=alpha*x*betai+(1-alpha)*(y-z);
	z = shrinkcpp(y-xbetah+u/rho-(2*tau-1)/rho,arma::ones<arma::vec>(n)/rho);

//updata u
	uold=u;
        newr=y-x*betai-z;
	u = u +rho*(newr);

//termination check
	rnorm = sqrt(accu(square(newr)));
	snorm = sqrt(accu(square(rho*x.t()*(z-zold))));

	comparev(0)=sqrt(accu(square(x*betai)));
	comparev(1)=sqrt(accu(square(z)));
	comparev(2)=sqrt(accu(square(y)));
		
	epspri = sqrt(n*1.0)*ABSTOL + RELTOL*arma::max(comparev);	
	epsdual = sqrt((p+1)*1.0)*ABSTOL + RELTOL*sqrt(accu(square(x.t()*u)));


	if (rnorm < epspri && snorm < epsdual) 
		{
			 iteration = maxit+1;
		} else {iteration = iteration + 1;}

	if (iteration>2&& arma::accu(arma::abs(betaold-betai))<toler&&rnorm < epspri)
	{ 

		iteration=maxit+1;
	}
}
arma::vec betanew=betai;
return betanew;
}

//CQRADMM
// [[Rcpp::export]]
arma::vec CQRADMMCPP (arma::mat xr,arma::vec yr,arma::vec betar,double to,int mr,arma::vec ta,double rhor,double pr){

double toler = (to);
int maxit=(mr);
arma::vec tau=(ta);
double rho=(rhor);
arma:: mat x=(xr);
arma:: vec signw,y=(yr),comparev=arma::zeros(3);
arma:: vec betaold,beta=(betar);

int n=x.n_rows;
int p=(pr);
double ABSTOL = 1e-4,RELTOL = 1e-2,alpha=1.4;
double rnorm,epspri,snorm,epsdual;
int iteration=1;
//toler=0.001;
int k=x.n_cols-p;

//im=inverse matrix (x^T x)
arma::mat m= (x.t()*x),xt;
arma::mat im=m.i()*x.t();
arma::vec betah=beta,gammah,uh,betai=betah,gammaold;
arma::vec u=arma::zeros(n),z=arma::zeros(n),zold,xbetah;
arma::vec yh=x*betah,uold,newr;
z=y-yh;
xt=x.t();

while (iteration<=maxit)
{
	betaold=betai;

//update beta
	betai = betah + im*(u/rho-z);
	//beta = betai.subvec(1,p);

// updata z
	zold = z;
	xbetah=alpha*x*betai+(1-alpha)*(y-z);
	z = shrinkcpp(y-xbetah+u/rho-(2*tau-1)/rho,arma::ones<arma::vec>(n)/rho);

//updata u
	uold=u;
        newr=y-x*betai-z;
	u = u +rho*(newr);

//termination check
	rnorm = sqrt(accu(square(newr)));
	snorm = sqrt(accu(square(rho*xt*(z-zold))));

	comparev(0)=sqrt(accu(square(x*betai)));
	comparev(1)=sqrt(accu(square(z)));
	comparev(2)=sqrt(accu(square(y)));
		
	epspri = sqrt(n*1.0)*ABSTOL + RELTOL*arma::max(comparev);	
	epsdual = sqrt((p+k)*1.0)*ABSTOL + RELTOL*sqrt(accu(square(xt*u)));


	if (rnorm < epspri && snorm < epsdual) 
		{
			 iteration = maxit+1;
		} else {iteration = iteration + 1;}

	if (iteration>2&& accu(arma::abs(betaold-betai))<toler)
	{ 
		iteration=maxit+1;
	}
}
arma::vec betanew=betai;
return betanew;
}


//QRPADMM depends on lselasso
arma::vec lselassocpp(arma::mat xr,arma::vec yr,int mr,double rhor,double lambdar)
{

int maxit=(mr);
double rho=(rhor),lambda=(lambdar);
arma:: mat x=(xr);
arma:: vec z,signw,y=(yr);
arma:: vec betaold;

int n=x.n_rows;
int p=x.n_cols;
x.insert_cols( 0, arma::ones(n) );
double ABSTOL = 1e-6,RELTOL = 1e-3,alpha=1.2;
double rnorm,epspri,snorm,epsdual;
int iteration=1;
arma::mat rhom=rho*arma::eye(p+1,p+1);
rhom(0,0)=0;

//im=inverse matrix (x^T x)
arma::mat m= (x.t()*x+rhom);
arma::mat im=m.i();
arma::vec betah=im*x.t()*y,gammah,uh,betai,gammaold;
arma::vec u=arma::zeros(p),gamma=arma::zeros(p),beta;
beta=betah;

while (iteration<=maxit)
{
	

//gammah=c(0,gamma), uh=c(0,u)
	gammah=gamma;
	gammah.insert_rows(0,1);
	uh=u;
	uh.insert_rows(0,1);

//update beta
	betaold=beta;
	betai = betah + im*(rho*gammah-uh);
	beta = betai.subvec(1,p);

// updata gamma
	gammaold = gamma;
        arma::vec hbeta=alpha*beta+(1-alpha)*gamma;
	gamma = shrinkcpp(hbeta+u/rho,lambda*arma::ones<arma::vec>(p)/rho);

//updata u
	u = u +rho*(beta-gamma);

//termination check
	rnorm = sqrt(accu(square(beta-gamma)));
	snorm = sqrt(accu(square(-rho*(gamma-gammaold))));
		
	epspri = sqrt(p*1.0)*ABSTOL + RELTOL*std::max(sqrt(accu(square(beta))),sqrt(accu(square(-gamma))));	
	epsdual = sqrt(p*1.0)*ABSTOL + RELTOL*sqrt(accu(square(u)));


	if (rnorm < epspri && snorm < epsdual) 
		{
			 iteration = maxit+1;
		} else {iteration = iteration + 1;}
}

arma::vec betanew=betai;
return (betanew);
}

// [[Rcpp::export]]
arma::vec QRPADMMCPP(arma::mat xr,arma::vec yr,arma::vec betar,
int m,double ta,double rhor,double lambdar){

int maxit=(m);
double tau=(ta),rho=(rhor),lambda=(lambdar);
arma:: mat x=(xr);
arma:: vec z,signw,y=(yr),r;
arma:: vec betaold,beta=(betar);

int n=x.n_rows;
int p=x.n_cols;
x.insert_cols( 0, arma::ones(n) );
double ABSTOL = 1e-6,RELTOL = 1e-3;
double rnorm,epspri,snorm,epsdual,alpha=1.4;
int iteration=1;
arma::mat rhom=rho*arma::eye(p+1,p+1);
rhom(0,0)=0;

arma::vec betah=beta,gammah,uh,gammaold;
arma::vec u=arma::zeros(n),gamma=arma::zeros(p),xbetah;
arma::vec comparev=arma::zeros(3);
r=y-x*betah;
beta=beta.subvec(1,p);


while (iteration<=maxit)
{
//updata r- residual
	xbetah=alpha*x*betah+(1-alpha)*(y-r);
	r = shrinkcpp(u/rho+y-xbetah-(2*tau-1)/rho,arma::ones<arma::vec>(n)/rho);

//update beta with intercept
	betaold=beta;
	betah = lselassocpp(x.cols(1,p),u/rho+y-r,maxit,rho,lambda/rho);
	beta=betah.subvec(1,p);

//updata u
	u = u +rho*(y-x*betah-r);

//termination check
	rnorm = sqrt(accu(square(y-x*betah-r)));
	snorm = sqrt(accu(square(rho*x.cols(1,p)*(beta-betaold) )));
	
        comparev(0)=sqrt(accu(square(x.cols(1,p)*beta)));
	comparev(1)=sqrt(accu(square(-r)));
	comparev(2)=sqrt(accu(square(y-betah(0))));	
	epspri = sqrt(n*1.0)*ABSTOL + RELTOL*arma::max(comparev);	
	epsdual = sqrt(n*1.0)*ABSTOL + RELTOL*sqrt(accu(square(u)));


	if (rnorm < epspri && snorm < epsdual) 
		{
			 iteration = maxit+1;
		} else {iteration = iteration + 1;}
}
arma::vec betanew=betah;
return betanew;
}

//CQRPADMM depends on cselasso
arma::vec cselassocpp(arma::mat xr,arma::vec yr,int mr,double rhor,double lambdar,int p)
{

int maxit=(mr);
double rho=(rhor),lambda=(lambdar);
arma:: mat x=(xr);
arma:: vec z,signw,y=(yr);
arma:: vec betaold,beta;

//int n=x.n_rows;
int k=x.n_cols-p;
double ABSTOL = 1e-6,RELTOL = 1e-3;
double rnorm,epspri,snorm,epsdual;
int iteration=1;
arma::mat rhom=rho*arma::eye(p+k,p+k);
for (int i=0;i<k;i++)
	{
		rhom(i,i)=0;
	}
//im=inverse matrix (x^T x)
arma::mat m= (x.t()*x+rhom);
arma::mat im=m.i();
arma::vec betah=im*x.t()*y,gammah,uh,betai,gammaold;
arma::vec u=arma::zeros(p),gamma=arma::zeros(p);


while (iteration<=maxit)
{

//gammah=c(rep(0,k),gamma), uh=c(rep(0,k),u)
	gammah=gamma;
	gammah.insert_rows(0,k);
	uh=u;
	uh.insert_rows(0,k);

//update beta
	betai = betah + im*(rho*gammah-uh);
	beta = betai.subvec(k,p+k-1);

// updata gamma
	gammaold = gamma;
	gamma = shrinkcpp(beta+u/rho,lambda*arma::ones<arma::vec>(p)/rho);

//updata u
	u = u +rho*(beta-gamma);

//termination check
	rnorm = sqrt(accu(square(beta-gamma)));
	snorm = sqrt(accu(square(-rho*(gamma-gammaold))));
		
	epspri = sqrt(p*1.0)*ABSTOL + RELTOL*std::max(sqrt(accu(square(beta))),sqrt(accu(square(-gamma))));	
	epsdual = sqrt(p*1.0)*ABSTOL + RELTOL*sqrt(accu(square(u)));


	if (rnorm < epspri && snorm < epsdual) 
		{
			 iteration = maxit+1;
		} else {iteration = iteration + 1;}
}

arma::vec betanew=betai;
return (betanew);
}

// [[Rcpp::export]]
arma::vec CQRPADMMCPP(arma::mat xr,arma::vec yr,arma::vec betar,
int mr,arma::vec ta,double rhor,double lambdar,int pr,int kr) {

int maxit=(mr),p=(pr);
double rho=(rhor),lambda=(lambdar);
arma:: mat x=(xr);
arma:: vec z,signw,y=(yr),tau=(ta),r;
arma:: vec betaold,beta=(betar);

int n=x.n_rows;
double ABSTOL = 1e-6,RELTOL = 1e-3;
double rnorm,epspri,snorm,epsdual,alpha=1.4;
int iteration=1;

int k=(kr);
//im=inverse matrix (x^T x)
arma::mat m= (x.t()*x);
arma::mat im=m.i(),b;
arma::vec betah=beta,gammah,uh,betai,gammaold;
arma::vec u=arma::zeros(n),gamma=arma::zeros(p),xbetah;
arma::vec comparev=arma::zeros(3);
r=y-x*betah;
beta=betah.subvec(k,p+k-1);

while (iteration<=maxit)
{
//updata r- residual

	xbetah=alpha*x*betah+(1-alpha)*(y-r);
	r = shrinkcpp(u/rho+y-xbetah-(2*tau-1)/rho,arma::ones<arma::vec>(n)/rho);


	betaold=beta;
	betah = cselassocpp(x,u/rho+y-r,maxit,rho,lambda/rho,p);
	beta=betah.subvec(k,p+k-1);

//updata u
	u = u +rho*(y-x*betah-r);

//termination check
	rnorm = sqrt(accu(square(y-x*betah-r)));
	snorm = sqrt(accu(square(rho*x.cols(k,p+k-1)*(beta-betaold) )));

        comparev(0)=sqrt(accu(square(x.cols(k,p+k-1)*beta)));
	comparev(1)=sqrt(accu(square(-r)));

	arma::vec bv=betah.subvec(0,k-1);
	arma::mat bm=arma::mat(bv);
	arma::mat bs=trans(arma::reshape(bm,k,n/k));

	b=reshape(bs,n,1);

	comparev(2)=sqrt(accu(square(y-b)));
	
	epspri = sqrt(n*1.0)*ABSTOL + RELTOL*arma::max(comparev);	
	epsdual = sqrt(n*1.0)*ABSTOL + RELTOL*sqrt(accu(square(u)));


	if (rnorm < epspri && snorm < epsdual) 
		{
			 iteration = maxit+1;
		} else {iteration = iteration + 1;}
}
arma::vec betanew=betah;
return betanew;
}




