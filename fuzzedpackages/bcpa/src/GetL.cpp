#include <Rcpp.h>
#include <math.h>
#include <Rmath.h> // use for normal distribution

using namespace Rcpp;

RNGScope scope;

NumericVector SubSet(NumericVector A, int start, int end)
{
	NumericVector B(A.begin() + start, A.begin() + end);
	return B;
}

// [[Rcpp::export]]
double GetL(NumericVector x, NumericVector t, double rho, bool tau = false)
{
	NumericVector Time(t);
	NumericVector X(x);
	double Rho = rho;
  bool Tau = tau;
	
	int N = Time.size()-1;
	
	double mux = mean(X); 
	double sigmax = sd(X);
	
	NumericVector dT(N);
	NumericVector RhoToTheDt(N); 

	NumericVector Xstart = SubSet(X, 0, N);
	NumericVector Xend = SubSet(X, 1, N+1);
	dT = diff(Time); 

  if(Tau == false)
	for(int i = 0; i < N; i++)
    RhoToTheDt[i] = pow(Rho, dT[i]);
    
  else
  for(int i = 0; i < N; i++)
    RhoToTheDt[i] = exp(-dT[i]/Rho);

	NumericVector Mu = RhoToTheDt*(Xstart-mux) + mux;
	NumericVector Sigma = sigmax*sqrt(1-RhoToTheDt*RhoToTheDt);
  	
	double LL = 0;
	
	for(int i = 0; i < N; i++)
   LL = LL + R::dnorm(Xend[i], Mu[i], Sigma[i], 1);

	return LL;
}
