
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP file4ef03015( SEXP A, SEXP B, SEXP C, SEXP D, SEXP E) ;
}

// definition

SEXP file4ef03015( SEXP A, SEXP B, SEXP C, SEXP D, SEXP E ){
BEGIN_RCPP




double currenth = as<double>(A);

double kerneltog = as<double>(E);



Rcpp::NumericVector probrank(B);

Rcpp::NumericVector theta(C);

Rcpp::NumericVector optresp(D);



int nex = probrank.size();

int ntheta = theta.size();



int i,j;



double arg, denom;



Rcpp::NumericVector num(nex);

Rcpp::NumericVector smoothed(ntheta);

Rcpp::NumericVector smoother(ntheta);

Rcpp::NumericVector sqrd(nex);

Rcpp::NumericVector str(ntheta);





for (i=0; i  < ntheta ; i++){



denom = 0.0;

smoothed[i] = 0.0;

smoother[i] = 0.0;

str[i] =0.0;



for(j=0; j < nex; j++){



arg=(theta[i] - probrank[j])/currenth;



if(kerneltog==1){

num[j]=exp(-1*(arg*arg)/2);

}



else if(kerneltog==2){

if(fabs(arg) <= 1){ 

num[j]=1-(arg*arg);

}

else{ num[j]= 0;}

}

else if(kerneltog==3){

if(fabs(arg) <= 1){ num[j]=1;}

else{num[j]=0;}



}



denom = denom + num[j];



}



for(j=0; j < nex; j++){

if(denom==0){
	num[j]=0;
	denom=1;
}


smoother[i] = smoother[i] + (num[j] / denom);

smoothed[i] = smoothed[i] + optresp[j] * (num[j] / denom);

str[i] = str[i] + (num[j]/denom)*(num[j]/denom)*smoothed[i]*(1-smoothed[i]);



}



str[i]=sqrt(str[i]);



 

}



return List::create(_["ICC"] = smoothed, _["stderr"] = str, _["weights"] = smoothed );








END_RCPP
}



