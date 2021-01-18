
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
SEXP file1292576f( SEXP A, SEXP B, SEXP C, SEXP D, SEXP E) ;
}

// definition

SEXP file1292576f( SEXP A, SEXP B, SEXP C, SEXP D, SEXP E ){
BEGIN_RCPP


	double currentband = as<double>(A);
	double kerneltog = as<double>(C);
	int leftout = as<int>(D);

	Rcpp::NumericVector prank(B);
	Rcpp::NumericMatrix itemanswers(E);

	int noptions = itemanswers.nrow();
	int nex = prank.size();
	int j,k,l,m;



	double arg, denom;


	NumericVector num(nex);
	NumericVector smoother(nex);
	NumericVector smoothed(noptions);


	double probleftout = prank[leftout-1];

	denom=0.0;

	for(j=0; j < nex; j++){
		
			arg = (probleftout - prank[j])/currentband;
 
			if(j==leftout-1){num[j]=0;}

			else if(kerneltog==1){
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


	for(k=0; k < nex; k++){
		smoother[k] = num[k] / denom;

	}


	
	for(l=0; l<noptions;l++){
		smoothed[l]=0.0;
		for(m=0; m<nex;m++){
			smoothed[l]=smoothed[l]+smoother[m]*itemanswers(l,m);
		}
		


	}


	return(smoothed);

	
END_RCPP
}



