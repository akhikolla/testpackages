 
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

Rcpp::DataFrame rcpp_completeness(Rcpp::StringVector Indiv, const arma::ivec& ArmanumSire, const arma::ivec& ArmanumDam, int maxd){
  int n, i, j, g, sDepth, dDepth;
  double sCont, dCont;
  int N = ArmanumSire.n_elem;
  int*     pedDepth = (int*)calloc(N,sizeof(int));
  double** pedCompl = (double**)calloc(N, sizeof(double*));  
  int*     numSire  = (int*)calloc(N, sizeof(int)); 
  int*     numDam   = (int*)calloc(N, sizeof(int)); 
  if(pedDepth == NULL){error_return("Memory allocation failed.");};	  
  if(pedCompl == NULL){error_return("Memory allocation failed.");};  
  if(numSire  == NULL){error_return("Memory allocation failed.");};
  if(numDam   == NULL){error_return("Memory allocation failed.");};
  for(i=0; i<N;i++){
    numSire[i] = ArmanumSire.at(i);
    numDam[i]  = ArmanumDam.at(i);
  }
  
  n = 0;
  for(i=0; i<N;i++){
    sDepth = ((numSire[i]>0)?(pedDepth[numSire[i]-1]):(0));
    dDepth = ((numDam[i]>0)?(pedDepth[numDam[i]-1]):(0));
    pedDepth[i] = 1 + ((sDepth>dDepth)?(sDepth):(dDepth));
    if(pedDepth[i]>(maxd+1)){pedDepth[i] = maxd+1;}
    n = n + pedDepth[i];
    pedCompl[i] = (double*)calloc(pedDepth[i],sizeof(double));
    if(pedCompl[i] == NULL){error_return("Memory allocation failed.");};
    pedCompl[i][0] = 1.0;
    for(g=1; g<pedDepth[i]; g++){
      sCont = ((sDepth>=g)?(pedCompl[numSire[i]-1][g-1]):(0.0));
      dCont = ((dDepth>=g)?(pedCompl[numDam[i]-1][ g-1]):(0.0));
      pedCompl[i][g] = 0.5*(sCont + dCont);
    }
  }
  
  Rcpp::NumericVector Completeness(n);
  Rcpp::IntegerVector Generation(n);
  Rcpp::StringVector  Individual(n);
  j=0;
  for(i=0;i<N;i++){
    for(g=0;g<pedDepth[i];g++){
      Completeness.at(j) = pedCompl[i][g];
      Generation.at(j)   = g;
      Individual.at(j)   = Indiv.at(i);
      j = j + 1;
    }
    free(pedCompl[i]);
  }
  free(pedDepth);
  free(pedCompl);
  free(numSire);
  free(numDam);
  
  return Rcpp::DataFrame::create(Named("Indiv")= Individual, Named("Generation") = Generation, Named("Completeness") = Completeness, Named("stringsAsFactors")=false);
}
