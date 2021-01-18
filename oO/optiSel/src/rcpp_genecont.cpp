
#include <RcppArmadillo.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]]

Rcpp::NumericMatrix rcpp_genecont(const arma::ivec& numSire, const arma::ivec& numDam, const arma::ivec&  numAnc, const arma::ivec& numKeep, const arma::ivec& ainKeep, const Rcpp::CharacterVector rNames, const Rcpp::CharacterVector cNames, const arma::ivec& anOff){
  int i, j, nSire, nDam, nP, ni, nj;
  double pCont;
  int N     = numSire.n_elem;
  int NAnc  = numAnc.n_elem;
  int NKeep = numKeep.n_elem;  

  
  double** GCont  = (double**)calloc(N, sizeof(double*));  
  int*  isAlloc   = (int*)calloc(N, sizeof(int));  
  int*  nOff      = (int*)calloc(N, sizeof(int));  
  int*  inKeep    = (int*)calloc(N, sizeof(int));  
  
  for(j=0; j<NAnc;j++){
    nj = numAnc.at(j);
    GCont[nj]   = (double*) calloc(NAnc, sizeof(double));
    isAlloc[nj] = 1;
  }

  
  for(i=0; i<N;i++){
    if(!isAlloc[i]){GCont[i] = (double*) calloc(NAnc, sizeof(double));}
    isAlloc[i] = 1;
    nOff[i]    = anOff.at(i);
    inKeep[i]  = ainKeep.at(i);
    if(i<NAnc){GCont[numAnc.at(i)][i] = 1.0;}
    nSire = numSire.at(i);
    nDam  = numDam.at(i);
    
    if(nSire>0){nOff[nSire-1] = nOff[nSire-1]-1;}
    if(nDam>0){ nOff[nDam-1]  = nOff[nDam-1]-1;}
    
    if((nSire>0) | (nDam>0)){
      nP = ((nSire<nDam)?nDam:nSire);
      for(j=0; j<NAnc; j++){
        if(numAnc.at(j)<=nP){
          pCont = ((nSire>0)?(GCont[nSire-1][j]):0.0) + ((nDam>0)?(GCont[nDam-1][j]):0.0);
          if(pCont>0){GCont[i][j] = 0.5*pCont;}
        }
      }
    }
    if(nSire>0){
      if((!nOff[nSire-1]) & (!inKeep[nSire-1])){ 
        free(GCont[nSire-1]);
        isAlloc[nSire-1] = 0;
      }
    }
    if(nDam>0){
      if((!nOff[nDam-1]) & (!inKeep[nDam-1])){ 
        free(GCont[nDam-1]);
        isAlloc[nDam-1] = 0;
      }
    }
  }

  Rcpp::NumericMatrix rGeneCont(NKeep, NAnc);
  arma::mat GeneCont(rGeneCont.begin(), rGeneCont.nrow(), rGeneCont.ncol(), false);
  GeneCont.zeros();
  
  for(i=0;i<NKeep;i++){
    ni = numKeep.at(i);
    for(j=0;j<NAnc;j++){
      GeneCont.at(i, j) = GCont[ni][j];
    }
  }
  
  
  for(i=0;i<N;i++){
    if(isAlloc[i]){
      free(GCont[i]);
    }
  }
  free(GCont);
  free(isAlloc);
  free(inKeep);
  
  rGeneCont.attr("dimnames") = Rcpp::List::create(rNames, cNames);
  
  return rGeneCont;
}
