#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


#include "Common.h"
#include "Normal.h"
#include "Logistic.h"
#include "Cauchy.h"
#include "Kmt.h"

//[[Rcpp::export]]
List KmtMain(arma::vec X, int bModified, arma::mat NormalMat, arma::mat LogisMat, arma::mat ReMat, arma::mat CauchyMat, String strDist, int bGraph, int nNum){
  
  int n = X.n_elem;
  Kmt kmt(X, n, NormalMat, LogisMat, ReMat, CauchyMat, strDist);
  
  
  kmt.SetGiMat(); 
  //arma::mat GMat = kmt.GetGiMat(); 
  
  
  double xVal = 0;
  double FVal = 0;
  
  double xValP = 0;
  double FValP = 0;
  double xValM = 0;
  double FValM = 0;
  
  
  
  
  if(bModified==0){
    kmt.SetT2();
    kmt.FindOptimal();
    
    xVal = kmt.GetOptX();
    FVal = kmt.GetOptFVal();
    
    
  }else{
    kmt.Modified_SetT2();
    kmt.Modified_FindOptimal();
    
    xValP = kmt.GetOptXP();
    FValP = kmt.GetOptFValP();
    
    xValM = kmt.GetOptXM();
    FValM = kmt.GetOptFValM();
    
    
  }
  
  int nGraph = (n+1)*nNum;
  arma::vec GraphVec(nGraph);
  arma::vec LineVec(nGraph);
  
  LineVec.zeros();
  LineVec = GetLineVec(X, nNum);
  
  GraphVec.zeros();
  
  
  double xi=0; 
  double tmp = 0;
  
  if(bGraph==1){
    
    
    for(int i=1;i<= nGraph;i++){
      
      if(i<= (nGraph-nNum) ){
        xi = LineVec[i-1];
        tmp = kmt.ObjVal(xi);
        GraphVec[i-1] = tmp;
        
      }else{
        if(i == (nGraph-nNum+1)){
          xi = LineVec[i-1];
          tmp = kmt.ObjVal(xi);
          GraphVec[i-1] = tmp;
        }else{
          GraphVec[i-1] = tmp;
        }
        
        
        
      }
      
      
    }
    
  }
  
  
  List lst(8);
  lst[0] = xVal;
  lst[1] = FVal;
  lst[2] = LineVec;
  lst[3] = GraphVec;
  
  lst[4] = xValP;
  lst[5] = FValP;
  lst[6] = xValM;
  lst[7] = FValM;
  
  return lst;
}






