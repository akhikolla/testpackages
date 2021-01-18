#include "subMat.h"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP subMat(SEXP x, SEXP nx, SEXP ny, SEXP nz, SEXP nper){

  Rcpp::NumericVector xx(x);
  
  int NX = Rcpp::as<int>(nx);
  int NY = Rcpp::as<int>(ny);
  int NZ = Rcpp::as<int>(nz);
  int N = NX + NY + NZ;
  int NPER = Rcpp::as<int>(nper);
      
  // Create the matrices
  mat X(N, N);
  mat XX(N, N);
  mat temp(N, N);
  mat ST(N, N);
  mat EQ(N, N);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
	X(j,i) = xx(i);
	XX(i,j) = xx(i);
    }
  }
  temp = X-XX;
           
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
        if(temp(i,j)>0.0) { ST(i,j) = 1.0;} else {ST(i,j) = 0.0;};
	if(temp(i,j)==0.0){ EQ(i,j) = 1.0;} else {EQ(i,j) = 0.0;};
    }
  }
  int ties=0;
  if(accu(EQ)>N) ties = 1;
      
// Create here the permutation matrices
  mat SST(N, N);
  mat SSTTemp(N, N);
  mat SEQ(N, N);
  mat SEQTemp(N, N);
  vec index(N);
  vec result(NPER);
  
  for(int i=0;i<N;i++) index(i)=i;
 
  if(NPER>1){
    if(ties==1){
      for(int pRun=0;pRun<NPER;pRun++){

	index = shuffle(index);
  	for(int i=0;i<N;i++) SSTTemp.row(i) = ST.row(index(i));
	for(int i=0;i<N;i++) SST.col(i) = SSTTemp.col(index(i));
	for(int i=0;i<N;i++) SEQTemp.row(i) = EQ.row(index(i));
	for(int i=0;i<N;i++) SEQ.col(i) = SEQTemp.col(index(i));
  
	result(pRun) = accu(SST.submat(0,NX,NX-1,NX+NY-1)*SST.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1)) + 
			0.5 * accu(SST.submat(0,NX,NX-1,NX+NY-1)*SEQ.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1)) +
			0.5 * accu(SEQ.submat(0,NX,NX-1,NX+NY-1)*SST.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1)) +
			accu(SEQ.submat(0,NX,NX-1,NX+NY-1)*SEQ.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1))/6;
      }      
    } else {
      for(int pRun=0;pRun<NPER;pRun++){

	index = shuffle(index);
  	for(int i=0;i<N;i++) SSTTemp.row(i) = ST.row(index(i));
	for(int i=0;i<N;i++) SST.col(i) = SSTTemp.col(index(i));
	result(pRun) = accu(SST.submat(0,NX,NX-1,NX+NY-1)*SST.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1));
	
      }
    }
  } else {
    result(0) = accu(ST.submat(0,NX,NX-1,NX+NY-1)*ST.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1))+ 
		0.5 * accu(ST.submat(0,NX,NX-1,NX+NY-1)*EQ.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1)) +
		0.5 * accu(EQ.submat(0,NX,NX-1,NX+NY-1)*ST.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1)) +
		accu(EQ.submat(0,NX,NX-1,NX+NY-1)*EQ.submat(NX,NX+NY,NX+NY-1,NX+NY+NZ-1))/6;;
  }
   
   
    List res;
    res["result"] =  result/(NX*NY*NZ);
    return res;
}
