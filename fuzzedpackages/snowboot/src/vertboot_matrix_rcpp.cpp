#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;
//[[Rcpp::export]]
IntegerMatrix vertboot_matrix_rcpp(IntegerMatrix m1,IntegerVector blist){
  int a;
  int b;
  int num=m1.nrow();
  IntegerMatrix x1(num,num);
  for(int k=0;k<num;k++){
    for(int q=0;q<num;q++){
      x1(k,q)=0;
    }
  }
  for(int i=0;i<num;i++){
    a=blist[i];
    for(int j=0;j<num;j++){
      b=blist[j];
      if(a!=b){
        x1(i,j)=m1(a,b);
      }
      else{
        a=round(R::runif(-0.49,num-0.5));
        b=round(R::runif(-0.49,num-0.5));
        while(a==b){
          b=round(R::runif(-0.49,num-0.5));
        }
        x1(i,j)=m1(a,b);
        x1(j,i)=m1(a,b);
      }
    }
  }
  return (x1);
}
