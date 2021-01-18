#include <Rcpp.h>
#include <cmath>  
using namespace Rcpp;

double psiC(double x,int id,double par){
  double h;
  h = std::abs(x);
  if (id==1){return (x*x); };   /* for Spearman-coeff.*/
  if (id==2){return (h); };  /* for Spearman's footrule */
  if (id==3){return std::pow(h,par); };  /* Potenzkoeff.*/
  if (id==4){if (h<par) {return 0.5*h*h; } else {return par*(h-0.5*par); }}  /* Huberfunktion*/
  return 0.0;
}


//' @title coeffpml
//' @description Computing whatever
//' @details Some details
//' @author Eckhard Liebscher
//' @aliases coeffpml
//' @param u1 Ranks of x-values in subregion
//' @param v1 Ranks of x-values in subregion
//' @param u2 ranks of y-values in subregion
//' @param v2 ranks of y-values in subregion
//' @param amin minimum fraction of sample items in a subregion
//' @param n total number of sample items
//' @param na number of data in subregion
//' @param mf int
//' @return NumericVector
//' @useDynLib depcoeff
//' @importFrom Rcpp evalCpp
//' @keywords internal
// [[Rcpp::export]]
NumericVector coeffpml(NumericVector u1, NumericVector v1, NumericVector u2, NumericVector v2, double amin, 
double parp, double parh, long int n, long int na, int mf)
{ long int i,j,jj,k,l,im,imm,id;
  double h,hh,f0,u1mr,u2mr;
  NumericVector m(8),psiqf(4),par(4);
  par[0] = 0.0;
  par[1] = 0.0;
  par[2] = parp;
  par[3] = parh;  
  psiqf[0] = 0.166666666667;
  psiqf[1] = 0.333333333333;
  psiqf[2] = 2.0/(parp*(parp+3.0)+2.0);
  psiqf[3] = parh*(2.0-parh)*(parh*(parh-2.0)+2.0)/12.0;

  k = u1.size();
  l = u2.size();
   
  NumericMatrix s(n-na+1,8);
  /* main loop */
  for (j=na;j<=n-na;++j){
    jj=n-j;
    if (mf==1) { f0=jj;} else {
		if (mf==2) {f0=jj+1;} else {f0=sqrt(jj*jj-1.0);}}
	for (id=0;id<4;++id) { s(jj,id)= 0.0; s(j,id+4)=0.0;}  /*initialization of s*/
    /* left-side part: for data with rank<=n-j*/ 
	u1mr = 0;   /* maximum rank of u1 */
    for (i=0;i<k;++i) {    
	  h = (u1[i]-v1[i])/f0;        /*normed difference of ranks of x and y */
	  for (id=0;id<4;++id) {
        s(jj,id) += psiC(h,id+1,par[id]); /* psi function */
      }
      if (u1[i] > u1mr) { im = i; u1mr = u1[i];}   /* highest rank of u1 at im */ 
    }
     
    /* right-side part: for data with rank >j*/
	u2mr = n;   /* minimum rank of u2 */
    for (i=0;i<l;++i) {
	  h = (u2[i]+v2[i]-1.0)/f0;
      for (id=0;id<4;++id) { 
        s(j,id+4) += psiC(-1.0+h,id+1,par[id]);  /* psi function  */
      }
      if (u2[i] < u2mr) { imm=i; u2mr = u2[i];}  /* smallest rank of u2 at imm */
    }
     
    if (j<n-na) {  
      k= k-1; hh= v1[im];  /*hh= value v1 at maximum u1-rank*/
      if (im<k){
      for (i=im;i<k;++i){  /*delete the pair of data with highest u1-rank*/
        u1[i]= u1[i+1];
        v1[i]= v1[i+1];
      }}
	 
      if (u1mr<k+1) {   /* handling of ties */
		for (i=0;i<k;++i){ if (u1[i]==u1mr) {u1[i] -= 0.5;} }
      }
      for (i=0;i<k;++i){   /* update the ranking for Y */
        if (v1[i]>hh) {v1[i] -= 1;} else {
			if (v1[i]==hh)  {v1[i] -= 0.5;}}   /* ties */
      }
	  
      l= l-1; hh= v2[imm];   /*hh= value v2 at minimum u2-rank*/
      for (i=imm;i<l;++i){  /*delete the pair of data with smallest u2-rank */
        u2[i]= u2[i+1];
        v2[i]= v2[i+1];
      }
      
      for (i=0;i<l;++i){  /* update the ranking for X, Y*/
        if ((u2mr>1.0)&&(u2[i]==u2mr)) {u2[i] -= 0.5;} else {u2[i] -=1;} /* ties in X */
        if (v2[i]>hh) {v2[i] -=1;} else {
			if (v2[i]==hh)  {v2[i] -= 0.5;}}  /* ties in Y*/
      }
    }
  }
  
  for (id=0;id<4;++id) {
    m[id]= -10000.0;
    for (j=na;j<=n-na;++j){
      s(j,id) = 1.0 - (s(j,id)+s(j,id+4))/(n*psiqf[id]); 
      if (s(j,id)>m[id]) {
        m[id]= s(j,id);
        m[id+4]= j;
      }
    }
  }
 
  return(m);
}  
