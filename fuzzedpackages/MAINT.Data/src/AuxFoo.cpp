#include <RcppArmadillo.h>
#include "AdMatAlgFoo.h"
#include "AuxFoo.h"

void outerprod(const int p,const vec& v1,const vec& v2,mat& res)
{
	for (int r=0;r<p;r++) for(int c=0;c<p;c++) res(r,c) = v1(r)*v2(c);
	return;
}

void outerprod(const int p,const vec& v,mat& res)
{
	for (int r=0;r<p;r++) for(int c=0;c<=r;c++) {
		res(r,c) = v(r)*v(c);
		if (c<r) res(c,r) = res(r,c);
	}
	return;
}

//mat RestCov(const int q, NumericVector::iterator xpos, const int Config, const bool FixedArrays)
mat RestCov(const int q, NumericVector::iterator xpos, const int Config, const bool FixedArrays, const bool Srpar)
{
	int p(2*q);
	/* static */	mat A,Sigma;
/*	int NullMatdim(0);
	mat A,Sigma,NullppMat;
	// if (NullMatdim.size()!=p) NullppMat.resize(p,p);
	if (NullMatdim!=p) { 
		SetZero(NullppMat,p,p,!FixedArrays);
 		NullMatdim = p;
	}
*/
 	switch (Config)  {
		case 1: 
//			if (A.size()!=p) A.resize(p,p);
//			A = NullppMat;  			
			SetZero(A,p,p,!FixedArrays);
			for (int i=0;i<p;i++) A(i,i) = *xpos++;
			for (int c=1;c<p;c++) for (int r=0;r<c;r++) A(r,c) = *xpos++;
//			return A.t() * A;
			if (Srpar) return A.t() * A;
			else {
			  for (int c=0;c<p-1;c++) for (int r=c+1;r<p;r++) A(r,c) = A(c,r); 
                          return A;
                        }
		case 2: 
//			if (A.size()!=p) A.resize(p,p);
//			A = NullppMat;  			
			SetZero(A,p,p,!FixedArrays);
 			for (int i=0;i<p;i++) A(i,i) = *xpos++;
  			for (int c=1;c<p;c++) for (int r=0;r<c;r++)  {
    		          if ( (r<q && c<q)  || (r>=q && c>=q) || c==r+q )  A(r,c) = *xpos++;
    			  else if (r>0)  {
			    double dbltmp = 0.;
			    for (int i=0;i<r;i++) dbltmp -= A(i,r)*A(i,c);   
      	                    A(r,c) = dbltmp/A(r,r);
			  }
  			}  
			return A.t() * A;
		case 3: 
			int qplsi;
			double A11,A12,A22;
//			if (Sigma.size()!=p) Sigma.resize(p,p);
//			Sigma = NullppMat;
			SetZero(Sigma,p,p,!FixedArrays);
			for (int i=0;i<q;i++)  {
			  A11 = *(xpos+i);
			  A12 = *(xpos+p+i);
			  A22 = *(xpos+(qplsi=q+i));
			  if (Srpar) {
			    Sigma(i,i) = A11*A11;
			    Sigma(i,qplsi) = Sigma(qplsi,i) = A11*A12;  
			    Sigma(qplsi,qplsi) = A12*A12 + A22*A22;
                          } else {
			    Sigma(i,i) = A11;
			    Sigma(i,qplsi) = Sigma(qplsi,i) = A12;  
			    Sigma(qplsi,qplsi) = A22;
                          }
			}
			return Sigma;  
		case 4: 
			if (A.size()!=q) {
			  A.set_size(q,q);
  			  for (int c=0;c<q-1;c++) for (int r=c+1;r<q;r++)  A(r,c) = 0.;
			}
//			if (Sigma.size()!=p) Sigma.resize(p,p);
//			Sigma = NullppMat;
			SetZero(Sigma,p,p,!FixedArrays);
 			for (int i=0;i<q;i++) A(i,i) = *(xpos+i);
  			for (int c=1,xind=p;c<q;c++) for (int r=0;r<c;r++,xind++)  A(r,c) = *(xpos+xind);
//			Sigma.block(0,0,q,q) = A.t() * A; 
			if (Srpar) Sigma.submat(0,0,q-1,q-1) = A.t() * A;
                        else Sigma.submat(0,0,q-1,q-1) = A;
 			for (int i=0;i<q;i++) A(i,i) = *(xpos+q+i);
  			for (int c=1,xind=p+q*(q-1)/2;c<q;c++) for (int r=0;r<c;r++,xind++)  A(r,c) = *(xpos+xind);
//			Sigma.block(q,q,q,q) = A.t() * A;
			if (Srpar) Sigma.submat(q,q,p-1,p-1) = A.t() * A;
                        else {
                          Sigma.submat(q,q,p-1,p-1) = A;
                          for (int c=0;c<p-1;c++) for (int r=c+1;r<p;r++) Sigma(r,c) = Sigma(c,r); 
                        } 
			return Sigma;
		case 5: 
			double x;
//			if (Sigma.size()!=p) Sigma.resize(p,p);
//			Sigma = NullppMat;
			SetZero(Sigma,p,p,!FixedArrays);
			for (int i=0;i<p;i++) {
                          x = *xpos++;
//                        Sigma(i,i) = x*x;
                          if (Srpar) Sigma(i,i) = x*x;
                          else Sigma(i,i) = x;
			}  
			return Sigma;
	}
        return Sigma;
}

