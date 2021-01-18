#ifndef _EigenRestCovGrad_cpp
#define _EigenRestCovGrad_cpp

#include "msnCP_Aux.h"
#include "AuxFoo.h"

int inline utdind(const int i) { return (i+1)*(i+2)/2 - 1; } 
int inline utind0(const int r,const int c) { return c*(c-1)/2 + r; } 
int inline utind1(const int r,const int c) { return c*(c+1)/2 + r; } 

template<class RCTMATTP>
void C1CPgrad(const int p,NumericVector::iterator xpos,RCTMATTP& Jacob)
{
  int lrind,lcind;

  for (int c=0,ind=0;c<p;c++,ind++)  {
    for (int r=0;r<c;r++,ind++)  for (int l=0;l<=r;l++)  {
      if (l<r) lrind = p+utind0(l,r); else lrind = r;
      lcind = p+utind0(l,c);
      Jacob(ind,lcind) = *(xpos+lrind);  
      Jacob(ind,lrind) = *(xpos+lcind);  
    }
    for (int l=0;l<=c-1;l++)  {
      lcind = p + utind0(l,c);
      Jacob(ind,lcind) = 2*(*(xpos+lcind));
    }  
    Jacob(ind,c) = 2*(*(xpos+c));
  }  
  return;  
}

template<class RCTMATTP,class SQMATTP,class EXTVCTTP>
void C2CPgrad(const int p,const int q,const int nvcovpar,
  NumericVector::iterator xpos,const bool FixedArrays,RCTMATTP& Jacob)
{
   static SQMATTP Sigma,SigmaSrU;
   static RCTMATTP JacobC1;
   static EXTVCTTP fullpar;
   static arma::ivec nullrows,validcols;
   static int nullrowsdim(0),validcolsdim(0);
   int cind,row,qtrdim(q*(q-1)/2),pplsqtrdim(p+qtrdim);
   int nnullrows(2*qtrdim),nvalidcols(nvcovpar-nnullrows);
   double Asum,Acc;

  if (nullrowsdim!=nnullrows && q>1)  {
    nullrows.set_size(nnullrows);
    for (int c=q,ind=0,ind1=q*(q+1)/2;c<p;c++) { 
      for (int r=0;r<q;r++) if (c!=q+r) nullrows(ind++) = ind1+r;
      ind1 += c+1;
    }
    nullrowsdim = nnullrows;
  }
  if (validcolsdim!=nvalidcols)  {
    validcols.set_size(nvalidcols);
    for (int ind1=0;ind1<pplsqtrdim;ind1++) validcols(ind1) = ind1; 
    for (int c=q,ind=pplsqtrdim,ind1=pplsqtrdim;c<p;c++) for (int r=0;r<c;r++,ind1++) 
    if (c==q+r || r>=q) validcols(ind++) = ind1;
    validcolsdim = nvalidcols;
  }
  if (!FixedArrays) {
    if (Sigma.n_rows!=p || Sigma.n_cols!=p) Sigma.set_size(p,p);
    if (SigmaSrU.n_rows!=p || SigmaSrU.n_cols!=p) SigmaSrU.set_size(p,p);
    if (fullpar.size()!=nvcovpar)  fullpar.set_size(nvcovpar);
  }
//  Sigma = RestCov(q,xpos,2,FixedArrays);
  Sigma = RestCov(q,xpos,2,FixedArrays,true);
  SigmaSrU = chol(Sigma);      
  for (int c=0,i=p;c<p;c++)  {		
    fullpar(c) = SigmaSrU(c,c);
    for (int r=0;r<c;r++,i++) fullpar(i) = SigmaSrU(r,c);
  }
  NumericVector fullparasNV(wrap(fullpar));
  SetZero(JacobC1,nvcovpar,nvcovpar,true); 
  C1CPgrad<RCTMATTP>(p,fullparasNV.begin(),JacobC1);
  if (q==1) return; 

  for (int r=q-2;r>=0;r--)  {
    for (int c=q;c<p;c++) for (int i=r+1;i<q;i++) if (c!=q+i)  {
      cind = p+utind0(r,c);
      for (int r1=0;r1<nvcovpar;r1++) 
        JacobC1(r1,cind) -=  *(xpos+p+utind0(r,i)) * JacobC1(r1,p+utind0(i,c))  / *(xpos+i); 
    }  
    for (int c=r+1;c<q;c++) for (int j=q;j<p;j++) if (j!=q+c)  {
      cind = p+utind0(r,c);
      for (int r1=0;r1<nvcovpar;r1++)
        JacobC1(r1,cind) -= *(xpos+p+utind0(r,j)) * JacobC1(r1,p+utind0(c,j))  / *(xpos+c); 
    }
  }  
  for (int c=q-1;c>0;c--) for (int j=q;j<p;j++) if (j!=q+c)  {
    Asum = 0.;
    for (int c0=0;c0<c;c0++)  Asum += *(xpos+p+utind0(c0,c)) * *(xpos+p+utind0(c0,j)) ; 
    for (int r1=0;r1<nvcovpar;r1++) {
      Acc = *(xpos+c);
      JacobC1(r1,c) += Asum*JacobC1(r1,p+utind0(c,j)) / (Acc*Acc);
    }
  } 
  for(int c=0;c<nvalidcols;c++) for(int r=0;r<nvalidcols;r++) {
    row = validcols(r);
    Jacob(row,c) = JacobC1(row,validcols(c));
  }
  for (int r=0;r<nnullrows;r++) for(int c=0;c<nvalidcols;c++) Jacob(nullrows(r),c) = 0.;

  return;
}

template<class RCTMATTP>
void C3CPgrad(const int p,const int q,NumericVector::iterator xpos,RCTMATTP& Jacob)
{
  int v12rind,v12cind;
  double xv,xv12cind;

  for (int v=0;v<q;v++) {
    Jacob(utdind(v),v) = 2*(xv=(*(xpos+v)));
    Jacob(utdind(q+v),q+v) = 2*(*(xpos+q+v));
    Jacob(v12rind=utind1(v,v+q),v) = xv12cind = *(xpos+(v12cind=p+v));
    Jacob(v12rind,v12cind) = xv;
    Jacob(utdind(q+v),v12cind) = 2*xv12cind;
  }
  return; 
}

template<class RCTMATTP,class SQMATTP,class EXTVCTTP>
void C4CPgrad(const int p,const int q,const int nvcovpar,
  NumericVector::iterator xpos,const bool FixedArrays,RCTMATTP& Jacob)
{
  static RCTMATTP JacobC1;
  static EXTVCTTP ppar;
  int npcovpar(q*(q+1)/2),qtrdim(q*(q-1)/2);

  if (!FixedArrays && ppar.size()!=npcovpar)  ppar.set_size(npcovpar);
  NumericVector pparasNV(wrap(ppar));
  SetZero(JacobC1,npcovpar,npcovpar,true); 

  for(int j=0;j<q;j++) pparasNV(j) = *(xpos+j);
  for(int j=q,j1=p;j<q+qtrdim;j++,j1++) pparasNV(j) = *(xpos+j1); 
  C1CPgrad<RCTMATTP>(q,pparasNV.begin(),JacobC1);  
  for (int i=0;i<npcovpar;i++) {
    for (int j=0;j<q;j++) Jacob(i,j) = JacobC1(i,j);
    for (int j=q,j1=p;j<q+qtrdim;j++,j1++) Jacob(i,j1) = JacobC1(i,j);
  }

  for (int j=0,j1=q;j<q;j++,j1++) pparasNV(j) = *(xpos+j1);
  for (int j=q,j1=p+qtrdim;j<q+qtrdim;j++,j1++) pparasNV(j) = *(xpos+j1);
  C1CPgrad<RCTMATTP>(q,pparasNV.begin(),JacobC1);
  for (int c=0,prind=0,frind=npcovpar;c<q;c++) {
    prind += c;
    frind += c+q;
    for (int r=0;r<=c;r++)  {
      for (int j=0,j1=q;j<q;j++,j1++) Jacob(frind+r,j1) = JacobC1(prind+r,j);
      for (int j=q,j1=p+qtrdim;j<q+qtrdim;j++,j1++) Jacob(frind+r,j1) = JacobC1(prind+r,j);
    }
  }
  return;
}

template<class RCTMATTP,class SQMATTP,class EXTVCTTP>
void RestCov_grad(const int p,const int q,const int nvcovpar,const int Config,
  NumericVector::iterator xpos,const bool FixedArrays,RCTMATTP& Jacob)
{  
  switch (Config) {
    case 1: C1CPgrad<RCTMATTP>(p,xpos,Jacob); break; 
    case 2: C2CPgrad<RCTMATTP,SQMATTP,EXTVCTTP>(p,q,nvcovpar,xpos,FixedArrays,Jacob);  break;
    case 3: C3CPgrad<RCTMATTP>(p,q,xpos,Jacob);  break;
    case 4: C4CPgrad<RCTMATTP,SQMATTP,EXTVCTTP>(p,q,nvcovpar,xpos,FixedArrays,Jacob);  break;
    case 5:  for (int i=0;i<p;i++) Jacob(utdind(i),i) = 2*(*(xpos+i));  break; 
  }
  return;
}


#endif

