#include <iostream>
#include <string>
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include "FBFsearch_H.h"

using namespace std;
using namespace arma;


//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------
//FUNZIONI MATEMATICHE GENERALI
//------------------------------------------------------------------------------------------------------------------------

//FUNZIONE LOGARITMO FATTORIALE

double lfactorial(int n)
{
  
 double x;

 x=lgamma(n+1);

 return(x);

} //end function


//FUNZIONE LOGARITMO COEFFICIENTE BINOMIALE

double lchoose(int n, int k)
{
  
  double x;
  
  x=lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);

  return(x);

} //end function


//FUNZIONE LOGARITMO DELLA SOMMA

double log_sum(colvec v)
{

  double lsum; 
 
  if(accu(v>=0)>0){
   if(max(v)<=700){lsum=log(accu(exp(v)));} else{lsum=max(v);}
  } else{
   if(min(v)>=-700){lsum=log(accu(exp(v)));} else {lsum=max(v);} 
  }

  return lsum;

} //end function


//FUNZIONE CHE RESTITUISCE UN SOTTOINSIEME DI UN VETTORE DATA UNA CONDIZIONE

vec sub_elem_eq(vec v, vec w, double x)
{
  uvec iw;

  iw=find(w==x);
  if(!iw.is_empty()){v=v.elem(iw);}else{v=datum::nan;}
   
  return v;

} //end function


//FUNZIONE CHE RESTITUISCE UNA SOTTOMATRICE CON COLONNE E RIGHE INDICATE IN UN VETTORE

mat sub_mat(mat M, vec vr, vec vc)
{
  int r, c, lvr, lvc;
  lvr=vr.n_elem;
  lvc=vc.n_elem;
  mat Q(lvr,lvc);

  for(c=0; c<lvc; c++){
   for(r=0; r<lvr; r++){ 
    Q(r,c)=M(vr(r),vc(c)); 
   }
  }

  return Q;

} //end function



//FUNZIONE CHE FA L'ELEVAMENTO A POTENZA ELEMENTO PER ELEMENTO

mat pow_vec(vec v, vec w)
{

   int k, lv; 
   vec vv;

   lv=v.n_elem;
   vv=zeros<vec>(lv);

   for(k=0; k<lv; k++){
    vv(k)=pow(v(k),w(k));
   }

 return vv;

}

//------------------------------------------------------------------------------------------------------------------------
//BINARY TREE
//------------------------------------------------------------------------------------------------------------------------

//FUNZIONE CHE AGGIUNGE IL MODELLO M ALL'ALBERO

field<mat> add_to_tree(vec M, double lM, int nM, mat tree, double ltree)
{

 int j;
 int k;

 field<mat> Res(2,1);

 if(nM==0){
    
  for(j=0; j<=lM; j++){
 
   if(M(j)==1){tree(j,0)=j+1;}
   else{tree(j,1)=j+1;}
 
  }
 
  ltree=lM; 
  tree.row(ltree)=tree.row(ltree)*0+nM;

 }

 int z;
 int h;
 int iM=datum::nan;

 if(nM>0){ //if1
   
  z=0;
  h=ltree+1;
  
  for(j=0; j<=lM; j++){ //for1  
  
   iM=1-M(j);
      
   if(!is_finite(tree(z,iM)) && (j<=lM) ){ //if2
     
    tree(z,iM)=h;
     
    for(k=(j+1); k<=lM; k++){ //for2
  
      if(M(k)==1){tree(h,0)=h+1;} else{tree(h,1)=h+1;}
      
      h=h+1;
    
    } //end for2    

    iM=1-M(lM);
    ltree=h-1;
    break;
   
   } //end if2
 
   if(j==lM){tree(z,iM)=nM; ltree=ltree+1; break;}
  
   if(tree(z,iM)>=0){z=tree(z,iM);}
  
  } //end for1
  

  tree(ltree,iM)=tree(ltree,iM)*0+nM;

 } //end if1


 Res(0,0)=tree;
 Res(1,0)=ltree;

 return Res;

} //end function






//FUNZIONE CHE RESTITUISCE I MOVIMENTI POSSIBILI DA M DATI I MODELLI PRECEDENTEMENTE VISITATI

vec mov_tree(mat tree, vec M, int lM, vec vlM, int max_lM)
{

 int q; int k; int z; int h; int iM2;
 double sumM;
 vec mov(lM+1); uvec imov; uvec umov; vec mov2; vec M2;
 
 
 mov.fill(-1);
 sumM=sum(M);
 q=0;

 for(k=0; k<=lM; k++){ //for1
  
  M2=M; 
  M2(k)=1-M(k);
  z=0;
  
  for(h=0; h<=lM; h++){ //for2
    
   iM2=1-M2(h);
   if(!is_finite(tree(z,iM2))){mov(q)=k; q=q+1; break;} else{z=tree(z,iM2);}
   
    } //end for2
  
 
 } //end for1

 imov=find(mov>-1);
 
 if(!imov.is_empty()){
 
  mov=mov.elem(imov);
  umov=conv_to<uvec>::from(mov);

  if(sumM>=max_lM){
  
   mov2=zeros<vec>(lM+1);
   mov2.elem(umov)=ones<vec>(mov.n_elem);
   mov=(mov2%M)%vlM;
   imov=find(mov>0);
   if(!imov.is_empty()){mov=mov.elem(imov); mov=mov-1;} else{mov=datum::nan;}
  }

 } else {mov=datum::nan;}
 

 return mov;

} //end function



//------------------------------------------------------------------------------------------------------------------------
//LEGATE AL PROBLEMA
//------------------------------------------------------------------------------------------------------------------------


double log_H_h_i(double mu, double sigma, int h, int i)
{

 double x;

 x=lfactorial(2*h)+i*log(sigma)-lfactorial(i)+(2*h-2*i)*log(abs(mu))-lfactorial(2*h-2*i);

 return x;

} //end function



double log_FBF_Ga_Gb(vec G_a, vec G_b, int edge, mat edges, mat YtY, int add, double n, double h)
{

  int e1, e2, i, iwi;
  double p, b, S2, mu, sigma, logS2, ilogS2, logHhi, ilog4, log_num1, log_den1, log_w_1, log_num0i0, log_den0i0, log_w_0, log_FBF_unpasso;
  vec V1, V2, G1, V11, pa1, pa0, betah, vv(1), z1;
  uvec iw, ipa1;
  mat e, yty, XtX, Xty, invXtX;   
  
  // std::ostream nullstream(0);
  // set_stream_err2(nullstream);

  e=edges.row(edge);
  e1=e(0);
  e2=e(1);

  V1=edges.col(0);
  V2=edges.col(1);
    
  if(add==1){G1=G_a;}else{G1=G_b;}   

  V11=(V1+1)%G1;
  iw=find(V2==e2); pa1=V11.elem(iw);
  iw=find(pa1>0); pa1=pa1.elem(iw); pa1=pa1-1;
 
  iw=find(pa1!=e1); if(!iw.is_empty()){pa0=pa1.elem(iw);}else{pa0=datum::nan;}
 
  p=pa1.n_elem;
  b=(p+2*h+1)/n;

  yty=YtY(e2,e2);
    
  //calcolo w1
   
  vv(0)=e2; Xty=sub_mat(YtY,pa1,vv);
  XtX=sub_mat(YtY,pa1,pa1);
  betah=solve(XtX,Xty);
  //invXtX=inv_sympd(XtX);
  //betah=invXtX*Xty;
  S2=conv_to<double>::from(yty-(trans(Xty)*betah));
  
  iw=find(pa1==e1); mu=conv_to<double>::from(betah.elem(iw));
  iwi=conv_to<int>::from(iw); 
  //_sigma=invXtX(iwi,iwi); sigma=conv_to<double>::from(_sigma);	
  z1=zeros<vec>(pa1.n_elem);
  z1(iwi)=1;
  z1=solve(XtX,z1);
  sigma=conv_to<double>::from(z1.elem(iw));
  
  if(S2>0){
   log_w_1=(-n*(1-b)/2)*log(datum::pi*b*S2);
   logS2=log(S2);
   log_num1=-datum::inf;
   log_den1=-datum::inf;

   for(i=0; i<=h; i++){
   
    ilogS2=i*logS2;
    logHhi=log_H_h_i(mu,sigma,h,i);
    ilog4=-i*log(4);
      
    log_num1=log_add(log_num1, (ilog4+logHhi+lgamma((n-p-2*i)/2)+ilogS2));
    log_den1=log_add(log_den1, (ilog4+logHhi+lgamma((n*b-p-2*i)/2)+ilogS2));
  
   }
    
   log_w_1=log_w_1+log_num1-log_den1;
  }else{
   log_w_1=datum::inf;
  }
     
  //calcolo w0

  if(!pa0.is_finite()){p=0;}else{p=pa0.n_elem;}
  
  log_num0i0=lgamma((n-p)/2);
  log_den0i0=lgamma((n*b-p)/2);

  if(p==0){S2=conv_to<double>::from(yty);}
  else{
   vv(0)=e2; Xty=sub_mat(YtY,pa0,vv);
   XtX=sub_mat(YtY,pa0,pa0);
   betah=solve(XtX,Xty);
   S2=conv_to<double>::from(yty-(trans(Xty)*betah));
  }

  if(S2>0){
   log_w_0=(-(n*(1-b)/2))*log(datum::pi*b*S2)+log_num0i0-log_den0i0;
  }else{
   log_w_0=datum::inf;
  }
  
  //calcolo FBF

  if(add==1){log_FBF_unpasso=log_w_1-log_w_0;}
  else{log_FBF_unpasso=log_w_0-log_w_1;} 
  
  if(!is_finite(log_FBF_unpasso)){
   log_FBF_unpasso=0;	  
  }  
 
  return log_FBF_unpasso;

} // end function




field<mat> FBF_heart(double nt, mat YtY, vec vG_base, double lcv, vec vlcv, mat edges, double n_tot_mod, double C, double maxne, double h)
{
 
   int t, add, edge, imq, limodR, s;
   double ltree, lM, sum_log_FBF, log_FBF_G, log_pi_G, log_num_MP_G, sum_log_RSMP, n_mod_r, log_FBF_t, log_FBF1;
   vec M_log_FBF, log_num_MP, log_sume, G, imod_R, M_log_RSMP, pRSMP, mov, vlM, qh, G_t, M_q, M_P;
   uvec iw;
   mat tree, SM, M_G; 
   field<mat> treeRes, Res(3,1);
   uword i_n_mod_r, imaxe;
 
   M_G=zeros<mat>(lcv,n_tot_mod);  
   M_P=zeros<vec>(n_tot_mod); 
   M_log_FBF=zeros<vec>(n_tot_mod); 
   log_num_MP=zeros<vec>(n_tot_mod); 
   M_q=zeros<vec>(lcv); 
   tree=zeros<mat>(n_tot_mod*lcv,2); tree.fill(datum::nan);
   ltree=datum::nan;
   lM=lcv-1; 
 
   sum_log_FBF=-datum::inf; 
   log_sume=zeros<vec>(lcv); log_sume.fill(-datum::inf);
  
   M_log_RSMP=zeros<vec>(n_tot_mod);
   sum_log_RSMP=-datum::inf;
   imod_R=zeros<vec>(n_tot_mod);

   lM=lcv-1;

   for(t=0; t<lcv; t++){ //for1
       
    G=vG_base;
    G(t)=1-vG_base(t);
    add=G(t);
    edge=t;

    log_FBF_G=log_FBF_Ga_Gb(G,vG_base,edge,edges,YtY,add,nt,h);
   
    M_G.col(t)=G;
    
    treeRes=add_to_tree(G,lM,t,tree,ltree);
    tree=treeRes(0,0);
    ltree=conv_to<double>::from(treeRes(1,0));
        
    M_log_FBF(t)=log_FBF_G;
    log_pi_G=-log(lcv+1)-lchoose(lcv,sum(G));
    log_num_MP_G=log_FBF_G+log_pi_G;
    log_num_MP(t)=log_num_MP_G;

    sum_log_FBF=log_add(sum_log_FBF, log_num_MP_G);
  
    for(imq=0; imq<lcv; imq++){
     if(G(imq)==1){log_sume(imq)=log_add(log_sume(imq), log_num_MP_G);}
    }
    
    M_q=exp(log_sume-sum_log_FBF);
   
    M_log_RSMP(t)=log_num_MP_G;
    sum_log_RSMP=log_add(sum_log_RSMP, log_num_MP_G);
   
    imod_R(t)=t;
 
   } //end for1

 
   limodR=t-1; 
   s=lcv;

    
   while(t<n_tot_mod){ //while1
   

    pRSMP=exp(M_log_RSMP.subvec(0,limodR)-sum_log_RSMP);
	pRSMP.max(i_n_mod_r);
    
    n_mod_r=imod_R(i_n_mod_r);
        
    G=M_G.col(n_mod_r);
    G_t=G;
    log_FBF_t=M_log_FBF(n_mod_r);
    
    
    vlM=vlcv+1;
    mov=mov_tree(tree,G,lM,vlM,maxne);
   
    if(!is_finite(mov)){ //if1
      
      imod_R(i_n_mod_r)=-1;
      iw=find(imod_R>-1); imod_R=imod_R.elem(iw);
      M_log_RSMP=M_log_RSMP.elem(iw);
          
      limodR=limodR-1;
      t=t-1;
        
    } else{
        
       qh=pow_vec((M_q+C)/(1-M_q+C), (2*(1-G))-1);
       qh=qh.elem(conv_to<uvec>::from(mov));
    
        
       if(mov.n_elem==1){ //if2
           
        imod_R(i_n_mod_r)=-1;
        iw=find(imod_R>-1); imod_R=imod_R.elem(iw);
        M_log_RSMP=M_log_RSMP.elem(iw);
      
         limodR=limodR-1;  
         edge=mov(0);
         
        } else{
           
           qh.max(imaxe);
           edge=mov(imaxe);
          
        } // end if2
    
        
        G(edge)=1-G(edge);
        add=G(edge);   
       
        log_FBF1=log_FBF_Ga_Gb(G,G_t,edge,edges,YtY,add,nt,h);
	    log_FBF_G=log_FBF1+log_FBF_t;
      
        M_G.col(t)=G;
        
        treeRes=add_to_tree(G,lM,t,tree,ltree);
        tree=treeRes(0,0);
        ltree=conv_to<double>::from(treeRes(1,0));
      
        M_log_FBF(t)=log_FBF_G;
        log_pi_G=-log(lcv+1)-lchoose(lcv,sum(G));
        log_num_MP_G=log_FBF_G+log_pi_G;
        log_num_MP(t)=log_num_MP_G;
    
        sum_log_FBF=log_add(sum_log_FBF, log_num_MP_G);
        
        for(imq=0; imq<lcv; imq++){
         if(G(imq)==1){log_sume(imq)=log_add(log_sume(imq), log_num_MP_G);}
        }
        M_q=exp(log_sume-sum_log_FBF);
	    
        limodR=limodR+1;
        imod_R(limodR)=t;
        M_log_RSMP(limodR)=log_num_MP_G; 
         
    } //end if1
    
    t=t+1;
    s=s+1;
    
  } //end while1 
   

  t=t-1; 
  s=s-1;
  
  M_P.subvec(0,t)=exp(log_num_MP.subvec(0,t)-sum_log_FBF);
  if(max(M_P.subvec(0,t))>0){
   M_P.subvec(0,t)=M_P.subvec(0,t)/sum(M_P.subvec(0,t));
  }else{
   M_P.subvec(0,t)=zeros<vec>(t+1);	  	  
  }
  
  M_G=M_G.submat(0,0,lcv-1,t);
  
  Res(0,0)=M_q;
  Res(1,0)=M_G;
  Res(2,0)=M_P;

  return Res;
   

} // end function




// FUNZIONE CHE RIEMPE LA MATRICE G_fin con gli elementi di M_q 

mat G_fin_fill(mat G, vec vr, int ic, vec x)
{
  int k, lvr;
  lvr=vr.n_elem;
 
  for(k=0; k<lvr; k++){
    G(vr(k),ic)=x(k); 
  }

  return G;

} //end function



//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------

//LOCAL SEARCH

RcppExport SEXP FBF_LS(SEXP Corr, SEXP nobs, SEXP G_base, SEXP h, SEXP C, SEXP n_tot_mod)
{

 int k, neq, rr; 
 double maxne, Mlogbin_sum, lcv, rrmax, q;
 vec V1, V2, vlcv, vG_base, M_q;
 mat edges, G_fin;
 field<mat> heartRes;
 
 Rcpp::NumericMatrix Corr_c1(Corr);
 Rcpp::NumericVector nobs_c1(nobs);
 Rcpp::NumericMatrix G_base_c1(G_base);
 Rcpp::NumericVector h_c1(h);
 Rcpp::NumericVector C_c1(C);
 Rcpp::NumericVector n_tot_mod_c1(n_tot_mod);

 q=Corr_c1.ncol();

 mat Corr_c(Corr_c1.begin(), q, q, false);
 mat G_base_c(G_base_c1.begin(), q, q, false);
 double nobs_c=nobs_c1[0];
 double h_c=h_c1[0];
 double C_c=C_c1[0];
 double n_tot_mod_c=n_tot_mod_c1[0];


 q=Corr_c.n_cols; 

 G_fin=zeros<mat>(q,q); 
 maxne=nobs_c-2*h_c-2;

 for(k=0; k<(q-1); k++){ //for1
   
   neq=k+1;

   V1=linspace<vec>(q-neq-1,0,q-neq);
   V2=zeros<vec>(q-neq); V2.fill(q-neq);
   edges=join_rows(V1,V2); 
   
   lcv=V1.n_elem;
   vlcv=linspace<vec>(0,lcv-1,lcv); 
  
   vG_base=flipud(G_base_c.submat(0,q-neq,q-neq-1,q-neq));
   
   rrmax=std::min(maxne,lcv); 
   Mlogbin_sum=0;
   
   for(rr=1; rr<=rrmax; rr++){ //for2
     Mlogbin_sum=log_add(Mlogbin_sum,lchoose(lcv,rr));   
   } //end for2
 
   n_tot_mod_c=std::min(Mlogbin_sum,log(n_tot_mod_c));
   n_tot_mod_c=round(exp(n_tot_mod_c)); 

   heartRes=FBF_heart(nobs_c, Corr_c*nobs_c, vG_base, lcv, vlcv, edges, n_tot_mod_c, C_c, maxne, h_c);
   M_q=heartRes(0,0);
   
   G_fin=G_fin_fill(G_fin,V1,V2[0],M_q);

  } //end for1

 
  return Rcpp::wrap(G_fin);

 
} // end function



// FUNZIONE FINALE REGRESSION SEARCH 

RcppExport SEXP FBF_RS(SEXP Corr, SEXP nobs, SEXP G_base, SEXP h, SEXP C, SEXP n_tot_mod, SEXP n_hpp)
{

 int neq, rr, j; 
 double maxne, Mlogbin_sum, lcv, rrmax, q;
 vec V1, V2, vlcv, vG_base, M_q, M_P, iM_P, M_P2;
 mat edges, G_fin, M_G, M_G2;
 field<mat> heartRes;
 
 Rcpp::NumericMatrix Corr_c1(Corr);
 Rcpp::NumericVector nobs_c1(nobs);
 Rcpp::NumericVector G_base_c1(G_base);
 Rcpp::NumericVector h_c1(h);
 Rcpp::NumericVector C_c1(C);
 Rcpp::NumericVector n_tot_mod_c1(n_tot_mod);
 Rcpp::NumericVector n_hpp_c1(n_hpp);

 q=Corr_c1.ncol();

 mat Corr_c(Corr_c1.begin(), q, q, false);
 vec G_base_c(G_base_c1.begin(), q-1, false);
 double nobs_c=nobs_c1[0];
 double h_c=h_c1[0];
 double C_c=C_c1[0];
 double n_tot_mod_c=n_tot_mod_c1[0];
 double n_hpp_c=n_hpp_c1[0];

 q=Corr_c.n_cols; 

 maxne=nobs_c-2*h_c-2;

 neq=1;

 V1=linspace<vec>(q-neq-1,0,q-neq);
 V2=zeros<vec>(q-neq); V2.fill(q-neq);
 edges=join_rows(V1,V2); 
   
 lcv=V1.n_elem;
 vlcv=linspace<vec>(0,lcv-1,lcv); 
  
 vG_base=flipud(G_base_c);
   
 rrmax=std::min(maxne,lcv); 
 Mlogbin_sum=0;
   
 for(rr=1; rr<=rrmax; rr++){
  Mlogbin_sum=log_add(Mlogbin_sum,lchoose(lcv,rr));   
 }
 
 n_tot_mod_c=std::min(Mlogbin_sum,log(n_tot_mod_c));
 n_tot_mod_c=round(exp(n_tot_mod_c)); 

 heartRes=FBF_heart(nobs_c, Corr_c*nobs_c, vG_base, lcv, vlcv, edges, n_tot_mod_c, C_c, maxne, h_c);
 
 M_q=heartRes(0,0);
 M_G=heartRes(1,0);
 M_P=heartRes(2,0);
 
 iM_P=conv_to<vec>::from(sort_index(conv_to<vec>::from(M_P), "descend"));

 M_P2=zeros<vec>(n_hpp_c);
 M_G2=zeros<mat>(lcv,n_hpp_c);

 for(j=0; j<n_hpp_c; j++){
  M_P2(j)=M_P(iM_P(j));
  M_G2.col(j)=M_G.col(iM_P(j));
 }

 
 return Rcpp::List::create(Rcpp::Named("M_q") = M_q,
                           Rcpp::Named("M_G") = M_G2,
                           Rcpp::Named("M_P") = M_P2);


} // end function



// FUNZIONE FINALE GLOBAL SEARCH 

RcppExport SEXP FBF_GS(SEXP Corr, SEXP nobs, SEXP G_base, SEXP h, SEXP C, SEXP n_tot_mod, SEXP n_hpp)
{

  int j, k, neq, rr; 
  double rrmax, Mlogbin_sum, maxne, lcv, nc_edges, nc_edges2, q;
  vec V1, V2, vlcv, vG_base, M_q, M_P, iM_P, M_P2, M_G_j;
  mat edges, G_fin, M_q2, M_G, M_G2, G_fin2, M_G2_j;
  field<mat> heartRes;
  
  Rcpp::NumericMatrix Corr_c1(Corr);
  Rcpp::NumericVector nobs_c1(nobs);
  Rcpp::NumericMatrix G_base_c1(G_base);
  Rcpp::NumericVector h_c1(h);
  Rcpp::NumericVector C_c1(C);
  Rcpp::NumericVector n_tot_mod_c1(n_tot_mod);
  Rcpp::NumericVector n_hpp_c1(n_hpp);

  q=Corr_c1.ncol();

  mat Corr_c(Corr_c1.begin(), q, q, false);
  mat G_base_c(G_base_c1.begin(), q, q, false);
  double nobs_c=nobs_c1[0];
  double h_c=h_c1[0];
  double C_c=C_c1[0];
  double n_tot_mod_c=n_tot_mod_c1[0];
  double n_hpp_c=n_hpp_c1[0]; 

  q=Corr_c.n_cols; 

  lcv=q*(q-1)/2;
  vlcv=linspace<vec>(0,lcv-1,lcv);

  maxne=nobs_c-2*h_c-2;
  
  edges=zeros<mat>(lcv,2);
  vG_base=zeros<vec>(lcv);

  nc_edges=0;

  for(k=0; k<(q-1); k++){ //for1
  
   neq=k+1;

   V1=linspace<vec>(q-neq-1,0,q-neq);
   V2=zeros<vec>(q-neq); V2.fill(q-neq);
   
   nc_edges2=nc_edges+(q-neq); 
  
   edges.rows(nc_edges,nc_edges2-1)=join_rows(V1,V2); 
   vG_base.subvec(nc_edges,nc_edges2-1)=flipud(G_base_c.submat(0,q-neq,q-neq-1,q-neq));

   nc_edges=nc_edges2;

  } //end for1

  rrmax=std::min(maxne,lcv);     
  Mlogbin_sum=0;
   
  for(rr=1; rr<=rrmax; rr++){
   Mlogbin_sum=log_add(Mlogbin_sum,lchoose(lcv,rr));   
  }
 
  n_tot_mod_c=std::min(Mlogbin_sum,log(n_tot_mod_c));
  n_tot_mod_c=round(exp(n_tot_mod_c)); 

  heartRes=FBF_heart(nobs_c, Corr_c*nobs_c, vG_base, lcv, vlcv, edges, n_tot_mod_c, C_c, maxne, h_c);
  
  M_q=conv_to<vec>::from(heartRes(0,0));
  M_G=heartRes(1,0);
  M_P=conv_to<vec>::from(heartRes(2,0));

 
  M_q2=zeros<mat>(q,q);
  for(j=0; j<lcv; j++){
    M_q2(edges(j,0),edges(j,1))=M_q(j);
  }

  iM_P=conv_to<vec>::from(sort_index(conv_to<vec>::from(M_P), "descend"));

  M_P2=zeros<vec>(n_hpp_c);
  M_G2_j=zeros<mat>(lcv,n_hpp_c);
  M_G2=zeros<mat>(q*lcv,q);
  
  for(j=0; j<n_hpp_c; j++){ //for2
  
   M_P2(j)=M_P(iM_P(j));
  
   M_G_j=M_G.col(iM_P(j));
   M_G2_j=zeros<mat>(q,q);
   
   for(k=0; k<lcv; k++){
    M_G2_j(edges(k,0),edges(k,1))=M_G_j(k);
   }

   M_G2.rows(j*q, (j+1)*q-1)=M_G2_j;

  } //end for2
  
 
  return Rcpp::List::create(Rcpp::Named("M_q2") = M_q2,
                            Rcpp::Named("M_G2") = M_G2,
                            Rcpp::Named("M_P2") = M_P2);

  
} // end function



