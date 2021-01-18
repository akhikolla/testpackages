#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "cpploss.h"
#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>  

using namespace Rcpp;

// [[Rcpp::export]]
SEXP  GCPM_cpploss(SEXP default_distr_a,SEXP link_function_a, SEXP S_a,SEXP Sigma_a, SEXP W_a, SEXP PD_a, SEXP PL_a, SEXP calc_rc_a, SEXP loss_thr_a, SEXP max_entries_a){
  
  NumericMatrix S(S_a), W(W_a),Sigma(Sigma_a);
  NumericVector PD(PD_a), PL(PL_a),max_entries(max_entries_a),default_distr(default_distr_a),link_function(link_function_a),calc_rc(calc_rc_a),loss_thr(loss_thr_a);
  List ret;
  
  bool full=false;
  int N=S.nrow(), NS=S.ncol();
  int NCP=W.nrow(),defaults=0;
  int print_every=(int)(floor((float) (N/100))),inflate=1000;
  double condPD=0,temp=0,temp2=0;
  std::vector<double> lossszenarios(0);
  std::vector<std::vector<double> > CPsimlosses(NCP, std::vector<double>(0));
  if(calc_rc(0)==1){
    lossszenarios.reserve(inflate);
    for(int i=0;i<NCP;i++)
      CPsimlosses[i].reserve(inflate);
  }
  double *CPlosses=new double[NCP];
  NumericVector simlosses(N),state(1);
  GetRNGstate();
  Progress p((int)floor((float)(N/print_every)), true);  
  if(link_function(0)==1){ //CRP
    for(int n=0;n<N;n++){
      if(n%print_every==0 && n>0){
        p.increment();
      }
      if(!Progress::check_abort()){
        for(int i=0;i<NCP;i++){
          condPD=W(i,0)*PD[i];
          for(int k=0;k<NS;k++)
            condPD+=W(i,k+1)*S(n,k)*PD[i];
          if(default_distr(i)==1){ //Bernoulli
            condPD=std::min(1.0,std::max(0.0000001,condPD));
            temp=unif_rand();
            if(temp<=condPD){
              CPlosses[i]=PL(i);
              simlosses(n)+=PL(i);
            }
            else
              CPlosses[i]=0;
          }
          else{ //Poisson
            condPD=std::max(0.0000001,condPD);
            
            defaults=(int)R::rpois(condPD);
            simlosses(n)+=PL(i)*defaults;
            CPlosses[i]=PL(i)*defaults;
          }
        }
        if(calc_rc(0)==1 && !full){
          if(simlosses(n)>=loss_thr(0)){
            if(lossszenarios.capacity()==lossszenarios.size())
              lossszenarios.reserve(lossszenarios.capacity()+inflate);
            if(CPsimlosses[0].capacity()==CPsimlosses[0].size()){
              for(int i=0;i<NCP;i++)
                CPsimlosses[i].reserve(CPsimlosses[0].capacity()+inflate);
            }
            lossszenarios.push_back(n);
            for(int i=0;i<NCP;i++)
                CPsimlosses[i].push_back(CPlosses[i]);
            if(lossszenarios.size()>=max_entries(0))
              full=true;
          }
        }
      }
      else{
        simlosses(n)=-1;
      }
    }
  }
  else if(link_function(0)==2){ //CM
    double *DD=new double[NCP];
    for(int i=0;i<NCP;i++)
      DD[i]=R::qnorm(PD[i],0,1,1,0);
    
    for(int n=0;n<N;n++){
      if(n%print_every==0 && n>0){
        p.increment();
      }
      if(!Progress::check_abort()){
        for(int i=0;i<NCP;i++){
          temp=DD[i];
          for(int k=0;k<NS;k++){
            temp-=W(i,k)*S(n,k);
          }
          temp2=0;
          for(int k=0;k<NS;k++){
            for(int l=0;l<NS;l++)
              temp2+=W(i,k)*Sigma(k,l)*W(i,l);
          }
          condPD=R::pnorm(temp/sqrt(1-temp2),0,1,1,0);
          if(default_distr(i)==1){ //Bernoulli
            condPD=std::min(1.0,std::max(0.0000001,condPD));
            temp=unif_rand();
            if(temp<=condPD){
              CPlosses[i]=PL(i);
              simlosses(n)+=PL(i);
            }
            else
              CPlosses[i]=0;
          }
          else{ //Poisson
            condPD=std::max(0.0000001,condPD);
            defaults=(int)R::rpois(condPD);
            simlosses(n)+=PL(i)*defaults;
            CPlosses[i]=PL(i)*defaults; 
          }
        }
        if(calc_rc(0)==1 && !full){
          if(simlosses(n)>=loss_thr(0)){
            if(lossszenarios.capacity()==lossszenarios.size())
              lossszenarios.reserve(lossszenarios.capacity()+inflate);
            if(CPsimlosses[0].capacity()==CPsimlosses[0].size()){
              for(int i=0;i<NCP;i++)
                CPsimlosses[i].reserve(CPsimlosses[0].capacity()+inflate);
            }
            lossszenarios.push_back(n);
            for(int i=0;i<NCP;i++)
                CPsimlosses[i].push_back(CPlosses[i]);
            if(lossszenarios.size()>=max_entries(0))
              full=true;
          }
        }
      }
      else{
        simlosses(n)=-1;
      }
    }
    delete[] DD;
  }
  PutRNGstate();
  p.increment();

  NumericMatrix CPsimlossesfinal(NCP,lossszenarios.size());
  NumericVector lossszenariosfinal(lossszenarios.size());
  if(calc_rc(0)==1){
    if(full){
      for(int i=0;i<(int)lossszenarios.size();i++){
        for(int j=0;j<NCP;j++)
          CPsimlossesfinal(j,i)=-1;
        lossszenariosfinal(i)=-1;
      }
    }
    else{
      for(int i=0;i<(int)lossszenarios.size();i++){
        for(int j=0;j<NCP;j++)
          CPsimlossesfinal(j,i)=CPsimlosses[j][i];
        lossszenariosfinal(i)=lossszenarios[i]+1;
      } 
    }  
  }
  ret["simlosses"]=simlosses;  
  ret["CPsimlosses"]=CPsimlossesfinal;
  ret["lossszenarios"]=lossszenariosfinal;
  
  delete[] CPlosses;
  return ret;
}

