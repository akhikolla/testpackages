#include "STCXEM.h"


void STCXEM::Estep(){
  m_tig = m_tig.each_col() / sum(m_tig,1);
  for (int g=0; g<m_model_p->m_G; g++){
    Cube<double> Normalise = cube(m_data_p->m_n, m_data_p->m_JJ, m_data_p->m_TT, fill::zeros);
    for (int k=0; k< m_model_p->m_K; k++)  Normalise = Normalise + m_sig[g][k];    
    if (any(m_tig.col(g) ==0)) Normalise.elem(find(Normalise==0)).ones();
    for (int k=0; k< m_model_p->m_K; k++){
      m_sig[g][k] /= Normalise;
      for (int t=0; t< m_data_p->m_TT; t++) m_sig[g][k].slice(t) = m_sig[g][k].slice(t).each_col() % m_tig.col(g);
    } 
  } 
}

void STCXEM::Output(S4 * reference_p){
  as<S4>(reference_p->slot("criteria")).slot("loglike") = wrap(max(m_loglikeSmall));
  as<S4>(reference_p->slot("criteria")).slot("degeneracy") = wrap(cpdegeneracy);
  as<S4>(reference_p->slot("param")).slot("proportions") = wrap(trans(m_paramCurrent_p->m_proportions));
  as<S4>(reference_p->slot("param")).slot("lambda") = wrap(m_paramCurrent_p->m_lambda);
  as<S4>(reference_p->slot("param")).slot("sigma") = wrap(m_paramCurrent_p->m_sigma);
  as<S4>(reference_p->slot("param")).slot("beta") = wrap(m_paramCurrent_p->m_beta);
  Estep();
  as<S4>(reference_p->slot("partitions")).slot("fuzzyind") = wrap(m_tig);
}

void STCXEM::OneEM(const int itermax, const double tol){
  m_nondegeneracy=1;
  double loglike = ComputeLogLike(), prec = log(0);
  int it=0;
  while ( (it<itermax) && ((loglike-prec)>tol) && (m_nondegeneracy==1) ){
    it ++;
    Estep();
    Mstep();
    if (m_nondegeneracy==1){
      prec = loglike;
      loglike = ComputeLogLike();
     // if (loglike< (prec-tol)) cout << "Error in EM algorithm (loglikelihood decreases)" << endl << "diff: "<< loglike - prec<< " loglike: " << loglike<< " prec: " << prec << " tol" << tol << " iteration " << it << endl;
    }else{
      loglike = log(0);
    }
  } 
}

void STCXEM::Run(){  
  cpdegeneracy=0;
  for (int ini=0; ini< m_tune_p->m_nbinitSmall ; ini++){
    SwitchParamCurrent(ini);
    OneEM(m_tune_p->m_nbiterSmall, m_tune_p->m_tol);
    m_loglikeSmall(ini) = ComputeLogLike();  
    cpdegeneracy = 1 - m_nondegeneracy;
  }
  //cout << "nb degenerate" << cpdegeneracy << endl;
  
  uvec indices = sort_index(m_loglikeSmall);
  double tmp1=0;
  cpdegeneracy=0;
  while (tmp1< m_tune_p->m_nbinitKept){ 
    //while (tmp1< min(m_tune_p->m_nbinitKept+cpdegeneracy, maxikeep)){
    SwitchParamCurrent(indices(m_tune_p->m_nbinitSmall - tmp1 - 1));
    OneEM(m_tune_p->m_nbiterKept, m_tune_p->m_tol);
    m_loglikeSmall(indices(m_tune_p->m_nbinitSmall - tmp1 - 1)) = ComputeLogLike(); 
    cpdegeneracy += 1 - m_nondegeneracy;
    tmp1++;
  }
  if (cpdegeneracy != m_tune_p->m_nbinitSmall){  
    indices = sort_index(m_loglikeSmall);
    SwitchParamCurrent(indices(m_tune_p->m_nbinitSmall  - 1));
    //OneEM(m_tune_p->m_nbiterKept, m_tune_p->m_tol);
  }
  cpdegeneracy = cpdegeneracy / tmp1;
}

