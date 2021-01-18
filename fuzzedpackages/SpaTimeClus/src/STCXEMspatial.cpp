#include "STCXEMspatial.h"


STCXEMspatial::STCXEMspatial(const S4 & input, const List & inputparam, const NumericMatrix & matT){
  m_data_p = new STCdata( as<S4>(input.slot("data")) );
  m_model_p = new STCmodel( as<S4>(input.slot("model")) );
  m_tune_p = new STCtune( as<S4>(input.slot("tune")) );
  m_paramCurrent_p = new STCparam( as<S4>(inputparam[0]) );
  m_paramlist.resize(m_tune_p->m_nbinitSmall);
  for (int ini=0; ini< m_tune_p->m_nbinitSmall ; ini++) m_paramlist[ini] = STCparam( as<S4>(inputparam[ini]) );
  m_loglikeSmall= ones<vec>(m_tune_p->m_nbinitSmall) * log(0);
  m_matT = as<mat>(matT); 
  m_poidspolynom= vec(m_model_p->m_K, fill::zeros);
  m_tig= mat(m_data_p->m_n, m_model_p->m_G, fill::zeros);
  m_sig.resize(m_model_p->m_G);
  for (int g=0; g<m_model_p->m_G; g++){
    m_sig[g].resize(m_model_p->m_K);
    for (int k=0; k<m_model_p->m_K; k++) m_sig[g][k] = cube( m_data_p->m_n, m_data_p->m_JJ, m_data_p->m_TT, fill::zeros);
  }
  
  m_hessian=mat(4*(m_model_p->m_K-1), 4*(m_model_p->m_K-1), fill::zeros);
  m_Mjte= cube(m_data_p->m_JJ, m_data_p->m_TT, 4, fill::ones); 
  m_Mjte.slice(1).each_col() =  m_data_p->m_map.col(0);
  m_Mjte.slice(2).each_col() =  m_data_p->m_map.col(1);
  m_Mjte.slice(3).each_row() =  trans(m_data_p->m_m);
}




double STCXEMspatial::ComputeLogLike(){
  m_tig.zeros();
  Mat<double> SumProbaRegression = zeros<mat>(m_data_p->m_n, m_model_p->m_K);
  Col<double> tmpsum= zeros<vec>(m_data_p->m_n);
  for (int g=0; g<m_model_p->m_G;g++){
    Mat<double> center = m_matT * trans(m_paramCurrent_p->m_beta[g]);
    for (int t=0; t<m_data_p->m_TT; t++){
      for (int j=0; j<m_data_p->m_JJ; j++){  
        m_poidspolynom = m_paramCurrent_p->m_lambda[g].col(0) + m_paramCurrent_p->m_lambda[g].cols(1,2) * trans(m_data_p->m_map.row(j)) + m_data_p->m_m(t) * m_paramCurrent_p->m_lambda[g].col(3);
        m_poidspolynom = exp(m_poidspolynom - max(m_poidspolynom));
        m_poidspolynom /= sum(m_poidspolynom);   
        m_poidspolynom = log(m_poidspolynom);
        for (int k=0; k< m_model_p->m_K; k++)  SumProbaRegression.col(k) = m_poidspolynom(k) - 0.5 * pow(m_data_p->m_x.slice(t).col(j) -  center(t, k), 2)/m_paramCurrent_p->m_sigma(g,k) - log(sqrt( 2*M_PI*m_paramCurrent_p->m_sigma(g,k)));
        tmpsum = max(SumProbaRegression, 1);
        SumProbaRegression = exp(SumProbaRegression.each_col() - tmpsum);
        for (int k=0; k< m_model_p->m_K; k++) m_sig[g][k].slice(t).col(j) = SumProbaRegression.col(k) ;
        m_tig.col(g) +=  tmpsum + log(sum(SumProbaRegression,1));
      } 
    }
    m_tig.col(g)  +=  log(m_paramCurrent_p->m_proportions(g));
  }
  tmpsum = max(m_tig, 1);
  m_tig = exp(m_tig.each_col() - tmpsum);
  double test=sum(tmpsum + log(sum(m_tig,1)));
  if ((test!=test) || (m_nondegeneracy==0)){
    test = log(0);
    m_nondegeneracy = 0 ;
  }
  return test;
}



void STCXEMspatial::NewtonLogitWeighted(const int g){
  Cube<double> u = cube(m_data_p->m_JJ, m_data_p->m_TT, m_model_p->m_K, fill::zeros); 
  Mat<double> totalsite = mat(m_data_p->m_JJ, m_data_p->m_TT, fill::zeros); 
  Cube<double> ratioexp = cube(m_data_p->m_JJ,  m_data_p->m_TT, m_model_p->m_K, fill::zeros);
  Col<double> gradient = vec(4*(m_model_p->m_K-1), fill::zeros) ;
  Col<double> currentparam=trans(vectorise(m_paramCurrent_p->m_lambda[g].rows(1, m_model_p->m_K-1),1));
  Col<double> tmp=trans(vectorise(m_paramCurrent_p->m_lambda[g].rows(1, m_model_p->m_K-1),1));
  for (int t=0; t< m_data_p->m_TT; t++){
    Mat<double> tmp = m_data_p->m_map  * trans(m_paramCurrent_p->m_lambda[g].cols(1,2)) ;
    tmp.each_row() +=  trans(m_paramCurrent_p->m_lambda[g].col(0)) + trans(m_paramCurrent_p->m_lambda[g].col(3) * m_data_p->m_m(t));
    tmp =  exp(tmp.each_col() - max(tmp, 1));
    tmp.each_col() /= sum(tmp,1);
    for (int k=0; k<m_model_p->m_K; k++){
      u.slice(k).col(t) += trans(sum(m_sig[g][k].slice(t)));
      totalsite.col(t) += u.slice(k).col(t);
      ratioexp.slice(k).col(t) = tmp.col(k);
    } 
  }
  for (int k1=0; k1< (m_model_p->m_K-1); k1++){
    for (int k2=k1; k2< (m_model_p->m_K-1); k2++){
      for (int h1=0; h1 < 4; h1++){
        for (int h2=0; h2 < 4; h2++){
          m_hessian(k1*4 + h1, k2*4 + h2) = accu(totalsite % m_Mjte.slice(h1) % m_Mjte.slice(h2) % ratioexp.slice(k1+1) % ratioexp.slice(k2+1));
          m_hessian(k2*4 + h2, k1*4 + h1) = m_hessian(k1*4 + h1, k2*4 + h2);
        }
      }
    }
    for (int h1=0; h1 < 4; h1++){
      for (int h2=h1; h2 < 4; h2++){
        m_hessian(k1*4 + h1, k1*4 + h2) -= accu(totalsite % m_Mjte.slice(h1) % m_Mjte.slice(h2) % ratioexp.slice(k1+1));
        m_hessian(k1*4 + h2, k1*4 + h1) = m_hessian(k1*4 + h1, k1*4 + h2);
      }
      gradient(k1*4+h1) = accu((u.slice(k1+1) - totalsite % ratioexp.slice(k1+1)) % m_Mjte.slice(h1) );
    }
  }   
  m_nondegeneracy = m_nondegeneracy * solve(tmp, m_hessian,  gradient, solve_opts::no_approx);
  if (m_nondegeneracy == 1){
    tmp = currentparam - tmp;
    for (int k=0; k< (m_model_p->m_K-1); k++) m_paramCurrent_p->m_lambda[g].row(k+1) = trans(tmp.subvec(k*4, (k*4)+3));
  }
}

void STCXEMspatial::Mstep(){
  m_paramCurrent_p->m_proportions = trans(sum(m_tig,0) / m_data_p->m_n);
  for (int g=0; g<m_model_p->m_G; g++){
    for (int k=0; k< m_model_p->m_K; k++){
      Col<double> weight = vec(m_data_p->m_TT, fill::ones);
      Col<double> weightX = vec(m_data_p->m_TT, fill::ones);
      for (int t=0; t<m_data_p->m_TT; t++){
        weight(t) = accu(m_sig[g][k].slice(t));
        weightX(t) = accu(m_data_p->m_x.slice(t) % m_sig[g][k].slice(t));
      }
      Col <double> inv=trans(m_paramCurrent_p->m_beta[g].row(k));
      m_nondegeneracy = m_nondegeneracy * solve(inv, trans(m_matT) * diagmat(weight) * m_matT, trans(m_matT) * (weightX), solve_opts::no_approx);
      if (m_nondegeneracy == 1) m_paramCurrent_p->m_beta[g].row(k) = trans(inv);     
      Row<double> tmpmean= trans(m_matT * trans(m_paramCurrent_p->m_beta[g].row(k)));
      for (int t=0; t<m_data_p->m_TT; t++){
        weight(t) = accu(pow(m_data_p->m_x.slice(t) - tmpmean(t), 2) % m_sig[g][k].slice(t));
        weightX(t) = accu(m_sig[g][k].slice(t));
      }
      m_paramCurrent_p->m_sigma(g,k) = accu(weight) / accu(weightX);
      if (m_paramCurrent_p->m_sigma(g,k) < 0.000001)m_nondegeneracy=0;
    }
    if (m_model_p->m_K>1) { NewtonLogitWeighted(g);}
  }
}
