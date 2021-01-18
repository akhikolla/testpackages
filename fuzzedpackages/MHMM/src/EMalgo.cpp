#include "EMalgo.h"
#include "Param.h"
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/tools/minima.hpp>

using namespace boost::math;

double EMalgo::FunctionGammaToOptimize(double a){
  return abs(digamma(a) - log(a) - m_cst);
}

double ComputeLogSumVec(const Col<double> & tmpVec){
  double mymax = tmpVec.max();
  return mymax + log( sum( exp(tmpVec - mymax)));
}

Col<double> ComputeLogRowMat(const Mat<double> & tmpMat){
  Col<double> norm = max(tmpMat, 1);
  return norm + log(sum(exp(tmpMat.each_col() - norm), 1));
}

Row<double> ComputeLogColMat(const Mat<double> & tmpMat){
  Row<double> norm = max(tmpMat, 0);
  return norm + log(sum(exp(tmpMat.each_row() - norm), 0));
}

Col<double> NormLogVec(const Col<double> & tmpVec){
  return exp(tmpVec - ComputeLogSumVec(tmpVec));
}


double ComputeLogSumMat(const Mat<double> & tmpMat){
  double mymax = tmpMat.max();
  return mymax + log( sum(sum( exp(tmpMat - mymax))));
}

Mat<double> NormLogMatFull(const Mat<double> & tmpMat){
  Mat<double> out =   exp(tmpMat - ComputeLogSumMat(tmpMat));
  return out;
}

Mat<double> NormLogMatRow(const Mat<double> & tmpMat){
  Mat<double> out=tmpMat * 0;
  for (int t=0; t<tmpMat.n_rows; t++){
    out.row(t) = trans(NormLogVec(trans(tmpMat.row(t))));
  }
  return out;
}

Mat<double> NormLogMatRow2(const Mat<double> & tmpMat){
  Col<double> norm = max(tmpMat, 1);
  norm = norm + log(sum(exp(tmpMat.each_col() - norm), 1));
  Mat<double> out=tmpMat.each_col() - norm;
  return out;
}


EMalgo::EMalgo(const Data * data, const List paramR, const double tolR, const int nbKeep, const int iterSmall){
  m_tol = tolR;
  m_nbKeep = nbKeep;
  m_nbiter = iterSmall;
  m_data_p = data;
  for (int i=0;i<paramR.size();i++) m_paramCand.push_back( Param(as<S4>(paramR[i]) ));
  m_allloglike = zeros(m_paramCand.size());
  SwitchParamCurrent(0);
  m_tau = zeros<mat>(m_data_p->m_nobs, m_param_current_p->m_K);
  m_maxtmplogproba = ones<vec>(m_data_p->m_nobs);
  m_rowsums = ones<mat>(m_data_p->m_nobs);
  m_weightTMP = zeros<vec>(m_data_p->m_nobs);
  m_logalpha.resize(m_param_current_p->m_K);
  m_logbeta.resize(m_param_current_p->m_K);
  m_gamma.resize(m_param_current_p->m_K);
  m_xi.resize(m_param_current_p->m_K);
  for (int k=0; k<m_param_current_p->m_K; k++){
    m_logalpha[k].resize(m_data_p->m_nobs);
    m_logbeta[k].resize(m_data_p->m_nobs);
    m_gamma[k].resize(m_data_p->m_nobs);
    m_xi[k].resize(m_data_p->m_nobs);
    for (int i=0; i < m_data_p->m_nobs; i++){
      m_logalpha[k][i].resize(m_data_p->nbseq(i));
      m_logbeta[k][i].resize(m_data_p->nbseq(i));
      m_gamma[k][i].resize(m_data_p->nbseq(i));
      m_xi[k][i].resize(m_data_p->nbseq(i));
      for (int b=0; b < m_data_p->nbseq(i); b++){
        m_logalpha[k][i][b] = zeros(m_data_p->m_nbtimeobs[i](b), m_param_current_p->m_M);
        m_logbeta[k][i][b] = zeros(m_data_p->m_nbtimeobs[i](b), m_param_current_p->m_M);
        m_gamma[k][i][b] = zeros(m_data_p->m_nbtimeobs[i](b)-1, m_param_current_p->m_M);
        m_xi[k][i][b].zeros(m_param_current_p->m_M, m_param_current_p->m_M, m_data_p->m_nbtimeobs[i](b)-1);
      }
    }
  }
  m_tihConta.resize(m_param_current_p->m_M);
  for (int h=0; h<m_param_current_p->m_M;h++){
    m_tihConta[h].resize(m_data_p->m_nobs);
    for (int i=0; i<m_data_p->m_nobs; i++){
      m_tihConta[h][i].resize(m_data_p->nbseq(i));
      for (int b=0; b < m_data_p->nbseq(i); b++){
        m_tihConta[h][i][b] = zeros(m_data_p->m_nbtimeobs[i](b), 2);
      }
    }
  }
  m_eta.resize(m_data_p->m_nobs);
  for (int i=0; i<m_data_p->m_nobs; i++){
    m_eta[i].resize(m_data_p->nbseq(i));
    for (int b=0; b < m_data_p->nbseq(i); b++){
      m_eta[i][b] = zeros(m_data_p->m_nbtimeobs[i](b), m_param_current_p->m_M);
    }
  }

}

void EMalgo::SwitchParamCurrent(int ini){m_param_current_p = &m_paramCand[ini];}

void EMalgo::Run(){
  // Partie Small EM
  int nbSmall = m_paramCand.size();
  for (int ini=0; ini<nbSmall; ini++){
    SwitchParamCurrent(ini);
    OneEM();
    m_allloglike(ini) = ComputeLogLike();
    // Filtre degenerescence
    if (m_allloglike(ini)!=m_allloglike(ini)) m_allloglike(ini)=-999999999999;
  }
  // On conserve les meilleurs initialisations
  uvec indices = sort_index(m_allloglike);
  m_nbiter = -log(0);
  for (int tmp1=0; tmp1<m_nbKeep; tmp1++){
    SwitchParamCurrent(indices(nbSmall - tmp1 - 1));
    OneEM();
    m_allloglike(indices(nbSmall - tmp1 - 1)) = ComputeLogLike();
    // Filtre degenerescence
    if (m_allloglike(indices(nbSmall - tmp1 - 1)) != m_allloglike(indices(nbSmall - tmp1 - 1))){
      m_allloglike(indices(nbSmall - tmp1 - 1)) = -999999999999;
    }
  }
  uword  index;
  double indicebest = (m_allloglike).max(index);
  SwitchParamCurrent(index);
  OneEM();
  m_loglikeoutput = ComputeLogLike();
  Estep();
  indices = sort_index(m_allloglike);
}

colvec EMalgo::FindZMAP(){
  Col<double> zMAP=ones<vec>(m_data_p->m_nobs);
  uword  index;
  double max_val=0;
  for (int i=0; i<m_data_p->m_nobs; i++){
    max_val = (m_tau.row(i)).max(index);
    zMAP(i)=index;
  }
  return zMAP;
}

double EMalgo::ComputeLogLike(){
  forwardbackward();
  double output=-99999999999999;
  ComputeTmpLogProba();
  m_maxtmplogproba = max(m_tau, 1);
  if (min(m_maxtmplogproba) == 0){
    output = log(0);
  }else{
    for (int k=0; k<m_param_current_p->m_K; k++) m_tau.col(k)-=m_maxtmplogproba;
    m_tau = exp(m_tau);
    m_rowsums = sum(m_tau,1);
    output = sum(m_maxtmplogproba) + sum(log(m_rowsums));
  }
  return output;
}


void EMalgo::forwardbackward(){
  for (int k=0; k<m_param_current_p->m_K; k++)
    m_param_current_p->m_A[k] = log(m_param_current_p->m_A[k]);
  for (int i=0; i<m_data_p->m_nobs; i++){
    for (int b=0; b < m_data_p->nbseq(i); b++){
      m_logdspe_oneind = m_param_current_p->m_lambda.dlogspecificAll(m_data_p->m_yi[i][b], m_param_current_p->m_M);
      for (int k=0; k<m_param_current_p->m_K; k++){

        m_logalpha[k][i][b].row(0) = log(m_param_current_p->m_pi.row(k)) + m_logdspe_oneind.row(0);
        for (int t=0; t < (m_data_p->m_nbtimeobs[i](b)-1); t++)
          m_logalpha[k][i][b].row(t+1) = ComputeLogColMat( m_param_current_p->m_A[k].each_col() + trans(m_logalpha[k][i][b].row(t)) ) + m_logdspe_oneind.row(t+1);


        m_logbeta[k][i][b].row(m_data_p->m_nbtimeobs[i](b)-1) = trans(zeros(m_param_current_p->m_M));
        for (int t=(m_data_p->m_nbtimeobs[i](b)-2); t >=0 ; t--)
          m_logbeta[k][i][b].row(t) = trans(ComputeLogRowMat(m_param_current_p->m_A[k].each_row() + (m_logbeta[k][i][b].row(t+1)+ m_logdspe_oneind.row(t+1))));

        m_gamma[k][i][b] = NormLogMatRow(m_logalpha[k][i][b] + m_logbeta[k][i][b]);

        for (int t=0; t < (m_data_p->m_nbtimeobs[i](b)-1); t++)
          m_xi[k][i][b].slice(t) = NormLogMatFull((m_param_current_p->m_A[k].each_col() + trans(m_logalpha[k][i][b].row(t))).each_row() + ( m_logbeta[k][i][b].row(t+1) + m_logdspe_oneind.row(t+1)));

      }
    }
  }
  for (int k=0; k<m_param_current_p->m_K; k++)
    m_param_current_p->m_A[k] = exp(m_param_current_p->m_A[k]);
}


void EMalgo::Estep(){
  for (int k=0; k<m_param_current_p->m_K; k++)
    m_tau.col(k) = m_tau.col(k)/m_rowsums;

  for (int i=0; i<m_data_p->m_nobs; i++){
    for (int b=0; b < m_data_p->nbseq(i); b++){
      m_eta[i][b] = zeros(m_data_p->m_nbtimeobs[i](b), m_param_current_p->m_M);
      for (int k=0; k< m_param_current_p->m_K; k++)
        m_eta[i][b] += m_gamma[k][i][b] * m_tau(i,k);
      for (int h=0; h< m_param_current_p->m_M; h++)
        m_tihConta[h][i][b] = m_param_current_p->m_lambda.dspecificProbaCond(m_data_p->m_yi[i][b], h);
    }
  }
}

void EMalgo::OneEM(){
  double loglike = ComputeLogLike(), prec = log(0);
   int it=0;
  while ( (it<m_nbiter) && ((loglike-prec)>m_tol))  {
    it ++;
    Estep();
    Mstep();
    prec = loglike;
    loglike = ComputeLogLike();
   }
}


void EMalgo::ComputeTmpLogProba(){
  double tmp = 0;
  for (int k=0; k<m_param_current_p->m_K; k++){
    for (int i=0; i<m_data_p->m_nobs; i++){
      tmp = 0;
      for (int b=0; b < m_data_p->nbseq(i); b++){
        tmp += ComputeLogSumVec(trans(m_logalpha[k][i][b].row(m_data_p->m_nbtimeobs[i](b)-1)));
      }
      m_tau(i,k) =  log(m_param_current_p->m_delta(k)) + tmp;
    }
  }
}

void EMalgo::Mstep(){
  MstepProp();
  MstepA();
  MstepSpecific();
}


void EMalgo::MstepProp(){
  m_param_current_p->m_delta = trans(sum(m_tau,0));
  m_param_current_p->m_delta = m_param_current_p->m_delta / sum(m_param_current_p->m_delta);
}

void EMalgo::MstepA(){
  Mat<double> num = zeros(m_param_current_p->m_M, m_param_current_p->m_M);
  for (int k=0; k < m_param_current_p->m_K; k++){
    num.zeros(m_param_current_p->m_M, m_param_current_p->m_M);
    for (int i=0; i<m_data_p->m_nobs; i++){
      for (int b=0; b < m_data_p->nbseq(i); b++){
        for (int t=0; t < (m_data_p->m_nbtimeobs[i](b)-1); t++)
          num = num + m_xi[k][i][b].slice(t) * m_tau(i,k);
      }
    }
    m_param_current_p->m_A[k] = num.each_col()/sum(num, 1);
  }
}



void EMalgo::MstepSpecific(){
  /*if (m_param_current_p->m_lambda.m_model ==  "GaussConta"){
    for (int h=0; h< m_param_current_p->m_M; h++){
      double upper_eps=0;
      double lower_eps=0;
      double upper_mean=0;
      double upper_sd=0;
      double lower_mean=0;
      for (int i=0; i<m_data_p->m_nobs; i++){
        for (int b=0; b < m_data_p->nbseq(i); b++){
          upper_eps += sum(m_tihConta[h][i][b].col(1) % m_eta[i][b].col(h));
          lower_eps += sum(m_eta[i][b].col(h));
          upper_mean += sum(m_data_p->m_yi[i][b] % m_tihConta[h][i][b].col(0) % m_eta[i][b].col(h));
          lower_mean += sum(m_tihConta[h][i][b].col(0) % m_eta[i][b].col(h));
        }
      }
      m_param_current_p->m_lambda.m_eps(h) =  (upper_eps / lower_eps);
      m_param_current_p->m_lambda.m_a(h) = upper_mean / lower_mean;
      for (int i=0; i<m_data_p->m_nobs; i++){
        for (int b=0; b < m_data_p->nbseq(i); b++){
          upper_sd += sum(pow(m_data_p->m_yi[i][b] - m_param_current_p->m_lambda.m_a(h), 2) % m_tihConta[h][i][b].col(0) % m_eta[i][b].col(h));
        }
      }
      m_param_current_p->m_lambda.m_b(h) = sqrt(upper_sd / lower_mean);
    }
  }*/
    for (int h=0; h< m_param_current_p->m_M; h++){
      double upper_eps=0;
      double lower_eps=0;
      double upper_mean=0;
      double upper_2=0;
      double lower_mean=0;
      Col<double> tmpCol;
      for (int i=0; i<m_data_p->m_nobs; i++){
        for (int b=0; b < m_data_p->nbseq(i); b++){
          tmpCol =  m_tihConta[h][i][b].col(0) % m_eta[i][b].col(h);
          upper_eps += sum(m_tihConta[h][i][b].col(1) % m_eta[i][b].col(h));
          lower_eps += sum(m_eta[i][b].col(h));
          uvec keepObs = find(tmpCol != 0);
          lower_mean += sum(tmpCol.elem(keepObs));
          upper_mean += sum(m_data_p->m_yi[i][b].elem(keepObs) % tmpCol.elem(keepObs));
          upper_2 += sum( log(m_data_p->m_yi[i][b].elem(keepObs)) % tmpCol.elem(keepObs));
        }
      }
      double cst = (log(upper_mean / lower_mean) - upper_2 / lower_mean);
      m_param_current_p->m_lambda.m_eps(h) =  (upper_eps / lower_eps);
      double tmp_prec_alpha=0;
      m_param_current_p->m_lambda.m_a(h) =  0.1;
      int cp = 0;
      while ( (abs(tmp_prec_alpha - m_param_current_p->m_lambda.m_a(h)) > m_tol) && (cp<50)){
        tmp_prec_alpha = m_param_current_p->m_lambda.m_a(h);
        m_param_current_p->m_lambda.m_a(h) = m_param_current_p->m_lambda.m_a(h) - (log(m_param_current_p->m_lambda.m_a(h)) - digamma(m_param_current_p->m_lambda.m_a(h)) - cst) / (1/m_param_current_p->m_lambda.m_a(h) - polygamma(1,m_param_current_p->m_lambda.m_a(h)) );
        cp++;
      }
      m_param_current_p->m_lambda.m_b(h) = m_param_current_p->m_lambda.m_a(h) * lower_mean / upper_mean;
    }
  
}
