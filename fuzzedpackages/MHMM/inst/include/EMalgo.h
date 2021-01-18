#ifndef EMalgo_H
#define EMalgo_H

#include "Data.h"
#include "Param.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

class EMalgo{
public:
  double m_tol, m_loglikeoutput,m_cst;
  int m_nbiter, m_nbKeep;
  const Data * m_data_p;
  Param * m_param_current_p;
  vector < Param > m_paramCand;
  vector < vector < vector < Mat <double> > > > m_logalpha, m_logbeta, m_gamma, m_tihConta;
  vector < vector < Mat <double> > > m_eta;
  vector < vector < vector < Cube <double> > > > m_xi;
  Mat<double> m_tau, m_dspe_oneind, m_logdspe_oneind;
  Col<double> m_maxtmplogproba, m_rowsums, m_weightTMP, m_allloglike;


  EMalgo(){};
  ~EMalgo(){};
  EMalgo(const Data *, const List, const double, const int, const int);

  void SwitchParamCurrent(int);
  void Run();
  void Estep();
  void MstepProp();
  void MstepA();
  void MstepSpecific();
  void Mstep();
  void OneEM();
  void forwardbackward();
  colvec FindZMAP();
  double ComputeLogLike();
  void ComputeTmpLogProba();
  void Output(S4 * reference_p){
    reference_p->slot("K") = wrap(m_param_current_p->m_K);
    reference_p->slot("M") = wrap(m_param_current_p->m_M);
    reference_p->slot("A") = wrap(m_param_current_p->m_A);
    reference_p->slot("delta") = wrap(m_param_current_p->m_delta);
    reference_p->slot("pi") = wrap(m_param_current_p->m_pi);
    as<List>(reference_p->slot("lambda"))[0] = wrap(m_param_current_p->m_lambda.m_eps);
    as<List>(reference_p->slot("lambda"))[1] = wrap(m_param_current_p->m_lambda.m_a);
    as<List>(reference_p->slot("lambda"))[2] = wrap(m_param_current_p->m_lambda.m_b);
  }
  double FunctionGammaToOptimize(double a);
};
#endif
