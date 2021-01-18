#ifndef XEMPEN_H
#define XEMPEN_H

#include "DataMixed.h"
#include "ParamMixed.h"

class XEMPen{
  public:
  DataMixed *data_p;
  S4 * strategy_p;
  Col<double> loglikepen, rowsums, maxtmplogproba, m_weightTMP, m_loglikenondis, munondisc, sdnondisc, lambdanondisc; 
  int nbSmall, iterSmall, nbKeep, iterKeep, iterCurrent, g, m_nbdegenere, degeneracy;
  double tolKeep, m_penalty;
  vector< Col<double> > omegaCand,  nbparamCand, alphanondisc;
  vector <ParamMixed> paramCand;
  Mat<double> tmplogproba ;
  Col<double> * omegaCurrent_p;
  ParamMixed * paramCurrent_p;
  Col<double> *  nbparamCurrent_p;

  XEMPen(){};
  XEMPen(const S4 *, const double);
  ~XEMPen(){};
  void Run();
  void Estep(){for (int k=0; k<g; k++) tmplogproba.col(k) = tmplogproba.col(k)/rowsums;};
  void OneEM();
  colvec FindZMAP();
  double ComputeLoglikepen();
  void Mstep();
  void ComputeTmpLogProba();
  void SwitchCurrent(int);
  void Output(S4 *);
};
#endif
