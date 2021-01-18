/*
Cette classes définie les paramètres pour des données continues

Ces éléments sont:
m_mu : matrice des moyennes
m_sd : matrice des ecarts-types
m_pi : proportions

*/
#ifndef ParamInteger_H
#define ParamInteger_H

#include "DataInteger.h"
#include "Param.h"

class ParamInteger : public Param{
  public:
  Mat<double> m_lambda;
    
  ParamInteger();
  ParamInteger(const ParamInteger & param);
  ParamInteger(const DataInteger *, const colvec & , const int &);
  ParamInteger(const DataInteger *, const colvec & , const int &, ivec);
  ~ParamInteger(){};
    void egalise(const DataInteger * data,const colvec );

};
#endif
