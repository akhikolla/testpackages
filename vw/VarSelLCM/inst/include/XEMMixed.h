/*
Cette classe est héritière d'XEM et permet d'obtenir le MLE pour des données continues

Ces élements sont:
ceux de la classe Algorithm
paramCand : vecteur contenant les paramètres issus des différentes initialisations
paramCurrent_p : pointeur vers les paramètres à partir desquels TOUS les calculs sont faits
data_p : pointeur vers les données
m_weightTMP : ulilisé pour le calcul des probabilités conditionelles

*/
#ifndef XEMMixed_H
#define XEMMixed_H
#include "DataMixed.h"
#include "ParamMixed.h"
#include "XEM.h"
#include "XEMContinuous.h"


class XEMMixed : public XEM{
  public:
    vector<ParamMixed> paramCand;
  ParamMixed * paramCurrent_p;
  const DataMixed * data_p;
  Col<double> m_weightTMP;
  uvec locationContinuous, locationInteger, locationCategorical;
  Col<double> omegaContinuous, omegaInteger, omegaCategorical;
  
  // Constructeurs et destructeurs par défaut (non utilisé)
  XEMMixed(){};
  ~XEMMixed(){};
  
  // Constructeurs avec et sans les paramètres de réglages
  XEMMixed(const DataMixed *, const colvec &, const int &);
  XEMMixed(const DataMixed *,  const S4 *);
  // Initialisation des paramètres spécifiques
  void InitSpecificParamXEMMixed(const DataMixed * datapasse);
  
  // Les trois fonction suivantes sont à redéfinir pour chaque classe héritiaire d'XEM
  // Etape M
  virtual void Mstep();
  // calcul des log-proba conditionnelles
  virtual void ComputeTmpLogProba();
  // change le pointeur des paramètres acutels (voir définietion spécifique des classes)
  virtual void SwitchParamCurrent(int);
  // Acutalise l'object S4 retourné sous R
  virtual void Output(S4 *);
};
#endif
