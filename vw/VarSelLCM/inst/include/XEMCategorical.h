/*
Cette classe est héritière d'XEM et permet d'obtenir le MLE pour des données qualitatives

Ces élements sont:
ceux de la classe Algorithm
paramCand : vecteur contenant les paramètres issus des différentes initialisations
paramCurrent_p : pointeur vers les paramètres à partir desquels TOUS les calculs sont faits
data_p : pointeur vers les données
tmpval : ulilisé pour le calcul des probabilités conditionelles

*/
#ifndef XEMCategorical_H
#define XEMCategorical_H


#include "DataCategorical.h"
#include "ParamCategorical.h"
#include "XEM.h"

class XEMCategorical : public XEM{
  public:
  vector<ParamCategorical> paramCand;
  ParamCategorical * paramCurrent_p;
  const DataCategorical * data_p;
  Col<double>tmpval;
  
  // Constructeurs et destructeurs par défaut (non utilisé)
  XEMCategorical(){};
  ~XEMCategorical(){};

  // Constructeurs avec et sans les paramètres de réglages
  XEMCategorical(const DataCategorical *, const colvec &, const int &);
  XEMCategorical(const DataCategorical *,  const S4 *);
  // Initialisation des paramètres spécifiques
  void InitSpecificParamXEMCategorical(const DataCategorical * datapasse);
  
  // Les trois fonction suivantes sont à redéfinir pour chaque classe héritiaire d'XEM
  // Etape M
  virtual void Mstep();
  // calcul des log-proba conditionnelles
  virtual void ComputeTmpLogProba();
  // calcul et renvoie de la logvraisemblance
  virtual double ComputeLogLike();
  // change le pointeur des paramètres acutels (voir définietion spécifique des classes)
  virtual void SwitchParamCurrent(int);
  // Acutalise l'object S4 retourné sous R
  virtual void Output(S4 *);
};
#endif
