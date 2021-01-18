/*
Cette classe virtuelle permet d'effectuer l'estimation par maximum de vraisemblance

Ces élements sont:
nbSmall: nombre de Small EM effectués
iterSmall:  nombre d'itération de chaque Small EM
nbKeep : nombre de chaines que l'on conserve
iterKeep: nombre maximum d'itération pour les chaines conservées
iterCurrent : nombre maximum d'itération pour l'algorithme en cours (Small EM ou EM)
g : nombre de composantes
tolKeep: arret si la différence entre deux log-vraisemblances successives est inférieure que tolKeep
loglikeoutput: valeur de la lo-vraisemblance associée au MLE
loglikeSmall: og-vraisemblance pour chaque chaine
omega: role des variables (1: discriminant, 0: non discriminant)
rowsums: pour calculer les vraisemblances
maxtmplogproba: maximum des log-proba conditionelles pour chaque individu
tmplogproba: matrices des log-proba conditionnelles ou des TIK en fonction de l'étape du EM
location: indices des variables discriminantes
paramEstim: booléen indiquant si on fait l'estimation
*/
#ifndef XEM_H
#define XEM_H

#include "DataContinuous.h"
#include "DataCategorical.h"

class XEM{
  public:
  int nbSmall, iterSmall, nbKeep, iterKeep, iterCurrent, g, m_nbdegenere, degeneracy;
  double tolKeep, loglikeoutput;
  Col<double> loglikeSmall, omega, rowsums, maxtmplogproba;
  Mat<double> tmplogproba;
  uvec location;
  bool paramEstim;
  
  // Constructeur et destructeur par défaut
  XEM(){};
  ~XEM(){};
  // C'est deux fontions initialisent les paramètres communs aux EMs
  void InitCommumParamXEM(const colvec &, const int &);
  void InitCommumParamXEM(const colvec &, const int &, const S4 &);
  
  // Les quatres fonctions suivantes sont communes pour les classes héritiaires d'XEM.  
  // Lance l'estimation
  void Run();
  // Etape E
  void Estep();
  // Effectue un EM
  void OneEM();
  // Renvoie la partition MAP
  colvec FindZMAP();
  // calcul et renvoie de la logvraisemblance (a modifier pour les données avec des poids comme les qualitatives)
  virtual double ComputeLogLike();
  
  // Les trois fonction suivantes sont à redéfinir pour chaque classe héritiaire d'XEM
  // Etape M
  virtual void Mstep() = 0;
  // calcul des log-proba conditionnelles
  virtual void ComputeTmpLogProba() = 0;
  // change le pointeur des paramètres acutels (voir définietion spécifique des classes)
  virtual void SwitchParamCurrent(int) = 0;
  // Acutalise l'object S4 retourné sous R
  virtual void Output(S4 *) = 0;
  
};
#endif
