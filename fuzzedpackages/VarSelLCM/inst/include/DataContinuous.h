/*
Cette classe contient les éléments relatifs aux données continues

Ces élements sont:
m_x : matrice des données
m_priors : matrice des hyper-paramètres (1 ligne par variable)
m_notNA : matrice binaire indiquant si l'observation est manquante ou non (attention dans m_x les valeurs
manquantes ont été artificellement remplacées par la valeur 0 pour pouvoir utiliser armadillo)
*/
#ifndef DataContinuous_H
#define DataContinuous_H

#include "Data.h"

class DataContinuous: public Data{
  public:
    Mat<double> m_x, m_priors, m_notNA;
  
  DataContinuous(){};
  DataContinuous(const S4 &);
  ~DataContinuous(){};
  
};
#endif
