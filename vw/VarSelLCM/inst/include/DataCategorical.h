/*
Cette classe contient les éléments relatifs aux données qualitatives

Ces élements sont:
m_profiles : matrice constituée des profils uniques
m_nmodalities: vecteurs indiquant le nombre de modalités pour chaque variable
m_w : poids des profils
m_nprofiles : nombre de profiles
m_whotakewhat : vector de vector où m_whotakewhate[j][h] liste les indices des profiles prenant
la modalité h pour la variable j. Ne sont pas présents dans m_whotakewhate[j], les individus 
ayant une valeur manquante pour la variable j
*/
#ifndef DataCategorical_H
#define DataCategorical_H

#include "Data.h"

class DataCategorical : public Data{
  public:
    Mat<double> m_profiles;
    rowvec m_nmodalities;
    colvec m_w;
    Col<double> m_dl;
    int m_nprofiles;
    vector < vector < uvec > > m_whotakewhat;
  
  DataCategorical(){};
  DataCategorical(const S4 &);
  ~DataCategorical(){};
};
#endif
