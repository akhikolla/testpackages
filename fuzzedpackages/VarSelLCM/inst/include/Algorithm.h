/*
Cette classe virtuelle permet d'effectuer le choix de modèle au sens du critère MICL

Ces élements sont:
m_zStarBest : partition associée au modèle actuel qui a obtenue la meilleure vraisemblance complétée intégrée
m_zStarCurrent :  partition actuelle du modèle actuel
m_zCandCurrent :  partition candidate associée au modèle actuelle (i.e. égale à m_zStarCurrent sauf pour un individu).
                  C'est sur cette partition que les vraisemblances complétées intégrées sont toujours calculée.
m_integralenondiscrim : vecteur correspondant à la vraisemblance intégrée de chaque variable (non discriminante).
m_omegainit : matrice des initialisations d'omega (un modèle correspond à une colonne).
m_miclCurrent : valeur actuelle de la vraisemblance complétée intégrée.
m_miclBest : meilleure valeur de la vraisemblance complétée intégrée.
m_g : nombre de classes
vbleSelec : booléen indiquant si on effectue la sélection de variables au sens de MICL
m_omegaCurrent : modèle actuel
m_omegaBest : meilleur modèle au sens de MICL

*/
#ifndef Algorithm_H
#define Algorithm_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Algorithm{ 
  public:
  colvec m_zStarBest, m_zStarCurrent, m_zCandCurrent, m_integralenondiscrim;
  Mat<double> omegainit;
  double m_miclCurrent, m_miclBest;
  int m_g;
  bool vbleSelec;
  Col<double> m_omegaCurrent, m_omegaBest;
  
  // Constructeur et destructeur par défaut
  Algorithm(){};
  ~Algorithm(){};
  
  // Les trois fonctions suivantes sont communes pour les classes héritiaires d'Algorithm.  
  // Permet d'inialiser les éléments de la classes
  void InitCommumParamAlgo(const int &, const int &, const int &, const int &) ;
  // Effectue la selection de variables (si vbleSelec=TRUE) et actualise les sorties
  void Run(S4 *);
  // Calcul MICL pour un modèle fixe
  void ComputeMICL(S4 *);
  // Optimisation de la partition m_zStarCurrent pour le modèle m_omegaCurrent
  void Optimize_partition();
  
  // Les trois fonction suivantes sont à redéfinir pour chaque classe héritiaire d'Algorithm
  // Calcul la vraisemblance complétée intégrée pour le modèle m_omegaCurrent et la partition m_zCandCurrent
  virtual double Integre_Complete_Like_Cand() = 0;
  // Définit m_omegaCurrent comme le modèle maximisant la vraisemblance complétée intégrée pour la partition m_zCandCurrent
  virtual void Optimize_model()  = 0;
  // Définit la partition m_zCandCurrent initiale associée au modèle initial
  virtual void zCandInit()  = 0;
};
#endif  
