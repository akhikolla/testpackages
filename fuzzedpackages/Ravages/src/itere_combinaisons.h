#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

#ifndef oz_comb
#define oz_comb
/*
  n = longueur des combinaisons
  n_cat = nombre de categories
  cur_disp = ce qui reste dispo dans chaque categorie compte tenu de la combinaison courante
  cur_comb = la combinaison courante
*/
class comb {
  public:
  int n, n_cat;
  std::vector<int> cur_disp, cur_comb;
  bool anything_left;

  // l'utilisateur passe n et disp = ce qui est disponible au total dans chaque categorie
  comb(int n_, std::vector<int> disp);
  std::vector<int> & current();
  std::vector<int> & disp();
  bool left() const;
  void itere();
  
  void print_current();
};
#endif
