#include <Rcpp.h>
#include <iostream>
#include "itere_combinaisons.h"
using namespace Rcpp;

/*
  n = longueur des combinaisons
  n_cat = nombre de categories
  cur_disp = ce qui reste dispo dans chaque categorie compte tenu de la combinaison courante
  cur_comb = la combinaison courante
*/

// l'utilisateur passe n et disp = ce qui est disponible au total dans chaque categorie
comb::comb(int n_, std::vector<int> disp) : n(n_), n_cat(disp.size()), cur_disp(disp) {
  cur_comb.resize(n);
  int j = 0;
  for(int i = 0; i < n; i++) {
    // on cherche une catégorie non vide
    while( cur_disp[j] == 0 && j < n_cat) j++;
    if(j == n_cat) {
      // on ne peut pas commencer
      anything_left = false;
      return;
    }
    cur_disp[j]--;
    cur_comb[i] = j;
  }
  anything_left = true;
}

std::vector<int> & comb::current() {
  return cur_comb;
}

std::vector<int> & comb::disp() {
  return cur_disp;
}
 
bool comb::left() const {
  return anything_left;
}

void comb::itere() {
  if(!anything_left) return;
  // en partant de la fin
  // on cherche une catégorie itérable
  int i = n-1;
  for(; i > -1; i--) {
    int c = cur_comb[i];
    cur_disp[c]++;  // on libere ça.
    // y a-t-il des dispos de c+1 à n_cat-1 ?
    bool found = false;
    for(int j = c+1; j < n_cat; j++) {
      if(cur_disp[j] > 0) { // ok
        cur_comb[i] = j; 
        cur_disp[j]--;
        found = true;
        break;
      }
    }
    if(found) { // il faut peupler de i à la fin
      int j = 0;
      for(int k = i+1; k < n; k++) {
        while( cur_disp[j] == 0 && j < n_cat) j++;
        if(j == n_cat) stop("y a un bug"); // ceci devrait donc ne jamais arriver
        cur_comb[k] = j;
        cur_disp[j]--;
      } 
      return;
    }
  }
  anything_left = false;
}

void comb::print_current() {
  for(int i = 0; i < n; i++)
    Rcout << cur_comb[i] << " ";
  Rcout << "\n";
  /* 
  Rcout << "[";
  for(int i = 0; i < n_cat; i++)
    Rcout << cur_disp[i] << " ";
  Rcout << "]\n";
  */
}
