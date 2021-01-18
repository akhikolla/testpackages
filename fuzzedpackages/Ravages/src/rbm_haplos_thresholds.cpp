#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "polytomic_haplotype_proba.h"

using namespace Rcpp;

// haplos = matrice d'haplotypes (dimensions u x p) (u haplotypes, p SNPs)
// burden = liste (longueur n) de vecteurs de burdens (length u)
// sd = l'écart type de la composante E
// thr1, thr2 = thresholds définissant les différents groupes (longueur K)
// size un vecteur de tailles (longueur K)
// reps = nombre de réplicats 
// [[Rcpp::export]]
XPtr<matrix4> rbm_haplos_thresholds(IntegerMatrix haplos, List burden, NumericVector sd,
                     NumericVector thr1, NumericVector thr2, NumericVector size, int reps) {
  int u = haplos.nrow(); // nb d'haplotypes
  int p = haplos.ncol(); // nb SNPs
  int K = thr1.size(); // nb groupes de pop
  if(burden.size() != K || thr2.size() != K || size.length() != K || sd.length() != K)
    stop("Dimensions mismatch");

  int N = sum(size);
  XPtr<matrix4> pA(new matrix4(reps*p, N));
 
  size_t n_diplos = (u * (u+1)) / 2;
  // pour ne calculer qu'une fois les fréquences diplotypiques dans chaque groupe,
  // première boucle sur les groupes
  int ind = 0; // premier individu du groupe en cours
  for(int k = 0; k < K; k++) {
    haplo_probs Probs(burden[k], sd[k], thr1[k], thr2[k]);
    NumericVector P(n_diplos); // toutes les freqs diplotypiques (à constante près)
    size_t l = 0;
    for(int i = 0; i < u; i++) {
      for(int j = 0; j <= i; j++) {
         P[l++] = Probs(i,j);
      }
    }
    // on tire d'un coup reps*size[k] diplotypes
    // les arguments booleens sont replace = true, one_based = false
    // La fonction sample commence par "normaliser" P
    // (cf code dans Rcpp/sugar/functions/sample.h)
    IntegerVector DIP = sample(n_diplos, reps*size[k], true, P, false);
    int d = 0; // numéro de diplotype
    for(int rep = 0; rep < reps; rep++) { // replicats
      for(int i = 0; i < size[k]; i++) { // individus
        int h1, h2;   // numero d'haplotypes
        one_to_pair(DIP[d++], h1, h2);
        for(int snp = 0; snp < p; snp++) {
          pA->set(p*rep + snp, ind + i, haplos(h1,snp) + haplos(h2,snp) );
        }
      }
    }
    // premier individu du prochain groupe
    ind += size[k];
  }
  return pA;
}


RcppExport SEXP rbm_haplos_thresholds(SEXP haplosSEXP, SEXP burdenSEXP, SEXP sdSEXP, SEXP thr1SEXP, SEXP thr2SEXP, SEXP sizeSEXP, SEXP repsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type haplos(haplosSEXP);
    Rcpp::traits::input_parameter< List >::type burden(burdenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thr1(thr1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thr2(thr2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type reps(repsSEXP);
    rcpp_result_gen = Rcpp::wrap(rbm_haplos_thresholds(haplos, burden, sd, thr1, thr2, size, reps));
    return rcpp_result_gen;
END_RCPP
}

