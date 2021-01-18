#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "statistics_class.h"
#include "allelecounter.h"
#include "weighted_left_product.h"

using namespace Rcpp;
using namespace RcppParallel;

class SKAT : public Stats {
  public:
  NumericVector full_p;      //  le vecteur de fréquences allélique (freq A2) [pour l'imputation données manquantes dans WLP]
  NumericMatrix Pi;      // matrice des probabilités d'appartenir à chaque groupe
  NumericMatrix Y_Pi;    // Y_Pi = Y - Pi
  NumericVector M1;      // les quatre premiers moments (non centrés), pour chaque SNPgroup
  NumericVector M2;
  NumericVector M3;
  NumericVector M4;
  std::vector<double> p;
  std::vector<double> full_W; // poids des SNPs (full data)
  std::vector<double> W; // poids des SNPs (données courantes)
  std::vector<int> nb_ind_in_group; // nb individu dans chaque groupe

  int iterates; // compteur (pour la mise à jour des moments)


// std::vector<double> debug;


  // pA = la matrice de genotypes
  // which_snps = les snp à considérer
  // SNP_group = facteur donnant le groupe pour tous les SNPs
  // ind_group = facteur donnant le groupe des individus
  // Pi = matrice donnant la proba a priori de chaque individu d'appartenir à chaque groupe	
  //      dimensions ncol (= nb individus)  x  nb_ind_groups 
  // W = vecteur des poids des SNPs
  SKAT(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup, IntegerVector ind_group, 
       NumericVector p_, NumericMatrix Pi_, NumericVector W_)
  : Stats(pA, which_snps, SNPgroup, ind_group), full_p(p_), Pi(Pi_), Y_Pi(ncol, nb_ind_groups), 
    M1(nb_snp_groups), M2(nb_snp_groups), M3(nb_snp_groups), M4(nb_snp_groups), 
    full_W(as<std::vector<double>>(W_)), nb_ind_in_group(nb_ind_groups), iterates(0) { 
    if(Pi.nrow() != ncol || Pi.ncol() != nb_ind_groups)
      stop("Pi dimensions mismatch");
    // pour compter combien d'individus dans chaque groupe
    for(int i : ind_group) nb_ind_in_group[i-1]++;
    // le constructeur Stats() n'a pas appelé la fonction redéfinie ici
    extra_update_snps();
  }


 // mise à jour du vecteur de poids et du vecteur p de fréquences alléliques
 void extra_update_snps() {
    p.resize(nb_snps);
    W.resize(nb_snps);
    size_t k = 0;
    for(size_t i = 0; i < full_nb_snps; i++) {
      if(which_snps[i]) {
        p[k] = full_p[i];
        W[k++] = full_W[i];
      }
    }
  }


  void compute_stats() {
    //Rcout << "nb_snps = " << nb_snps << "\n";
    //Rcout << "nb_snp_groups = " << nb_snp_groups << "\n";
    if(nb_snps == 0 || nb_snp_groups == 0) {
      return;
    }

    // Fabriquer la matrice des Y - Pi
    for(int j = 0; j < nb_ind_groups; j++) { // classes d'individus
      for(int i = 0; i < ncol; i++) { // individus
        if(ind_group[i] == j+1) 
          Y_Pi(i,j) = 1 - Pi(i,j);
        else
          Y_Pi(i,j) = -Pi(i,j);
      }
    }

    // Produit (Y - Pi) * G * W [ dimensions nb_snps x nb_ind_groups : elle est transposée ]
    NumericMatrix Z = WLP(&data[0], &p[0], nb_snps, ncol, true_ncol, W, Y_Pi);
    
    // Vecteur de stats
    for(int i = 0; i < nb_snp_groups; i++) stats[i] = 0;

    for(int j = 0; j < nb_ind_groups; j++) {
      for(int i = 0; i < nb_snps; i++) {
        stats[ snp_group[i] - 1 ] += Z(i,j)*Z(i,j) / nb_ind_in_group[j];
      }
    }

    // Mise à jour moments;
    if(iterates > 0) { // iterates == 0 correspond à la stat observée
      for(int i = 0; i < nb_snp_groups; i++) {
        if(nb_snp_in_group[i] == 0) 
          continue;
        double s = stats[i];
        double s2 = s * s;
        double s3 = s2 * s;
        double s4 = s3 * s;
        M1[i] += (s  - M1[i])/iterates;
        M2[i] += (s2 - M2[i])/iterates;
        M3[i] += (s3 - M3[i])/iterates;
        M4[i] += (s4 - M4[i])/iterates;
      }
    }
    iterates++;
  }

};

//[[Rcpp::export]]
List skat(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector region, IntegerVector group, 
          NumericVector p, NumericMatrix Pi, NumericVector weights, int A_target, int B_max) {

  SKAT B(p_A, which_snps, region, group, p, Pi, weights);
  if(B_max > 0) {
    List L = B.permute_stats(A_target, B_max);
    L["M1"] = B.M1;
    L["M2"] = B.M2;
    L["M3"] = B.M3;
    L["M4"] = B.M4;
    return L;
  } else {
    B.compute_stats();
    List L;
    L["statistic"] = B.stats;
    return L;
  }
}


RcppExport SEXP skat(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP pSEXP, SEXP PiSEXP, SEXP weightsSEXP, SEXP A_targetSEXP, SEXP B_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type A_target(A_targetSEXP);
    Rcpp::traits::input_parameter< int >::type B_max(B_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(skat(p_A, which_snps, region, group, p, Pi, weights, A_target, B_max));
    return rcpp_result_gen;
END_RCPP
}

