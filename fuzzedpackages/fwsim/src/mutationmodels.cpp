#include "common.h"
#include "mutationmodels.h"

/*
##############################################################################
# MutationModel
##############################################################################
*/
MutationModel::MutationModel() {
}

MutationModel::MutationModel(const Rcpp::NumericMatrix _mutpars) {
  this->mutpars = _mutpars;
  this->loci = _mutpars.ncol();
}

std::vector<double> MutationModel::prob_mut_dw(const std::vector<int> profile) {
  assert((int)(profile.size()) == this->loci);
  
  std::vector<double> mut_dw(this->loci);
  
  for (int locus = 0; locus < this->loci; locus++) {
    mut_dw[locus] = this->prob_mut_dw(profile[locus], locus);
    //Rcpp::Rcout << "P(" << profile[locus] << " -> " << profile[locus] - 1 << ") = " << mut_dw[locus] << std::endl;
  }
  
  return mut_dw;
}

std::vector<double> MutationModel::prob_mut_up(const std::vector<int> profile) {
  assert((int)(profile.size()) == this->loci);

  std::vector<double> mut_up(this->loci);
  
  for (int locus = 0; locus < this->loci; locus++) {
    mut_up[locus] = this->prob_mut_up(profile[locus], locus);
    //Rcpp::Rcout << "P(" << profile[locus] << " -> " << profile[locus] + 1 << ") = " << mut_up[locus] << std::endl;
  }
  
  return mut_up;
}

std::vector<double> MutationModel::prob_not_mut(const std::vector<int> profile) {
  assert((int)(profile.size()) == this->loci);
  
  std::vector<double> mut_dw = prob_mut_dw(profile);
  std::vector<double> mut_up = prob_mut_up(profile);
  
  std::vector<double> not_mut(this->loci);
  
  for (int locus = 0; locus < this->loci; locus++) {
    not_mut[locus] = 1.0 - mut_dw[locus] - mut_up[locus];
    //Rcpp::Rcout << "P(" << profile[locus] << " -> " << profile[locus] << ") = " << not_mut[locus] << std::endl;
  }
  
  return not_mut;
}

/*
Rcpp::NumericMatrix MutationModel::mutation_table(const std::vector<int> profile) {
  assert((int)(profile.size()) == this->loci);
  
  std::vector<double> mut_dw = prob_mut_dw(profile);
  std::vector<double> mut_up = prob_mut_up(profile);
  
  Rcpp::NumericMatrix tab(3, this->loci);
  
  for (int locus = 0; locus < this->loci; locus++) {
    tab(0, locus) = mut_dw[locus];
    tab(1, locus) = mut_up[locus];
    tab(2, locus) = 1.0 - mut_dw[locus] - mut_up[locus];
  }
  
  return tab;
}
*/


void MutationModel::mutation_table(const int allele, const int locus, double* vec) { 
  vec[0] = this->prob_mut_dw(allele, locus);
  vec[1] = this->prob_mut_up(allele, locus);
  vec[2] = 1.0 - vec[0] - vec[1];
}

/*
##############################################################################
# Traditionel SMM
##############################################################################
# 1: mu_d
# 2: mu_u
#
# P(i -> i-1) = mu_d
# P(i -> i+1) = mu_u
# P(i -> i)   = 1 - mu_d - mu_u
##############################################################################
*/
SMM::SMM() {
}

SMM::SMM(const Rcpp::NumericMatrix _mutpars) : MutationModel(_mutpars) {
  //assert(this->mutpars.nrow() == 2);
  if (this->mutpars.nrow() != 2) {
    Rcpp::stop("The mutation parameter matrix must have 2 rows.");
  }
}

double SMM::prob_mut_dw(const int allele, const int locus) {
  return this->mutpars(0, locus);
}

double SMM::prob_mut_up(const int allele, const int locus) {
  return this->mutpars(1, locus);
}


/*
##############################################################################
# Jochens et al. 2011
##############################################################################
# 1: gamma_d
# 2: alpha_d
# 3: beta_d
# 4: gamma_u
# 5: alpha_u
# 6: beta_u
#
# P(i -> i-1) = gamma_d / (1 + exp(alpha_d*(beta_d - i)))
# P(i -> i+1) = gamma_u / (1 + exp(alpha_u*(beta_u - i)))
# P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)
##############################################################################
*/
LMM::LMM() {
}

LMM::LMM(const Rcpp::NumericMatrix _mutpars) : MutationModel(_mutpars) {
  //assert(this->mutpars.nrow() == 6);
  if (this->mutpars.nrow() != 6) {
    Rcpp::stop("The mutation parameter matrix must have 6 rows.");
  }
}

double LMM::prob_mut_dw(const int allele, const int locus) {
  return this->mutpars(0, locus) / (1.0 + exp(this->mutpars(1, locus)*(this->mutpars(2, locus) - (double)allele)));
}

double LMM::prob_mut_up(const int allele, const int locus) {
  return this->mutpars(3, locus) / (1.0 + exp(this->mutpars(4, locus)*(this->mutpars(5, locus) - (double)allele)));
}

/*
##############################################################################
# Svantes
##############################################################################
# 1: a
# 2: b
# 3: alpha
# 4: beta
#
# P(i -> i-1) = 1/((1 + exp(a + b*i))*(1 + exp(alpha + beta*i)))
# P(i -> i+1) = exp(alpha + beta*i)/((1 + exp(a + b*i))*(1 + exp(alpha + beta*i)))
# P(i -> i)   = exp(a + b*i)/(1 + exp(a + b*i))
##############################################################################
*/
EMM::EMM() {
}

EMM::EMM(const Rcpp::NumericMatrix _mutpars) : MutationModel(_mutpars) {
  //assert(this->mutpars.nrow() == 4);
  if (this->mutpars.nrow() != 4) {
    Rcpp::stop("The mutation parameter matrix must have 4 rows.");
  }
}

double EMM::prob_mut_dw(const int allele, const int locus) {
  double h = (double)allele;
  return 1.0 / ( (1.0 + exp(this->mutpars(0, locus) + this->mutpars(1, locus)*h)) * (1.0 + exp(this->mutpars(2, locus) + this->mutpars(3, locus)*h)) );
}

double EMM::prob_mut_up(const int allele, const int locus) {
  double h = (double)allele;
  return exp(this->mutpars(2, locus) + this->mutpars(3, locus)*h)/((1.0 + exp(this->mutpars(0, locus) + this->mutpars(1, locus)*h))*(1.0 + exp(this->mutpars(2, locus) + this->mutpars(3, locus)*h)));
}

