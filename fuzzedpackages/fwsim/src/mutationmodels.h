#ifndef _haptools_MUTATIONMODELS_H
#define _haptools_MUTATIONMODELS_H

#include "common.h"

class MutationModel {
  protected:
    Rcpp::NumericMatrix mutpars;
    int loci;
    //Rcpp::NumericMatrix
    
    virtual double prob_mut_dw(const int allele, const int locus) = 0;
    virtual double prob_mut_up(const int allele, const int locus) = 0;
    //double prob_not_mut(const int allele, const int locus);
    
    std::vector<double> prob_mut_dw(const std::vector<int> profile);
    std::vector<double> prob_mut_up(const std::vector<int> profile);
    std::vector<double> prob_not_mut(const std::vector<int> profile);
    
  public:
    MutationModel();
    MutationModel(const Rcpp::NumericMatrix _mutpars);    
    //Rcpp::NumericMatrix mutation_table(const std::vector<int> profile);
    void mutation_table(const int allele, const int locus, double* vec);
};

class SMM : public MutationModel {
  protected:  
    virtual double prob_mut_dw(const int allele, const int locus);
    virtual double prob_mut_up(const int allele, const int locus);

  public:
    SMM();
    SMM(const Rcpp::NumericMatrix _mutpars);
};

class LMM : public MutationModel {
  protected:  
    virtual double prob_mut_dw(const int allele, const int locus);
    virtual double prob_mut_up(const int allele, const int locus);

  public:
    LMM();
    LMM(const Rcpp::NumericMatrix _mutpars);
};

class EMM : public MutationModel {
  protected:  
    virtual double prob_mut_dw(const int allele, const int locus);
    virtual double prob_mut_up(const int allele, const int locus);

  public:
    EMM();
    EMM(const Rcpp::NumericMatrix _mutpars);
  
};

#endif

