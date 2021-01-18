//=================================
// include guard
#pragma once

//=================================
// forward declared dependencies

//=================================
// included dependencies
#include "vHMM.h"

//=================================
// the actual class
class HMMpoisson : public vHMM  // Parent object
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of HMMpoisson.    
    HMMpoisson(unsigned short int  numberStates);
    HMMpoisson(Rcpp::CharacterVector stateNames);
    HMMpoisson(Rcpp::CharacterVector stateNames, Rcpp::NumericMatrix A, Rcpp::NumericVector B, Rcpp::NumericVector Pi);

    //! Destructor of cHHMMpoissonMM.
    virtual ~HMMpoisson(void);

    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:    
    //Getters
    Rcpp::NumericVector getB(void) const;
    // Setters
    // hidden states -> rows    
    void setB(Rcpp::NumericVector B);
    void setParameters(Rcpp::NumericMatrix A, Rcpp::NumericVector B, Rcpp::NumericVector  Pi);

    // Evaluation methods 
    double evaluation(Rcpp::IntegerVector sequence, char method = 'f');

    // Decoding methods
    Rcpp::CharacterVector viterbi(Rcpp::IntegerVector sequence);
    Rcpp::CharacterVector forwardBackward(Rcpp::IntegerVector sequence);

    // Learning methods 
    double loglikelihood(Rcpp::IntegerMatrix sequence);
    void learnBW(Rcpp::IntegerVector sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);
    void learnEM(Rcpp::IntegerMatrix sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);
    //void learnEMParallel(Rcpp::IntegerMatrix sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);

    // Simulation
    Rcpp::List generateObservations( unsigned short int length);

    // Miscellaneous
    Rcpp::List toList(void) const;
    //std::ostream& print(std::ostream& out) const;        
    //--------------------------------------------------------------------------
    // PROTECTED METHODS:
    //--------------------------------------------------------------------------

protected:
    // Miscellaneous
    void randomInit(double min, double max);   

    //--------------------------------------------------------------------------
    // PRIVATE METHODS:
    //--------------------------------------------------------------------------

private:
    // Evaluation methods
    void  forwardMatrix(Rcpp::IntegerVector sequence, unsigned int length , scaledMatrix & forward);
    void  backwardMatrix(Rcpp::IntegerVector sequence, unsigned int length , scaledMatrix & backward);

    // Decoding methods
    Rcpp::NumericMatrix forwardBackwardGamma(Rcpp::IntegerVector sequence);
    void forwardBackwardGamma(Rcpp::IntegerVector index, scaledMatrix & forward, scaledMatrix & backward,  Rcpp::NumericVector & scaledf, Rcpp::NumericVector & scaledb, Rcpp::NumericMatrix & matrix, unsigned int length);

    // Learning methods  
    void BaumWelch( Rcpp::IntegerVector sequences, unsigned int pseudo);  
    void expectationMaximization(Rcpp::IntegerMatrix sequences, unsigned int pseudo);
    //void expectationMaximizationParallel( Rcpp::IntegerMatrix sequences, unsigned int pseudo);

    //--------------------------------------------------------------------------
    // PRIVATE MEMBERS:
    //--------------------------------------------------------------------------

private:
    Rcpp::NumericVector m_B;  // Emission matrix         

};
