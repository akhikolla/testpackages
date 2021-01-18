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
class HMM : public vHMM  // Derived class
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:
    //! Constructor of HMM.
    HMM(unsigned short int  numberStates, unsigned short int  numberEmissions);
    HMM(Rcpp::CharacterVector stateNames, Rcpp::CharacterVector emissionNames);
    HMM(Rcpp::CharacterVector stateNames, Rcpp::CharacterVector emissionNames, Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::NumericVector Pi);

    //! Destructor of HMM.
    virtual ~HMM(void);

    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:
    // Getters       
    Rcpp::CharacterVector getEmissionNames(void) const;
    Rcpp::NumericMatrix getB(void) const;

    // Setters
    // Hidden states -> rows, Possible emissions -> columns
    void setEmissionNames(Rcpp::CharacterVector emissionNames);
    void setB(Rcpp::NumericMatrix B);
    void setParameters(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::NumericVector  Pi);

    // Evaluation methods 
    double evaluation(Rcpp::CharacterVector sequence, char method = 'f');

    // Decoding methods
    Rcpp::CharacterVector viterbi(Rcpp::CharacterVector sequence);
    Rcpp::CharacterVector forwardBackward(Rcpp::CharacterVector sequence);

    // Learning methods 
    double loglikelihood(Rcpp::CharacterMatrix sequence);    
    void learnBW(Rcpp::CharacterVector sequences, unsigned short int iter, double delta, unsigned char pseudo, bool print );
    void learnEM(Rcpp::CharacterMatrix sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);
    //void learnEMParallel(Rcpp::CharacterMatrix sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);

    // Simulation
    Rcpp::List generateObservations( unsigned short int length);

    //  Miscellaneous
    Rcpp::List toList(void) const;
    //std::ostream& print(std::ostream& out) const;      

    //--------------------------------------------------------------------------
    // PROTECTED METHODS:
    //--------------------------------------------------------------------------

protected:
    //  Miscellaneous
    void randomInit();    
    Rcpp::CharacterVector toName( Rcpp::IntegerVector index, char vectorName);
    Rcpp::IntegerVector toIndex( Rcpp::CharacterVector observations);  
        
    //--------------------------------------------------------------------------
    // PRIVATE METHODS:
    //--------------------------------------------------------------------------

private:
    // Evaluation methods
    void  forwardMatrix(Rcpp::IntegerVector sequence, unsigned int length , scaledMatrix & forward);
    void  backwardMatrix( Rcpp::IntegerVector sequence, unsigned int length , scaledMatrix & backward);

    // Decoding methods
    Rcpp::NumericMatrix forwardBackwardGamma(Rcpp::CharacterVector sequence);
    void forwardBackwardGamma(Rcpp::IntegerVector index, scaledMatrix & forward, scaledMatrix & backward,  Rcpp::NumericVector & scaledf, Rcpp::NumericVector & scaledb, Rcpp::NumericMatrix & matrix, unsigned int length);

    // Learning methods 
    
    void BaumWelch( Rcpp::CharacterVector sequence, unsigned int pseudo);
    void expectationMaximization(Rcpp::CharacterMatrix sequences, unsigned int pseudo);
    //void expectationMaximizationParallel( Rcpp::CharacterMatrix sequences, unsigned int pseudo);    
  
    //--------------------------------------------------------------------------
    // PRIVATE MEMBERS:
    //--------------------------------------------------------------------------

private:  
    // Possible observations
    unsigned short int m_M; 
    Rcpp::CharacterVector m_ObservationNames;  

    Rcpp::NumericMatrix m_B;  // Emission matrix

};
