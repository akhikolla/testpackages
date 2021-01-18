//=================================
// include guard
#pragma once

//=================================
// included dependencies
//=================================
#include <vector>
#include <string>
#include "mvnorm.h"

#define EPSILON 1.0e-5
#define MINSD 0.399  // (1/sqrt(2*pi))

//=================================
// new structure definition

struct scaledMatrixM 
{
    arma::rowvec scaling;
    arma::mat matrix;
};

//=================================
// the actual class
class MultiGHMM
{
    public:
        //! Constructor of MultiGHMM.
        MultiGHMM(unsigned short int numberStates, unsigned short int numberObservations);    
        MultiGHMM(Rcpp::CharacterVector stateNames, arma::mat A, arma::mat mu, arma::cube sigma, arma::rowvec Pi);    

        //! Destructor of HMMcont.
        virtual ~MultiGHMM(void);
    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

    public:
        // Getters   
        // The first const (before the float type statement) means it will return a const value
        // The second const (after the parameters) means the function will not modify any member variables of the class    
        Rcpp::CharacterVector getStateNames(void) const;    
        arma::mat getA(void) const; 
        arma::mat getMu(void) const;
        arma::cube getSigma(void) const;   
        arma::rowvec getPi(void) const;

        // Setters
        // Hidden states -> rows        
        void setStateNames(Rcpp::CharacterVector stateNames);
        void setA(arma::mat A);        
        void setMu(arma::mat mu);  
        void setSigma(arma::cube sigma);      
        void setPi(arma::rowvec Pi);    
        void setParameters(arma::mat A, arma::mat  mu, arma::cube sigma, arma::rowvec Pi);
        
        // Evaluation methods 
        double evaluation(arma::mat sequences, char method = 'f');

        // Decoding methods
        Rcpp::CharacterVector viterbi(arma::mat sequences);
        Rcpp::CharacterVector forwardBackward(arma::mat sequences);

        // Learning methods 
        double loglikelihood(arma::cube sequences);
        void learnBW(arma::mat sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);
        void learnEM(arma::cube sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);
        //void learnEMParallel(arma::cube sequences, unsigned short int iter = 100, double delta = EPSILON, unsigned char pseudo = 0, bool print = true);

        // Simulation
        Rcpp::List generateObservations(unsigned short int length);

        //  Miscellaneous
        //  Returns all the model parameterrs as a list
        Rcpp::List toList(void) const;
        //  All the data expressed as discrete {0,1,2,...} is returned to categorical
        Rcpp::CharacterVector toName( Rcpp::IntegerVector index, char vectorName);
    
        //--------------------------------------------------------------------------
        // PROTECTED METHODS:
        //--------------------------------------------------------------------------

    protected:
        //  Initialize the model with random parameters
        void randomInit(double min, double max);   
        //  Verify if it is a stochastic matrix or vector (row elements sum equals to 1)
        bool verifyMatrix(arma::mat matrix);
        bool verifyVector(arma::rowvec vector);
    
        //--------------------------------------------------------------------------
        // PRIVATE METHODS:
        //--------------------------------------------------------------------------

    private:
        // Evaluation methods
        void  forwardMatrix(arma::mat sequences, unsigned int length , scaledMatrixM & forward);
        void  backwardMatrix(arma::mat sequences, unsigned int length , scaledMatrixM & backward);

        // Decoding methods
        arma::mat forwardBackwardGamma(arma::mat sequences);
        void forwardBackwardGamma(arma::mat index, scaledMatrixM & forward, scaledMatrixM & backward,  arma::rowvec & scaledf, arma::rowvec & scaledb, arma::mat & matrix, unsigned int length);

        // Learning methods  
        bool BaumWelch( arma::mat sequences, unsigned int pseudo);  
        bool expectationMaximization(arma::cube sequences, unsigned int pseudo);
        //void expectationMaximizationParallel( arma::cube sequences, unsigned int pseudo);
    
        //--------------------------------------------------------------------------
        // PROTECTED MEMBERS:
        //--------------------------------------------------------------------------

    protected:
        // Hidden states
        unsigned short int m_N; 
        Rcpp::CharacterVector m_StateNames; 

        //  Number of observations per each time
        unsigned short int m_M; 

        // Parameters    
        arma::mat m_A;        // Transition matrix
        arma::mat m_mu;       // Emission matrix         
        arma::cube m_sigma;    // Emission matrix 
        arma::rowvec m_Pi;       // Initial state probability vector              
};
