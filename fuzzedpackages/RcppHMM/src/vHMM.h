//=================================
// include guard
#pragma once

//=================================
// forward declared dependencies

//=================================
// included dependencies
#include <vector>
#include <string>
#include "RcppArmadillo.h"

#define EPSILON 1.0e-5
#define MINSD 0.399  // (1/sqrt(2*pi))

//=================================
// new structure definition

struct scaledMatrix 
{
    Rcpp::NumericVector scaling;
    Rcpp::NumericMatrix matrix;
};

//=================================
// the actual class
class vHMM  // Base class
{

    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:
    // Getters   
    // The first const (before the float type statement) means it will return a const value
    // The second const (after the parameters) means the function will not modify any member variables of the class    
    virtual Rcpp::CharacterVector getStateNames(void) const;    
    Rcpp::NumericMatrix getA(void) const;    
    Rcpp::NumericVector getPi(void) const;

    // Setters
    // Hidden states -> rows, Possible emissions -> columns
    void setB(Rcpp::NumericMatrix B);
    void setParameters(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B, Rcpp::NumericVector  Pi);
    void setStateNames(Rcpp::CharacterVector stateNames);
    void setA(Rcpp::NumericMatrix A);    
    void setPi(Rcpp::NumericVector  Pi);    

    //  Miscellaneous
    //  Returns all the model parameterrs as a list
    virtual Rcpp::List toList(void) const;
    //  All the data expressed as discrete {0,1,2,...} is returned to categorical
    Rcpp::CharacterVector toName( Rcpp::IntegerVector index, char vectorName);
    //  (C++ only) functions used to overload the << operator and print to console the model parameters
    //virtual std::ostream& print(std::ostream& out) const;
    //friend std::ostream& operator<< (std::ostream & out, const vHMM& data);
        
    //--------------------------------------------------------------------------
    // PROTECTED METHODS:
    //--------------------------------------------------------------------------

protected:
    //  Miscellaneous

    //  Initialize the model with random parameters
    void randomInit();
    //  Verify if it is a stochastic matrix or vector (row elements sum equals to 1)
    bool verifyMatrix(Rcpp::NumericMatrix matrix);
    bool verifyVector(Rcpp::NumericVector vector);
  
    //--------------------------------------------------------------------------
    // PROTECTED MEMBERS:
    //--------------------------------------------------------------------------

protected:
    // Hidden states
    unsigned short int m_N; 
    Rcpp::CharacterVector m_StateNames; 

    // Parameters    
    Rcpp::NumericMatrix m_A;  // Transition matrix
    Rcpp::NumericVector m_Pi; // Initial state probability vector

    //--------------------------------------------------------------------------
    // PRIVATE METHODS:
    //--------------------------------------------------------------------------

private:
    Rcpp::NumericMatrix getB(void) const;

    //--------------------------------------------------------------------------
    // PRIVATE MEMBERS:
    //--------------------------------------------------------------------------

private:
    Rcpp::NumericMatrix m_B;  // Emission matrix    
};
