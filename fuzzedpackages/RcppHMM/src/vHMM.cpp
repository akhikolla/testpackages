// included dependencies
#include "vHMM.h"

//  System libraries
#include <iostream>
#include <stdexcept>    /* invalid_argument */
#include <random>

//  namespaces
using namespace std;
using namespace Rcpp;

//--------------------------------------------------------------------------
// PUBLIC
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// GETTERS
//--------------------------------------------------------------------------

CharacterVector vHMM::getStateNames(void) const
{
    return m_StateNames;
}

NumericMatrix vHMM::getA(void) const
{    
    return m_A;
}

NumericMatrix vHMM::getB(void) const
{
    return m_B;
}

NumericVector vHMM::getPi(void) const
{
    return m_Pi;
}

//--------------------------------------------------------------------------
//  SETTERS
//
//  Every method has an error handler
//--------------------------------------------------------------------------

void vHMM::setStateNames(CharacterVector stateNames)
{
    if( stateNames.size() != m_N)
        Rf_error("The number of state names does not coincide with the one declared.");
    m_StateNames = CharacterVector(clone(stateNames));
}

void vHMM::setA(NumericMatrix A)
{    
    if(A.ncol() != m_N || A.nrow() != m_N)
        Rf_error("The transition matrix size is wrong");
    if(verifyMatrix(A) == false)
        Rf_error("The transition matrix is not normalized");
    m_A = A;       
}

void vHMM::setPi(NumericVector Pi)
{
    if(Pi.size() != m_N)
        Rf_error("The initial probability vector size is wrong");
    if(verifyVector(Pi) == false)        
        Rf_error("The initial probability vector is not normalized");
    m_Pi = Pi;
}

void vHMM::setParameters(NumericMatrix A, NumericMatrix B, NumericVector Pi)
{
    if(Pi.size() != m_N)
        Rf_error("The initial probability vector size is wrong");
    if(verifyVector(Pi) == false)
        Rf_error("The initial probability vector is not normalized");
    
    if(A.ncol() != m_N || A.nrow() != m_N)
        Rf_error("The transition matrix size is wrong");
    if(verifyMatrix(A) == false)
        Rf_error("The transition matrix is not normalized");  
      
    m_Pi = NumericVector(clone(Pi));
    m_A = NumericMatrix(clone(A));            
}

//  Not set. To be overriden by the inherited classes
void vHMM::setB(NumericMatrix B){}  

//--------------------------------------------------------------------------
//  MISCELLANEOUS
//--------------------------------------------------------------------------

//  Returns all the model parameterrs as a list
Rcpp::List vHMM::toList(void) const
{
    return List::create(
            Named("Model", "vHMM"),
            Named("StateNames", getStateNames() ),        
            Named("A", getA()),
            Named("B", getB()),
            Named("Pi", getPi() )
    );    
}

//  All the data expressed as discrete {0,1,2,...} is returned to categorical
//  S: State, O: Observation
CharacterVector vHMM::toName( IntegerVector index, char vectorName)
{
    unsigned int length, i; 
    length = index.size();
    CharacterVector names(length);

    switch(vectorName)
    {
        case 'S':
            for(i = 0; i < length; i++)
                names[i] = m_StateNames[index[i]];
        break;
    }

    return names;
}

//  Function used to print to console the model parameters
/*
ostream& vHMM::print(ostream& out) const
{
   //Nombres
    CharacterVector nombre = m_StateNames;
    unsigned int N = m_N;    

    for(unsigned int i = 0; i < N; i++ )
    {
        out << "Estado : " << nombre[i] << ", ";
    }
    out << "\n";

    //Matrices
    NumericVector Pi = m_Pi;
    NumericMatrix A = m_A;    

   for(unsigned int i = 0; i < N; i++ )
    {
        out << "Pi : " << Pi[i] << ",";
    }
    out << "\n";

    for(unsigned int i = 0; i < N; i++ )
    {
        for(unsigned int j = 0; j < N; j++ )
        {
            out << "A : " << A(i,j)<< ", ";
        }
        out << "\n";
    }

    return out ;
}

//  Function used to overload the << operator 
ostream& operator<< (ostream & out, const vHMM& data) 
{
    return data.print(out);
}
*/

//--------------------------------------------------------------------------
//  PROTECTED
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
//  MISCELLANEOUS
//--------------------------------------------------------------------------

//  Initialize the model with random parameters
//  Not set. To be overriden by the inherited classes
void vHMM::randomInit(){}

//  Verify if it is a stochastic matrix (row elements sum equals to 1)
bool vHMM::verifyMatrix(NumericMatrix matrix)
{
    bool flag = true;  
    double normalized;
    int i;
    
    //  Sum of each row must be 1
    for(i = 0 ; i < matrix.nrow() ; i++)
    {
        normalized = sum(matrix.row(i)); 
        if(normalized < 1 - EPSILON || normalized > 1 + EPSILON)
        {
            flag = false;
            break;
        } 
    }
    return flag ;    
}

//  Verify if the sum of all its elements sum 1
bool vHMM::verifyVector(NumericVector vector)
{
    bool flag = true;       

    double normalized = sum(vector); 
    if(normalized < 1 - EPSILON || normalized > 1 + EPSILON)
            flag = false;
 
    return flag;
}
