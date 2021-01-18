/**
 * Operations.h
 *   Purpose: Calculates matrix and vectors multiplications.
 * @authors Hmamouche Youssef
 * @date 2016
 **/

#ifndef OPERATEURS_H
#define OPERATEURS_H

#include <Rcpp.h>

#include<iostream>
#include "struct.h"
#include "exception.h"


using namespace std;
using namespace Struct;

template<typename T>
void show (const T & a)
{
    cout << a;
}

template<typename T>
void show (const std::vector<T> & vect)
{
    cout << "[";
    for (auto & val : vect)
    {
        show (val);
        if (val == vect. back ())
             cout << "]\n";
        else
            cout << " ";
    }

}


/*****************************/
template<typename T>
void show (const std::vector<std::vector<T>> & matrix)
{
    cout << "[";
    // Print the matrix
    for (auto & row : matrix)
    {
           show (row);
           if (row == matrix. back ())
               cout << "] \n";
           else
                cout << "\n";
    }
}
/**
    Compute a vector-vector multiplication.

    @param A the first vector.
    @param B the second vector.
    @param Res the vector where tu put the result.
*/
void MultCVDouble (const CVDouble & A, const CVDouble & B, CVDouble & Res);

/**
    Compute a matrix-vector multiplication.

    @param A the matrix.
    @param B the  vector.
    @param Res the vector where tu put the result.
*/
void MultCVDouble (const CMatDouble & A, const CVDouble & B, CVDouble & Res);


/**
    Compute a matrix-matrix multiplication.

    @param A the first matrix.
    @param B the second matrix.
    @param Res the matrix where tu put the result.
*/
void MultCVDouble (const CMatDouble & A, const CMatDouble & B, CMatDouble & Res);

#endif // OPERATEURS_H
