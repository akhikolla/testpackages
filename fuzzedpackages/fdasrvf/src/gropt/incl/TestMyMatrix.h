/*
This is the test file to test the wrapper functions, of lapack and blas, in MyMatrix.h and MyMatrix.cpp.

---- WH
*/

#ifndef TESTMYMATRIX_H
#define TESTMYMATRIX_H

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "Others/randgen.h"

/*Computational time*/
#include <ctime>

/*The global head file*/
#include "Others/def.h"

/*Test MyMatrix class*/
#include "Others/MyMatrix.h"

using namespace ROPTLIB;

void testEigenSymmetricM(void);
void testExpSymmetricM(void);
void testLogSymmetricM(void);
void testDSYL(void);

#endif // end of TESTSTIEBROCKETT_H
