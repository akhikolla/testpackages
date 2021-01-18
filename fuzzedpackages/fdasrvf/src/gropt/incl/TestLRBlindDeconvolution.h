/*
This is the test file to run the problem defined in LRBlindDeconvolution.h and LRBlindDeconvolution.cpp.

---- WH
*/

#ifndef TESTLRBLINDDECONVOLUTION_H
#define TESTLRBLINDDECONVOLUTION_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/LRBlindDeconvolution/LRBlindDeconvolution.h"
#include "Problems/EucBlindDeconvolution/EucBlindDeconvolution.h"
#include "Problems/SphereTxRQ/SphereTxRQ.h"
#include "Manifolds/LowRank/LowRank.h"
#include "Manifolds/LowRank/LowRankVariable.h"
#include "Manifolds/SphereTx/SphereTx.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

#include "Others/fftw/fftw3.h"

#ifdef ROPTLIB_WITH_FFTW

using namespace ROPTLIB;

//void testfft();
void testLRBlindDeconvolution(void);
void testLRBlindDeconvolutionSparse(void);

#endif
#endif
