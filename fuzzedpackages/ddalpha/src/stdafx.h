/*
  File:             stdafx.h
  Created by:       Oleksii Pokotylo
  First published:  28.02.2013
  Last revised:     28.02.2013
  
  Defines the Includes needed.
*/

#pragma once

#define BOOST_UBLAS_NO_STD_CERR

#include <time.h>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <vector>
#include <set>
#include <stdlib.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#ifndef _MSC_VER
#include <Rcpp.h> 
using namespace Rcpp;
#endif

using namespace std;

#include "DataStructures.h"
#include "Common.h"
#include "AlphaProcedure.h"
#include "TukeyDepth.h"
#include "HD.h"
#include "ZonoidDepth.h"
#include "Mahalanobis.h"
#include "SimplicialDepth.h"
#include "OjaDepth.h"
#include "Knn.h"
#include "Polynomial.h"
#include "PotentialDepth.h"
#include "ProjectionDepth.h"
#include "DKnn.h"
#include "LensDepth.h"
#include "BandDepth.h"

// global rEngine is defined in ddalpha.cpp, extern rEngine defined in stdafx.h
#define ran(x) rEngine()%x
#define setseed(x) rEngine.seed(x)

int random(int x);

