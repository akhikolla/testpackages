// Copyright (C) 2014 Mohammad H. Ferdosi
//
// HSPhase is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// HSPhase program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


/*
 * File:   hsp.h
 * Author: Mohammad H. Ferdosi
 *
 * Created on 13 October 2012, 6:14 PM
 */

#ifndef HSP_H
#define	HSP_H

//#include "R.h"
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <cmath>

#include "TypeConversion.h"
#define uint unsigned int

int c2rBlocks2(const int* matrix,int *nrow,const int *ncol, int* result);
int c2rBlocks(int const * matrix, int *nrow, int *ncol, int* result, int* MaxBlock);
int c2rStrandF(const int* matBlock, const int* matGenotype, const int* nrow, const int* ncol, double* result);
int frequencyVector(vector<int>::iterator strat, int seachNumber, vector<int>::iterator end);
int reverseConvert(int *col, int *nrow);
int reverseConvert(int *col,const int *nrow);
int memMaker(int *lastMemory, int* newCol, int *nrow);
int strandOrigin(int *col, int *nrow);
int strandOrigin(int *col, const int *nrow);
int switchDetector(int *Memory, int *tempCol2, int *nrow);
int switchDetector(int *Memory, int *tempCol2, const int *nrow);
int switchDetector(int *Memory, int *tempCol2, vector<int> &switches, int *nrow);
int switchDetector(int *Memory, int *tempCol2, vector<int> &switches, const int *nrow);
int phaseFunction(int const * genotypeMat, int const *nrow, int const *ncol, int const* blockMat, int const* sirePhasedMat, int* result);
int c2rRecombinations(const uint* matrix, const uint* nrow, const uint* ncol, const uint* method, double* result);
int recombinationFun(int const *matrix,int const * nrow,int const * ncol,int * result);
int phaseFunctionNoGenotype(int const *nrow, int const *ncol, int const* blockMat, int const* sirePhasedMat, int* result);

#endif	/* HSP_H */

