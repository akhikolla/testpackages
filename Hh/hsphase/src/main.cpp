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
 * File:   main.cpp
 * Author: Mohammad H. Ferdosi
 *
 * Created on 13 October 2012, 6:14 PM
 */

#include <cstdlib>
#include "block.h"
#include "block4Phase.h"
#include "hsp.h"
#include "swDetect.h"
#include "TypeConversion.h"
#include "diag.h"

using namespace std;


extern "C"
{

void hbphased(uint* matrix, uint *nrow, uint *ncol, uint* result,
		uint* siregenotype, uint *str)
{

	block4Phase haplotypeblockphased(matrix, nrow, ncol, result, siregenotype, str);

}
}
extern "C"
{

void hblock(int* matrix, int *nrow, int *ncol, int* result, int* MaxBlock)
{

	c2rBlocks(matrix, nrow, ncol, result, MaxBlock);

}
}

extern "C"
{

void ssp(const int* matBlock, const int* matGenotype, int * hetsite,
		const int* nrow, const int* ncol, double* result)
{

	c2rStrandF(matBlock, matGenotype, nrow, ncol, result);

}
}
extern "C"
{

void phase(int const * genotypeMat, int const *nrow, int const *ncol,
		int const* blockMat, int const* sirePhasedMat, int* result)
{

	phaseFunction(genotypeMat, nrow, ncol, blockMat, sirePhasedMat, result);

}
}
extern "C"
{

void phaseNogenotype(int const *nrow, int const *ncol,
		int const* blockMat, int const* sirePhasedMat, int* result)
{

	phaseFunctionNoGenotype(nrow, ncol, blockMat, sirePhasedMat, result);

}
}
extern "C"
{
void sw(int* matrix1, int *matrix2, int *nrow, int *ncol, int *result)
{

	swFun(matrix1, matrix2, nrow, ncol, result);

}
}
extern "C"
{

void bmh(int const * matrix, int* zeroFrq, int* oneFrq, int* twoFrq, int *nrow,
		int *ncol, int* result, int* sire, int* forwardVectorSize, bool *FP,
		int *nsap)
{

	vector<int> sireVec;

	for (int i = 0; i < (*ncol); i++)
	{
		if (sire[i] == 1)
			sireVec.push_back(i);
	}

	c2rphaseOPT(matrix, zeroFrq, oneFrq, twoFrq, nrow, ncol, result, sireVec,
			forwardVectorSize, FP, nsap);
}

}
extern "C"
{

void bmhr(int const * matrix, int* zeroFrq, int* oneFrq, int* twoFrq, int *nrow,
		int *ncol, int* result, int* sire, int* forwardVectorSize, bool *FP,
		int *nsap, int* recombinationMatrix)
{

	vector<int> sireVec;

	for (int i = 0; i < (*ncol); i++)
	{
		if (sire[i] == 1)
			sireVec.push_back(i);
	}

	c2rphaseOPTRecombination(matrix, zeroFrq, oneFrq, twoFrq, nrow, ncol, result, sireVec,
			forwardVectorSize, FP, nsap, recombinationMatrix);
}

}
extern "C"
{

void pm(uint* matrix, uint *nrow, uint *ncol, uint *method,double* result)
{

	c2rRecombinations(matrix, nrow, ncol,method , result);

}
}

extern "C"
{
void recombinations(int const * matrix, int const *nrow, int const *ncol,
		int * result)
{

	recombinationFun(matrix, nrow, ncol, result);

}

extern "C"
{
void phaseDiag(int const * matrix, int *nrow, int *ncol, int * result,
		int* sire,int* zeroFrq, int* twoFrq)
{
	vector<int> sireVec;

	for (int i = 0; i < (*ncol); i++)
	{
		if (sire[i] == 1)
			sireVec.push_back(i);
	}

	diagnostic(matrix, nrow, ncol, result, sireVec,zeroFrq,twoFrq);
}
}
}

