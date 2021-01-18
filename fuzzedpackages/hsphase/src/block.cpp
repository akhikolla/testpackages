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
 * File:   block.cpp
 * Author: Mohammad H. Ferdosi
 *
 * Created on 21 October 2012, 2:17 PM
 */
#include "block.h"
#include "hsp.h"
#include "memory.h"
#include "swDetect.h"

#include <R.h>

using namespace std;

block::block(int const **pRowsMat, int ** result, int * tbackwardmemory, vector<int> switchLocus,
		vector<int>::iterator iElement, int *nrow, vector<int> &hetsite)
{
	itspRowsMat_ = pRowsMat;
	itsResult_ = result;

	itsiElement_ = iElement;
	itsNrow_ = nrow;
	itsHetsite_ = hetsite;
	itsHetEnd_ = hetsite.end();
	itsBackwardMemoryMain_ = tbackwardmemory;
	itstbackwardmemoryBlock_ = new int[*nrow];
	itsValidateMem_ = new int[*nrow];
	itsValidateConst_ = new int[*nrow];

	fill_n(itstbackwardmemoryBlock_, *nrow, 0);
	fill_n(itsValidateMem_, *nrow, 0);
	fill_n(itsValidateConst_, *nrow, 0);

	for (int i = 0; i < *nrow; i++)
	{
		if (tbackwardmemory[i] == 3)
		{
			itstbackwardmemoryBlock_[i] = 0;
		}
		if (tbackwardmemory[i] == 4)
		{
			itstbackwardmemoryBlock_[i] = 2;
		}
		if (tbackwardmemory[i] == 9)
		{
			itstbackwardmemoryBlock_[i] = 1;
		}

	}

}

block::~block()
{
	delete[] itstbackwardmemoryBlock_;
	delete[] itsValidateMem_;
	delete[] itsValidateConst_;
}

int block::recombinationDetector(int ur, int ul, int dr, int dl)
{
	if (ur == 1 || ul == 1 || dr == 1 || dl == 1)
	{
		return (0);
	}
	int sum = ur + ul + dr + dl;

	return (((sum == 2) || (sum == 6)) ? -1 : 1);
}
/**
 *
 * @param windowsWidth hetsite in the sire in advance for validating the recombinations
 * @return
 */
int block::validitySwitch(int windowsWidth)
{
	vector<int>::iterator iElemnt;

	int *tempcol2 = new int[*itsNrow_];
	fill_n(tempcol2, *itsNrow_, 0);
	vector<int> SwitchLocusA;
	vector<int> SwitchLocusB;

	for (iElemnt = itsiElement_; (iElemnt != (itsiElement_ + windowsWidth)) && (iElemnt != itsHetEnd_); ++iElemnt)
	{
		for (int i = 0; i < *itsNrow_; i++)
		{
			*(*(itsResult_ + *iElemnt) + i) = 0;
			tempcol2[i] = *(*(itspRowsMat_ + *iElemnt) + i);
		}
		strandOrigin(tempcol2, itsNrow_);
		switchDetector(itsBackwardMemoryMain_, tempcol2, SwitchLocusA, itsNrow_);
		reverseConvert(tempcol2, itsNrow_);
		switchDetector(itsBackwardMemoryMain_, tempcol2, SwitchLocusB, itsNrow_);

		if ((SwitchLocusA.size() == SwitchLocusB.size()) && (SwitchLocusB.size() != 0 && iElemnt == itsiElement_))
		{

		}

		if (SwitchLocusA.size() < SwitchLocusB.size())
		{
			vector<int>::iterator iElement2;
			for (iElement2 = SwitchLocusA.begin(); iElement2 != SwitchLocusA.end(); ++iElement2)
			{
				*(*(itsResult_ + *iElemnt) + *iElement2) = -2;
			}
		}
		else if (SwitchLocusA.size() > SwitchLocusB.size())
		{
			vector<int>::iterator iElement2;
			for (iElement2 = SwitchLocusB.begin(); iElement2 != SwitchLocusB.end(); ++iElement2)
			{
				*(*(itsResult_ + *iElemnt) + *iElement2) = -2;
			}
		}
		else if ((SwitchLocusA.size() == SwitchLocusB.size()) && SwitchLocusB.size() != 0)
		{
			for (int k = 0; k < *itsNrow_; k++)
				*(*(itsResult_ + *iElemnt) + k) = -2;
		}

		for (int k = 0; k < *itsNrow_; k++)
		{
			if (*(*(itsResult_ + *iElemnt) + k) < -1 && *(itsValidateConst_ + k) == 0)
			{
				*(itsValidateMem_ + k) = *(itsValidateMem_ + k) + 1;
			}
			else if (*(*(itsResult_ + *iElemnt) + k) != -1)
			{
				*(itsValidateConst_ + k) = 1;
			}
		}

	}

	delete[] tempcol2;

	return (0);

}

int block::makeMemory(int * nsap)
{

	for (int k = 0; k < *itsNrow_; k++)
	{
		//K must change to catch SNP error (0 Not valid)
		if (*(itsValidateMem_ + k) > *nsap)
		{
			if (*(itsBackwardMemoryMain_ + k) == 3)
				*(itsBackwardMemoryMain_ + k) = 5;
			if (*(itsBackwardMemoryMain_ + k) == 4)
				*(itsBackwardMemoryMain_ + k) = 3;
			if (*(itsBackwardMemoryMain_ + k) == 5)
				*(itsBackwardMemoryMain_ + k) = 4;
		}

	}
	for (int i = 0; i < *itsNrow_; i++)
	{
		if (*(*(itspRowsMat_ + *(itsiElement_)) + i) != 1)
		{
			*(*(itsResult_ + *(itsiElement_)) + i) = *(itsBackwardMemoryMain_ + i);
		}
		else
		{

			*(*(itsResult_ + *(itsiElement_)) + i) = 1;
		}
	}

	return (0);
}

/**
 *
 * @param matrix Genotype matrix 0, 1, 2 -> aa, ab, bb
 * @param zeroFrq Vector of frequency of zero
 * @param oneFrq Vector of frequency of one
 * @param twoFrq Vector of frequency of two
 * @param nrow Number of individuals
 * @param ncol Number of SNPs
 * @param result A block matrix with gaps
 * @param hetSite Heterozygote sites of sire
 * @param forwardvectorsize The number of het to search in advance in the case of possible recombinations
 * @param excludeFP
 * @param nsap Number of SNPs per block
 * @return
 */
int c2rphaseOPT(int const * matrix, int const *zeroFrq, int const * oneFrq, int const* twoFrq, int* nrow, int* ncol,
		int* result, vector<int> &hetSite, int *forwardvectorsize /* = 0 */, bool *excludeFP, int *nsap)
{
	for (int i = 0; i < (*nrow) * (*ncol); i++)
	{
		result[i] = matrix[i];
	}

//For matrix
	int const**pRowsMat = new int const*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pRowsMat[i] = matrix + (*nrow) * i;
	}

// For result
	int **pRowsRes = new int*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pRowsRes[i] = result + (*nrow) * i;
	}

	int *tempCol1 = new int[*nrow];
	int *tempCol2 = new int[*nrow];
	int *tempCol2R = new int[*nrow];
	int *keeptempCol2 = new int[*nrow];
	vector<int> switchA, switchB;
	memoryCLS tbackwardMemory(nrow);

	for (int i = 0; i < *nrow; i++)
	{
		tempCol1[i] = pRowsMat[hetSite[0]][i];
	}

	strandOrigin(tempCol1, nrow);
	for (int j = 0; j < *nrow; j++)
		pRowsRes[hetSite[0]][j] = tempCol1[j];

// Initialise backward memory
	tbackwardMemory.memoryMaker(tempCol1);

// Main loop
	int n = 1;
	vector<int>::iterator iElement;

	for (iElement = hetSite.begin() + 1; iElement != hetSite.end(); ++iElement)
	{

		double zero = (*(twoFrq + *iElement) < *(zeroFrq + *iElement)) ? *(zeroFrq + *iElement) : *(twoFrq + *iElement);
		double two = (*(twoFrq + *iElement) > *(zeroFrq + *iElement)) ? *(zeroFrq + *iElement) : *(twoFrq + *iElement);

		if (*excludeFP)
		{

			if (((double) (two / zero)) < .4 || ((double) (zero / (*nrow)) > .8))
			{
				continue;
			}
		}

		for (int i = 0; i < *nrow; i++)
		{
			*(keeptempCol2 + i) = *(*(pRowsMat + *iElement) + i);
			*(tempCol2 + i) = *(*(pRowsMat + *iElement) + i);
			*(tempCol2R + i) = *(*(pRowsMat + *iElement) + i);
		}

		strandOrigin(tempCol2, nrow);
		strandOrigin(tempCol2R, nrow);
		reverseConvert(tempCol2R, nrow);
		switchDetector(tbackwardMemory.pMemory, tempCol2, switchA, nrow);
		switchDetector(tbackwardMemory.pMemory, tempCol2R, switchB, nrow);

		if (switchA.size() < switchB.size() && switchA.size() == 0)
		{
			for (int j = 0; j < (*nrow); j++)
			{
				if (tempCol2[j] == 3 || tempCol2[j] == 4)
				{
					*(*(pRowsRes + *iElement) + j) = *(tempCol2 + j);
				}
			}
			tbackwardMemory.memoryMaker(tempCol2);
		}
		else if (switchA.size() > switchB.size() && switchB.size() == 0)
		{
			for (int j = 0; j < (*nrow); j++)
			{
				if (tempCol2[j] == 3 || tempCol2[j] == 4)
				{
					*(*(pRowsRes + *iElement) + j) = *(tempCol2R + j);
				}
			}
			tbackwardMemory.memoryMaker(tempCol2R);
		}
		else
		{

			block mainBlock(pRowsMat, pRowsRes, tbackwardMemory.pMemory, switchA, iElement, nrow, hetSite);

			mainBlock.validitySwitch(*forwardvectorsize);

			mainBlock.makeMemory(nsap);
			int *backward = mainBlock.getItsBackwardMemoryMain_();
			for (int i = 0; i < *nrow; i++)
			{
				tbackwardMemory.pMemory[i] = *(backward + i);
			}

		}

		n = n + 1;
		switchA.clear();
		switchB.clear();
	}

	delete[] pRowsMat;
	delete[] pRowsRes;
	delete[] tempCol1;
	delete[] tempCol2;
	delete[] tempCol2R;
	delete[] keeptempCol2;

	return (0);
}

int c2rphaseOPTRecombination(int const * matrix, int const *zeroFrq, int const * oneFrq, int const* twoFrq, int* nrow,
		int* ncol, int* result, vector<int> &hetSite, int *forwardvectorsize /* = 0 */, bool *excludeFP, int *nsap,
		int* recombinationMatrix)
{
	//Rprintf("Start\n");
	for (int i = 0; i < (*nrow) * (*ncol); i++)
	{
		result[i] = matrix[i];
	}

//For matrix (genotype)
	int const**pRowsMat = new int const*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pRowsMat[i] = matrix + (*nrow) * i;
	}

// For result (blocks)
	int **pRowsRes = new int*[*ncol];

	for (int i = 0; i < *ncol; i++)
	{
		pRowsRes[i] = result + (*nrow) * i;
	}

	// For result (recombinations)
	int **pRowsRecombionations = new int*[*ncol - 1];

	for (int i = 0; i < *ncol - 1; i++)
	{
		pRowsRecombionations[i] = recombinationMatrix + (*nrow) * i;
	}

	/*for (int i = 0; i < *nrow; i++)
	 {
	 for (int j = 0; j < (*ncol) - 1; j++)
	 {
	 Rprintf("%d ", pRowsRecombionations[j][i]);
	 }
	 Rprintf("\n");
	 }
	 */
	int *tempCol1 = new int[*nrow];
	int *tempCol2 = new int[*nrow];
	int *tempCol2R = new int[*nrow];
	int *keeptempCol2 = new int[*nrow];
	vector<int> tempCrossOver;
	vector<int> switchA, switchB;
	memoryCLS tbackwardMemory(nrow);

	for (int i = 0; i < *nrow; i++)
	{
		tempCol1[i] = pRowsMat[hetSite[0]][i];
	}

	strandOrigin(tempCol1, nrow);
	for (int j = 0; j < *nrow; j++)
		pRowsRes[hetSite[0]][j] = tempCol1[j];

// Initialise backward memory
	tbackwardMemory.memoryMaker(tempCol1);

// Main loop
	int n = 1;
	vector<int>::iterator iElement;

	for (iElement = hetSite.begin() + 1; iElement != hetSite.end(); ++iElement)
	{

		double zero = (*(twoFrq + *iElement) < *(zeroFrq + *iElement)) ? *(zeroFrq + *iElement) : *(twoFrq + *iElement);
		double two = (*(twoFrq + *iElement) > *(zeroFrq + *iElement)) ? *(zeroFrq + *iElement) : *(twoFrq + *iElement);

		if (*excludeFP)
		{

			if (((double) (two / zero)) < .4 || ((double) (zero / (*nrow)) > .8))
			{
				continue;
			}
		}

		for (int i = 0; i < *nrow; i++)
		{
			*(keeptempCol2 + i) = *(*(pRowsMat + *iElement) + i);
			*(tempCol2 + i) = *(*(pRowsMat + *iElement) + i);
			*(tempCol2R + i) = *(*(pRowsMat + *iElement) + i);
		}

		strandOrigin(tempCol2, nrow);
		strandOrigin(tempCol2R, nrow);
		reverseConvert(tempCol2R, nrow);
		switchDetector(tbackwardMemory.pMemory, tempCol2, switchA, nrow);
		switchDetector(tbackwardMemory.pMemory, tempCol2R, switchB, nrow);

		if (switchA.size() < switchB.size() && switchA.size() == 0)
		{
			for (int j = 0; j < (*nrow); j++)
			{
				if (tempCol2[j] == 3 || tempCol2[j] == 4)
				{
					*(*(pRowsRes + *iElement) + j) = *(tempCol2 + j);
				}
			}
			tbackwardMemory.memoryMaker(tempCol2);

		}
		else if (switchA.size() > switchB.size() && switchB.size() == 0)
		{
			for (int j = 0; j < (*nrow); j++)
			{
				if (tempCol2[j] == 3 || tempCol2[j] == 4)
				{
					*(*(pRowsRes + *iElement) + j) = *(tempCol2R + j);
				}
			}
			tbackwardMemory.memoryMaker(tempCol2R);

		}
		else if (switchA.size() > 3 && switchB.size() > 3)
		{

			tempCrossOver.clear();
			for (int l = 0; l < *nrow; l++)
			{
				tempCrossOver.push_back(pRowsRecombionations[(*iElement) - 1][l]);
			}
			vector<int>::iterator lElement1;
			double aSum = 0, bSum = 0;
			for (lElement1 = switchA.begin(); lElement1 != switchA.end(); ++lElement1)
			{
				aSum += tempCrossOver[*lElement1];
			}
			aSum = aSum / switchA.size();
			vector<int>::iterator lElement2;
			for (lElement2 = switchB.begin(); lElement2 != switchB.end(); ++lElement2)
			{
				bSum += tempCrossOver[*lElement2];
			}
			bSum = bSum / switchB.size();

			if (bSum < aSum)
			{
				for (int j = 0; j < (*nrow); j++)
				{
					if (tempCol2[j] == 3 || tempCol2[j] == 4)
					{
						*(*(pRowsRes + *iElement) + j) = *(tempCol2 + j);
					}
				}
				tbackwardMemory.memoryMaker(tempCol2);

			}
			else if (aSum > bSum)
			{
				for (int j = 0; j < (*nrow); j++)
				{
					if (tempCol2[j] == 3 || tempCol2[j] == 4)
					{
						*(*(pRowsRes + *iElement) + j) = *(tempCol2R + j);
					}
				}
				tbackwardMemory.memoryMaker(tempCol2R);

			}


		}

		n = n + 1;
		switchA.clear();
		switchB.clear();
	}

	delete[] pRowsMat;
	delete[] pRowsRes;
	delete[] tempCol1;
	delete[] tempCol2;
	delete[] tempCol2R;
	delete[] keeptempCol2;

	return (0);
}

